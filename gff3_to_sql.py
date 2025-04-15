import sqlite3
import sys
import time
import ray
import threading
import queue

def process_line(line, tablename):
    """
    Process a single GTF file line and return an SQL INSERT statement.
    
    The function splits the line on tab characters,
    processes the first 8 fields and subtracts 1 from the Start field.
    Then it parses the attribute field (the ninth column) by splitting on ';'
    and further splitting each attribute on the first space.
    It finally builds an INSERT statement inserting the collected key/value pairs.
    """
    # dont process comment lines i.e lines starting with '#'
    if line.startswith('#'):
        return None
    line = line.rstrip("\n")
    # Strip the trailing newline and split by tab:
    fields = line.strip().split('\t')
    if len(fields) < 9:
        raise ValueError("Error: Line has less than 9 tab-separated fields.")
    
    # Prepare the first columns in order:
    keys_list = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"]
    data = {}
    for i, key in enumerate(keys_list):
        data[key] = fields[i]
    
    # Adjust the Start field (subtract 1)
    try:
        data["Start"] = str(int(data["Start"]) - 1)
    except Exception as e:
        raise ValueError(f"Error processing 'Start' value: {data['Start']}") from e
    
    # Process the attribute field (9th field)
    attributes = fields[8]
    for attribute in attributes.split(';'):
        attribute = attribute.lstrip()  # Remove leading whitespace
        if not attribute:
            continue  # Skip empty tokens
        # Split attribute into key and value at first space
        parts = attribute.split('=', 1)
        if len(parts) == 2:
            key, value = parts
            # Remove any quotation marks and extra whitespace from the value
            value = value.replace('"', '').strip()
            data[key] = value

    # Build the INSERT statement.
    # Note: This builds the SQL command as a string, similar to the C++ code.
    # In production code you would generally want to use parameterized queries.
    columns = ", ".join(data.keys())
    # Escape single quotes in values by doubling them up:
    values = ", ".join("'" + v.replace("'", "''") + "'" for v in data.values())
    sql = f"INSERT INTO {tablename} ({columns}) VALUES ({values});"
    return sql

@ray.remote
def process_batch(lines_batch, table_name):
    # Process each line in the batch
    sql_statements = ""
    for line in lines_batch:
        # You can call your original process_line function here.
        # Or, directly implement the processing logic.
        sql_statement = process_line(line, table_name)
        if sql_statement:
            sql_statements += sql_statement + "\n"
    # Return the accumulated SQL statements for this batch
    return sql_statements

def run_gff3_to_sql(db_name, table_name, input_file, batch_size=20000000, chunk_size=4000):
    """
    Main function to process a GTF file and write its data into an SQLite database.
    
    Parameters:
        db_name    : Name (path) of the SQLite database file.
        table_name : The name of the table ("human" or "mouse").
        input_file : GTF file to process.
        batch_size : Number of bytes/characters to read from the file at once. Default is 1000000.
        chunk_size : Number of lines to process in each batch. Default is 8000.
    """
    # Connect to the SQLite database.
    # Set check_same_thread=False so that the connection can be used by multiple threads.
    start_time = time.time()
    conn = sqlite3.connect(db_name, check_same_thread=False)
    cur = conn.cursor()

    # Drop the table if it exists.
    drop_query = f"DROP TABLE IF EXISTS {table_name};"
    cur.execute(drop_query)
    
    # Create table based on the table name.
    if table_name == "human":
        create_query = f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            "Chromosome" TEXT,
            "Source" TEXT,
            "Feature" TEXT,
            "Start" INTEGER,
            "End" INTEGER,
            "Score" TEXT,
            "Strand" TEXT,
            "Frame" TEXT,
            "gene_id" TEXT,
            "gene_version" TEXT,
            "gene_source" TEXT,
            "gene_biotype" TEXT,
            "transcript_id" TEXT,
            "transcript_version" TEXT,
            "transcript_source" TEXT,
            "transcript_biotype" TEXT,
            "tag" TEXT,
            "transcript_support_level" TEXT,
            "exon_number" TEXT,
            "exon_id" TEXT,
            "exon_version" TEXT,
            "gene_name" TEXT,
            "transcript_name" TEXT,
            "protein_id" TEXT,
            "protein_version" TEXT,
            "ccds_id" TEXT
        );
        '''
    elif table_name == "mouse":
        create_query = f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            "Chromosome" TEXT,
            "Source" TEXT,
            "Feature" TEXT,
            "Start" INTEGER,
            "End" INTEGER,
            "Score" TEXT,
            "Strand" TEXT,
            "Frame" TEXT,
            "gene_id" TEXT,
            "gene_type" TEXT,
            "gene_name" TEXT,
            "level" TEXT,
            "mgi_id" TEXT,
            "havana_gene" TEXT,
            "transcript_id" TEXT,
            "transcript_type" TEXT,
            "transcript_name" TEXT,
            "transcript_support_level" TEXT,
            "tag" TEXT,
            "havana_transcript" TEXT,
            "exon_number" TEXT,
            "exon_id" TEXT,
            "protein_id" TEXT,
            "ccdsid" TEXT,
            "ont" TEXT
        );
        '''
    elif table_name == "arabidopsis":
        create_query = f'''
            CREATE TABLE IF NOT EXISTS {table_name} (
                "Chromosome" TEXT,
                "Source" TEXT,
                "Feature" TEXT,
                "Start" INTEGER,
                "End" INTEGER,
                "Score" TEXT,
                "Strand" TEXT,
                "Frame" TEXT,
                "ID" TEXT,
                "Alias" TEXT,
                "Name" TEXT,
                "biotype" TEXT,
                "description" TEXT,
                "gene_id" TEXT,
                "logic_name" TEXT,
                "Parent" TEXT,
                "tag" TEXT,
                "transcript_id" TEXT,
                "constitutive" TEXT,
                "ensembl_end_phase" TEXT,
                "ensembl_phase" TEXT,
                "exon_id" TEXT,
                "rank" TEXT,
                "protein_id" TEXT,
                "Is_circular" TEXT
            );
        '''
    else:
        raise ValueError(f"Unsupported table name: {table_name}")

    cur.execute(create_query)
    
    # Set PRAGMA options.
    cur.executescript("""
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = 0;
        PRAGMA cache_size = 1000000;
        PRAGMA locking_mode = EXCLUSIVE;
        PRAGMA temp_store = MEMORY;
    """)
    conn.commit()
    elapsed = time.time() - start_time
    print(f"Time taken to create table and set PRAGMA options: {elapsed*1000:.0f}ms")

    # Create a shared queue for SQL batches and an event to signal when reading is done.
    sql_queue = queue.Queue(maxsize=1)
    read_done = threading.Event()

    # Writer thread function: waits for SQL batches and executes them.
    def writer_thread_func():
        while not (read_done.is_set() and sql_queue.empty()):
            try:
                batch_sql_statements = sql_queue.get(timeout=1)
            except queue.Empty:
                continue
            start_time = time.time()
            cur.executescript(batch_sql_statements)
            conn.commit()
            elapsed = time.time() - start_time
            print(f"Pushed a batch to SQLite in {elapsed*1000:.0f}ms")

    # Reader thread function: reads the file and puts processed batches into the queue.
    def reader_thread_func():
        with open(input_file, "r") as f:
            while True:
                # Read a chunk of lines from the file
                start_time = time.time()
                lines = f.readlines(batch_size)
                if not lines:
                    break
                elapsed = time.time() - start_time
                print(f"Read {len(lines)} lines from file in {elapsed*1000:.0f}ms")
                
                # Process the lines in batches using Ray
                start_time = time.time()
                futures = []
                for i in range(0, len(lines), chunk_size):
                    # Each batch processes chunk_size lines at once.
                    lines_batch = lines[i:i+chunk_size]
                    futures.append(process_batch.remote(lines_batch, table_name))
                elapsed = time.time() - start_time
                print(f"Submitted {len(lines)} lines to ray in {elapsed*1000:.0f}ms")
                
                # Wait for all futures to complete and collect the results.
                start_time = time.time()
                sql_statements = ray.get(futures)
                elapsed = time.time() - start_time
                print(f"Processed {len(lines)} lines in {elapsed*1000:.0f}ms")
                
                # Flatten the list of SQL statements into a single string.
                start_time = time.time()
                batch_sql_statements = "".join(sql_statements)
                elapsed = time.time() - start_time
                print(f"Flattened SQL statements in {elapsed*1000:.0f}ms")
                
                # Push the batch of SQL statements to the queue
                sql_queue.put(batch_sql_statements)
        read_done.set()
        print("Reader thread finished reading file")

    # Create and start the reader and writer threads.
    reader_thread = threading.Thread(target=reader_thread_func)
    writer_thread = threading.Thread(target=writer_thread_func)
    reader_thread.start()
    writer_thread.start()

    # Wait for both threads to finish
    reader_thread.join()
    writer_thread.join()

    # Create indices and time the operation.
    start_time = time.time()
    index_query = f'CREATE INDEX "ix_{table_name}_index" ON "{table_name}" ("index");'
    cur.execute(index_query)
    index_query2 = f'CREATE INDEX "ix_{table_name}_Chromosome_Strand" ON "{table_name}" ("Chromosome", "Strand");'
    cur.execute(index_query2)
    conn.commit()
    elapsed = time.time() - start_time
    print(f"Time taken to create indices: {elapsed*1000:.0f}ms")

    # Close the connection only if the database is file-based.
    if db_name.endswith(".sqlite3"):
        conn.close()

def main():
    if len(sys.argv) != 4:
        print("Usage: python gtf_to_sql.py <db_name> <table_name> <input_file>")
        sys.exit(1)

    db_name = sys.argv[1]
    table_name = sys.argv[2]
    input_file = sys.argv[3]
    run_gff3_to_sql(db_name, table_name, input_file)

if __name__ == "__main__":
    start_time = time.time()
    ray.init()
    elapsed = time.time() - start_time
    print(f"Ray initialized in {elapsed*1000:.0f}ms")
    # Run the main function
    start_time = time.time()
    main()
    elapsed = time.time() - start_time
    print(f"Main function completed in {elapsed*1000:.0f}ms")
    # Shutdown Ray
    start_time = time.time()
    ray.shutdown()
    elapsed = time.time() - start_time
    print(f"Ray shutdown in {elapsed*1000:.0f}ms")
