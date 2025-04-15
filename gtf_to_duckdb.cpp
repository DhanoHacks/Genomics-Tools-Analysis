#include <duckdb.hpp>               // DuckDB C++ API header
#include <ray/api.h>                // Experimental Ray C++ API header
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// -----------------------------------------------------------------------------
// Thread-safe queue
template<typename T>
class TSQueue {
public:
    void push(const T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        queue_.push(item);
        cond_.notify_one();
    }
    
    bool try_pop(T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (queue_.empty())
            return false;
        item = queue_.front();
        queue_.pop();
        return true;
    }
    
    bool empty() {
        std::unique_lock<std::mutex> lock(mutex_);
        return queue_.empty();
    }
    
private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_;
};

// -----------------------------------------------------------------------------
// Process a single line from the GTF file and return an SQL INSERT statement.
std::string process_line(const std::string &line, const std::string &tablename) {
    if (line.empty() || line[0] == '#')
        return "";
    
    // Remove trailing newline (if any) and tokenize by tab.
    std::istringstream iss(line);
    std::vector<std::string> fields;
    std::string token;
    while (std::getline(iss, token, '\t')) {
        fields.push_back(token);
    }
    if (fields.size() < 9)
        throw std::runtime_error("Error: Line has less than 9 tab-separated fields.");
    
    // Process first eight fields.
    std::vector<std::pair<std::string, std::string>> data;
    std::vector<std::string> keys_list = {"Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"};
    for (size_t i = 0; i < keys_list.size(); i++) {
        data.push_back({keys_list[i], fields[i]});
    }
    
    // Adjust the Start field (subtract 1)
    try {
        int start_val = std::stoi(data[3].second);
        data[3].second = std::to_string(start_val - 1);
    } catch (...) {
        throw std::runtime_error("Error processing 'Start' value: " + data[3].second);
    }
    
    // Process attribute field (9th field)
    std::string attributes = fields[8];
    std::istringstream attr_stream(attributes);
    std::string attribute;
    while (std::getline(attr_stream, attribute, ';')) {
        // Trim leading white-space:
        size_t first = attribute.find_first_not_of(" \t\n\r");
        if (first == std::string::npos)
            continue;
        attribute = attribute.substr(first);
        // Split on first space
        size_t pos = attribute.find(' ');
        if (pos != std::string::npos) {
            std::string key = attribute.substr(0, pos);
            std::string value = attribute.substr(pos + 1);
            // Remove any quotation marks
            value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
            // Trim spaces from value.
            size_t val_start = value.find_first_not_of(" \t\n\r");
            size_t val_end = value.find_last_not_of(" \t\n\r");
            if (val_start != std::string::npos && val_end != std::string::npos)
                value = value.substr(val_start, val_end - val_start + 1);
            data.push_back({key, value});
        }
    }
    
    // Build the INSERT statement. Column names are double-quoted to avoid conflicts with reserved words.
    std::ostringstream oss;
    oss << "INSERT INTO " << tablename << " (";
    bool first_col = true;
    for (auto &p : data) {
        if (!first_col)
            oss << ", ";
        first_col = false;
        oss << "\"" << p.first << "\"";
    }
    oss << ") VALUES (";
    bool first_val = true;
    for (auto &p : data) {
        if (!first_val)
            oss << ", ";
        first_val = false;
        // Escape single quotes: replace each ' with ''
        std::string val = p.second;
        size_t pos = 0;
        while ((pos = val.find("'", pos)) != std::string::npos) {
            val.insert(pos, "'");
            pos += 2;
        }
        oss << "'" << val << "'";
    }
    oss << ");";
    return oss.str();
}

// -----------------------------------------------------------------------------
// Define a Ray remote function to process a batch of lines.
// The function takes a vector of lines and the table name, and returns a single string
// with the concatenated SQL INSERT statements.
std::string ProcessBatch(const std::vector<std::string> &lines_batch) {
    std::string table_name = "human"; // Default table name
    std::ostringstream batch_sql;
    for (const auto &line : lines_batch) {
        std::string stmt = process_line(line, table_name);
        if (!stmt.empty())
            batch_sql << stmt << "\n";
    }
    return batch_sql.str();
}
RAY_REMOTE(ProcessBatch);

// -----------------------------------------------------------------------------
// Main function
int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <db_name> <table_name> <input_file>" << std::endl;
        return 1;
    }
    std::string db_name = argv[1];
    std::string table_name = argv[2];
    std::string input_file = argv[3];

    // Initialize Ray.
    ray::Init();

    auto start = std::chrono::high_resolution_clock::now();

    // Connect to DuckDB.
    // For a file-based database, pass the filename; otherwise, use nullptr for an in-memory DB.
    duckdb::DuckDB db(db_name.c_str());
    duckdb::Connection conn(db);

    // Drop existing table.
    conn.Query("DROP TABLE IF EXISTS " + table_name + ";");

    // Create table.
    std::string create_query;
    if (table_name == "human") {
        create_query = "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                       "\"Chromosome\" TEXT, "
                       "\"Source\" TEXT, "
                       "\"Feature\" TEXT, "
                       "\"Start\" INTEGER, "
                       "\"End\" INTEGER, "
                       "\"Score\" TEXT, "
                       "\"Strand\" TEXT, "
                       "\"Frame\" TEXT, "
                       "\"gene_id\" TEXT, "
                       "\"gene_version\" TEXT, "
                       "\"gene_source\" TEXT, "
                       "\"gene_biotype\" TEXT, "
                       "\"transcript_id\" TEXT, "
                       "\"transcript_version\" TEXT, "
                       "\"transcript_source\" TEXT, "
                       "\"transcript_biotype\" TEXT, "
                       "\"tag\" TEXT, "
                       "\"transcript_support_level\" TEXT, "
                       "\"exon_number\" TEXT, "
                       "\"exon_id\" TEXT, "
                       "\"exon_version\" TEXT, "
                       "\"gene_name\" TEXT, "
                       "\"transcript_name\" TEXT, "
                       "\"protein_id\" TEXT, "
                       "\"protein_version\" TEXT, "
                       "\"ccds_id\" TEXT"
                       ");";
    } else if (table_name == "mouse") {
        create_query = "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                       "\"Chromosome\" TEXT, "
                       "\"Source\" TEXT, "
                       "\"Feature\" TEXT, "
                       "\"Start\" INTEGER, "
                       "\"End\" INTEGER, "
                       "\"Score\" TEXT, "
                       "\"Strand\" TEXT, "
                       "\"Frame\" TEXT, "
                       "\"gene_id\" TEXT, "
                       "\"gene_type\" TEXT, "
                       "\"gene_name\" TEXT, "
                       "\"level\" TEXT, "
                       "\"mgi_id\" TEXT, "
                       "\"havana_gene\" TEXT, "
                       "\"transcript_id\" TEXT, "
                       "\"transcript_type\" TEXT, "
                       "\"transcript_name\" TEXT, "
                       "\"transcript_support_level\" TEXT, "
                       "\"tag\" TEXT, "
                       "\"havana_transcript\" TEXT, "
                       "\"exon_number\" TEXT, "
                       "\"exon_id\" TEXT, "
                       "\"protein_id\" TEXT, "
                       "\"ccdsid\" TEXT, "
                       "\"ont\" TEXT"
                       ");";
    } else {
        std::cerr << "Unsupported table name: " << table_name << std::endl;
        return 1;
    }
    conn.Query(create_query);
    
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    std::cout << "Time taken to create table: " << elapsed.count() << "ms" << std::endl;

    // Set up a thread-safe queue for SQL batches and a flag for read completion.
    TSQueue<std::string> sql_queue;
    std::atomic<bool> read_done(false);

    // Writer thread function: pop SQL batches from the queue and execute each statement.
    auto writer_thread_func = [&]() {
        while (!read_done.load() || !sql_queue.empty()) {
            std::string batch_sql;
            if (sql_queue.try_pop(batch_sql)) {
                auto start_time = std::chrono::high_resolution_clock::now();
                std::istringstream iss(batch_sql);
                std::string stmt;
                while (std::getline(iss, stmt)) {
                    if (stmt.find_first_not_of(" \t") != std::string::npos) {
                        conn.Query(stmt);
                    }
                }
                auto commit_time = std::chrono::high_resolution_clock::now();
                auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(commit_time - start_time);
                std::cout << "Pushed a batch to DuckDB in " << diff.count() << "ms" << std::endl;
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    };

    // Reader thread function: read the input file, divide into batches, and process via Ray tasks.
    auto reader_thread_func = [&]() {
        std::ifstream infile(input_file);
        if (!infile) {
            std::cerr << "Failed to open input file: " << input_file << std::endl;
            return;
        }
        // Although the Python code uses a byte-based batch size,
        // here we simply read line-by-line and group into chunks.
        const size_t chunk_size = 8000; // lines per batch
        std::vector<std::string> all_lines;
        std::string line;
        while (std::getline(infile, line)) {
            all_lines.push_back(line);
            if (all_lines.size() >= chunk_size) {
                // Submit a Ray remote task for this batch.
                std::vector<ray::ObjectRef<std::string>> futures;
                futures.push_back(ray::Task(ProcessBatch).Remote(all_lines));
                std::string batch_sql;
                for (auto &future : futures) {
                    batch_sql += *future.Get();
                }
                sql_queue.push(batch_sql);
                all_lines.clear();
            }
        }
        // Process any remaining lines.
        if (!all_lines.empty()) {
            auto future = ray::Task(ProcessBatch).Remote(all_lines);
            std::string batch_sql = *future.Get();
            sql_queue.push(batch_sql);
        }
        read_done.store(true);
        std::cout << "Reader thread finished reading file" << std::endl;
    };

    // Start the reader and writer threads.
    std::thread reader_thread(reader_thread_func);
    std::thread writer_thread(writer_thread_func);
    reader_thread.join();
    writer_thread.join();

    // Create indices.
    start = std::chrono::high_resolution_clock::now();
    try {
        std::string index_query = "CREATE INDEX \"ix_" + table_name + "_Chromosome_Strand\" ON \"" + table_name + "\" (\"Chromosome\", \"Strand\");";
        conn.Query(index_query);
    } catch (std::exception &e) {
        std::cout << "Warning: Skipped index creation due to error: " << e.what() << std::endl;
    }
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    std::cout << "Time taken to create indices: " << elapsed.count() << "ms" << std::endl;

    // For file-based databases (db_name containing ".ddb"), the connection will be saved when closed.
    // DuckDB cleans up when objects go out of scope.
    
    // Shutdown Ray.
    ray::Shutdown();
    return 0;
}
