#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <duckdb.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
using namespace std;

mutex gtf_mtx, sql_mutex;
condition_variable cv_read, cv_write;
int num_finished = 0, num_threads;
bool writer_should_exit = false;
int batch_len = 200000;

string processLine(const string &line, string tablename) {
    istringstream lineStream(line);
    string token;
    vector<string> fields;

    while (getline(lineStream, token, '\t')) {
        fields.push_back(token);
    }

    if (fields.size() < 9) {
        throw runtime_error("Error: Line has less than 9 tab-separated fields.");
    }

    unordered_map<string, string> data;
    vector<string> keys = {"Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"};

    for (size_t i = 0; i < 8; ++i) {
        data[keys[i]] = fields[i];
    }
    data["Start"] = to_string(stoi(data["Start"]) - 1);

    string attributes = fields[8];
    istringstream attrStream(attributes);
    string attribute;

    while (getline(attrStream, attribute, ';')) {
        attribute = string(find_if(attribute.begin(), attribute.end(), [](char c) { return !isspace(c); }), attribute.end());

        size_t pos = attribute.find(' ');
        if (pos != string::npos) {
            string key = attribute.substr(0, pos);
            string value = attribute.substr(pos + 1);
            value.erase(remove(value.begin(), value.end(), '"'), value.end());
            data[key] = value;
        }
    }
    string sql = "INSERT INTO " + tablename + " (";
    string keysStr, valuesStr;

    for (const auto &pair : data) {
        if (!keysStr.empty()) {
            keysStr += ", ";
            valuesStr += ", ";
        }
        keysStr += pair.first;
        valuesStr += "'" + pair.second + "'";
    }

    sql += keysStr + ") VALUES (" + valuesStr + ");";
    return sql;
}

void myreadthread(ifstream *inputFile, string *commonlinedata, int *common_num_rows, string tablename) {
    string line, linedata = "";
    vector<string> lines;
    int num_rows = 0;
    int linecount = 0;
    bool eof = false;

    while (!eof) {
        gtf_mtx.lock();
        auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < 10000; i++) {
            if (!getline(*inputFile, line)) {
                eof = true;
                break;
            }
            if (line[0] == '#') {
                continue;
            }
            lines.push_back(line);
            num_rows++;
        }
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time taken to read data: " << duration.count() << "ms" << endl;
        gtf_mtx.unlock();

        linecount += lines.size();
        if (linecount % 100000 == 0) {
            cout << "Processed " << linecount << " lines" << endl;
        }
        
        start = chrono::high_resolution_clock::now();
        for (auto &line : lines) {
            linedata += processLine(line, tablename);
        }
        lines.clear();
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time taken to process data: " << duration.count() << "ms" << endl;

        unique_lock<mutex> lock(sql_mutex);
        cv_read.wait(lock, [&]() { return *common_num_rows < batch_len || writer_should_exit; });

        *commonlinedata += linedata;
        *common_num_rows += num_rows;
        num_rows = 0;
        linedata.clear();

        if (*common_num_rows >= batch_len) {
            cv_write.notify_one();
        }
    }

    unique_lock<mutex> lock(sql_mutex);
    if (!linedata.empty()) {
        *commonlinedata += linedata;
    }
    num_finished++;
    if (num_finished == num_threads) {
        writer_should_exit = true;
        cv_write.notify_one();
    }
}

void mywritethread(duckdb_connection *conn, string *commonlinedata, int *common_num_rows) {
    while (true) {
        unique_lock<mutex> lock(sql_mutex);
        cv_write.wait(lock, [&]() { return *common_num_rows >= batch_len || writer_should_exit; });

        if (!commonlinedata->empty()) {
            auto start = chrono::high_resolution_clock::now();
            string data_to_write = move(*commonlinedata);
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            cout << "Time taken to move data: " << duration.count() << "ms" << endl;
            int num_rows = *common_num_rows;
            commonlinedata->clear();
            *common_num_rows = 0;
            cv_read.notify_all();

            start = chrono::high_resolution_clock::now();
            duckdb_query(*conn, data_to_write.c_str(), nullptr);
            end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            cout << "Pushed " << num_rows << " lines to db in time " << duration.count() << "ms" << endl;
        }

        if (writer_should_exit && commonlinedata->empty()) {
            break;
        }
    }
}

int main(int argv, char **argc) {
    duckdb_database db;
    duckdb_connection conn;
    
    // Open database
    if (duckdb_open(argc[1], &db) != DuckDBSuccess) {
        cerr << "Failed to open database" << endl;
        return 1;
    }
    if (duckdb_connect(db, &conn) != DuckDBSuccess) {
        cerr << "Failed to connect to database" << endl;
        return 1;
    }

    string tablename = argc[2];
    string query = "DROP TABLE IF EXISTS " + tablename + ";";
    duckdb_query(conn, query.c_str(), nullptr);
    
    if(tablename == "human"){
        query = "CREATE TABLE IF NOT EXISTS " + tablename + " (\"Chromosome\" TEXT, \"Source\" TEXT, \"Feature\" TEXT, \"Start\" INTEGER, \"End\" INTEGER, \"Score\" TEXT, \"Strand\" TEXT, \"Frame\" TEXT, \"gene_id\" TEXT, \"gene_version\" TEXT, \"gene_source\" TEXT, \"gene_biotype\" TEXT, \"transcript_id\" TEXT, \"transcript_version\" TEXT, \"transcript_source\" TEXT, \"transcript_biotype\" TEXT, \"tag\" TEXT, \"transcript_support_level\" TEXT, \"exon_number\" TEXT, \"exon_id\" TEXT, \"exon_version\" TEXT, \"gene_name\" TEXT, \"transcript_name\" TEXT, \"protein_id\" TEXT, \"protein_version\" TEXT, \"ccds_id\" TEXT);";
    }
    else if(tablename == "mouse"){
        query = "CREATE TABLE IF NOT EXISTS " + tablename + " (\"Chromosome\" TEXT, \"Source\" TEXT, \"Feature\" TEXT, \"Start\" INTEGER, \"End\" INTEGER, \"Score\" TEXT, \"Strand\" TEXT, \"Frame\" TEXT, \"gene_id\" TEXT, \"gene_type\" TEXT, \"gene_name\" TEXT, \"level\" TEXT, \"mgi_id\" TEXT, \"havana_gene\" TEXT, \"transcript_id\" TEXT, \"transcript_type\" TEXT, \"transcript_name\" TEXT, \"transcript_support_level\" TEXT, \"tag\" TEXT, \"havana_transcript\" TEXT, \"exon_number\" TEXT, \"exon_id\" TEXT, \"protein_id\" TEXT, \"ccdsid\" TEXT, \"ont\" TEXT);";
    }
    duckdb_query(conn, query.c_str(), nullptr);

    ifstream inputFile(argc[3]);
    string *commonlinedata = new string();
    int *common_num_rows = new int(0);
    vector<thread> readers;

    num_threads = stoi(argc[4]);
    for (int i = 0; i < num_threads; ++i) {
        readers.emplace_back(myreadthread, &inputFile, commonlinedata, common_num_rows, tablename);
    }

    thread writer(mywritethread, &conn, commonlinedata, common_num_rows);

    for (auto &t : readers) {
        t.join();
    }

    writer.join();

    auto start = chrono::high_resolution_clock::now();
    query = "CREATE INDEX \"ix_" + tablename + "_index\" ON \"" + tablename + "\" (\"index\");";
    duckdb_query(conn, query.c_str(), nullptr);
    query = "CREATE INDEX \"ix_" + tablename + "_Chromosome_Strand\" ON \"" + tablename + "\" (\"Chromosome\", \"Strand\");";
    duckdb_query(conn, query.c_str(), nullptr);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Time taken to create indices: " << duration.count() << "ms" << endl;

    duckdb_disconnect(&conn);
    duckdb_close(&db);
}

extern "C" {
    void run_gtf_to_sql(const char* db_name, const char* table_name, const char* input_file, int num_threads) {
        const char* argv[] = {"gtf_to_sql", db_name, table_name, input_file, to_string(num_threads).c_str()};
        int argc = 5;
        main(argc, const_cast<char**>(argv));
    }
}