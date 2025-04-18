#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sqlite3.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
// add library for timing purposes (ms)
#include <chrono>
using namespace std;

mutex gtf_mtx, sql_mutex;
condition_variable cv_write;
int num_finished = 0, num_threads;
bool writer_should_exit = false;
int batch_len = 200000;

std::string processLine(const std::string &line, string tablename) {
    std::istringstream lineStream(line);
    std::string token;
    std::vector<std::string> fields;

    while (std::getline(lineStream, token, '\t')) {
        fields.push_back(token);
    }

    if (fields.size() < 9) {
        throw std::runtime_error("Error: Line has less than 9 tab-separated fields.");
    }

    std::unordered_map<std::string, std::string> data;
    std::vector<std::string> keys = {"Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"};

    for (size_t i = 0; i < 8; ++i) {
        data[keys[i]] = fields[i];
    }
    // reduce End by 1 keeping the datatype as string
    data["Start"] = std::to_string(std::stoi(data["Start"]) - 1);

    std::string attributes = fields[8];
    std::istringstream attrStream(attributes);
    std::string attribute;

    while (std::getline(attrStream, attribute, ';')) {
        attribute = std::string(std::find_if(attribute.begin(), attribute.end(), [](char c) { return !std::isspace(c); }), attribute.end());

        size_t pos = attribute.find(' ');
        if (pos != std::string::npos) {
            std::string key = attribute.substr(0, pos);
            std::string value = attribute.substr(pos + 1);
            value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
            data[key] = value;
        }
    }
    std::string sql = "INSERT INTO " + tablename + " (";
    std::string keysStr, valuesStr;

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
        // store time taken to read data
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
        
        // store time taken to process data
        start = chrono::high_resolution_clock::now();
        for (auto &line : lines) {
            linedata += processLine(line, tablename);
        }
        lines.clear();
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time taken to process data: " << duration.count() << "ms" << endl;

        // unique_lock<mutex> lock(sql_mutex);
        // cv_read.wait(lock, [&]() { return *common_num_rows < batch_len || writer_should_exit; });
        sql_mutex.lock();

        *commonlinedata += linedata;
        *common_num_rows += num_rows;
        sql_mutex.unlock();
        num_rows = 0;
        linedata.clear();

        if (*common_num_rows >= batch_len) {
            cv_write.notify_one();
        }
    }

    // unique_lock<mutex> lock(sql_mutex);
    sql_mutex.lock();
    if (!linedata.empty()) {
        *commonlinedata += linedata;
    }
    num_finished++;
    if (num_finished == num_threads) {
        writer_should_exit = true;
        sql_mutex.unlock();
        cv_write.notify_one();
    }
    else {
        sql_mutex.unlock();
    }
}

void mywritethread(sqlite3 **db, string *commonlinedata, int *common_num_rows) {
    while (true) {
        unique_lock<mutex> lock(sql_mutex);
        cv_write.wait(lock, [&]() { return *common_num_rows >= batch_len || writer_should_exit; });

        if (!commonlinedata->empty()) {
            // store time taken to move data
            auto start = chrono::high_resolution_clock::now();
            string data_to_write = move(*commonlinedata);
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            cout << "Time taken to move data: " << duration.count() << "ms" << endl;
            int num_rows = *common_num_rows;
            commonlinedata->clear();
            *common_num_rows = 0;
            lock.unlock();

            // store time taken to push data to db
            start = chrono::high_resolution_clock::now();
            sqlite3_exec(*db, data_to_write.c_str(), NULL, NULL, NULL);
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
    sqlite3 *db;
    // database name is first argument
    sqlite3_open(argc[1], &db);
    // table name is second argument
    string tablename = argc[2];
    string query = "DROP TABLE IF EXISTS " + tablename + ";";
    sqlite3_exec(db, query.c_str(), NULL, NULL, NULL);
    // check if tablename is equal to "human" or "mouse"
    if(tablename == "human"){
        query = "CREATE TABLE IF NOT EXISTS " + tablename + " (\"Chromosome\" TEXT, \"Source\" TEXT, \"Feature\" TEXT, \"Start\" INTEGER, \"End\" INTEGER, \"Score\" TEXT, \"Strand\" TEXT, \"Frame\" TEXT, \"gene_id\" TEXT, \"gene_version\" TEXT, \"gene_source\" TEXT, \"gene_biotype\" TEXT, \"transcript_id\" TEXT, \"transcript_version\" TEXT, \"transcript_source\" TEXT, \"transcript_biotype\" TEXT, \"tag\" TEXT, \"transcript_support_level\" TEXT, \"exon_number\" TEXT, \"exon_id\" TEXT, \"exon_version\" TEXT, \"gene_name\" TEXT, \"transcript_name\" TEXT, \"protein_id\" TEXT, \"protein_version\" TEXT, \"ccds_id\" TEXT);";
    }
    else if(tablename == "mouse"){
        query = "CREATE TABLE IF NOT EXISTS " + tablename + " (\"Chromosome\" TEXT, \"Source\" TEXT, \"Feature\" TEXT, \"Start\" INTEGER, \"End\" INTEGER, \"Score\" TEXT, \"Strand\" TEXT, \"Frame\" TEXT, \"gene_id\" TEXT, \"gene_type\" TEXT, \"gene_name\" TEXT, \"level\" TEXT, \"mgi_id\" TEXT, \"havana_gene\" TEXT, \"transcript_id\" TEXT, \"transcript_type\" TEXT, \"transcript_name\" TEXT, \"transcript_support_level\" TEXT, \"tag\" TEXT, \"havana_transcript\" TEXT, \"exon_number\" TEXT, \"exon_id\" TEXT, \"protein_id\" TEXT, \"ccdsid\" TEXT, \"ont\" TEXT);";
    }
    sqlite3_exec(db, query.c_str(), NULL, NULL, NULL);
    sqlite3_exec(db, "PRAGMA journal_mode = OFF; PRAGMA synchronous = 0; PRAGMA cache_size = 1000000; PRAGMA locking_mode = EXCLUSIVE; PRAGMA temp_store = MEMORY;", NULL, NULL, NULL);

    // input file is third argument
    ifstream inputFile(argc[3]);
    string *commonlinedata = new string();
    int *common_num_rows = new int(0);
    vector<thread> readers;

    // number of threads is fourth argument
    num_threads = stoi(argc[4]);
    for (int i = 0; i < num_threads; ++i) {
        readers.emplace_back(myreadthread, &inputFile, commonlinedata, common_num_rows, tablename);
    }

    thread writer(mywritethread, &db, commonlinedata, common_num_rows);

    for (auto &t : readers) {
        t.join();
        // cout << "Reader thread finished execution" << endl;
    }

    writer.join();
    // cout << "Write Thread finished execution" << endl;

    // store time taken to create indices
    auto start = chrono::high_resolution_clock::now();
    query = "CREATE INDEX \"ix_" + tablename + "_index\" ON \"" + tablename + "\" (\"index\");";
    sqlite3_exec(db, query.c_str(), NULL, NULL, NULL);
    // create index on chromosome,strand
    query = "CREATE INDEX \"ix_" + tablename + "_Chromosome_Strand\" ON \"" + tablename + "\" (\"Chromosome\", \"Strand\");";
    sqlite3_exec(db, query.c_str(), NULL, NULL, NULL);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Time taken to create indices: " << duration.count() << "ms" << endl;
    // only close if database location is not in memory, i.e. if it ends with .sqlite3 then close it
    // check if argc[1] ends with .sqlite3
    string db_name = argc[1];
    if(db_name.substr(db_name.find_last_of(".") + 1) == "sqlite3"){
        sqlite3_close(db);
    }
}

extern "C" {
    void run_gtf_to_sql(const char* db_name, const char* table_name, const char* input_file, int num_threads) {
        // Convert arguments to argc/argv format
        const char* argv[] = {"gtf_to_sql", db_name, table_name, input_file, to_string(num_threads).c_str()};
        int argc = 5;

        // Call the main function
        main(argc, const_cast<char**>(argv));
    }
}