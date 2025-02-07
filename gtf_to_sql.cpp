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
using namespace std;

mutex gtf_mtx, sql_mutex;
condition_variable cv_read, cv_write;
int num_finished = 0, num_threads = 8;
bool writer_should_exit = false;

std::string processLine(const std::string &line) {
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

    std::string sql = "INSERT INTO human2(";
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

void myreadthread(ifstream *inputFile, string *commonlinedata) {
    string line, linedata = "";
    vector<string> lines;
    int linecount = 0;
    bool eof = false;

    while (!eof) {
        gtf_mtx.lock();
        for (int i = 0; i < 10; i++) {
            if (!getline(*inputFile, line)) {
                eof = true;
                break;
            }
            if (line[0] == '#') {
                continue;
            }
            lines.push_back(line);
        }
        gtf_mtx.unlock();

        linecount += lines.size();
        if (linecount % 100000 == 0) {
            // cout << "Processed " << linecount << " lines" << endl;
        }

        for (auto &line : lines) {
            linedata += processLine(line);
        }
        lines.clear();

        unique_lock<mutex> lock(sql_mutex);
        cv_read.wait(lock, [&]() { return commonlinedata->length() < 80 || writer_should_exit; });

        *commonlinedata += linedata;
        linedata.clear();

        if (commonlinedata->length() >= 80) {
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

void mywritethread(sqlite3 **db, string *commonlinedata) {
    while (true) {
        unique_lock<mutex> lock(sql_mutex);
        cv_write.wait(lock, [&]() { return commonlinedata->length() >= 80 || writer_should_exit; });

        if (!commonlinedata->empty()) {
            string data_to_write = move(*commonlinedata);
            commonlinedata->clear();
            lock.unlock();

            sqlite3_exec(*db, data_to_write.c_str(), NULL, NULL, NULL);
            // cout << "Pushed " << data_to_write.length() << " bytes to db" << endl;

            lock.lock();
            cv_read.notify_all();
        }

        if (writer_should_exit && commonlinedata->empty()) {
            break;
        }
    }
}

int main() {
    sqlite3 *db;
    sqlite3_open(":memory:", &db);
    sqlite3_exec(db, "DROP TABLE IF EXISTS human2;", NULL, NULL, NULL);
    sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS \"human2\" (\"Chromosome\" TEXT, \"Source\" TEXT, \"Feature\" TEXT, \"Start\" INTEGER, \"End\" INTEGER, \"Score\" TEXT, \"Strand\" TEXT, \"Frame\" TEXT, \"gene_id\" TEXT, \"gene_version\" TEXT, \"gene_source\" TEXT, \"gene_biotype\" TEXT, \"transcript_id\" TEXT, \"transcript_version\" TEXT, \"transcript_source\" TEXT, \"transcript_biotype\" TEXT, \"tag\" TEXT, \"transcript_support_level\" TEXT, \"exon_number\" TEXT, \"exon_id\" TEXT, \"exon_version\" TEXT, \"gene_name\" TEXT, \"transcript_name\" TEXT, \"protein_id\" TEXT, \"protein_version\" TEXT, \"ccds_id\" TEXT);", NULL, NULL, NULL);
    sqlite3_exec(db, "PRAGMA journal_mode = OFF; PRAGMA synchronous = 0; PRAGMA cache_size = 1000000; PRAGMA locking_mode = EXCLUSIVE; PRAGMA temp_store = MEMORY;", NULL, NULL, NULL);

    ifstream inputFile("Homo_sapiens.GRCh38.112.chr.gtf");
    string *commonlinedata = new string();

    vector<thread> readers;
    for (int i = 0; i < num_threads; ++i) {
        readers.emplace_back(myreadthread, &inputFile, commonlinedata);
    }

    thread writer(mywritethread, &db, commonlinedata);

    for (auto &t : readers) {
        t.join();
        // cout << "Reader thread finished execution" << endl;
    }

    writer.join();
    // cout << "Write Thread finished execution" << endl;

    sqlite3_exec(db, "CREATE INDEX \"ix_human_index\" ON \"human\" (\"index\");", NULL, NULL, NULL);
    sqlite3_close(db);
}
