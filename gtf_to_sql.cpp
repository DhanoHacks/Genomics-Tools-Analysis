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
#include <semaphore>
using namespace std;

mutex gtf_mtx, sql_mutex;
int num_finished = 0, num_threads = 8;

std::string processLine(const std::string &line) {
    // Split the line by tabs
    std::istringstream lineStream(line);
    std::string token;
    std::vector<std::string> fields;

    while (std::getline(lineStream, token, '\t')) {
        fields.push_back(token);
    }

    if (fields.size() < 9) {
        throw std::runtime_error("Error: Line has less than 9 tab-separated fields.");
    }

    // Map the first 8 fields to keys
    std::unordered_map<std::string, std::string> data;
    std::vector<std::string> keys = {"Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"};

    for (size_t i = 0; i < 8; ++i) {
        data[keys[i]] = fields[i];
    }

    // Parse the additional key-value pairs from the 9th field
    std::string attributes = fields[8];
    std::istringstream attrStream(attributes);
    std::string attribute;

    while (std::getline(attrStream, attribute, ';')) {
        attribute = std::string(std::find_if(attribute.begin(), attribute.end(), [](char c) { return !std::isspace(c); }), attribute.end());

        size_t pos = attribute.find(' ');
        if (pos != std::string::npos) {
            std::string key = attribute.substr(0, pos);
            std::string value = attribute.substr(pos + 1);

            // Remove quotes from value if present
            value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
            data[key] = value;
        }
    }

    // Construct the SQL INSERT statement
    std::string sql = "INSERT INTO human2(";
    std::string keysStr;
    std::string valuesStr;

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

void myreadthread(ifstream *inputFile, string *commonlinedata){
    string line, linedata="";
    vector<string> lines;
    int linecount = 0;
    bool eof = false;
    while (!eof) {
        gtf_mtx.lock();
        for(int i=0; i<10; i++){
            if(!getline(*inputFile, line)){
                eof = true;
                break;
            }
            // ignore first few lines beginning with #
            if(line[0] == '#'){
                continue;
            }
            lines.push_back(line);
        }
        gtf_mtx.unlock();

        linecount+=10;
        if (linecount%100000 == 0){
            cout<<"Processed "<<linecount<<" lines"<<endl;
        }

        for(auto line: lines){
            linedata += processLine(line);
        }
        lines.clear();

        if(linecount%50 == 0){
            if(sql_mutex.try_lock()){
                *commonlinedata += linedata;
                sql_mutex.unlock();
                linedata = "";
            }
        }

        // if(linecount%50 == 0){
        //     sql_mutex.lock();
        //     sqlite3_exec(*db, linedata.c_str(), NULL, NULL, NULL);
        //     sql_mutex.unlock();
        //     linedata = "";
        // }
    }
    sql_mutex.lock();
    if(linedata != ""){
        // sqlite3_exec(*db, linedata.c_str(), NULL, NULL, NULL);
        *commonlinedata += linedata;
    }
    num_finished++;
    sql_mutex.unlock();
}

void mywritethread(sqlite3 **db, string *commonlinedata){
    while(true){
        if(sql_mutex.try_lock()){
            if(commonlinedata->length() > 0){
                string *copy_commonlinedata = new string();
                copy_commonlinedata->append(*commonlinedata);
                delete commonlinedata;
                commonlinedata = new string();
                sql_mutex.unlock();
                sqlite3_exec(*db, copy_commonlinedata->c_str(), NULL, NULL, NULL);
                delete copy_commonlinedata;
            }
            else{
                if(num_finished == num_threads){
                    sql_mutex.unlock();
                    return;
                }
                sql_mutex.unlock();
            }
        }
    }
}

int main(){
    sqlite3 *db;
    int rc;
    rc = sqlite3_open("db.sqlite3", &db);
    rc = sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS \"human2\" (\"Chromosome\" TEXT, \"Source\" TEXT,\"Feature\" TEXT,\"Start\" INTEGER,\"End\" INTEGER,\"Score\" TEXT,\"Strand\" TEXT,\"Frame\" TEXT,\"gene_id\" TEXT,\"gene_version\" TEXT,\"gene_source\" TEXT,\"gene_biotype\" TEXT,\"transcript_id\" TEXT,\"transcript_version\" TEXT,\"transcript_source\" TEXT,\"transcript_biotype\" TEXT,\"tag\" TEXT,\"transcript_support_level\" TEXT,\"exon_number\" TEXT,\"exon_id\" TEXT,\"exon_version\" TEXT,\"gene_name\" TEXT,\"transcript_name\" TEXT,\"protein_id\" TEXT,\"protein_version\" TEXT,\"ccds_id\" TEXT);" , NULL, NULL, NULL);
    rc = sqlite3_exec(db, "PRAGMA journal_mode = OFF; PRAGMA synchronous = 0; PRAGMA cache_size = 1000000; PRAGMA locking_mode = EXCLUSIVE; PRAGMA temp_store = MEMORY;", NULL, NULL, NULL);
    // open file "Homo_sapiens.GRCh38.112.chr.gtf"
    ifstream inputFile("Homo_sapiens.GRCh38.112.chr.gtf");
    string *commonlinedata = new string();
    // string line;
    // int linecount = 0;
    // while (getline(inputFile, line)) {
    //     linecount++;
    //     if (linecount%100000 == 0){
    //         cout<<"Processed "<<linecount<<" lines"<<endl;
    //     }
    //     // ignore first few lines beginning with #
    //     if(line[0] == '#'){
    //         continue;
    //     }
    //     // string linedata = processLine(line);
    //     // rc = sqlite3_exec(db, linedata.c_str(), NULL, NULL, NULL);
    // }
    thread t1(myreadthread, &inputFile, commonlinedata);
    thread t2(myreadthread, &inputFile, commonlinedata);
    thread t3(myreadthread, &inputFile, commonlinedata);
    thread t4(myreadthread, &inputFile, commonlinedata);
    thread t5(myreadthread, &inputFile, commonlinedata);
    thread t6(myreadthread, &inputFile, commonlinedata);
    thread t7(myreadthread, &inputFile, commonlinedata);
    thread t8(myreadthread, &inputFile, commonlinedata);
    thread w(mywritethread, &db, commonlinedata);

    t1.join(); cout<<"Thread1 finished execution"<<endl;
    t2.join(); cout<<"Thread2 finished execution"<<endl;
    t3.join(); cout<<"Thread3 finished execution"<<endl;
    t4.join(); cout<<"Thread4 finished execution"<<endl;
    t5.join(); cout<<"Thread5 finished execution"<<endl;
    t6.join(); cout<<"Thread6 finished execution"<<endl;
    t7.join(); cout<<"Thread7 finished execution"<<endl;
    t8.join(); cout<<"Thread8 finished execution"<<endl;
    w.join(); cout<<"Write Thread finished execution"<<endl;
    rc = sqlite3_exec(db, "CREATE INDEX \"ix_human_index\"ON \"human\" (\"index\");", NULL, NULL, NULL);
    sqlite3_close(db);
}