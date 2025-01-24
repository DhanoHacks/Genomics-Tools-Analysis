#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sqlite3.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>

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
    std::string sql = "INSERT INTO human(";
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

using namespace std;
int main(){
    sqlite3 *db;
    int rc;
    rc = sqlite3_open(":memory:", &db);
    rc = sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS \"human\" (\"index\" INTEGER,\"Chromosome\" TEXT,Source\" TEXT,\"Feature\" TEXT,\"Start\" INTEGER,\"End\" INTEGER,\"Score\" TEXT,\"Strand\" TEXT,\"Frame\" TEXT,\"gene_id\" TEXT,\"gene_version\" TEXT,\"gene_source\" TEXT,\"gene_biotype\" TEXT,\"transcript_id\" TEXT,\"transcript_version\" TEXT,\"transcript_source\" TEXT,\"transcript_biotype\" TEXT,\"tag\" TEXT,\"transcript_support_level\" TEXT,\"exon_number\" TEXT,\"exon_id\" TEXT,\"exon_version\" TEXT,\"gene_name\" TEXT,\"transcript_name\" TEXT,\"protein_id\" TEXT,\"protein_version\" TEXT,\"ccds_id\" TEXT);" , NULL, NULL, NULL);
    rc = sqlite3_exec(db, "PRAGMA journal_mode = OFF; PRAGMA synchronous = 0; PRAGMA cache_size = 1000000; PRAGMA locking_mode = EXCLUSIVE; PRAGMA temp_store = MEMORY;", NULL, NULL, NULL);
    // open file "Homo_sapiens.GRCh38.112.chr.gtf"
    ifstream inputFile("Homo_sapiens.GRCh38.112.chr.gtf");
    string line;
    int linecount = 0;
    while (getline(inputFile, line)) {
        linecount++;
        if (linecount%100000 == 0){
            cout<<"Processed "<<linecount<<" lines"<<endl;
        }
        // ignore first few lines beginning with #
        if(line[0] == '#'){
            continue;
        }
        string linedata = processLine(line);
        // Output the parsed data
        // vector<string> keys={"Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"};
        // sqlite3_stmt *ppStmt;
        // rc = sqlite3_prepare_v2(db,linedata.c_str(),-1,&ppStmt,NULL);
        // rc = sqlite3_step(ppStmt);
        rc = sqlite3_exec(db, linedata.c_str(), NULL, NULL, NULL);
    }
    sqlite3_close(db);
}