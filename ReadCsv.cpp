#include "ReadCsv.h"

#include <fstream>
// #include <iostream>
// #include <cstdio>
using namespace std;

ReadCsv::ReadCsv(){};
ReadCsv::~ReadCsv(){};

bool ReadCsv::getInputs(std::string fileName){
  string line, word;
  int supportCount = 0;
  int loadCount = 0;
  int match = 0;

  ifstream inputFile(fileName, ios::in);
  if (!inputFile) {
    cout << "Error while reading input file PREPROCESSOR.csv!" << endl;
    return false;
  }

  while ( getline(inputFile, line)){
    stringstream s(line);
    // If line in the file starts with # symbol it means this is a comment
    // so the line is ignored
    if (line.rfind("#", 0 ) == 0 ) {
      continue;
    // If line in the file starts with $GEOMETRY support count is read from this line (4th column)
    // Next n lines (where n is support count) will be read and stored to supports matrix
    // Same applies for $LOADS
    } else if (line.rfind("$GEOMETRY", 0 ) == 0 ){
      int i=0;
      while (getline(s, word, ';')){
        if (i == 3) {
          try {
            supportCount = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - you need to provide valid support count in the 4th column!" << endl;
            return false;
          }
          break;
        }
        i++;
      }
      supports = MatrixXd::Zero(supportCount,2);
      continue;
    } else if (line.rfind("$LOADS", 0 ) == 0 ){
      int i=0;
      while (getline(s, word, ';')){
        if (i == 3) {
          try {
            loadCount = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - you need to provide valid load count in the 4th column!" << endl;
            return false;
          }
          break;
        }
        i++;
      }
      continue;
    }

    // while support count is not 0  csv file lines are read and stored in supports matrix
    if (supportCount > 0){
      int i=0;
      int coordX = 0;
      int supportType = 0;
      while (getline(s, word, ';')){
        if (i == 1) {
          try {
          coordX = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - valid coordinates for the supports!" << endl;
            return false;
          }
        } else if (i == 2) {
          try {
            supportType = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - valid support type!" << endl;
            return false;
          }
          break;
        }
        i++;
      }
      supports(supports.rows() - supportCount,0) = coordX;
      supports(supports.rows() - supportCount,1) = supportType;
      supportCount--;
    }

    // while point Load count is not 0 csv file lines are read and stored in pointLoads matrix
    if (loadCount > 0){
      int i=0;
      int loadType = 0;
      int coordX_0 = 0;
      int coordX_1 = 0;
      int loadValue_0 = 0;
      int loadValue_1 = 0;
      while (getline(s, word, ';')){
        if (i == 1) {
          try {
            loadType = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - load type is not valid!" << endl;
            return false;
          }
        } else if (i == 2) {
          try {
            coordX_0 = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - coordX_0 value is not valid!" << endl;
            return false;
          }
        } else if (i == 3) {
          try {
            loadValue_0 = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - load_0 value is not valid!" << endl;
            return false;
          }
        } else if (i == 4) {
          try {
            coordX_1 = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - coordX_1 value is not valid!" << endl;
            return false;
          }
        } else if (i == 5) {
          try {
            loadValue_1 = stoi(word);
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - load_1 value is not valid!" << endl;
            return false;
          }
        }
        i++;
      }
      Load f1(loadType, coordX_0, loadValue_0, coordX_1, loadValue_1);
      loads.push_back(f1);
      loadCount--;
    }

    // while loop to get additional beam properties - currently only beam modulos of elasticity is supported
    // and beam width and beam height (rectangular sections only)
    while (getline(s, word, ';'))
    { 
      if (word == "E"){
        match = 1;
        continue;
      } else if (word == "B"){
        match = 2;
        continue;
      } else if (word == "H"){
        match = 3;
        continue;
      } 
 
      switch (match) {
        case 1:
          try {
            data.push_back(stoi(word));
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - E value is not valid!" << endl;
            return false;
          }
          match = 0;
          break;
        case 2:
          try {
            data.push_back(stoi(word));
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - B value is not valid!" << endl;
            return false;
          }
          match = 0;
          break;
        case 3:
          try {
            data.push_back(stoi(word));
          }
          catch(std::invalid_argument& e){
            cout << "PREPROCESSOR failure - H value is not valid!" << endl;
            return false;
          }
          match = 0;
          break;
        default:
          break;
      }
    }
  }
  return true;
}

std::vector<int> ReadCsv::getData(){
    return data;
};

std::vector<Load> ReadCsv::getLoads(){
    return loads;
};

MatrixXd ReadCsv::getSupports(){
    return supports;
};

