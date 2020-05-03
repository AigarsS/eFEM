#include "ReadCsv.h"

#include <fstream>
// #include <iostream>
// #include <cstdio>
using namespace std;

ReadCsv::ReadCsv(std::string fileName){
  string line, word;
  int supportCount = 0;
  int loadCount = 0;

  int match = 0;

  ifstream inputFile(fileName, ios::in);

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
          supportCount = stoi(word);
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
          loadCount = stoi(word);
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
          coordX = stoi(word);
        } else if (i == 2) {
          supportType = stoi(word);
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
          loadType = stoi(word); 
        } else if (i == 2) {
          coordX_0 = stoi(word);
        } else if (i == 3) {
          loadValue_0 = stoi(word);
        } else if (i == 4) {
          coordX_1 = stoi(word);
        } else if (i == 5) {
          loadValue_1 = stoi(word);
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
          data.push_back(stoi(word));
          match = 0;
          break;
        case 2:
         data.push_back(stoi(word));
          match = 0;
          break;
        case 3:
          data.push_back(stoi(word));
          match = 0;
          break;
        default:
          break;
      }
    }

  }

};

ReadCsv::~ReadCsv(){};

std::vector<int> ReadCsv::getData(){
    return data;
};

std::vector<Load> ReadCsv::getLoads(){
    return loads;
};

MatrixXd ReadCsv::getSupports(){
    return supports;
};

