#ifndef READCSV_H
#define READCSV_H

#include <iostream>
#include <vector> 
#include <Eigen/Dense>
#include "Load.h"


using namespace Eigen;


class ReadCsv 
{
    private:
        std::vector<int> data;
        std::vector<Load> loads;
        MatrixXd supports;

    public:
        ReadCsv();
        ~ReadCsv();
        std::vector<int> getData();
        std::vector<Load> getLoads();
        MatrixXd getSupports();
        bool getInputs(std::string fileName);
};

#endif