#ifndef READCSV_H
#define READCSV_H

#include <iostream>
#include <vector> 
#include <Eigen/Dense>
#include "PointLoad.h"


using namespace Eigen;


class ReadCsv 
{
    private:
        std::vector<int> data;
        std::vector<PointLoad> pointLoads;
        MatrixXd supports;

    public:
        ReadCsv(std::string fileName);
        ~ReadCsv();
        std::vector<int> getData();
        std::vector<PointLoad> getPointLoads();
        MatrixXd getSupports();
};

#endif