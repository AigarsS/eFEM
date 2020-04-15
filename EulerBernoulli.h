#ifndef EULERBERNOULLI_H
#define EULERBERNOULLI_H

#include <Eigen/Dense>
 
using namespace Eigen;

class EulerBernoulli
{
public:
    EulerBernoulli();
    ~EulerBernoulli();
    static Matrix4d getStiffnessMatrix(double l_e);
};

#endif