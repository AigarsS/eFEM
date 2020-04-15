#include "EulerBernoulli.h"
#include <Eigen/Dense>

using namespace Eigen;

EulerBernoulli::EulerBernoulli(){};

EulerBernoulli::~EulerBernoulli(){};

//Method takes element length as parameter and return stiffenss matrix
Matrix4d EulerBernoulli::getStiffnessMatrix(double l_e){
  //Stiffness matrix for Euler-Bernoulli beam
  Matrix4d M;
  M <<  12.0,       6.0 * l_e,       -12.0,        6.0 * l_e,
        6.0 * l_e,  4.0 * l_e * l_e,  -6.0 * l_e,   2.0 * l_e * l_e,
        -12.0,      -6.0 * l_e,        12.0,        -6.0 * l_e,
        6.0 * l_e,  2.0 * l_e * l_e,   -6.0 * l_e,   4.0 * l_e * l_e;
  return M;
};


