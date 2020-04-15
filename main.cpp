#include <iostream>
#include <Eigen/Dense>
#include "EulerBernoulli.h"
 
using namespace Eigen;
using namespace std;
 
int main()

{
  int spanCount=1;
  int spanLenght=5000; 
  int current = 0;
  bool flag = false;
  double stiffness = 0;

  //Jungs modulus of the beam in N/mm2 - equivalent for wood
  int E = 11000;

  //parameters of the beam section
  int b = 50;
  int h = 200;

  // Second moment of area
  double I = (b * pow(h, 3))/12;
  cout << "This is second moment of area for the beam: " << I << "mm4" << endl;

  //Coordinate where concentrated force is applied
  int forcePoint = 2500;

  VectorXd dx(spanCount*10 + 1);

  //foreach loop - supported only for eigen 3.4 version
  //while iterating should use refernce - so we could set values
  for (auto &v : dx) {
    v = current;
    current += spanLenght/10;
  } 

  //resizing vector and adding node coordinate for the force applied (2500mm in this case) to the end of the vector
  //node is added to vector only if value is not already found in the vector
  for (int i=0; i<dx.size(); i++){
    if (dx(i) == forcePoint){
      flag = true;
    }
  }

  if (!flag){
    dx.resize(dx.size()+1);
    dx(dx.size()-1) = forcePoint;
  }
  
  //Sorting vector - so that added node values are in ascending order
  sort(dx.begin(), dx.end());

  //Global stiffnes matrix
  MatrixXd K = MatrixXd::Zero(dx.size()*2, dx.size()*2 );
  MatrixXd stiffK = MatrixXd::Zero(K.rows(), K.cols());

  int delta = 0;
  for(int i=0; i < (dx.size()-1); i++){
    double L = dx(i+1) - dx(i);

    stiffness = E*I/pow(L,3);

    MatrixXd k_e = EulerBernoulli::getStiffnessMatrix(L);
    // cout << "++++++++++++++++++++++++++++++++++++" << endl;
    // cout << k_e << endl;
    for(int n=0; n<4; n++){
      for(int m=0; m<4; m++){
        K(n+delta, m+delta ) = K(n+delta, m+delta ) + stiffness * k_e(n, m);
        stiffK(n+delta, m+delta ) = stiffK(n+delta, m+delta ) + stiffness * k_e(n, m);
      }
    }
    delta = delta +2;
  }

  //Inverse of the matrix
  MatrixXd invK = K.inverse();

  // double elemLength = 2;
  // cout << EulerBernoulli::getStiffnessMatrix(elemLength) << endl;
  cout << "The size of the Global stiffness matrix is " << K.rows() << "x" << K.cols() << endl;
  cout << K << endl;
}
