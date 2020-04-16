#include <iostream>
#include <Eigen/Dense>
#include "EulerBernoulli.h"
#include "PointLoad.h"
#include <fstream>
#include <vector> 
 
using namespace Eigen;
using namespace std;

//Do not need this at the moment, but could come in handy
enum Support {HINGE, RIGID};
 
int main()

{
  int spanCount=1;
  // int spanLenght=5000; 


  std::vector<PointLoad> pointLoads;
  PointLoad f1("F1", 2500, 5000);
  pointLoads.push_back(f1);
  

  //Matrix for mapping support coordinates to support type enum (Hinge/RIGID)
  MatrixXd supports(spanCount+1,2);
  
  //For testing purposes hardcoded values are set
  supports << 0,    0,
              5000, 1;
              // 7500, 1;

  
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

  //Vector for storing X coordinates - element end points
  VectorXd dx(spanCount*10 + pointLoads.size());

  //Number of elements per beam
  int elemsPerBeam = 10;

  //span is divided into 10 parts - but maybe this can be avoided?
  dx(0) = 0;
  int current = 0;
  for (int i=0; i<supports.rows()-1; i++){
    int spanLength = supports(i+1,0) - supports(i,0);
    for (int j=1; j<=elemsPerBeam; j++){
      dx(i*10+j)= dx(i*10+j-1) + spanLength/elemsPerBeam;
    }
  }

  //resizing vector and adding node coordinate for the force applied (2500mm in this case) to the end of the vector
  //node is added to vector only if value is not already found in the vector
  for(int i=0; i<pointLoads.size(); i++){
      for (int j=0; j<dx.size(); j++){
        if (dx(j) == pointLoads[i].getCoordX()){
          flag = true;
        }
      }
      if (!flag){
        dx.resize(dx.size()+1);
        dx(dx.size()-1) = pointLoads[i].getCoordX();
      }
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
    for(int n=0; n<4; n++){
      for(int m=0; m<4; m++){
        K(n+delta, m+delta ) = K(n+delta, m+delta ) + stiffness * k_e(n, m);
        stiffK(n+delta, m+delta ) = stiffK(n+delta, m+delta ) + stiffness * k_e(n, m);
      }
    }
    delta = delta +2;
  }

  //vector for applied boundary conditions (maybe not necessary)
  VectorXd FBcs = VectorXd::Zero(K.rows());

  for(int i=0; i<supports.rows(); i++){
    int supportType = supports(i,1);
    for(int j=0; j< dx.size(); j++){
      if ( supports(i, 0) == dx(j) ){
          switch (supportType) {
                //This corresponds to case when support is hinge
                case 0:
                  K.row(2*j).setZero();
                  K(2*j, 2*j) = 1;
                  FBcs(2*j) = 1;
                  break;
                //This corresponds to case when support is rigid
                case 1:
                  K.row(2*j).setZero();
                  K.row(2*j+1).setZero();
                  K(2*j, 2*j) = 1;
                  K(2*j+1, 2*j+1) = 1;  
                  FBcs(2*j) = 1;
                  break;

                default:
                  break;
            }
      }
    }
  }  

    
  
  //Inverse of the matrix
  // MatrixXd invK = K.inverse();


  //Load Vector
  VectorXd F = VectorXd::Zero(K.rows());

  for(int i=0; i<pointLoads.size(); i++){
    for(int j=0; j< dx.size(); j++){
        if(pointLoads[i].getCoordX() == dx(j)){
          F(2*j) = pointLoads[i].getLoadValue();  
        } 
    }  
  }

  MatrixXd result = K.inverse() * F;

  ofstream file("matrix.txt");

  if (file.is_open())
  {
    // file << "m" << '\n' <<  K << '\n';
    file << "m" << '\n' <<  result << '\n';
  }

  // double elemLength = 2;
  // cout << EulerBernoulli::getStiffnessMatrix(elemLength) << endl;
  // cout << "The size of the Global stiffness matrix is " << K.rows() << "x" << K.cols() << endl;
  // cout << K << endl;
}
