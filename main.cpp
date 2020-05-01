#include <iostream>
#include <cstdio>
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
  
  std::vector<PointLoad> pointLoads;
  PointLoad f1("F1", 2500, 5000);
  pointLoads.push_back(f1);
  

  //Matrix for mapping support coordinates to support type enum (Hinge/RIGID)
  //For testing purposes hardcoded values are set
  int spanCount=1;
  MatrixXd supports(spanCount+1,2);
  supports << 0,    0,
              5000, 0;
              // 7500, 1;

  //Jungs modulus of the beam in N/mm2 - equivalent for wood
  int E = 11000;

  //parameters of the beam section
  // Second moment of area
  int b = 50;
  int h = 200;
  double I = (b * pow(h, 3))/12;


  

  //Number of elements per beam
  // int elemsPerBeam = 10;

  //Vector for storing X coordinates - element end points
  // VectorXd dx(spanCount*elemsPerBeam + pointLoads.size());
   VectorXd dx(spanCount + pointLoads.size());

  //span is divided into 10 parts - but maybe this can be avoided?
  // dx(0) = 0;
  // int current = 0;
  // for (int i=0; i<supports.rows()-1; i++){
  //   int spanLength = supports(i+1,0) - supports(i,0);
  //   for (int j=1; j<=elemsPerBeam; j++){
  //     dx(i*10+j)= dx(i*10+j-1) + spanLength/elemsPerBeam;
  //   }
  // }

  for (int i=0; i<supports.rows(); i++){
    dx(i) = supports(1, i);
  }

  //resizing vector and adding node coordinate for the force applied (2500mm in this case) to the end of the vector
  //node is added to vector only if value is not already found in the vector
  bool flag = false;
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
  int delta = 0;
  for(int i=0; i < (dx.size()-1); i++){
    double L = dx(i+1) - dx(i);
    double stiffness = E*I/pow(L,3);
    MatrixXd k_e = EulerBernoulli::getStiffnessMatrix(L);
    for(int n=0; n<4; n++){
      for(int m=0; m<4; m++){
        K(n+delta, m+delta ) = K(n+delta, m+delta ) + stiffness * k_e(n, m);
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

  //Load Vector
  VectorXd F = VectorXd::Zero(K.rows());
  MatrixXd result = MatrixXd::Zero(pointLoads.size(),K.rows());



  for(int i=0; i<pointLoads.size(); i++){
    for(int j=0; j< dx.size(); j++){
        if(pointLoads[i].getCoordX() == dx(j)){
          F(2*j) -= pointLoads[i].getLoadValue();  
        } 
    }

    for(int j=0; j < F.size(); j++){
      if(FBcs(j) == 1) {
        F(j) = 0;
      }
    }

    // Solution matrix is calculated for each loading case and stored in matrix
    // Later with superposition principle this can be added together for final results
    MatrixXd tempResult = K.inverse() * F;  

    for(int j=0; j < tempResult.size(); j++){
        result(i, j) = tempResult(j);
    }

  }


  //Result postprocessing

  // The size of the mesh elements (not the count)
  // Mesh is needed for results displaying
  int mesh = 500;
  
  //Displacement result vector
  MatrixXd u = MatrixXd::Zero(pointLoads.size(),1);
  VectorXd uX = VectorXd::Zero(1);

  //counter
  int m =0;


  for(int i=1; i < (dx.size()); i++){
    double L = dx(i) - dx(i-1);
    // cout << L << endl;

    for(int j=0; j < L; j+=mesh ) {
          // cout << dx(i-1) + j << endl;
          double n1 = 1 / pow(L, 3) * (2 * pow(j, 3) - 3 * pow(j, 2) * L + pow(L, 3));
          double n2 = 1 / pow(L, 3) * (pow(j, 3) * L - 2 * pow(j, 2) * pow(L, 2) + j * pow(L, 3));
          double n3 = 1 / pow(L, 3) * (-2 * pow(j, 3) + 3 * pow(j, 2) * L);
          double n4 = 1 / pow(L, 3) * (pow(j, 3) * L - pow(j, 2) * pow(L, 2));
          for (int k=0; k<pointLoads.size(); k++){
            u(k, m) = n1 * result(k, 2 * (i+1) - 4) + n2 * result(k, 2 * (i+1) - 3) + n3 * result(k, 2 * (i+1) - 2) + n4 * result(k, 2 * (i+1) - 1);
          }

          uX(m) = dx(i-1) + j;

          m++;
          u.conservativeResize(pointLoads.size(),m+1); //resize(pointLoads.size(),m+1);
          uX.conservativeResize(m+1);
    }

    if ( i==(dx.size()-1)){
      uX(m) = dx(i);
    }
    
 
  }




  //Temporary solution for outputting files to txt document

  ofstream file("matrix.txt");
  // file.open("matrix.txt", ios::out | ios::app);

  if (file.is_open())
  {
    // file << "m" << '\n' <<  K << '\n';
    for (int i=0; i<u.cols(); i++){
      file << uX(i) << "\t";
      file << u(0, i) << endl;
    }
  }

  // plotting graph
  FILE* pipe = _popen("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe", "w");
    if (pipe != NULL)
    {
        fprintf(pipe, "set term win\n");
        fprintf(pipe, "set style line 1 lt 1 lw 3 pt 3 linecolor rgb \"red\"\n");
        
        fprintf(pipe, "plot \"matrix.txt\" with linespoints linestyle 1\n");
        fprintf(pipe, "set term pngcairo\n");
        fprintf(pipe, "set output \"myFile.png\"\n" );
        fprintf(pipe, "replot\n");
        fprintf(pipe, "set term win\n");
        fflush(pipe);
    }
    else puts("Could not open the file\n");
    _pclose(pipe);

  double elemLength = 2;
  cout << EulerBernoulli::getStiffnessMatrix(elemLength) << endl;
  cout << "The size of the Global stiffness matrix is " << K.rows() << "x" << K.cols() << endl;
  cout << K << endl;
}
