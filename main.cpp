#include <iostream>
#include <cstdio>
#include <Eigen/Dense>
#include "EulerBernoulli.h"
#include "PointLoad.h"
#include "ReadCsv.h"
#include <fstream>
#include <vector> 
 
using namespace Eigen;
using namespace std;

//Do not need this at the moment, but could come in handy
enum Support {HINGE, RIGID};
 
int main()
{


  /*
    --------------------------------------------------------------------------------------------------------------
    |                                    PREPROCESSING PHASE                                                     |
    --------------------------------------------------------------------------------------------------------------
  */

  // Parameters for the beam section
  int E = 0; 
  int b = 0;
  int h = 0;

  ReadCsv preprocessorData("PREPROCESSOR.csv");
  std::vector<int> beamData = preprocessorData.getData();
  MatrixXd supports = preprocessorData.getSupports();
  std::vector<PointLoad> pointLoads = preprocessorData.getPointLoads();

  // Beam element properties are retrieved from beamData vector and assignec accordingly
  // Correct sequence in beamData is expected
  int n = 0;
  for ( auto &i : beamData ) {
    if ( n == 0){
      E = i;
    } else if ( n == 1){
      b = i;
    } else if ( n == 2){
      h = i;
    }
    n++;
  }
  
  // Second moment of area
  double I = (b * pow(h, 3))/12;


  /*
    --------------------------------------------------------------------------------------------------------------
    |                                    SOLUTION PHASE                                                          |
    --------------------------------------------------------------------------------------------------------------
  */
 
  //Vector for storing X coordinates - element end points
  // VectorXd dx(spanCount*elemsPerBeam + pointLoads.size());
  VectorXd dx(supports.rows());
  for (int i=0; i<supports.rows(); i++){
    dx(i) = supports(i, 0);
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
        dx.conservativeResize(dx.size()+1);
        dx(dx.size()-1) = pointLoads[i].getCoordX();
      }
  }
  
  //Sorting vector - so that added node values are in ascending order
  sort(dx.begin(), dx.end());

  // Global stiffnes matrix is created in this loop according
  // to Euler-Bernoulli beam bending 
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

  // Boundary conditions are set for Gloabal stiffness matrix
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

  // Results matrix - displacements and rotation at d(x) points
  // Results are calculated for each load case - so that
  // principle of supperposition could be applied later for final displacements  
  MatrixXd result = MatrixXd::Zero(pointLoads.size(),K.rows());
  for(int i=0; i<pointLoads.size(); i++){
    VectorXd F = VectorXd::Zero(K.rows());
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
    MatrixXd tempResult = K.inverse() * F;  
    for(int j=0; j < tempResult.size(); j++){
        result(i, j) = tempResult(j);
    }

  }


  /*
    --------------------------------------------------------------------------------------------------------------
    |                                    POSTPROCESSING PHASE                                                    |
    --------------------------------------------------------------------------------------------------------------
  */

  // The length of the mesh elements
  // Mesh is needed for displaying results - to refine calculation points, set finer mesh (shorter elements)
  int mesh = 500;
  
  // Resultant displacement vector
  MatrixXd u = MatrixXd::Zero(pointLoads.size(),1);

  // Coordinates for mesh points and corresponding result values
  VectorXd uX = VectorXd::Zero(1);

  // This is result recovery from calculated displacements matrix - results (displacements) are retrieved for
  // each mesh point
  int m =0;
  for(int i=1; i < (dx.size()); i++){
    double L = dx(i) - dx(i-1);
    for(int j=0; j < L; j+=mesh ) {
      double n1 = 1 / pow(L, 3) * (2 * pow(j, 3) - 3 * pow(j, 2) * L + pow(L, 3));
      double n2 = 1 / pow(L, 3) * (pow(j, 3) * L - 2 * pow(j, 2) * pow(L, 2) + j * pow(L, 3));
      double n3 = 1 / pow(L, 3) * (-2 * pow(j, 3) + 3 * pow(j, 2) * L);
      double n4 = 1 / pow(L, 3) * (pow(j, 3) * L - pow(j, 2) * pow(L, 2));
      for (int k=0; k<pointLoads.size(); k++){
        u(k, m) = n1 * result(k, 2 * (i+1) - 4) + n2 * result(k, 2 * (i+1) - 3) + n3 * result(k, 2 * (i+1) - 2) + n4 * result(k, 2 * (i+1) - 1);
      }
      uX(m) = dx(i-1) + j;
      m++;
      u.conservativeResize(pointLoads.size(),m+1); 
      uX.conservativeResize(m+1);
    }
    if ( i==(dx.size()-1)){
      uX(m) = dx(i);
    }
  }

  //Temporary solution for outputting files to txt document

  ofstream file("matrix.txt");
  if (file.is_open())
  {
    for (int i=0; i<u.cols(); i++){
      file << uX(i) << "\t";
      file << u(1, i) << "\t";
      file << u(0, i) << endl;
    }
  }

  // plotting graph
  // FILE* pipe = _popen("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe", "w");
  //   if (pipe != NULL)
  //   {
  //       fprintf(pipe, "set term win\n");
  //       fprintf(pipe, "set style line 1 lt 1 lw 3 pt 3 linecolor rgb \"red\"\n");
        
  //       fprintf(pipe, "plot \"matrix.txt\" with linespoints linestyle 1\n");
  //       fprintf(pipe, "set term pngcairo\n");
  //       fprintf(pipe, "set output \"myFile.png\"\n" );
  //       fprintf(pipe, "replot\n");
  //       fprintf(pipe, "set term win\n");
  //       fflush(pipe);
  //   }
  //   else puts("Could not open the file\n");
  //   _pclose(pipe);
}
