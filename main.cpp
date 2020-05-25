#include <iostream>
#include <cstdio>
#include <Eigen/Dense>
#include "EulerBernoulli.h"
#include "Load.h"
#include "ReadCsv.h"
#include <fstream>
#include <vector> 
#include <cmath>
 
using namespace Eigen;
using namespace std;

//Do not need this at the moment, but could come in handy
enum Support {HINGE, RIGID};
 
int main()
{
  double maxDeflection;
  double maxV;
  double maxM;

  /*
    --------------------------------------------------------------------------------------------------------------
    |                                    PREPROCESSING PHASE                                                     |
    --------------------------------------------------------------------------------------------------------------
  */

  // Parameters for the beam section
  int E = 0; 
  int b = 0;
  int h = 0;

  ReadCsv preprocessorData;
  if(!preprocessorData.getInputs("./PREPROCESSOR.csv")){
    cout << "Correct inputs and restart calculation!" << endl;
    exit (EXIT_FAILURE);
  };

  std::vector<int> beamData = preprocessorData.getData();
  MatrixXd supports = preprocessorData.getSupports();
  std::vector<Load> loads = preprocessorData.getLoads();

  // Beam element properties are retrieved from beamData vector and assignec accordingly
  // Correct sequence in beamData is expected

  if(beamData.size() != 3){
    cout << "There is not sufficient $MATERIAL OR $SECTION data in PREPROCESSOR FILE" << endl;
    exit (EXIT_FAILURE);
  };

  for (int i = 0; i < beamData.size(); i++){
    if (i == 0) {
      E = beamData[i];
    } else if (i ==1){
      b = beamData[i];
    } else if (i ==2){
      h = beamData[i];
    }
  }

  // Second moment of area
  double I = (b * pow(h, 3))/12;


  /*
    --------------------------------------------------------------------------------------------------------------
    |                                    SOLUTION PHASE                                                          |
    --------------------------------------------------------------------------------------------------------------
  */

 //Number of elements per beam (for some reason result recovery function does not give precise results)
  int elemsPerBeam = 50;
 
  //Vector for storing X coordinates - element end points
  VectorXd dx((supports.rows()-1)*elemsPerBeam + 1);
  
  //span is divided into 10 parts - but maybe this can be avoided?
  dx(0) = 0;
  for (int i=0; i<supports.rows()-1; i++){
    int spanLength = supports(i+1,0) - supports(i,0);
    for (int j=1; j<=elemsPerBeam; j++){
      dx(i*elemsPerBeam+j)= dx(i*elemsPerBeam+j-1) + spanLength/elemsPerBeam;
    }
  }

  //resizing vector and adding node coordinate for the force applied (2500mm in this case) to the end of the vector
  //node is added to vector only if value is not already found in the vector
  for(int i=0; i<loads.size(); i++){
      bool flag = false;
      for (int j=0; j<dx.size(); j++){
        if (dx(j) == loads[i].getCoordX0()){
          flag = true;
        } 
      }
      if (!flag){
        dx.conservativeResize(dx.size()+1);
        dx(dx.size()-1) = loads[i].getCoordX0();
      }
      flag = false;
      if (loads[i].getLoadType() == 1){
        for (int j=0; j<dx.size(); j++){
          if (dx(j) == loads[i].getCoordX1()){
            flag = true;
          } 
        }
        if (!flag){
          dx.conservativeResize(dx.size()+1);
          dx(dx.size()-1) = loads[i].getCoordX1();
        }
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
        K(n+delta, m+delta ) += stiffness * k_e(n, m);
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
          case 1:
            K.row(2*j).setZero();
            K(2*j, 2*j) = 1;
            FBcs(2*j) = 1;
            break;
          //This corresponds to case when support is rigid
          case 2:
            K.row(2*j).setZero();
            K.row(2*j+1).setZero();
            K(2*j, 2*j) = 1;
            K(2*j+1, 2*j+1) = 1;  
            FBcs(2*j) = 1;
            FBcs(2*j+1) = 1;
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
  MatrixXd result = MatrixXd::Zero(loads.size(),K.rows());
  MatrixXd f0 = MatrixXd::Zero(loads.size(),(dx.size()-1)*2);
  MatrixXd m0 = MatrixXd::Zero(loads.size(),(dx.size()-1)*2);
  for(int i=0; i<loads.size(); i++){
    VectorXd F = VectorXd::Zero(K.rows());
    double loadValue = loads[i].getLoadValue();
    switch (loads[i].getLoadType()) {
      case 0:
        for(int j=0; j< dx.size(); j++){
          if(loads[i].getCoordX0() == dx(j)){
            F(2*j) -= loadValue; 
          } 
        }
        break;
      case 1:
        for(int j=0; j< (dx.size()-1); j++){
          if(dx(j) >= loads[i].getCoordX0() && dx(j) < loads[i].getCoordX1()){
            int dLength = dx(j+1) - dx(j);
            F(2*j) -= loadValue * dLength/2;  
            F(2*j+1) -= loadValue * pow(dLength,2)/12;
            F(2*j+2) -= loadValue * dLength/2;  
            F(2*j+3) += loadValue * pow(dLength,2)/12;

            f0(i,j*2) =  -loadValue * dLength/2; 
            f0(i,j*2+1) =  -loadValue * dLength/2;
            m0(i,j*2) =  -loadValue * pow(dLength,2)/12;
            m0(i,j*2+1) =  loadValue * pow(dLength,2)/12;    
          } 
        }
        break;
      default:
        break;
    }
    for(int j=0; j < F.size(); j++){
      if(FBcs(j) == 1) {
        F(j) = 0;
      }
    }
    MatrixXd tempResult = K.inverse() * F;
    for(int j=0; j < FBcs.size(); j++){
      if(FBcs(j) == 1) {
        tempResult(j) = 0;
      }
    }
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
  int mesh = 100;
  
  // Results vectors - displacements, shear forces, bending moment
  MatrixXd u = MatrixXd::Zero(loads.size(),1);
  MatrixXd v_x = MatrixXd::Zero(loads.size(),(dx.size()-1)*2);
  MatrixXd m_x = MatrixXd::Zero(loads.size(),(dx.size()-1)*2);

  // Coordinates for result values
  VectorXd uX = VectorXd::Zero(1);
  VectorXd forcesCoord = VectorXd::Zero((dx.size()-1)*2);
  forcesCoord(0) = 0;

  // This is result recovery from calculated displacements matrix - 
  // results (displacements) are retrieved for each mesh point
  int m =0;
  int n = 0;
  for(int i=1; i < dx.size(); i++){
    double L = dx(i) - dx(i-1);
    for(int j=0; j < L; j+=mesh ) {
      double n1 = (1 / pow(L, 3) )* (2 * pow(j, 3) - 3 * pow(j, 2) * L + pow(L, 3));
      double n2 = 1 / pow(L, 3) * (pow(j, 3) * L - 2 * pow(j, 2) * pow(L, 2) + j * pow(L, 3));
      double n3 = 1 / pow(L, 3) * (-2 * pow(j, 3) + 3 * pow(j, 2) * L);
      double n4 = 1 / pow(L, 3) * (pow(j, 3) * L - pow(j, 2) * pow(L, 2));
      for (int k=0; k<loads.size(); k++){
        u(k, m) = n1 * result(k, 2 * (i+1) - 4) + n2 * result(k, 2 * (i+1) - 3) + n3 * result(k, 2 * (i+1) - 2) + n4 * result(k, 2 * (i+1) - 1);
      }
      uX(m) = dx(i-1) + j;
      m++;
      u.conservativeResize(loads.size(),m+1); 
      uX.conservativeResize(m+1);
    }
    if ( i==(dx.size()-1)){
      uX(m) = dx(i);
      for (int k=0; k<loads.size(); k++){
        u(k, m) = result(k, 2 * (i+2) - 4);
      }
    }

    
    //Recovery of the internal forces
    if (i != 1){
      forcesCoord(2*(i-1)) = forcesCoord(2*(i-1)-1);
    }
    forcesCoord(2*(i-1)+1) = forcesCoord(2*(i-1)) + L;
    double stiffness = E*I/pow(L,3);
    MatrixXd k_e = EulerBernoulli::getStiffnessMatrix(L);
    for (int k=0; k<loads.size(); k++){
      v_x(k, 2*n)   +=  stiffness*(k_e(0, 0)*result(k,2*n) + k_e(0, 1)*result(k,2*n+1) + k_e(0, 2)*result(k,2*n+2) +
                        k_e(0, 3)*result(k,2*n+3)) - f0(k,2*n);
      v_x(k, 2*n+1) +=  -(stiffness*(k_e(2, 0)*result(k,2*n) + k_e(2, 1)*result(k,2*n+1) + k_e(2, 2)*result(k,2*n+2) +
                        k_e(2, 3)*result(k,2*n+3)) - f0(k,2*n+1));
      m_x(k, 2*n)   +=  stiffness*(k_e(1, 0)*result(k,2*n) + k_e(1, 1)*result(k,2*n+1) + k_e(1, 2)*result(k,2*n+2) +
                        k_e(1, 3)*result(k,2*n+3)) + m0(k,2*n);
      m_x(k, 2*n+1) +=  -stiffness*(k_e(3, 0)*result(k,2*n) + k_e(3, 1)*result(k,2*n+1) + k_e(3, 2)*result(k,2*n+2) +
                        k_e(3, 3)*result(k,2*n+3)) - m0(k,2*n+1);
    }
    n++;
  }

  //Result output to txt files
  ofstream fileDeflection("./results/results-deflection.txt");
  if (fileDeflection.is_open())
  {
    for (int i=0; i<u.cols(); i++){
      fileDeflection << uX(i) << "\t";
      double deflection = 0;
      for(int j=0; j<u.rows(); j++){
          deflection += u(j, i);
      }
      maxDeflection = abs(deflection)>abs(maxDeflection) ? deflection : maxDeflection;
      fileDeflection << deflection << endl;
    }
  }
  fileDeflection.close();

  ofstream fileVx("./results/results-Vx.txt");
  if (fileVx.is_open())
  {
    for (int i=0; i<v_x.cols(); i++){
      double totalVx= 0;
      for(int j=0; j<v_x.rows(); j++){
          totalVx += v_x(j, i)/1000;
      }
      maxV = abs(totalVx)>abs(maxV) ? totalVx : maxV;
      fileVx << forcesCoord(i) << "\t" << totalVx << endl;
    }
  }
  fileVx.close();

  ofstream fileMx("./results/results-Mx.txt");
  if (fileMx.is_open())
  {
    for (int i=0; i<m_x.cols(); i++){
      double totalMx= 0;
      for(int j=0; j<m_x.rows(); j++){
          totalMx += m_x(j, i)/1000000;
      }
      maxM = abs(totalMx)>abs(maxM) ? totalMx : maxM;
      fileMx << forcesCoord(i) << "\t" << totalMx << endl;
    }
  }
  fileMx.close();

  char answer;

  cout << "Maximum deflection of the beam is u_max = " << maxDeflection << "mm" << endl;
  cout << "Maximum shear force in the beam is V_max = " << maxV << "kN" << endl;
  cout << "Maximum bending moment in the beam is M_max = " << maxM << "kNm" << endl;

  cout << "Do you want to plot Vx, Mx and deflection graphs (gnuplot has to be set up on the system)? (Y/N)" << endl;
  cin >> answer;

  if (answer == 'Y' || answer == 'y'){
    //plotting graph
    #ifdef _WIN32
    FILE* pipe = _popen("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe", "w");
      if (pipe != NULL)
      {
          fprintf(pipe, "set term win\n");
          fprintf(pipe, "plot \"./results/results-deflection.txt\" with filledcurves y1=0 fs transparent lt rgb \"red\" notitle\n");
          fprintf(pipe, "set term pngcairo\n");
          fprintf(pipe, "set output \"./results/deflection.png\"\n" );
          fprintf(pipe, "replot\n");
          fprintf(pipe, "set term win\n");
          fflush(pipe);

          fprintf(pipe, "plot \"./results/results-Vx.txt\" with filledcurves y1=0 fs transparent lt rgb \"blue\" notitle \n");
          fprintf(pipe, "set term pngcairo\n");
          fprintf(pipe, "set output \"./results/Vx.png\"\n" );
          fprintf(pipe, "replot\n");
          fprintf(pipe, "set term win\n");
          fflush(pipe);

          fprintf(pipe, "plot \"./results/results-Mx.txt\" with filledcurves y1=0 fs transparent lt rgb \"green\" notitle\n");
          fprintf(pipe, "set term pngcairo\n");
          fprintf(pipe, "set output \"./results/Mx.png\"\n" );
          fprintf(pipe, "replot\n");
          fprintf(pipe, "set term win\n");
          fflush(pipe);
      }
      else puts("Could not open the file\n");
      _pclose(pipe);
    #endif

    #ifdef linux
    FILE* pipe = popen("gnuplot", "w");
      if (pipe != NULL)
      {
          fprintf(pipe, "set term wxt\n");
          fprintf(pipe, "plot \"./results/results-deflection.txt\" with filledcurves y1=0 fs transparent lt rgb \"red\" notitle\n");
          fprintf(pipe, "set term pngcairo\n");
          fprintf(pipe, "set output \"./results/deflection.png\"\n" );
          fprintf(pipe, "replot\n");
          fprintf(pipe, "set term wxt\n");
          fflush(pipe);

          fprintf(pipe, "plot \"./results/results-Vx.txt\" with filledcurves y1=0 fs transparent lt rgb \"blue\" notitle \n");
          fprintf(pipe, "set term pngcairo\n");
          fprintf(pipe, "set output \"./results/Vx.png\"\n" );
          fprintf(pipe, "replot\n");
          fprintf(pipe, "set term wxt\n");
          fflush(pipe);

          fprintf(pipe, "plot \"./results/results-Mx.txt\" with filledcurves y1=0 fs transparent lt rgb \"green\" notitle\n");
          fprintf(pipe, "set term pngcairo\n");
          fprintf(pipe, "set output \"./results/Mx.png\"\n" );
          fprintf(pipe, "replot\n");
          fprintf(pipe, "set term wxt\n");
          fflush(pipe);
      }
      else puts("Could not open the file\n");
          pclose(pipe);
    #endif
  }


}
