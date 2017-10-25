#include <iostream>
#include "harvib.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
using namespace std;
const static double masses[]={0.0000,1.00782503223, 3.0160293201,  6.0151228874,9.012183065,10.01293695,
12.000000,14.003074004,15.99491461957,18.99840316273,20.1797,22.9897,24.305,26.9815,28.0855,30.9738,32.065,35.453};

//constructors 
harvib::harvib()
{
strcpy(title, "nomolecule");
strcpy(title2, "nomolecule");
natm = 5;
charges = new double[natm];
xcoord = new double[natm];
ycoord = new double[natm];
zcoord = new double[natm];
int row =  3*pow(natm,2);
hessianinput = new double *[row];
for (int i=0; i<row; i++)
{
hessianinput[i] = new double[3];
}
}

harvib::harvib(const char *m_title, const char *m_title2)
{
strcpy(title, m_title);
ifstream coord(title);

coord >> natm;

charges = new double[natm];
xcoord = new double[natm];
ycoord = new double[natm];
zcoord = new double[natm];

for (int i=0; i<natm; i++)
coord >> charges[i] >>  xcoord[i] >> ycoord[i] >> zcoord[i];

coord.close();

strcpy(title2, m_title2);
ifstream hmatrix(title2);

hmatrix >> natm;

int row =  3*pow(natm,2);
hessianinput = new double *[row];
for (int i=0; i<row; i++)
{
hessianinput[i] = new double[3];
}

for (int i=0; i<row; i++)
hmatrix>>hessianinput[i][0]>>hessianinput[i][1]>>hessianinput[i][2];
hmatrix.close(); 
}

//destructor 
harvib::~harvib()
{
int row =  3*pow(natm,2);
delete [] zcoord;
delete [] ycoord;
delete [] xcoord;
delete [] charges;
for (int i=0; i<row; i++)
delete [] hessianinput[i];
delete [] hessianinput;
}

void harvib::showinput()
{
cout << "Number of atoms:" << natm << endl;
cout << "Input atomic charges and cartesian coordinates:" << endl;

for (int i=0;i<natm;i++)
   printf("%20.12f %20.12f %20.12f %20.12f\n", charges[i],xcoord[i],ycoord[i],zcoord[i]);

cout << "Hessian input matrix"<<endl;

int row =  3*pow(natm,2);
for (int i=0; i<row; i++)
{
    printf("%20.12f %20.12f %20.12f\n", hessianinput[i][0], hessianinput[i][1], hessianinput[i][2]);
}
}

void harvib::frequency()
{
double **hessian = new double * [natm*3];
for (int i=0; i<natm*3; i++)
hessian[i] = new double [natm*3];
int count=0;
for (int i=0; i< 3*natm ; i++)
{
for (int k=0; k<3*natm; k++)

{hessian[i][k] = 0.0;
hessian[i][k] = hessianinput[(int) count/3][k%3]; 

//cout<< "i  " << i << "  k  " << k << "  count/3  " << count/3 << "  k%3  " << k%3 <<endl ;
count+=1;	

}

}


//for (int i=0; i<natm*3;i++)
//cout << hessian[i][i]<<endl;


//cout << "test" << hessian[1][1]<<endl;
for (int i =0; i< 3*natm; i++)
{
for (int j=0; j <3*natm; j++)
{
hessian[i][j] = hessian[i][j] / pow(masses[(int) charges[(i/3)]]*masses[(int) charges[(j/3)]],0.5);
//cout<< "i  " << i << "  j  " << j << "  charges[(i/3)  " << (int) charges[(i/3)] << "  charges[(j/3)  " << (int)charges[(j/3)] <<endl ;

}
}

//for (int i=0; i<natm*3;i++)
//cout << hessian[i][i]<<endl;


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
Matrix I (3*natm,3*natm);
for (int i=0; i<3*natm; i++)
{
for (int j=0; j<3*natm; j++)
{I(i,j) = 0.0;
I(i,j) = hessian[i][j];
}}

for (int i =0; i<3*natm; i++)
delete [] hessian[i];
delete [] hessian;


//cout << I<< endl;
cout << "Here is a list of frequencies of harmonic vibrational modes (cm-1):" <<endl;
Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();

// the eigenvalues are in the units of hartree/(bohr)^2/(amu); to get cm-1, we need to convert amu to mass of electrons
// first and then take a square root of hartree^2 to then convert hartree to cm-1. 1amu = 1822.89 electron rest mass

for (int i=0; i<3*natm; i++)
{evals(i) = pow(evals(i)/1822.89,0.5) * 219474;
cout << evals(i) << endl;
}
}
