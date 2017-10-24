#include <iostream>
#include "molecule.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#define PI 3.14159265
#define c 3.0e8
#define h 6.626e-34

using namespace std;

const double molecule::masses[] = {0.00, 1.00782503223,3.0160293201, 6.0151228874,9.012183065,10.01293695,12.0000,14.003074004,15.99491461957,18.99840316273};
//constructors
molecule::molecule()
{
strcpy(title,"nomolecule");
natm=5;
charges = new int[natm];
xcoord = new double[natm];
ycoord = new double[natm];
zcoord = new double[natm];
unitx = new double * [natm];
unity = new double * [natm];
unitz = new double * [natm];
for (int i=0; i<natm; i++)
{    unitx[i] = new double [natm];
    unity[i] = new double [natm];
    unitz[i] = new double [natm];
}
length = new double * [natm];
for (int i=0; i <natm; i++)
    length[i] = new double[natm];


bangles = new double ** [natm];
for (int i = 0; i<natm; i++)
{
bangles[i] = new double * [natm];
for (int j = 0; j<natm; j++)
bangles[i][j] = new double [natm];
}
}

molecule::molecule(const char *m_title)
{
strcpy(title,m_title);
ifstream input(title);

input >> natm;

charges = new int[natm];
xcoord = new double[natm];
ycoord = new double[natm];
zcoord = new double[natm];
unitx = new double * [natm];
unity = new double * [natm];
unitz = new double * [natm];
for (int i=0; i<natm; i++)
{    unitx[i] = new double [natm];
    unity[i] = new double [natm];
    unitz[i] = new double [natm];
}
length = new double * [natm];
for (int i=0; i <natm; i++)
    length[i] = new double[natm];

bangles = new double ** [natm];
for (int i = 0; i<natm; i++)
{
bangles[i] = new double * [natm];
for (int j = 0; j<natm; j++)
bangles[i][j] = new double [natm];
}

for (int i=0;i<natm;i++)
    input >> charges[i] >> xcoord[i] >> ycoord[i] >> zcoord[i];

input.close();
}

//destructor
molecule::~molecule()
{
delete [] zcoord;
delete [] ycoord;
delete [] xcoord;
delete [] charges;
for (int i=0; i<natm; i++)
{ delete [] unitx[i];
  delete [] unity[i];
  delete [] unitz[i];
}
 delete [] unitx; 
 delete [] unity;
 delete [] unitz;
for (int i=0; i<natm;i++)
delete [] length[i];
delete [] length;
for (int i = 0; i<natm; i++)
{
for (int j = 0; j<natm; j++)
{
delete [] bangles[i][j];
}
delete [] bangles[i];
}
delete [] bangles;

}

void molecule::showinput()
{
cout << "Number of atoms:" << natm << endl;
cout << "Input atomic charges and cartesian coordinates:" << endl;

for (int i=0;i<natm;i++)
   printf("%d %20.12f %20.12f %20.12f\n", (int) charges[i],xcoord[i],ycoord[i],zcoord[i]);
}

void molecule::bondlength()
{
for (int i=0; i <natm; i++)
{
    for (int j=0; j <natm; j++)
    length[i][j] = sqrt(pow((xcoord[i]-xcoord[j]),2.0)+pow((ycoord[i]-ycoord[j]),2.0)+pow((zcoord[i]-zcoord[j]),2.0));
}
cout << "Bond length between atoms: atom number 1, atom number 2, bond length" << endl;
for (int i=0; i<natm; i++)
{
for (int j=0; j<i; j++)
{    printf("%d %d %8.5f\n", i, j, length[i][j]);
}}
}

void molecule::unitvector()
{
for (int i=1; i<natm; i++)
{
for (int j=0; j<i; j++)
{
unitx[i][j] = -(xcoord[i]-xcoord[j])/length[i][j];
unity[i][j] = -(ycoord[i]-ycoord[j])/length[i][j];
unitz[i][j] = -(zcoord[i]-zcoord[j])/length[i][j];
unitx[j][i] = -unitx[i][j] ;
unity[j][i] = -unity[i][j] ;
unitz[j][i] = -unitz[i][j] ;
}
}
}

void molecule::bondangles()
{
unitvector();

for (int i = 0; i<natm; i++)
{
for (int j = 0; j<natm; j++)
{
for (int k = 0; k<natm; k++)
{
if (i!=j && k!=i && k!=j)
bangles[i][j][k] = acos(-unitx[i][j]*unitx[j][k]-unity[i][j]*unity[j][k]-unitz[i][j]*unitz[j][k])*180.0/PI;
}
}
}

cout<<"here is a list of bond angles" << endl;
for (int i = 2; i<natm; i++)
{
for (int j = 1; j<natm; j++)
{
for (int k = 0; k<j && k<i; k++)
{
if (j !=i && length[i][j] <4.0 && length[j][k] <4.0)
{
printf("%d %d %d %8.5f\n",i,j,k,bangles[i][j][k]);
}
}
}}
}

void molecule::outofplaneangles()
{
unitvector();
double theta=0.0;
double ejkl_x=0.0;
double ejkl_y=0.0;
double ejkl_z=0.0;

cout << "Here is a list of out of plane angles:"<< endl;
for (int i=0; i<natm; i++)
{
for (int j=0; j<natm; j++)
{
for (int k=0; k<natm;k++)
{
for (int l=0; l<j; l++)
{
if (j != i && k!=i && k!=j && l!=i && l!=k && l!=j){
if (length[i][k]<4.0 && length[j][k]<4.0 && length[k][l]<4.0 )
{
ejkl_x = (unity[k][j]*unitz[k][l] - unitz[k][j]*unity[k][l]);
ejkl_y = (unitz[k][j]*unitx[k][l] - unitx[k][j]*unitz[k][l]);
ejkl_z = (unitx[k][j]*unity[k][l] - unity[k][j]*unitx[k][l]);
theta = (ejkl_x*unitx[k][i]+ejkl_y*unity[k][i]+ejkl_z*unitz[k][i])/sin(bangles[j][k][l]/180.0*PI) ;
if (theta>1.0) {theta = asin(1.0)*180/PI;}
else if (theta<-1.0) {theta = asin(-1.0)*180/PI;}
else {theta = asin(theta)*180/PI;}
printf("%d %d %d %d %8.5f\n",i,j,k,l,theta);
}
}}
}
}
}
}

void molecule::torsionangles()
{
double cross1x=0.0;
double cross1y=0.0;
double cross1z=0.0;
double cross2x=0.0;
double cross2y=0.0;
double cross2z=0.0;
double theta=0.0;
unitvector();

cout<< "Here is a list of torsion/dihedral angles ijkl"<<endl;

for (int i=0; i<natm; i++)
{
for (int j=0; j<i; j++)
{
for (int k=0; k<j;k++)
{
for (int l=0; l<k; l++)
{
if (j != i && k!=i && k!=j && l!=i && l!=k && l!=j){
if (length[i][j]<4.0 && length[j][k]<4.0 && length[k][l]<4.0 )
{
cross1x= (unity[i][j]*unitz[j][k] - unitz[i][j]*unity[j][k]);
cross1y= (unitz[i][j]*unitx[j][k] - unitx[i][j]*unitz[j][k]);
cross1z= (unitx[i][j]*unity[j][k] - unity[i][j]*unitx[j][k]);
cross2x= (unity[j][k]*unitz[k][l] - unitz[j][k]*unity[k][l]);
cross2y= (unitz[j][k]*unitx[k][l] - unitx[j][k]*unitz[k][l]);
cross2z= (unitx[j][k]*unity[k][l] - unity[j][k]*unitx[k][l]);
theta = (cross1x*cross2x+cross1y*cross2y+cross1z*cross2z)/(sin(bangles[i][j][k]/180*PI)*sin(bangles[j][k][l]/180*PI));
if (theta>1.0) {theta = acos(1.0)*180/PI;}
else if (theta<-1.0) {theta = acos(-1.0)*180/PI;}
else {theta = acos(theta)*180/PI;}
printf("%d %d %d %d %8.5f\n",i,j,k,l,theta);
}}}}}}}

void molecule::centerofmass()
{
double masssum=0.0;
double comx=0.0;
double comy=0.0;
double comz=0.0;
int atomicnum=0;
for (int i=0; i<natm; i++)
{
atomicnum=charges[i];
masssum += masses[atomicnum]; 
comx+=masses[atomicnum] * xcoord[i]; 
comy+=masses[atomicnum] * ycoord[i];
comz+=masses[atomicnum] * zcoord[i];
}
comx=comx/masssum;
comy=comy/masssum;
comz=comz/masssum;
cout<<"Center of mass: " << comx <<"  " << comy <<"  " << comz<< endl;
translation(-comx,-comy,-comz);
}

void molecule::translation(double m_x, double m_y, double m_z)
{
for (int i =0; i<natm; i++)
{xcoord[i]+= m_x;
ycoord[i]+= m_y;
zcoord[i]+= m_z;
}
}

void molecule::momentofinertia()
{
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

Matrix I(3,3);
I(0,0) = 0;
I(1,1) = 0;
I(2,2) = 0;
I(0,1) = 0;
I(0,2) = 0;
I(1,2) = 0;
I(1,0) = 0;
I(2,0) = 0;
I(2,1) = 0;


int atomicnum=0;

for(int i=0; i<natm; i++)
{
atomicnum = charges[i];
I(0,0) += masses[atomicnum]*(pow(ycoord[i],2)+pow(zcoord[i],2));
I(1,1) += masses[atomicnum]*(pow(xcoord[i],2)+pow(zcoord[i],2));
I(2,2) += masses[atomicnum]*(pow(ycoord[i],2)+pow(xcoord[i],2));
I(0,1) += masses[atomicnum]*ycoord[i]*xcoord[i];
I(0,2) += masses[atomicnum]*zcoord[i]*xcoord[i];
I(1,2) += masses[atomicnum]*zcoord[i]*ycoord[i];
}
I(1,0) = I(0,1);
I(2,0) = I(0,2);
I(2,1) = I(1,2);
cout << "Moment of inertia tensor (amu bohr^2): " << endl;
cout << I <<endl;

Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
Matrix evecs = solver.eigenvectors();
Matrix evals = solver.eigenvalues();
cout << "Principal moments of inertia tensor (amu bohr^2): " << endl;
for(int i=0; i<3; i++)
cout << evals(i)<<"  ";

cout << "\n Principal moments of inertia tensor (amu A^2): " << endl;
for(int i=0; i<3; i++)
cout << evals(i)*0.529*0.529 << "  ";

cout << "\n Principal moments of inertia tensor (g cm^2): " << endl;
for(int i=0; i<3; i++)
cout << evals(i)*(5.29e-9)*(5.29e-9)*1.66054e-24 << "  ";

if (evals(0)<1.0 && abs(evals(2) - evals(1))<1.0e-4) {cout << "\n the molecule is a linear rotor \n" ;}
else if (abs(evals(0) - evals(2))<1.0e-4 && abs(evals(2) - evals(1))<1.0e-4) {cout << "\n the molecule is a spherical rotor \n" ;}
else if (abs(evals(0) - evals(1))<1.0e-4 && evals(1) < evals(2))  {cout << "\n the molecule is an oblate symmetric rotor \n" ;}
else if (abs(evals(2) - evals(1))<1.0e-4 && evals(0) < evals(1))  {cout << "\n the molecule is an prolate symmetric rotor \n" ;}
else {cout << "\n the molecule is an asymmetric rotor \n" ;}

Ia = evals(0);
Ib = evals(1);
Ic = evals(2);
}

void molecule::roconstants()
{

double A=0;
double B=0;
double C=0;

A=h*100.0/(8*pow(PI,2)*c*Ia*(5.29e-9)*(5.29e-9)*1.66054e-27);
B=h*100.0/(8*pow(PI,2)*c*Ib*(5.29e-9)*(5.29e-9)*1.66054e-27);
C=h*100.0/(8*pow(PI,2)*c*Ic*(5.29e-9)*(5.29e-9)*1.66054e-27);
cout << "\n rotational constants (cm-1): " << A << "  "<< B << "  " << C <<endl;
cout << "\n rotational constants (MHz): " << A*29979.2458 << "  "<< B*29979.2458 << "  " << C*29979.2458 <<endl;
}
