#include <iostream>
#include "molecule.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#define PI 3.14159265

using namespace std;

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
cout<<"destructor is called";
delete [] zcoord;
delete [] ycoord;
delete [] xcoord;
delete [] charges;
for (int i=0; i<natm; i++)
{ delete [] unitx[natm];
  delete [] unity[natm];
  delete [] unitz[natm];
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
double theta;
double ejkl_x;
double ejkl_y;
double ejkl_z;

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
