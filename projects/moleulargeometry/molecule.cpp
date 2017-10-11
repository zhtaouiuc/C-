#include <iostream>
#include "molecule.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>

using namespace std;

int *charges;
double *xcoord;
double *ycoord;
double *zcoord;


//constructors
molecule::molecule()
{
strcpy(title,"nomolecule");
natm=5;
charges = new int[natm];
xcoord = new double[natm];
ycoord = new double[natm];
zcoord = new double[natm];
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
length = new double * [natm];
for (int i=0; i <natm; i++)
    length[i] = new double[natm];
for (int i=0; i <natm; i++)
{
    for (int j=0; j <natm; j++)
    length[i][j] = sqrt(pow((xcoord[i]-xcoord[j]),2.0)+pow((ycoord[i]-ycoord[j]),2.0)+pow((zcoord[i]-zcoord[j]),2.0));
}
cout << "Bond length between atoms: atom number 1, atom number 2, bond length" << endl;
for (int i=0; i<natm; i++)
{
for (int j=0; j<i; j++)
    printf("%d %d %8.5f\n", i, j, length[i][j]);
}
}

