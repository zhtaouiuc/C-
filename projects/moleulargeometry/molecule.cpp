#include <iostream>
#include "molecule.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>

using namespace std;

//constructors
molecule::molecule()
{
strcpy(title,"nomolecule");
natm=0;
int *charges = new int[natm];
double *xcoord = new double[natm];
double *ycoord = new double[natm];
double *zcoord = new double[natm];
}

molecule::molecule(const char *m_title)
{
strcpy(title,m_title);
ifstream input(title);

input >> natm;

int *charges = new int[natm];
double *xcoord = new double[natm];
double *ycoord = new double[natm];
double *zcoord = new double[natm];

for (int i=0;i<natm;i++)
    input >> charges[i] >> xcoord[i] >> ycoord[i] >> zcoord[i];

input.close();

//for (int i=0;i<natm;i++)
//    printf("%d %20.12f %20.12f %20.12f\n", (int) charges[i],xcoord[i],ycoord[i],zcoord[i]);
}

//destructor
molecule::~molecule()
{
//delete [] charges;
//delete [] xcoord;
//delete [] ycoord;
//delete [] zcoord;
}

void molecule::showinput()
{
cout << "Number of atoms:" << natm << endl;
cout << "Input atomic charges and cartesian coordinates:" << endl;

cout << charges[0];

//for (int i=0;i<natm;i++)
 //   printf("%d %20.12f %20.12f %20.12f\n", (int) charges[i],xcoord[i],ycoord[i],zcoord[i]);
//    cout << charges[i]<<xcoord[i]<<ycoord[i]<<zcoord[i];
}


