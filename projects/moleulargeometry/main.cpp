#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <array>
#include <sstream>
using namespace std;

// this program calculates the internal coordinates, moment of inertia,
// and rotational constants of a polyatomic molecule. 

int main()
{
ifstream input("geo.dat");

int natm;
input >> natm;

int *zcharges = new int[natm];
double *x = new double[natm];
double *y = new double[natm];
double *z = new double[natm];



for (int i=0;i<natm;i++)

    input >> zcharges[i] >> x[i] >>y[i]>>z[i];
 
cout <<natm<<endl;
cout << x[0] <<x[1]<<x[2]<<x[3];
}
