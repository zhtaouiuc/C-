#include <iostream>
#include "hartreefock.h"
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

//constructors
hartreefock::hartreefock()
{
strcpy(enucdat,"nodat")
strcpy(overlapdat,"nodat")
strcpy(kineticdat,"nodat")
strcpy(nucleardat,"nodat")
num = 5;
enuc = 0;
overlap = new double*[num];
kinetic = new double*[num];
nuclear = new double*[num];
core = new double*[num];
for (int i=0; i<natm; i++)
{
overlap[i] = new double [num];
kinetic[i] = new double [num];
nuclear[i] = new double [num];
core[i]= new double [num];
}
}

hartreefock::hartreefock(const char *m_enucdat, const char *m_overlapdat, const char *m_kineticdat,
                         const char *m_nucleardat);
{
strcpy(enucdat,"m_enucdat")
ifstream enucinput(enucdat);
enucinput >> enuc;

strcpy(overlapdat,"m_overlapdat")
ifstream overlapinput(overlapdat);

string line;
num =0;
for (int i=0; getline(overlapintput,line);++i)
num =i ;
cout<< num;

int a=0;
int b=0;
for (int i=0; i<num; i++)
{cin << a << b << overlap[a-1][b-1];
overlap[b-1][a-1] = overlap[a-1][b-1];
}

strcpy(kineticdat,"m_kineticdat")
ifstream kineticinput(kineticdat);
for (int i=0; i<num; i++)
{cin << a << b << kinetic[a-1][b-1];
kinetic[b-1][a-1] = kinetic[a-1][b-1];
}

strcpy(nucleardat,"m_nucleardat)")
ifstream nuclearinput(nucleardat);
for (int i=0; i<num; i++)
{cin << a << b << nuclear[a-1][b-1];
nuclear[b-1][a-1] = nuclear[a-1][b-1];
}
}



