#include <iostream>
#include "hartreefock.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include "../../resources/Eigen/Dense"
#include "../../resources/Eigen/Eigenvalues"
#include "../../resources/Eigen/Core"
#include <string>
using namespace std;

//constructors
hartreefock::hartreefock()
{
strcpy(enucdat,"nodat");
strcpy(overlapdat,"nodat");
strcpy(kineticdat,"nodat");
strcpy(nucleardat,"nodat");
num = 5;
enuc = 0.0;
overlap = new double*[num];
kinetic = new double*[num];
nuclear = new double*[num];
core = new double*[num];
for (int i=0; i<num; i++)
{
overlap[i] = new double [num];
kinetic[i] = new double [num];
nuclear[i] = new double [num];
core[i]= new double [num];
}
}

hartreefock::hartreefock(const char *m_enucdat, const char *m_overlapdat, const char *m_kineticdat,
                         const char *m_nucleardat)
{
strcpy(enucdat, m_enucdat);

ifstream enucinput(enucdat);
enucinput >> enuc;
//cout << enuc;
enucinput.close();


strcpy(overlapdat,m_overlapdat);
string *line = new string[30];
int nline=0;
ifstream overlapinput(overlapdat);
if (overlapinput.is_open())
{ while (getline(overlapinput, line[nline]) && line[nline].size()>0)
 {nline++;}
}


double *values = new double [nline];
int * first = new int[nline];
int *second = new int[nline];

overlapinput.clear();
overlapinput.seekg(0,ios::beg);

for (int i=0; i<nline; i++)
overlapinput >> first[i]>> second[i] >> values[i];

num = first[nline-1];

overlap = new double*[num];
kinetic = new double*[num];
nuclear = new double*[num];
core = new double*[num];
for (int i=0; i<num; i++)
{
overlap[i] = new double [num];
kinetic[i] = new double [num];
nuclear[i] = new double [num];
core[i]= new double [num];
}


for (int i=0; i<num; i++)
{overlap[first[i]][second[i]]= values[i];
overlap[second[i]][first[i]] = values[i];
}
overlapinput.close();

int a=0;
int b=0;
strcpy(kineticdat,m_kineticdat);
ifstream kineticinput(kineticdat);
for (int i=0; i<nline; i++)
{kineticinput>> a >> b >> kinetic[a-1][b-1];
kinetic[b-1][a-1] = kinetic[a-1][b-1];
}

kineticinput.close();

strcpy(nucleardat,m_nucleardat);
ifstream nuclearinput(nucleardat);
for (int i=0; i<nline; i++)
{nuclearinput >> a >> b >> nuclear[a-1][b-1];
nuclear[b-1][a-1] = nuclear[a-1][b-1];
//cout << " a " << a << " b " << b << endl;
}

nuclearinput.close();
}

hartreefock::~hartreefock()
{
for (int i=0; i<num; i++)
{
delete [] overlap[i];
delete [] kinetic[i];
delete [] nuclear[i];
delete [] core[i];
}
delete [] overlap;
delete [] kinetic;
delete [] nuclear;
delete [] core;
}


void hartreefock::print1e()
{
cout << "\n this is the core matrix:\n"; 
for (int i=0; i<num; i++)
{
for (int j=0; j<num; j++)
{
core[i][j] = kinetic[i][j] + nuclear[i][j];
cout << "  " << core[i][j] << "  ";
}
cout << endl;
}

//for (int i=0; i<num; i++)
//cout << core[i][i];
}









