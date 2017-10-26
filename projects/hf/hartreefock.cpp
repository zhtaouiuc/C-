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
cout << "what:";
strcpy(enucdat, m_enucdat);
cout << enucdat;

ifstream enucinput(enucdat);
enucinput >> enuc;
cout << enuc;
enucinput.close();

strcpy(overlapdat,m_overlapdat);
ifstream overlapinput(overlapdat);

string line;
num =0;

while(getline(overlapinput, line))
{num +=1;
cout<< num;}

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

int a=0;
int b=0;
for (int i=0; i<num; i++)
{overlapinput >> a >> b >>overlap[a-1][b-1];
overlap[b-1][a-1] = overlap[a-1][b-1];
}

overlapinput.close();

strcpy(kineticdat,m_kineticdat);
ifstream kineticinput(kineticdat);
for (int i=0; i<num; i++)
{kineticinput>> a >> b >> kinetic[a-1][b-1];
kinetic[b-1][a-1] = kinetic[a-1][b-1];
}

kineticinput.close();

strcpy(nucleardat,m_nucleardat);
ifstream nuclearinput(nucleardat);
for (int i=0; i<num; i++)
{nuclearinput >> a >> b >> nuclear[a-1][b-1];
nuclear[b-1][a-1] = nuclear[a-1][b-1];
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
cout << core[i][j];
}
cout << endl;
}
}








