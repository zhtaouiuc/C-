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

const double energyconv = 1.0e-5;
const double densityconv = 1.0e-5;


//constructors
hartreefock::hartreefock()
{
strcpy(enucdat,"nodat");
strcpy(overlapdat,"nodat");
strcpy(kineticdat,"nodat");
strcpy(nucleardat,"nodat");
num = 5;
enuc = 0.0;
Eelec1 = 0.0;
Eelec2 = 0.0;
iter = 0;
super = 0;
overlap = new double*[num];
kinetic = new double*[num];
nuclear = new double*[num];
core = new double*[num];
densa = new double*[num];
densb = new double *[num];
fock = new double *[num];
index = new int [super];
ortho = new double *[num];
for (int i=0; i<num; i++)
{
overlap[i] = new double [num];
kinetic[i] = new double [num];
nuclear[i] = new double [num];
core[i]= new double [num];
densa[i] = new double [num];
densb[i] = new double [num];
fock [i] = new double [num];
ortho[i] = new double [num];
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

iter = 0;
double *values = new double [nline];
int * first = new int[nline];
int *second = new int[nline];
overlapinput.clear();
overlapinput.seekg(0,ios::beg);

for (int i=0; i<nline; i++)
overlapinput >> first[i]>> second[i] >> values[i];

num = first[nline-1];
super = num*(num+1)/2+num;
index = new int [super];

ortho = new double *[num];
overlap = new double*[num];
kinetic = new double*[num];
nuclear = new double*[num];
core = new double*[num];
densa = new double*[num];
densb = new double *[num];
fock = new double *[num];
for (int i=0; i<num; i++)
{
overlap[i] = new double [num];
kinetic[i] = new double [num];
nuclear[i] = new double [num];
core[i]= new double [num];
densa[i] = new double [num];
densb[i] = new double [num];
fock [i] = new double [num];
ortho[i] = new double [num];
}


for (int i=0; i<nline-1; i++)
{overlap[first[i]][second[i]]= values[i];
overlap[second[i]][first[i]] = values[i];
cout << "first[i] " << first[i] << "second[i] " << second[i] << "values[i] " << values[i] << endl;
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

delete [] first;
delete [] second;
}

hartreefock::~hartreefock()
{
for (int i=0; i<num; i++)
{
delete [] overlap[i];
delete [] kinetic[i];
delete [] nuclear[i];
delete [] core[i];
delete [] densa[i];
delete [] densb[i];
delete [] fock[i];
delete [] ortho[i];
}
delete [] overlap;
delete [] kinetic;
delete [] nuclear;
delete [] core;
delete [] densa;
delete [] densb;
delete [] fock;
delete [] index;
delete [] ortho;
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

cout << "\n this is the overlap matrix:\n";
for (int i=0; i<num; i++)
{
for (int j=0; j<num; j++)
{
cout << "  " << overlap[i][j] << "  ";
}
cout << endl;
}

}

void hartreefock::store2e(const char *m_2edat)
{
FILE *input;
double val;
int i =0;
int j=0;
int k=0;
int l=0;
int ij =0;
int kl =0;
int ijkl =0;
TEI = new double [super*(super+1)/2+super];

index[0] = 0;
for (int i=1; i< super; i++)
index[i] = index[i-1] +i ;
 
input = fopen(m_2edat,"r");
while (fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF)
{
ij = (i>j) ? index[i]+j : index[j]+i;
kl = (k>l) ? index[k]+l : index[l]+l;
ijkl = (ij>kl) ? index[ij]+kl : index[kl]+ij;
//cout << " i " <<i<< " j " <<j<<" ij "<< ij <<" k " << k << " l " <<l<< " kl "<<kl<<" ijkl " <<ijkl << endl;
TEI[ijkl] = val;
//cout<< " ijkl " <<ijkl << " val " << val<<endl;
}
}


void hartreefock::matrixmu(double **ma, double **mb, double **mc,int rownuma, int colnuma, int rownumb, int colnumb)
{
if (colnuma != rownumb) {cout << "error in the multiplication."<<endl;}
for (int i=0; i<rownuma; i++)
{
for (int j=0; j<colnumb; j++)
{mc[i][j] = 0.0;
{for (int k=0; k<colnuma; k++)
mc[i][j] += ma[i][k] * mb[k][j];
}
}
}
}
void hartreefock::orthomatrix()
{

double **eigenvalues = new double *[num];
double **eigenvectors = new double *[num];

for (int i = 0; i< num ; i++)
{eigenvalues[i]= new double [num];
eigenvectors[i] = new double [num];
}

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
Matrix I (num, num);
for (int i =0; i<num; i++)
{
for (int j=0; j<num; j++)
I(i,j) = overlap[i][j];
}
Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();

cout << evals<<endl;
for (int i = 0 ; i<num; i++)
{
eigenvalues[i][i] = pow(evals(i),-0.5);
cout << "evals(i) " << evals(i) << " eigenvalues[i][i] " << eigenvalues[i][i]<<endl;
for (int j=0; j<num; j++)
eigenvectors[i][j]= evecs(i,j);
}


double **orthotemp = new double *[num];
for (int i = 0; i< num ; i++)
orthotemp[i] = new double [num];

matrixmu(eigenvalues,eigenvectors,orthotemp, num, num, num, num);
matrixmu(eigenvectors, orthotemp, ortho, num,num,num, num);

for (int i = 0; i< num ; i++)
{
for (int j=0; j<num; j++)
cout << "  " << eigenvalues[i][j] << "  ";
cout << endl;
}
for (int i = 0; i< num ; i++)
{
delete [] eigenvalues[i];
delete [] eigenvectors[i];
delete [] orthotemp[i];
}
delete [] eigenvalues;
delete [] eigenvectors; 
delete [] orthotemp;
}

void hartreefock::fockbuild()
{

int ij =0;
int kl =0;
int ijkl =0;
int ik =0;
int jl=0;
int ikjl = 0;

if (iter ==0)
{for (int i =0; i<num; i++)
{
for (int j=0; j<num; j++)
{
fock[i][j] = 0.0;
fock[i][j] = core[i][j];
}}}

else
{
for (int i =0; i<num; i++)
{
for (int j=0; j<num; j++)
{fock[i][j] = core[i][j];
for (int k=0; k<num; k++)
{
for (int l=0; l<num; l++)
{
ij = (i>j) ? index[i]+j : index[j]+i;
kl = (k>l) ? index[k]+l : index[l]+l;
ijkl = (ij>kl) ? index[ij]+kl : index[kl]+ij;
cout << " i " <<i<< " j " <<j<<" ij "<< ij <<" k " << k << " l " <<l<< " kl "<<kl<<" ijkl " <<ijkl << endl;
ik = (i>k) ? index[i]+k : index[k]+i;
jl = (j>l) ? index[j]+l : index[l]+j;
ikjl = (ik>jl) ? index[ik]+jl : index[jl]+ik;
fock[i][j]+= densa[k][l]*2*(TEI[ijkl]-TEI[ikjl]);
}}}}
}
}

void hartreefock::den()
{

orthomatrix();
fockbuild();
double ** tempfock = new double *[num];
double ** tempfock2 = new double *[num];
for (int i =0; i<num; i++)
{
tempfock[i] = new double [num];
tempfock2[i] = new double [num];
}

matrixmu(fock,ortho,tempfock,num,num,num,num);
matrixmu(ortho,tempfock,tempfock2,num,num,num,num);

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
Matrix F0 (num,num);

for (int i =0; i<num; i++)
{
for (int j=0; j<num; j++)
F0(i,j) = tempfock2[i][j];
}

Eigen::SelfAdjointEigenSolver<Matrix> solver(F0);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();

double **eigenvectors = new double *[num];
double **temp = new double *[num];
for (int i = 0; i< num ; i++)
{
eigenvectors[i] = new double [num];
temp[i] = new double [num];
for (int j=0; j<num; j++)
{
eigenvectors[i][j] = evecs(i,j);
}}

matrixmu(ortho,eigenvectors,temp,num,num,num,num);

for (int i=0; i< num; i++)
{
for (int j=0; j< num; j++)
{
for (int k=0; k< nocc; k++)
densb[i][j]+= temp[i][k]*temp[k][j];
}
}

for (int i = 0; i< num ; i++)
{
delete [] tempfock[i];
delete [] tempfock2[i];
delete [] eigenvectors[i];
delete [] temp[i];
}
delete [] tempfock;
delete [] tempfock2;
delete [] eigenvectors;
delete [] temp;

convergence();
}

void hartreefock::convergence()
{
double dendiff=0;

for (int i=0; i< num; i++)
{
for (int j=0; j< num; j++)
{Eelec2 += densb[i][j]  * (core[i][j]+fock[i][j]); 
dendiff += pow(pow(densb[i][j]-densa[i][j],2),0.5);
}}

while (dendiff>densityconv && abs(Eelec2-Eelec1)>energyconv &&iter<100)
{
iter+=1;
for (int i=0; i< num; i++)
{
for (int j =0; j< num; j++)
densa[i][j] = densb[i][j];
}
Eelec1 = Eelec2; 
den();
}

if (iter >100) {cout << " the calculation is not converged."<< endl;}
else
{cout<< "after # iteration = " << iter<< endl;
cout << "the total energy of your system is " << Eelec2+enuc << endl;
}
}

