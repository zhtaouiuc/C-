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

typede::Matrix<double, Eigen::Dynamicgen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;, Eigen::Dynamic, Eigen::RowMajor> Matrix;

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
}

void hartreefock::store2e(const char *m_2edat)
{
FILE *input;
double val;
int super = num*(num+1)/2+num;
int *index = new int [super];
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

Matrix I (num, num);
for (int i =0; i<num; i++)
{
for (int j=0; j<num; j++)
I(i,j) = overlap[i][j];
}
Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();
for (int i = 0 ; i<num; i++)
{
eigenvalues[i][i] = pow(evals(i),-0.5);
for (int j=0; j<num; j++)
eigenvectors[i][j]= evecs(i,j);
}

ortho = new double *[num];
for (int i = 0; i< num ; i++)
ortho[i] = new double [num];

matrixmu(eigenvalues,eigenvectors,ortho, num, num, num, num);
matrixmu(eigenvectors, ortho, ortho, num,num,num, num);

}

void hartreefock::deninitial();
{
if (iter == 0)
{
double **fockinitial = new double *[num];
for (int i=0; i< num; i++)
fockinitial[i] = new double [num];

matrixmu(core,ortho,fockinitial,num,num,num,num);
matrixmu(ortho,fockinitial,num,num,num,num);

Matrix F0 (num,num);

for (int i =0; i<num; i++)
{
for (int j=0; j<num; j++)
F0(i,j) = fockinitial[i][j];
}

Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();

double **eigenvectorsini = new double *[num];
double **temp = new double *[temp];
for (int i = 0; i< num ; i++)
{
eigenvectorsini[i] = new double [num];
temp[i] = new double [num];
for (int j=0; j<num; j++)
{
eigenvectorsini[i][j] = evecs(i,j);
}}

matrixmu(ortho,eigenvectorsini,temp,num,num,num,num);


double **densa = new double *[num];
for (int i=0; i< num; i++)
densa[i] = new double [num];

for (int i=0; i< num; i++)
{
for (int j=0; j< num; j++)
{
for (int k=0; k< nocc; k++)
densa[i][j]+= temp[i][k]*temp[k][j];
}
}
}

else 
{

}
}
