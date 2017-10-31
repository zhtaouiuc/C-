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

const double energyconv = 1.0e-7;
const double densityconv = 1.0e-7;


//constructors
hartreefock::hartreefock()
{
strcpy(enucdat,"nodat");
strcpy(overlapdat,"nodat");
strcpy(kineticdat,"nodat");
strcpy(nucleardat,"nodat");
num = 5;
nocc = 0;
enuc = 0.0;
Eelec1 = 0.0;
Eelec2 = 0.0;
iter = 0;
super = 0;
TEI = new double [super*(super+1)/2+super];

overlap = new double*[num];
kinetic = new double*[num];
nuclear = new double*[num];
core = new double*[num];
densa = new double*[num];
densb = new double *[num];
fock = new double *[num];
index = new int [super];
ortho = new double *[num];
transortho = new double *[num];
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
transortho[i] = new double [num];
}
}

hartreefock::hartreefock(const char *m_enucdat, const char *m_overlapdat, const char *m_kineticdat,
                         const char *m_nucleardat, int m_nocc)
{
nocc = m_nocc;
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
super = num*(num+1)/2+num+1;
index = new int [super];
TEI = new double [super*(super+1)/2+super];

ortho = new double *[num];
transortho = new double *[num];
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
transortho[i] = new double [num];
}



for (int i=0; i<nline; i++)
{overlap[first[i]-1][second[i]-1]= values[i];
overlap[second[i]-1][first[i]-1] = values[i];

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
delete [] transortho[i];
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
delete [] transortho;
delete [] TEI;
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

index[0] = 0;
for (int i=1; i< super; i++)
index[i] = index[i-1] +i ;
 
input = fopen(m_2edat,"r");
while (fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF)
{
ij = (i>j) ? index[i]+j : index[j]+i;
kl = (k>l) ? index[k]+l : index[l]+k;
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

void hartreefock::transpose(double **ma, double **mb, int row, int col)
//mb gives the transposed ma
{
for (int i=0; i< row; i++)
{
for (int j=0; j< col; j++)
{
mb[i][j] = ma[j][i];
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

//cout << evals<<endl;
for (int i = 0 ; i<num; i++)
{
eigenvalues[i][i] = pow(evals(i),-0.5);
//cout << "evals(i) " << evals(i) << " eigenvalues[i][i] " << eigenvalues[i][i]<<endl;
for (int j=0; j<num; j++)
eigenvectors[i][j]= evecs(i,j);
}


double **orthotemp = new double *[num];
double **eigenvectrans = new double *[num];
for (int i = 0; i< num ; i++)
{
orthotemp[i] = new double [num];
eigenvectrans[i] = new double [num];
}
transpose(eigenvectors,eigenvectrans,num,num);
matrixmu(eigenvalues,eigenvectrans,orthotemp, num, num, num, num);
matrixmu(eigenvectors, orthotemp, ortho, num,num,num, num);

transpose(ortho,transortho,num,num);
//for (int i = 0; i< num ; i++)
//{
//for (int j=0; j<num; j++)
//cout << "  " << ortho[i][j] << "  ";
//cout << endl;
//}
for (int i = 0; i< num ; i++)
{
delete [] eigenvalues[i];
delete [] eigenvectors[i];
delete [] orthotemp[i];
delete [] eigenvectrans[i];
}
delete [] eigenvalues;
delete [] eigenvectors; 
delete [] orthotemp;
delete [] eigenvectrans;
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
{
fock[i][j] = 0.0;
fock[i][j] = core[i][j];
for (int k=0; k<num; k++)
{
for (int l=0; l<num; l++)
{
ij = (i>j) ? index[i+1]+j+1 : index[j+1]+i+1;
kl = (k>l) ? index[1+k]+l+1 : index[l+1]+k+1;
ijkl = (ij>kl) ? index[ij]+kl : index[kl]+ij;
//cout << " i " <<i<< " j " <<j<<" ij "<< ij <<" k " << k << " l " <<l<< " kl "<<kl<<" ijkl " <<ijkl << endl;
ik = (i>k) ? index[i+1]+k+1 : index[k+1]+i+1;
jl = (j>l) ? index[j+1]+l+1 : index[l+1]+j+1;
ikjl = (ik>jl) ? index[ik]+jl : index[jl]+ik;
fock[i][j]+= densa[k][l]*(2.0*TEI[ijkl]-TEI[ikjl]);
//cout << " densa[k][l] " << densa[k][l] << " TEI[ijkl] " << TEI[ijkl] << " TEI[ikjl] " << TEI[ikjl] << " fock[i][j] " << fock[i][j] <<endl;
}}}}
}
}




void hartreefock::den()
{

fockbuild();
double ** tempfock = new double *[num];
double ** tempfock2 = new double *[num];

for (int i =0; i<num; i++)
{
tempfock[i] = new double [num];
tempfock2[i] = new double [num];
}

matrixmu(fock,ortho,tempfock,num,num,num,num);
matrixmu(transortho,tempfock,tempfock2,num,num,num,num);
//
//cout << "fock" << endl;
//for (int i = 0; i< num ; i++)
//{
//for (int j=0; j<num; j++)
//cout << "  " << fock[i][j] << "  ";
//cout << endl;
//}
//
//cout << "ortho" << endl;
//for (int i = 0; i< num ; i++)
//{
//for (int j=0; j<num; j++)
//cout << "  " << ortho[i][j] << "  ";
//cout << endl;
//}


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
densb[i][j] = 0;
for (int k=0; k< nocc; k++)
{densb[i][j]+= temp[i][k]*temp[j][k];
//cout << " temp[i][k] "<< temp[i][k] << " temp[k][j] " << temp[k][j]<< " i "<< i << "j " << j << " densb[i][j] " << densb[i][j]<< endl;
}}
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

Eelec2 = 0.0;
for (int i=0; i< num; i++)
{
for (int j=0; j< num; j++)
{
Eelec2 += densb[i][j]  * (core[i][j]+fock[i][j]); 
dendiff += pow(pow(densb[i][j]-densa[i][j],2),0.5);
}}

printf("%d %20.12f %20.12f %20.12f %20.12f\n", iter, Eelec2, Eelec2+enuc,  abs(Eelec2-Eelec1),dendiff);

//cout << "densitya from last iteration" << endl;
//for (int i = 0; i< num ; i++)
//{
//for (int j=0; j<num; j++)
//cout << "  " << densa[i][j] << "  ";
//cout << endl;
//}

if ((dendiff > densityconv) || (abs(Eelec2-Eelec1)>energyconv))
{
if (iter >100) {cout << " the calculation is not converged."<< endl;
exit (EXIT_FAILURE);}

iter+=1;
for (int i=0; i< num; i++)
{
for (int j =0; j< num; j++)
{densa [i][j] = 0.0;
densa[i][j] = densb[i][j];
}}

//cout << "densityb" << endl;
//for (int i = 0; i< num ; i++)
//{
//for (int j=0; j<num; j++)
//cout << "  " << densb[i][j] << "  ";
//cout << endl;
//}

Eelec1 = 0.0;
Eelec1 = Eelec2; 
den();
}

}


void hartreefock::hfenergy()
{
orthomatrix();
printf("%s %10s %20s %20s %20s\n", "Iter", "E(elec)", "E(tot)", "Ediff", "Dendiff");
den();
}
