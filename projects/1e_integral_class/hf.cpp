#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "hf.h"

using namespace std;
//global terms declaration
long double pi = 3.1415926;

//constructors
hf::hf() //default constructor
{
//exponent[100] = {0.0};
//coeff[100] = {0.0};
//overlap[100] = {0.0};
//coord[100] = {0.0};
natm= 0;
num = 0;
}

hf::hf(const long double *m_exponent, const long double *m_coeff, const long double *m_coord, int m_natm, int m_num)
{
natm=m_natm;
num=m_num;
for (int i =0; i<natm; i++)
{
    coord[i]=m_coord[i];
}
for(int i=0;i<num*natm;i++)
{
    exponent[i] = m_exponent[i];
    coeff[i] = m_coeff[i];
}
overlapcalc(exponent,coeff,coord,natm,num);
}

//class desctructor 
hf::~hf()
{
}

//other methods

void hf::overlapcalc(const long double *m_exponent, const long double *m_coeff, const long double *m_coord, int natm, int num)
{
//this creates 2d arrays for exponent (temp) and coeff(temp2)
int count=0;
int count2=0;
long double norma;
long double normb;
long double norm;
long double distance;
long double ** temp = new long double *[natm];
long double ** temp2 = new long double *[natm];
long double ** temp3 = new long double *[natm];
for (int i = 0; i<natm; i++)
{    temp[i]= new long double[num];
    temp2[i]= new long double[num];
    temp3[i]= new long double[natm];
    for (int j=0;j<num;j++)
{
        temp[i][j]=exponent[count];
        temp2[i][j]=coeff[count];
        count+=1;
}}
    for (int i = 0; i<natm; i++)
{   
    for (int j = 0; j<natm; j++)
{   
    for (int k = 0; k<num; k++)
{   
    for (int l = 0; l<num; l++)
{   
    distance = abs(coord[i]-coord[j]);
    norm= pow(pi/(temp[i][k] + temp[j][l]),1.5);
    norma=pow(2.0 * (temp[i][k]) /(pi),0.75) * temp2[i][k];
    normb=pow(2.0 * (temp[j][l]) /(pi),0.75) * temp2[j][l];
    temp3[i][j]+=norm* (norma) * (normb) * exp (-temp[i][k] * temp[j][l]/(temp[i][k] + temp[j][l])* pow(distance,2));
}}
    overlap[count2]=temp3[i][j];
    count2+=1;
}}
//now deletes the temporary matrices
for (int i = 0; i<natm; i++)
{   delete [] temp[i];
    delete [] temp2[i];
    delete [] temp3[i];
}
    delete [] temp;
    delete [] temp2;
    delete [] temp3;
}
    
void hf::show()
{ 
    cout << "print overlap matrix"<<endl;
    for (int i= 0; i<4; i++)
    cout <<overlap[i]<<endl;
}


