#include <iostream>
#include <cmath>
#include <stdlib.h>
//this program calculates the HF energy of HeH+ using STO-3G minimal basis set

using namespace std;

//function prototype declaration
long double overlapcalc(const long double *expon, const long double *coord,
                        const long double *expon2, const long double *coord2, 
                        const long double *norma,const long double *normb);
long double kineticcalc(const long double *expon, const long double *coord,
                        const long double *expon2, const long double *coord2,
                        const long double *norma,const long double *normb);
long double expotentialcalc(const long double *expon, const long double *coord,
                            const long double *expon2, const long double *coord2,
                            const long double *norma,const long double *normb,
                            const long double *charge,const long double *nuclei);
//global terms declaration 
long double pi = 3.1415926;

int main()
{    
    const long double zeta[2] = {2.0925,1.24}; //slater exponents of He and H
    const long double charge[2] = {2.0,1.0};   // atomic charges
    const long double coord[2] = {0.0,1.4632}; // spherical coordinates 
    const long double exponsd[3] = {2.22766,0.405771,0.109818} ; //STO-3G exponents for zeta=1
//    const long double expon[2][3] = {{6.36242139,1.15892300,0.31364979},{3.42525091,0.62391373,0.16885540}};
    const long double coeff[2][3] = {{0.15432897,0.53532814,0.44463454},{0.15432897,0.53532814,0.44463454}};

// initialize matrices     
    long double overlap[2][2] = {0};
    long double expon[2][3] = {0};
    long double kinetic[2][2] = {0};
    long double expotential1[2][2] = {0}; 
    long double expotential2[2][2] = {0};    
    long double H1E[2][2] = {0};
// initialize variables
    long double norma=0.0;   //normalization coefficient of the first atomic orbital
    long double normb=0.0;   //normalization coefficient of the second atomic orbital
   
//  first scale the exponents by zeta 
    for (int i = 0; i<2; i++)
{
    for (int j = 0; j<3; j++)
    expon[i][j]=exponsd[j]*(pow(zeta[i],2));
}    

//  calculates overlap      
    for (int i = 0; i<2; i++)
{   
    for (int j = 0; j<2; j++)
{
    for (int k = 0; k<3; k++)
{
    for (int l = 0; l<3; l++)
{
    norma = pow(2.0 * (expon[i][k]) /(pi),0.75) * coeff[i][k];
    normb = pow(2.0 * (expon[j][l]) /(pi),0.75) * coeff[j][l];
    overlap[i][j]+=overlapcalc(&expon[i][k],&coord[i],&expon[j][l],&coord[j],&norma,&normb);
    kinetic[i][j]+=kineticcalc(&expon[i][k],&coord[i],&expon[j][l],&coord[j],&norma,&normb);
    expotential1[i][j]+=expotentialcalc(&expon[i][k],&coord[i],&expon[j][l],&coord[j],&norma,
                                     &normb,&charge[1],&coord[1]);
    expotential2[i][j]+=expotentialcalc(&expon[i][k],&coord[i],&expon[j][l],&coord[j],&norma,
                                     &normb,&charge[2],&coord[2]);
}}}}

// build 1-e Hamiltonian
    for (int i = 0; i<2; i++)
{
    for (int j = 0; j<2; j++)
        H1E[i][j]=kinetic[i][j]+expotential1[i][j]+expotential2[i][j];
}
    cout << expotential2[0][0] <<endl;
    cout << expotential1[0][0]<<endl;
    cout << H1E[1][1]<<endl;
    return 0;
}

long double overlapcalc(const long double *expon, const long double *coord,
                        const long double *expon2, const long double *coord2,
                        const long double *norma,const long double *normb)
{   
    long double norm;
    long double distance;
    long double overlap;
    norm = pow(pi/(*expon + *expon2),1.5);
    distance= abs(*coord - *coord2);
    overlap = norm* (*norma) * (*normb) * exp (-(*expon) * (*expon2)/(*expon + *expon2)* pow(distance,2));
    return overlap;
}

long double kineticcalc(const long double *expon, const long double *coord,
                        const long double *expon2, const long double *coord2,
                        const long double *norma,const long double *normb)
{ 
    long double distance;
    long double kinetic;
    long double norm;
    distance = abs(*coord - *coord2);
    norm = (*expon) * (*expon2)/(*expon + *expon2)*(3.0-2.0*(*expon) * (*expon2)/(*expon + *expon2)*pow(distance,2))
           * pow(pi/(*expon + *expon2),1.5);
    kinetic =norm*(*norma) * (*normb) * exp (-(*expon) * (*expon2)/(*expon + *expon2)* pow(distance,2)); 
    return kinetic;
} 

long double expotentialcalc(const long double *expon, const long double *coord,
                            const long double *expon2, const long double *coord2,
                            const long double *norma,const long double *normb,
                            const long double *charge,const long double *nuclei)
{
    long double distance;
    long double distance2;
    long double center;
    long double expotential;
    distance = abs(*coord - *coord2);
    center = ((*expon)*(*coord) + (*expon2)*(*coord2))/(*expon + *expon2);
    distance2 = abs(center-(*nuclei));
    expotential = -2.0*pi/(*expon + *expon2)*(*charge) *exp (-(*expon) * (*expon2)/(*expon + *expon2)* pow(distance,2))
                   * erf ((*expon + *expon2)*pow(distance2,2));
    return expotential;
}

