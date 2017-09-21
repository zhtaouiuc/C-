#include <iostream>
#include <cmath>
#include <stdlib.h>
//this program calculates the HF energy of HeH+ using STO-3G minimal basis set
using namespace std;
long double overlapcalc(const long double *expon, const long double *coeff, const long double *coord,
const long double *expon2, const long double *coeff2, const long double *coord2);
long double sum(const long double *expon, const long double *coeff);
long double pi = 3.1415926;

int main()
{    
    const long double zeta[2] = {2.0925,1.24}; //slater exponents of He and H
    const long double charge[2] = {2.0,1.0};   // atomic charges
    const long double coord[2] = {0.0,1.4632}; // spherical coordinates 
    const long double expon[2][3] = {{6.36242139,1.15892300,0.31364979},{3.42525091,0.62391373,0.16885540}};
    const long double coeff[2][3] = {{0.15432897,0.53532814,0.44463454},{0.15432897,0.53532814,0.44463454}};
    
    long double overlap[2][2] = {0};

//    cout<< overlapcalc(&expon[1][1],&coeff[1][1],&coord[1],&expon[1][1],&coeff[1][1],&coord[1])<< endl;
//    cout<< overlapcalc(&expon[1][1],&(coeff[1][1]),&coord[1],&(expon[1][1]),&(coeff[1][1]),&coord[2])<< endl;
    
//    cout << sum(&expon[0][0], &coeff[0][0])<<endl;
//    abort() ;
    
    for (int i = 0; i<2; i++)
{   
    for (int j = 0; j<2; j++)
{
    for (int k = 0; k<3; k++)
{
    for (int l = 0; l<3; l++)
    overlap[i][j]+=overlapcalc(&expon[i][k],&coeff[i][k],&coord[i],&expon[j][l],&coeff[j][l],&coord[j]);
}}}    
    cout << overlap[0][0] <<endl;
    cout << overlap[0][1]<<endl;
    cout << overlap[1][1]<<endl;

}

long double sum(const long double *expon, const long double *coeff)

{
   long double sum;

   sum = *expon + *coeff;
   return sum;
}

long double overlapcalc(const long double *expon, const long double *coeff, const long double *coord,
const long double *expon2, const long double *coeff2, const long double *coord2)
{   
    long double norma;
    long double normb;
    long double norm;
    long double distance;
    long double overlap;
    norma = pow(2.0 * (*expon) /(pi),0.75);
    normb = pow(2.0 * (*expon2) /(pi),0.75);
    norm = pow(pi/(*expon + *expon2),1.5);
    distance= abs(*expon - *expon2);
    overlap = norm* norma * normb * exp (-(*expon) * (*expon2)/(*expon + *expon2)* pow(distance,2));
    return overlap;
}
