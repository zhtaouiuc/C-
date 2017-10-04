#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "hf.h"

int main()
{
const long double exponsd[3] = {2.22766,0.405771,0.109818} ; //STO-3G exponents for zeta=1
const long double coeff[6] = {0.15432897,0.53532814,0.44463454,0.15432897,0.53532814,0.44463454};
const long double coord[2] = {0.0,1.4632}; // spherical coordinates
const long double zeta[2] = {2.0925,1.24}; 
long double expon[6] = {0.0};
int count=0;
    for (int i = 0; i<2; i++)
{   
    for (int j = 0; j<3; j++)
{    expon[count]=exponsd[j]*(pow(zeta[i],2));
    count+=1;
}}   

hf heh(expon,coeff,coord,2,3);
heh.show();
}
