#include <iostream>
#include <cmath>
#include <stdlib.h>
//this program tests the concept of passing array by reference
//to do this, we store 2 2*2 matrix into 2 1d array and we pass
// them to a function to add 2 matrix to return a 2*2 matrix stored
// in 2 1d array.

using namespace std;

void sum(long double *c,const long double *a, const long double *b,int row,int column);
 
int main()
{
    long double a[4]={1.0,2.0,3.0,4.0};
    long double b[4]={1.0,3.0,2.0,1.0};
    long double c[4]={0};
    sum(c,a,b,2,2);
    for(int i=0; i<4;i++)
    cout << c[i] << endl;
}

void sum(long double *c,const long double *a, const long double *b,int row,int column)
{    
    for (int i=0; i<row+column; i++)
        c[i]=a[i]+b[i];
}    
