#include <iostream>

//This program investigates how for loop is different between i++ and ++i for a 2by3 dynamic array
using namespace std;

int main()
{
    int count=0;
    int count2=0;
    long double matrix[6]= {1.0,2.0,3.0,4.0,5.0,6.0};
    long double ** temp = new long double *[2];
    long double ** temp2 = new long double *[2];
    for (int i=0; i<2; i++)
{        temp[i] = new long double [3];
        for (int j=0; j<3; j++)
{
            temp[i][j] = matrix [count];
            count+=1;
}}
    for (int i=0; i<2; ++i)
{        temp2[i] = new long double [3];
        for (int j=0; j<3; ++j)
{           
            temp2[i][j] = matrix [count];
            count2+=1;
}}
    cout << "i++ matrix"<<endl;
    cout << temp[0][0]<<endl;
    cout << temp[0][1]<<endl;
    cout << temp[0][2]<<endl;
    cout << temp[1][0]<<endl;
    cout << temp[1][1]<<endl;
    cout << temp[1][2]<<endl;
    cout << "++i matrix"<<endl;
    cout << temp2[0][0]<<endl;
    cout << temp2[0][1]<<endl;
    cout << temp2[0][2]<<endl;
    cout << temp2[1][0]<<endl;
    cout << temp2[1][1]<<endl;
    cout << temp2[1][2]<<endl;
}
