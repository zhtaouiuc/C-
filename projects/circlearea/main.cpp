#include <iostream>
#include <cmath>
//this program calculates the area of circle with an input radius
int main()
{
    using namespace std;
    long double r;
    long double pi=3.1415926;
    long double area;
    cout << "Enter the radius of your circle \n";
    cin >> r;
    if (r>=0)
    {
    area = pi * pow(r,2);
    cout << setprecision (8) << "the area of the circle is \n" << area << endl;}
    else 
       cout << "You cannot have a negative radius \n";
    return 0;
}
