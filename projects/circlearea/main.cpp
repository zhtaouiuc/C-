#include <iostream>
#include <cmath>
//this program calculates the area of circle with an input radius
int main()
{
    using namespace std;
    double r;
    double pi=3.1415926;
    double area;
    cout << "Enter the radius of your circle \n";
    cin >> r;
    area = pi * pow(r,2);
    cout << "the area of the circle is \n" << area << endl;
    return 0;
}
