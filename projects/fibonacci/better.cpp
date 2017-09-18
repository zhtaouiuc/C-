#include <iostream>
#include <cmath>
//this program takes in an integer n and prints out the nth number in the fibonacci sequence

int fibonacci(int);

int main()
{
    using namespace std;
    int n;
    cout << "Enter the nth number you hope to print out in the fibonacci sequence \n";
    cin >> n;
    if (n<1)
        cout << "n has to integers greater than 0 \n";
    else 
        cout << "the" << n << "th number in the fibonacci sequence is" << fibonacci(n) << endl; 
}

int fibonacci(int number)
{   
    int sequence[3]={};
    int first =1 ;
    int second =1 ;
    
    for(int temp=1; temp<=number-2; temp++){
    sequence[0]=first;
    sequence[1]=second;
    sequence[2]=first+second;
    first=sequence[0]+sequence[1];
    second= sequence[1]+sequence[2];
}
    return first+second;
}
