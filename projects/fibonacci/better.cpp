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
    {
    int sequence[3] = {0};
    int first = 1;
    int second = 1;
    for(int temp=1; temp<=n-2; temp++){
    sequence[0]=first;
    sequence[1]=second;
    sequence[2]=first+second;
    first=sequence[1];
    second= sequence[2];        
    }
    cout << "the" << n << "th number in the fibonacci sequence is" << sequence[2] << endl; 
    }
}
