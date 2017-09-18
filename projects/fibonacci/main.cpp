#include <iostream>
#include <cmath>
//this program takes in an integer n and prints out the nth number in the fibonacci sequence
// only works up to 46th number in the fibonacci sequence
int fibonacci(int);

int main()
{
    using namespace std;
    int n;
    cout << "Enter the nth number you hope to print out in the fibonacci sequence \n";
    cin >> n;
    if (n<1)
        cout << "n has to integers greater than 0 \n";
    else if (n>46)
        cout << "too bad...this little program cannot give that many..";
    else 
        cout << "the" << n << "th number in the fibonacci sequence is" << fibonacci(n) << endl; 
}

int fibonacci(int number)
{   
    int *sequence;
    sequence= new int[number-1];
    sequence[0]=1;
    sequence[1]=1;    

    for(int temp=1; temp<=number-2; temp++){
    sequence[1+temp]=sequence[1+temp-1]+sequence[1+temp-2];
    }
    return sequence[number-1];
    delete [] sequence;
}
