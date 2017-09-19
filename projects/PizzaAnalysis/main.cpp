// this program stores data for a pizza-analyssi service
#include<iostream>

struct pizza{
    char name[20];
    double diameter;
    double weight;
    double dollar;
};

int main()
{ 
    using namespace std;
    int number;
    cout<< "Enter how many companies data you want to store: \n";
    cin >> number;
    
    pizza data[number];
    for(int i=0; i<=number-1; i++)
    {
    cout<< "Name of the company: \n";
    cin >> data[i].name;
    cout << "diameter of the pizza: \n";
    cin >> data[i].diameter;
    cout << "weight of the pizza:\n";
    cin >> data[i].weight;
    cout << "price of the pizza in $ \n";
    cin >> data[i].dollar;
    cout << "End of this data set. \n";
    }

    for(int j=0; j<=number-1; j++)
    {
    cout << data[j].name << " company provides " << data[j].diameter <<"-inches " << data[j].weight
    << "-lb pizza for " << data[j].dollar << " dollars. \n";
    }

}     
