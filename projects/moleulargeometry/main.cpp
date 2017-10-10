#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <array>
#include <sstream>
#include "molecule.h"
using namespace std;

// this program calculates the internal coordinates, moment of inertia,
// and rotational constants of a polyatomic molecule. 

int main()
{
molecule one("geo.dat");
one.showinput();
return 0;
}
