#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <array>
#include <sstream>
#include "harvib.h"
using namespace std;

int main()
{
harvib h2o("h2ogeo.dat","h2ohessian.dat");
h2o.showinput();
h2o.frequency();
}
