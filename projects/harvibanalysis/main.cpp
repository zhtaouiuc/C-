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

harvib c6h6("benzenegeo.dat","benzenehessian.dat");
c6h6.showinput();
c6h6.frequency();

harvib clc4h4("clbutenegeo.dat", "clbutenehessian.dat");
clc4h4.showinput();
clc4h4.frequency();

}
