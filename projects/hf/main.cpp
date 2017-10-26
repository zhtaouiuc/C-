#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <array>
#include <sstream>
#include "hartreefock.h"
using namespace std;

int main()
{
hartreefock h2o("h2o_enuc.dat", "h2o_overlap.dat", "h2o_kinetic.dat", "h2o_nuclear.dat");
h2o.print1e();
}
