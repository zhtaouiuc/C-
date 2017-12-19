#include <iostream>
#include "string.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <string>
using namespace std;

int main()
{
stringread tomato("h2o_enuc.dat");
stringread cucumber("h2o_enuc_DZ.dat");
tomato.print();
cucumber.print();
tomato.plus();
cucumber.plus();
}
