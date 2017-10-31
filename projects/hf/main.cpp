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
hartreefock h2o("h2o_enuc.dat", "h2o_overlap.dat", "h2o_kinetic.dat", "h2o_nuclear.dat",5);
h2o.print1e();
h2o.store2e("h2o_2e.dat");
h2o.hfenergy();

hartreefock h2odz("h2o_enuc_DZ.dat", "h2o_s_DZ.dat", "h2o_k_DZ.dat", "h2o_v_DZ.dat",5);
h2odz.print1e();
h2odz.store2e("h2o_2e.dat");
h2odz.hfenergy();

}
