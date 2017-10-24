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
one.bondlength();
one.bondangles();
one.outofplaneangles();
one.torsionangles();
one.centerofmass();
one.momentofinertia();
one.roconstants();      

molecule acetaldehyde("acetaldehyde.dat");
acetaldehyde.showinput();
acetaldehyde.bondlength();
acetaldehyde.bondangles();
acetaldehyde.outofplaneangles();
acetaldehyde.torsionangles();   
acetaldehyde.centerofmass();
acetaldehyde.momentofinertia();  
acetaldehyde.roconstants();      

molecule benzene("benzene.dat");
benzene.showinput();
benzene.bondlength();
benzene.bondangles();
benzene.outofplaneangles();
benzene.torsionangles();
benzene.centerofmass();
benzene.momentofinertia();
benzene.roconstants(); 

molecule allene("allene.dat");
allene.showinput();            
allene.bondlength();
allene.bondangles();
allene.outofplaneangles();
allene.torsionangles();
allene.centerofmass();
allene.momentofinertia();
allene.roconstants();     

return 0;
} 
