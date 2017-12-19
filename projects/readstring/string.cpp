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

stringread::stringread()
{
len = 4;
name = new char[4];
strcpy(name,"nostring");
enuc=0.0;
}

stringread::stringread(const char * m_name)
{
len = strlen(m_name);
name = new char[len+1];
strcpy(name,m_name);
ifstream enucinput(name);
if (enucinput.is_open())
{enucinput >> enuc;
enucinput.close();
}
}

stringread::~stringread()
{
cout<<"calling destructor";
delete [] name;
}

void stringread::print()
{
cout<< name<< endl;
cout<<enuc<< endl;
}

void stringread::plus()
{
enuc = enuc+1;
cout<<enuc<<endl;
conv();
}

void stringread::conv()
{
if (enuc<10)
{
plus();
}
}
