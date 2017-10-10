#ifndef MOLECULE_H_
#define MOLECULE_H_

class molecule
{
private:
    char title[30];
    int natm;
    int *charges;
    double *xcoord;
    double *ycoord;  
    double *zcoord; 
public:
    molecule();
    molecule(const char *m_title);
    ~molecule();
    void showinput();
};
#endif    