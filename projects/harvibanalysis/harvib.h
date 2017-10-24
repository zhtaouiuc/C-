#ifndef HARVIB_H_
#define HARVIB_H_

class harvib
{
private:
    char title[30];
    char title2[30];
    int natm;
    double *charges;
    double *xcoord;
    double *ycoord;
    double *zcoord;
    double **hessianinput;
public:
    harvib();
    harvib(const char *m_title, const char *m_title2);
    ~harvib();
    void showinput();
    void frequency();
};
#endif 
