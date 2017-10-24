#ifndef MOLECULE_H_
#define MOLECULE_H_

class molecule
{
private:
    char title[50];
    int natm;
    int *charges;
    double Ia;
    double Ib;
    double Ic;
    double *xcoord;
    double *ycoord;  
    double *zcoord; 
    double **length;
    double **unitx;
    double **unity;
    double **unitz;
    double ***bangles;
    const static double masses[];
    void unitvector();
    void translation(double x,double y, double z);
public:
    molecule();
    molecule(const char *m_title);
    ~molecule();
    void showinput();
    void bondlength();
    void bondangles();
    void outofplaneangles();
    void torsionangles();
    void centerofmass();
    void momentofinertia();
    void roconstants();
};
#endif    
