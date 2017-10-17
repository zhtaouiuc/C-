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
    double **length;
//    double **unitx;
//    double **unity;
//    double **unitz;
//    double ***bangles;
//    void unitvector();
public:
    molecule();
    molecule(const char *m_title);
    ~molecule();
    void showinput();
    void bondlength();
    void unitvector();
//    void bondangles();
};
#endif    
