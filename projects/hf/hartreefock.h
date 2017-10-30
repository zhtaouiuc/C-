#ifndef HARTREEFOCK_H_
#define HARTREEFOCK_H_

class hartreefock
{
private:
// build Hcore
    char overlapdat[100];
    double **overlap;
    char kineticdat[30];
    double **kinetic;
    char nucleardat[30];
    double **nuclear;
    double **core;    
    int num; //number of AO basis set

// read in nuclear energy
    char enucdat[30];
    double enuc;

// read in two electron integral
    double *TEI;
    int *index;   
    int super;
 
// build orthogonal matrix 
    double **ortho;
    void orthomatrix();

// build fock matrix
    double **fock; 
    void fockbuild();

// build density matrix
    int nocc; //number of occupied orbitals
    int iter; //number of iteration 
    double **densa;
    double **densb;
    
// calculate energy/convergence test
    double Eelec1;
    double Eelec2;
    void convergence();

// matrix multiplication
    void matrixmu(double **ma, double **mb, double **mc,int rownuma, int colnuma, int rownumb, int colnumb);
      
public:
    hartreefock();
    hartreefock(const char *m_enucdat, const char *m_overlapdat, const char *m_kineticdat,
                const char *m_nucleardat);
    ~hartreefock();
    void print1e();  //build Hcore
    void store2e(const char *m_2edat);  //store 2e integral
    void den();
};
#endif
