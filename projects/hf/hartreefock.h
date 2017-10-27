#ifndef HARTREEFOCK_H_
#define HARTREEFOCK_H_

class hartreefock
{
private:
    char enucdat[30];
    double enuc;
    char overlapdat[100];
    double **overlap;
    char kineticdat[30];
    double **kinetic;
    char nucleardat[30];
    double **nuclear;
    double **core;
    int num; //number of AO basis set
public:
    hartreefock();
    hartreefock(const char *m_enucdat, const char *m_overlapdat, const char *m_kineticdat,
                const char *m_nucleardat);
    ~hartreefock();
    void print1e();
};
#endif
