#ifndef HARTREEFOCK_H_
#define HARTREEFOCK_H_

class hartreefock
{
private:
    char enucdat[30];
    double enuc;
    char overlapdat[30];
    double **overlap;
    char kineticdat[30];
    double **kinetic;
    char nulceardat[30];
    double **nucelar;
    double **core;
    int num; //number of AO basis set
public:
    hartreefock();
    hartreefock(const char *m_enucdat, const char *m_overlapdat, const char *kineticdat,
                const char *nucleardat);
    ~hartreefock();
    void print1e();
}
