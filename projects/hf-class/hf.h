#ifndef HF_H_
#define HF_H_

class hf
{
private:
    long double exponent[100];
    long double coeff[100];
    long double overlap[100];
    long double coord[100];
    int natm;
    int num; //how many gaussian basis functions for each atomic orbital
    void overlapcalc(const long double *m_exponent, const long double *m_coeff, const long double *m_coord, int natm, int num); 
    
public:
    hf();  //default constructor
    hf(const long double *m_exponent, const long double *m_coeff, const long double *m_coord, int m_natm=0, int m_num=0);
    ~hf(); 
    void show();
};
#endif
