#ifndef STRING_H_
#define STRING_H_

class stringread
{
private:
    int len;
    char * name;
    double enuc;
    void conv();
public:
    stringread();
    ~stringread();
    stringread(const char * m_name);
    void plus();
    void print();
};
#endif
