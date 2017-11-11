#ifndef LSQ_H
#define LSQ_H

#include <vector>
#include <matvec.h>

using namespace std;

typedef GNU_gama::Index Index;
typedef GNU_gama::Mat<double> Mat;
typedef GNU_gama::Vec<double> Vec;

typedef Vec::iterator Vit;

class Lsq
{
public:
    Lsq(Mat* A, Vec* b)
    {
        m_A = A;
        m_b = b;

        m_rows = A->rows();
        m_cols = A->cols();
    }

    void solve()
    {
        Mat N = trans(*m_A ) * *m_A ;
        Vec n = trans(*m_A ) * *m_b;

        m_x = inv(N) * n;
        m_v = *m_A * m_x - *m_b;
    }

    void print_solution()
    {
        cout << "LSQ: " << endl;
        cout << "X: " << endl;
        cout << trans(m_x) << endl;

        cout << "v: " << endl;
        cout << trans(m_v) << endl;
    }

    Vec* get_x(){ return &m_x;}
    Vec* get_v(){ return &m_v;}

private:
    Mat* m_A;
    Vec* m_b;

    Vec m_x;
    Vec m_v;

    Index m_rows;
    Index m_cols;
};

#endif // LSQ_H
