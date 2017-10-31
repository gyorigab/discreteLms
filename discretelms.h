#ifndef DISCRETELMS_H
#define DISCRETELMS_H

#include <matvec.h>
#include <vector>

using namespace std;

typedef GNU_gama::Index Index;
typedef GNU_gama::Mat<double> Mat;
typedef GNU_gama::Vec<double> Vec;

typedef Vec::iterator Vit;

class DiscreteLms
{
public:
    DiscreteLms(Mat* A, Vec* b);

    int generateRandomSampleSize();
    int generateRandomSample();


    int solve();
    int deleteMatVecRows();

private:

    Mat* m_A;
    Vec* m_b;

    Mat m_A_cpy;
    Vec m_b_cpy;

    Vec m_x;
    Vec m_v;

    Index m_rows;
    Index m_cols;
    Index m_sample_size;

    vector<int> m_distribution;
    int m_K;

    vector<int> m_deleted_rows;

};

#endif // DISCRETELMS_H
