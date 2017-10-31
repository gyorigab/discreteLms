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

    struct Pair
    {
        double  val;
        int index;
    };

    struct byVal
    {
        bool operator()(const Pair &left,const Pair &right)
        {
            return left.val < right.val;
        }
    };

    DiscreteLms(Mat* A, Vec* b);

    int generateRandomSampleSize();
    int generateRandomSample();

    int solve();
    int deleteMatVecRows();

    double cmpMedian(std::vector<Pair> &vsort);
    int reweightDistribution(const std::vector<Pair> &vsort);

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

    unsigned int m_K;
    vector<int> m_distribution;
    vector<int> m_deleted_rows;

    double m_median;

};

#endif // DISCRETELMS_H
