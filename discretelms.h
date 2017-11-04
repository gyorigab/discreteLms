#ifndef DISCRETELMS_H
#define DISCRETELMS_H

#include <matvec.h>
#include <vector>
#include <set>

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
    int generateRandomSampleSet();
    int generateRandomSample();

    int solve();
    int deleteMatVecRows();
    int deleteMatVecRowsSet();

    bool exists(unsigned int index);

    double cmpMedian(std::vector<Pair> &vsort);
    int reweightDistribution(const std::vector<Pair> &vsort);

    void initDistribution();

    void cmpFinalEstimation();

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
    vector<double> m_distribution;
    vector<unsigned int> m_deleted_rows;
    set<unsigned int> m_deleted_rows_set;

    vector<unsigned int> m_sample;
    set<unsigned int> m_sample_set;

    vector<Vec> m_results;
    vector< set<unsigned int> > m_best_sample_sets;

    unsigned int m_outlieres_count;

    unsigned int m_experiments_count;
    unsigned int m_end_experiment;
    unsigned int m_max_iter;

    double m_max_estimation_diff;

    double m_median;

};

#endif // DISCRETELMS_H
