#ifndef DISCRETELMS_H
#define DISCRETELMS_H

#include <matvec.h>
#include <vector>
#include <set>
#include <functional>

using namespace std;

typedef GNU_gama::Index Index;
typedef GNU_gama::Mat<double> Mat;
typedef GNU_gama::Vec<double> Vec;

typedef Vec::iterator Vit;
typedef std::vector<double>::const_iterator CIT;

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
    int generateRandomSampleSet(CIT beg, CIT end);
    int generateRandomSampleSet();

    //int solve(void (*rewieghtFunction)(const std::vector<Pair> &, DiscreteLms *));
    int solve(std::function< std::vector<double>*(const std::vector<Pair> &vsort)> reweight);

    int deleteMatVecRowsSet(CIT beg, CIT end);

    bool exists(unsigned int index);

    double cmpMedian(std::vector<Pair> &vsort);
/*  static void reweightLinear(const std::vector<Pair> &vsort, DiscreteLms *pdlms);
    static void reweightExponential(const std::vector<Pair> &vsort, DiscreteLms *pdlms);
    static void probabilityGroups(const std::vector<Pair> &vsort, DiscreteLms *pdlms);
    static void randomSamples(const std::vector<Pair> &vsort, DiscreteLms *pdlms);*/

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
    vector<unsigned int> m_deleted_rows;
    set<unsigned int> m_deleted_rows_set;

    vector<unsigned int> m_sample;
    set<unsigned int> m_sample_set;

    unsigned int m_max_iter;

    double m_max_estimation_diff;

    double m_median;


};

#endif // DISCRETELMS_H
