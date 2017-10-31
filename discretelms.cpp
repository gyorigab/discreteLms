#include "discretelms.h"

#include <random>
#include <algorithm>

DiscreteLms::DiscreteLms(Mat *A, Vec *b): m_distribution(A->rows(),100.0)
{
    m_A = A;
    m_b = b;

    m_A_cpy = *A;
    m_b_cpy = *b;

    m_rows = A->rows();
    m_cols = A->cols();

    m_K = generateRandomSampleSize();

    m_median = 10e8;
}

int DiscreteLms::solve()
{
    // TODO pocitat SVD pseudoinverziu

    //m_K = generateRandomSampleSize();

    generateRandomSample();
    deleteMatVecRows();


    Vec x(m_cols);
    Vec v(m_rows);

    // vytvaram normalne rovnice
    Mat N = trans(m_A_cpy) * m_A_cpy;
    Vec n = trans(m_A_cpy) * m_b_cpy;

    // riesim sustavu x = A^(-1)*b
    x = inv(N) * n;

    // pocitam vyrovnane opravy z povodnej sustavy
    v = *m_A * x - *m_b;

    cout << x;

    std::vector<Pair> vsort;
    Pair pair;

    for ( unsigned int i = 1; i <= m_rows; i++ )
    {
        // pocitam stvorec vyrovnanych oprav
        v(i) *= v(i);

        // aby som nestratil informaciu o indexoch po triedeni
        pair.val   = v(i);
        pair.index = i;

        vsort.push_back(pair);
    }

    double median = cmpMedian(vsort);

    if( median < m_median)
    {
        m_median = median;
        m_x = x;
        m_v = v;

        cout << m_x;
    }

    return 0;
}

int DiscreteLms::deleteMatVecRows()
{
    generateRandomSample();

    m_A_cpy = *m_A;
    m_b_cpy = *m_b;

    int vecIndex = 0;
    double zero = 0.0;

    for(unsigned int i = 1; i <= m_A->rows(); i++)
    {
        if( i == m_deleted_rows[vecIndex] )
        {
            for(unsigned int j = 1; j<= m_A->cols(); j++)
            {
                m_A_cpy(i,j) = zero;
            }

            m_b_cpy(i) = zero;
            vecIndex++;
        }
    }

    return 0;
}

double DiscreteLms::cmpMedian( std::vector<Pair> &vsort)
{
    int mid = vsort.size()/2;
    double median = 0.0;

    std::sort(vsort.begin(), vsort.end(), byVal());

    if(vsort.size()%2 == 0) { median = (vsort[mid-1].val + vsort[mid].val)/2; }
    else                    { median =  vsort[mid].val; }

    return median;
}


int DiscreteLms::reweightDistribution(const std::vector<Pair> &vsort)
{
    // TODO implementation of the linear distribution

    int index = 0;

    for(unsigned int i = vsort.size() - 1; i > m_rows/2; i--)
    {
        index = vsort[i].index - 1;
        m_distribution[index] -= 1;
    }

    return 0;
}

// generujem pocet odstranenych riadkov 1 <= K <= N/2-1
int DiscreteLms::generateRandomSampleSize()
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<Index> d(1, m_rows/2-1);

    return d(generator);
}

// generujem nahodny vzor na zaklade diskretneho rozdelenia
int DiscreteLms::generateRandomSample()
{
    std::random_device rd;
    std::mt19937 generator(rd());

    std::discrete_distribution<> d(m_distribution.begin(), m_distribution.end());

    m_deleted_rows.clear();

    // vygenerujem mnozinu nahodnych riadkov ktore budu vynulovane
    for(unsigned int i=1; i <= m_K; i++)
    {
        m_deleted_rows.push_back(d(generator));
    }

    return 0;
}
