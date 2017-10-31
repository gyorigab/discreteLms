#include "discretelms.h"
#include <random>

DiscreteLms::DiscreteLms(Mat *A, Vec *b): m_distribution(A->rows(),100.0)
{
    m_A = A;
    m_b = b;

    m_A_cpy = *A;
    m_b_cpy = *b;

    m_rows = A->rows();
    m_cols = A->cols();

    m_K = generateRandomSampleSize();
}

int DiscreteLms::solve()
{

    m_K = generateRandomSampleSize();

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
