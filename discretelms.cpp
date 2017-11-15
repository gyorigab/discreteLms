#include "discretelms.h"

#include <random>
#include <algorithm>
#include <matvec/pinv.h>
#include <functors.h>
#include <functional>

DiscreteLms::DiscreteLms(Mat *A, Vec *b)
{
    m_A = A;
    m_b = b;

    m_A_cpy = *A;
    m_b_cpy = *b;

    m_rows = A->rows();
    m_cols = A->cols();

    m_K = generateRandomSampleSize();

    m_median = 10e8;
    m_max_iter = 1000;
}

int DiscreteLms::solve(std::function< std::vector<double>*(const std::vector<DiscreteLms::Pair> &vsort)> reweight)
{
    std::vector<double> *distribution;
    m_K = generateRandomSampleSize();
    generateRandomSampleSet();

    for(unsigned int expIter = 0; expIter < m_max_iter; expIter++)
    {
        Vec x(m_cols);
        Vec v(m_rows);

        // vytvaram normalne rovnice
        Mat N = trans(m_A_cpy) * m_A_cpy;
        Vec n = trans(m_A_cpy) * m_b_cpy;

        // riesim sustavu x = A^(-1)*b SVD - pseudoinverze
        x = pinv(N) * n;

        // pocitam vyrovnane opravy z povodnej sustavy
        v = *m_A * x - *m_b;

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
            // prerozdelim pravdpodobnost vynulovania jednotlivych riadkov

            distribution = reweight(vsort);

            m_median = median;
            m_x = x;
            m_v = v;

            m_sample_set = m_deleted_rows_set;

            cout << " << Iterace: "           << expIter
                 << " << Hodnota medianu: "   << median << endl;
        }

        m_K = generateRandomSampleSize();

        // ziskam nahodny vzor sustavy
        deleteMatVecRowsSet(distribution->begin(), distribution->end());
    }

    cout << "\nOdstranene riadky najlepsieho odhadu: ";

    for(auto a: m_sample_set)
    {
        cout << a << " ";
    }

    cout << "\nOdhad s min medianom: \n";
    cout << m_x;

    return 0;
}

int DiscreteLms::deleteMatVecRowsSet(CIT beg, CIT end)
{
    // vygenerujem nahodne riadky vektoru
    generateRandomSampleSet(beg, end);

    // kopirujem povodnu maticu planu a vektora pravej strany
    // do pracovnych objektov
    m_A_cpy = *m_A;
    m_b_cpy = *m_b;

    set<unsigned int>::const_iterator it = m_deleted_rows_set.begin();
    double zero = 0.0;

    for(unsigned int i = 1; i <= m_A->rows(); i++)
    {
        // pokial je riadkovy index v zozname na odstranenie
        if( i-1 == *it && it != m_deleted_rows_set.end() )
        {
            for(unsigned int j = 1; j<= m_A->cols(); j++)
            {
                // vynulujem riadok matice planu
                m_A_cpy(i,j) = zero;
            }

            // vynulujem pvok vektora pravej strany
            m_b_cpy(i) = zero;
            it++;
        }
    }
    return 0;
}

double DiscreteLms::cmpMedian( std::vector<Pair> &vsort)
{
    int mid = vsort.size()/2;
    double median = 0.0;

    // zoradim prvky vektora vyrovnanych oprav
    std::sort(vsort.begin(), vsort.end(), byVal());

    // urcim median
    if(vsort.size()%2 == 0) { median = (vsort[mid-1].val + vsort[mid].val)/2; }
    else                    { median =  vsort[mid].val; }

    return median;
}


// generujem pocet odstranenych riadkov 1 <= K <= N/2-1
int DiscreteLms::generateRandomSampleSize()
{
    std::random_device rd;
    std::mt19937 generator(rd());

    int upperBound = m_rows/2-1 > m_cols ? m_rows/2-1 : m_cols-1;
    std::uniform_int_distribution<Index> d(1, upperBound);

    return d(generator);
}

// generujem nahodny vzor na zaklade diskretneho rozdelenia
int DiscreteLms::generateRandomSampleSet(CIT beg, CIT end)
{
    std::random_device rd;
    std::mt19937 generator(rd());

    std::discrete_distribution<> d(beg, end);

    m_deleted_rows_set.clear();

    // vygenerujem mnozinu nahodnych riadkov ktore budu vynulovane
    while(m_deleted_rows_set.size() <= m_K)
    {
        m_deleted_rows_set.insert(d(generator));
    }

    return 0;
}

// generujem nahodny vzor na zaklade diskretneho rozdelenia
int DiscreteLms::generateRandomSampleSet()
{
    std::random_device rd;
    std::mt19937 generator(rd());

    std::uniform_int_distribution<Index> d(1, m_rows);

    m_deleted_rows_set.clear();

    // vygenerujem mnozinu nahodnych riadkov ktore budu vynulovane
    while(m_deleted_rows_set.size() <= m_K)
    {
        m_deleted_rows_set.insert(d(generator));
    }

    return 0;
}
