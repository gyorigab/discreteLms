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

    //m_K = generateRandomSampleSize();
    m_K = 1;

    m_median = 10e8;

    m_experiments_count = m_rows/2;
    m_end_experiment = 100;
    m_max_iter = 1000;
}

void DiscreteLms::initDistribution()
{
    m_distribution.assign(m_A->rows(),100.0);
}

int DiscreteLms::solve()
{
    // TODO pocitat SVD pseudoinverziu (momentalne regularny priklad takze OK)

    for(unsigned int expNum = 0; expNum < m_experiments_count; expNum++ )
    {
        // INFO Ales: Tu je asi zbytocne generovat m_K
        // mozem ho postupne znizovat a vyhotovit 1000 experimentov pre
        // K v intervale 0..N/2

        // generujem nahodny pocet odstranenych riadkov
        //m_K = generateRandomSampleSize();


        // INFO Ales: so zvacsujucim poctom riadkov stupa pravdepodobnost
        // ze trafim odlahle merania. To znamena ze pocet odstranenych riadkov
        // bude vzdy o dost vacsi ako je skutocny pocet odlahlych merani

        // INFO Ales: pozitivne ale je ze LMS da casto spravny vysledok aj pri pouziti
        // mensieho poctu odstranenych riadkov

        // pocet odstranenych riadkov
        m_K++;

        // cout << "Pocet odstranenych riadkov: " << m_K << endl;

        initDistribution();

        for(unsigned int expIter = 0, xconvergence = 0; expIter < m_max_iter; expIter++)
        {
            // ziskam nahodny vzor sustavy
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

            // kontrola indexov a sortovani
            /*for(auto vs : vsort)
            {
                cout << vs.index << ":" << vs.val << endl;
            }*/

            // INFO Ales: Dohodli sme sa ze budeme menit rozdeleni len pri
            // poklesu medianu ale tak to "konverguje" rychlejsie

            // prerozdelim pravdpodobnost vynulovania jednotlivych riadkov
            reweightDistribution(vsort);

            if( median < m_median)
            {
                m_median = median;
                m_x = x;
                m_v = v;

                m_sample = m_deleted_rows;

                cout << endl;

                cout << " Experiment No.: "       << expNum
                     << " << Iterace: "           << expIter
                     << " << Hodnota medianu: "   << median << endl;
            }
            else
            {
                //INFO Ales: toto je asi irelevantne

                // pokud median nekonverguje zvol novu velkost vzoru
                // a experiment opakuj
                // if(xconvergence++ == m_end_experiment) break;
            }
        }
    }

    // Pre tieto vynulovane riadky je median extremny

    cout << "\n Odstranene riadky: ";

    for(auto a: m_sample)
    {
        cout << a << " ";
    }

    cout << m_x;

    return 0;
}

int DiscreteLms::deleteMatVecRows()
{
    // vygenerujem nahodne riadky vektoru
    generateRandomSample();

    // kopirujem povodnu maticu planu a vektora pravej strany
    // do pracovnych objektov
    m_A_cpy = *m_A;
    m_b_cpy = *m_b;

    int vecIndex = 0;
    double zero = 0.0;

    for(unsigned int i = 1; i <= m_A->rows(); i++)
    {
        // pokial je riadkovy index v zozname na odstranenie
        if( i-1 == m_deleted_rows[vecIndex] )
        {
            for(unsigned int j = 1; j<= m_A->cols(); j++)
            {
                // vynulujem riadok matice planu
                m_A_cpy(i,j) = zero;
            }

            // vynulujem pvok vektora pravej strany
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

    // zoradim prvky vektora vyrovnanych oprav
    std::sort(vsort.begin(), vsort.end(), byVal());

    // urcim median
    if(vsort.size()%2 == 0) { median = (vsort[mid-1].val + vsort[mid].val)/2; }
    else                    { median =  vsort[mid].val; }

    return median;
}

// linearne rozdelenie pravdepodobnosti vyberu riadku v zavislosti na
// velkosti vyrovnanej opravy
int DiscreteLms::reweightDistribution(const std::vector<Pair> &vsort)
{
    int index = 0;

    // pociatocna zmena pravdepodobnosti rozdelenia
    double changedDistribution  = 10.0;
    // linearny prirastok prirastok
    double diff = changedDistribution/(m_rows/2);

    // najvacsia oprava ma najvacsi prirastok pravdepodobnosti
    // ze bude odstranena z nahodneho vzoru v dalsom kroku
    for(unsigned int i = vsort.size() - 1; i > m_rows/2; i--)
    {
        index = vsort[i].index - 1;

        if(m_distribution[index] < 1000.0 )
        {
            m_distribution[index] += changedDistribution;
        }

        changedDistribution -= diff;
    }

    changedDistribution = 10.0;

    // najmensia oprava ma najvacsi ubytok pravdepodobnosti
    // ze bude odstranena z nahodneho vzoru v dalsom kroku
    for(unsigned int i = 0; i < m_rows/2; i++)
    {
         index = vsort[i].index - 1;

         if(m_distribution[index] > 10.0)
         {
            m_distribution[index] -= changedDistribution;
         }
         changedDistribution -= diff;
    }

    /*cout << "Vahy: ";
    for(auto a : m_distribution )
    {
        cout << a << " \n";
    }
    cout << "\n" << endl;*/

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
