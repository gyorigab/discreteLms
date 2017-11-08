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
    m_max_iter = 100;

    m_max_estimation_diff = 10e-4; // [mm]
    m_outlieres_count = 0;
}

void DiscreteLms::initDistribution()
{
    m_distribution.assign(m_A->rows(),100.0);
}

void DiscreteLms::cmpFinalEstimation()
{
    bool found = false;

    // hladam stejne dobry odhad ako nejlepsi odhad ale s minimalnym
    // poctom vylucenych merani
    for(int i = m_results.size()-1; i >= 0; i--)
    {
        for(unsigned int j = 1; j <= m_cols; j++)
        {
            // testujem ci nieje nejaky podobny odhad ako vysledny s mensim
            // poctom odstranenych merani
            if(abs(m_x(j) - m_results[i](j)) > m_max_estimation_diff)
            {
                // pokial som tu tak to znamena ze tento odhad uz nieje spravny
                // a ulozim do vysledneho odhadu ten predchdzi posledny spravny
                if( i+1 < m_results.size() )
                {
                    m_x = m_results[i+1];
                    m_sample_set = m_best_sample_sets[i+1];
                }
                found = true;
                break;
            }
        }
        if(found) break;
    }

    // pokial to nenaslo ani jeden horsi odhad od odhadu
    // s najlepsim medianom tak to vyzera ze v merani bude max
    // jedno odlahle meranie
    if(!found)
    {
        m_x = m_results[0];
        m_sample_set = m_best_sample_sets[0];
    }

    cout << "\nOdstranene riadky vysledneho odhadu: ";

    for(auto a: m_sample_set)
    {
        cout << a << " ";
    }

    cout << "\nVysledny odhad: \n";
    cout << m_x;
}

int DiscreteLms::solve()
{
    // TODO pocitat SVD pseudoinverziu (momentalne regularny priklad takze OK)

    initDistribution();

    for(unsigned int expIter = 0; expIter < m_max_iter; expIter++)
    {
        // generujem nahodny pocet odstranenych riadkov
        m_K = generateRandomSampleSize();

        // ziskam nahodny vzor sustavy
        deleteMatVecRowsSet();

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

        if( median < m_median)
        {
            // prerozdelim pravdpodobnost vynulovania jednotlivych riadkov
            reweightDistribution(vsort);

            m_median = median;
            m_x = x;
            m_v = v;

            m_sample_set = m_deleted_rows_set;

            cout << endl;

            cout << " << Iterace: "           << expIter
                 << " << Hodnota medianu: "   << median << endl;
        }

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

bool DiscreteLms::exists(unsigned int index)
{
    for(auto v: m_deleted_rows)
    {
        if(v==index)
            return true;
    }
    return false;
}

int DiscreteLms::deleteMatVecRows()
{
    // vygenerujem nahodne riadky vektoru
    generateRandomSample();

    // kopirujem povodnu maticu planu a vektora pravej strany
    // do pracovnych objektov
    m_A_cpy = *m_A;
    m_b_cpy = *m_b;

    double zero = 0.0;

    for(unsigned int i = 1; i <= m_A->rows(); i++)
    {
        // INFO Ales tu bola chyba. Nevynuloval som niektore vybrane riadky

        // pokial je riadkovy index v zozname na odstranenie
        if(exists(i-1))
        {
            for(unsigned int j = 1; j<= m_A->cols(); j++)
            {
                // vynulujem riadok matice planu
                m_A_cpy(i,j) = zero;
            }

            // vynulujem pvok vektora pravej strany
            m_b_cpy(i) = zero;
        }
    }

    return 0;
}

int DiscreteLms::deleteMatVecRowsSet()
{
    // vygenerujem nahodne riadky vektoru
    generateRandomSampleSet();

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

// generujem nahodny vzor na zaklade diskretneho rozdelenia
int DiscreteLms::generateRandomSampleSet()
{
    std::random_device rd;
    std::mt19937 generator(rd());

    std::discrete_distribution<> d(m_distribution.begin(), m_distribution.end());

    m_deleted_rows_set.clear();

    // vygenerujem mnozinu nahodnych riadkov ktore budu vynulovane
    while(m_deleted_rows_set.size() <= m_K)
    {
        m_deleted_rows_set.insert(d(generator));
    }

    return 0;
}
