#ifndef NETWORKGENERATOR_H
#define NETWORKGENERATOR_H

#include <matvec.h>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

typedef GNU_gama::Index Index;
typedef GNU_gama::Mat<double> Mat;
typedef GNU_gama::Vec<double> Vec;

typedef Vec::iterator Vit;

class NetworkGenerator
{
public:

    NetworkGenerator()
    {
        m_outliers_count = 0;
        m_outliers_indexes.clear();
        m_set_by_user = false;
    }

    void reset()
    {
        m_outliers_count = 0;
        m_outliers_indexes.clear();
        m_set_by_user = false;
    }

    int generate()
    {
        std::random_device rd;
        std::mt19937 generator(rd());

        // standardne chyby z normalneho rozdelenia
        std::normal_distribution<>           norm(0,0.01);
        // hrube chyby s normalneho rozdelenia
        std::normal_distribution<>           normo(10,100);
        // Koeficienty matice planu
        std::uniform_real_distribution<>     rdist(-0.5,0.5);
        // Nahodna velkost sustavy
        std::uniform_int_distribution<Index> idist(8,100);

        Index M = idist(generator);

        // Pocet neznamych parametrov
        std::uniform_int_distribution<Index> idist_ns(2,Index(M/2)-1);

        Index k = idist_ns(generator);
        Index N = Index(M/2)-k;

        // Pocet hrubych chyb v N prvych meraniach
        std::uniform_int_distribution<Index> outliers_count_dist(1,Index(M/2)-1);

        if ( M == N ) {  return -1; }
        if ( M < N ) std::swap(M,N);

        if( !m_set_by_user )
            m_outliers_count = outliers_count_dist(generator);

        Mat A(M,N);
        Vec X(N);
        Vec e(M), E(M);

        E.set_zero();
        e.set_zero();

        // generate first design matrix
        for (Index i=1; i<=M; i++)
            for (Index j=1; j<=N; j++)
                A(i,j) = rdist(generator);

        for (Index i=1; i<=N; i++) X(i) = i;

        m_b = A*X;
        m_B = m_b;

        for (Index i=1; i<=M; i++)
        {
            e(i) = norm(generator);
        }

         m_outliers_indexes.clear();

        for(Index i=1; i <= m_outliers_count; i++ )
        {
            E(i) = normo(generator);
            m_outliers_indexes.push_back(i);
        }

        m_b += e; // add error to meassurement

        m_o = m_b;
        m_O = m_B;

        m_o += E; // add error to meassurement with erros
        m_O += E; // add error to meassurement without erros

        m_A = A;
        m_e = e;
        m_E = E;

        return 0;
    }

    Mat* get_A(){ return &m_A; }
    Vec* get_b(){ return &m_b; }
    Vec* get_B(){ return &m_B; }
    Vec* get_O(){ return &m_O; }
    Vec* get_o(){ return &m_o; }
    Vec* get_e(){ return &m_e; }
    Vec* get_E(){ return &m_E; }

    Index get_outliers_count(){ return m_outliers_count; }
    Index get_cols(){ return m_cols; }
    Index get_rows(){ return m_rows; }

    void set_outliers_count( Index outliers_count )
    {
        m_outliers_count = outliers_count;
        m_set_by_user = true;
    }

    void print_outliers_indexes()
    {
        for(unsigned int i = 0; i < m_outliers_indexes.size(); i++ )
        {
            cout << m_outliers_indexes[i] << " ";
        }
        cout <<"\n";
    }

private:

    Mat m_A;
    Vec m_b;
    Vec m_B;
    Vec m_o;
    Vec m_O;
    Vec m_e;
    Vec m_E;

    Index m_rows;
    Index m_cols;
    Index m_outliers_count;
    std::vector<Index> m_outliers_indexes;

    bool m_set_by_user;
};

#endif // NETWORKGENERATOR_H
