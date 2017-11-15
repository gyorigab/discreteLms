#ifndef FUNCTORS_H
#define FUNCTORS_H

#include <vector>
#include <discretelms.h>

class Linear
{
public:

    Linear(unsigned int size) : m_size(size), m_distribution(m_size,1.0/m_size){}

    std::vector<double>* operator()(const std::vector<DiscreteLms::Pair> &vsort)
    {
        int index = 0;

        // pociatocna zmena pravdepodobnosti rozdelenia
        double changedDistribution  = 0.1;
        // linearny prirastok prirastok
        double diff = changedDistribution/(m_size/2);

        // najvacsia oprava ma najvacsi prirastok pravdepodobnosti
        // ze bude odstranena z nahodneho vzoru v dalsom kroku
        for(unsigned int i = vsort.size() - 1; i > m_size/2; i--)
        {
            index = vsort[i].index - 1;

            if(m_distribution[index] < 1.0)
            {
                m_distribution[index] += changedDistribution;
            }

            changedDistribution -= diff;
        }

        changedDistribution = 0.1;

        // najmensia oprava ma najvacsi ubytok pravdepodobnosti
        // ze bude odstranena z nahodneho vzoru v dalsom kroku
        for(unsigned int i = 0; i < m_size/2; i++)
        {
             index = vsort[i].index - 1;

             if(m_distribution[index] > 0.05)
             {
                m_distribution[index] -= changedDistribution;
             }
             changedDistribution -= diff;
        }
        return &m_distribution;
    }

private:
    unsigned int m_size;
    std::vector<double> m_distribution;
};

class Eexponential
{
public:

    Eexponential(unsigned int size) : m_size(size), m_distribution(m_size,1.0/m_size){}

    std::vector<double>* operator()(const std::vector<DiscreteLms::Pair> &vsort)
    {
        int index = 0;

       // pociatocna zmena pravdepodobnosti rozdelenia
       double changedDistribution  = 0.1;
       double diff = (changedDistribution/2);

       // najvacsia oprava ma najvacsi prirastok pravdepodobnosti
       // ze bude odstranena z nahodneho vzoru v dalsom kroku
       for(unsigned int i = vsort.size() - 1; i > m_size/2; i--)
       {
           index = vsort[i].index - 1;

           if(m_distribution[index] < 1.0)
           {
               m_distribution[index] += changedDistribution;
           }

           changedDistribution -= diff;
           diff /= 2;
       }

       changedDistribution = 0.1;
       diff = (changedDistribution/2);

       // najmensia oprava ma najvacsi ubytok pravdepodobnosti
       // ze bude odstranena z nahodneho vzoru v dalsom kroku
       for(unsigned int i = 0; i < m_size/2; i++)
       {
            index = vsort[i].index - 1;

            if(m_distribution[index] > 0.05)
            {
               m_distribution[index] -= changedDistribution;
            }
            changedDistribution -= diff;
            diff /= 2;
       }
       return &m_distribution;
    }

private:
    unsigned int m_size;
    std::vector<double> m_distribution;
};

class ProbabilityGroups
{
public:

    ProbabilityGroups(unsigned int size) : m_size(size), m_distribution(m_size,1.0/m_size){}

    std::vector<double>* operator()(const std::vector<DiscreteLms::Pair> &vsort)
    {
        int index = 0;

        // pociatocna zmena pravdepodobnosti rozdelenia
        double changedDistribution  = 0.1;
        double diff = 0.01;

        // 10 skupin
        int groupSize = (m_size/2) / 10;
        int j = 0;

        // najvacsia oprava ma najvacsi prirastok pravdepodobnosti
        // ze bude odstranena z nahodneho vzoru v dalsom kroku
        for(unsigned int i = vsort.size() - 1; i > m_size/2; i--)
        {
            index = vsort[i].index - 1;

            if(m_distribution[index] < 1.0)
            {
                m_distribution[index] += changedDistribution;
            }

            if(j < groupSize) j++;
            else
            {
                changedDistribution -= diff;
                j=0;
            }
        }

        changedDistribution = 0.1;
        j=0;

        // najmensia oprava ma najvacsi ubytok pravdepodobnosti
        // ze bude odstranena z nahodneho vzoru v dalsom kroku
        for(unsigned int i = 0; i < m_size/2; i++)
        {
             index = vsort[i].index - 1;

             if(m_distribution[index] > 0.05)
             {
                m_distribution[index] -= changedDistribution;
             }

             if(j < groupSize) j++;
             else
             {
                 changedDistribution -= diff;
                 j=0;
             }
        }
        return &m_distribution;
    }

private:
    unsigned int m_size;
    std::vector<double> m_distribution;
};

class Random
{
public:

    Random(unsigned int size) : m_size(size), m_distribution(m_size,1.0/m_size){}

    std::vector<double>* operator()(const std::vector<DiscreteLms::Pair> &vsort)
    {
        return &m_distribution;
    }

private:
    unsigned int m_size;
    std::vector<double> m_distribution;
};

#endif // FUNCTORS_H
