#include <iostream>
#include <fstream>
#include <matvec.h>
#include <discretelms.h>
#include <lsq.h>
#include <networkgenerator.h>

using namespace std;

typedef GNU_gama::Index Index;
typedef GNU_gama::Mat<double> Mat;
typedef GNU_gama::Vec<double> Vec;

int main(int argc, char *argv[])
{
    cout << "\n#####################################################\n ";
    cout << "                  Input Values                         " ;
    cout << "\n#####################################################\n\n";

    Mat *A;
    Vec *b;

    NetworkGenerator ng;

    if( argc > 2 )
    {
        int c = atoi(argv[2]);
        ng.set_outliers_count(c);
    }

    ng.generate();

    A = ng.get_A();

    cout << "Design Matrix:" << endl;
    cout << *A << endl;

    cout << "Generated standard deviations:\n";
    cout << trans(*(ng.get_e())) << endl;

    cout << "Generated outlieres:\n";
    cout << trans(*(ng.get_E())) << endl;

    cout << "Right side ";
    if( argc > 1 )
    {
         int c = atoi(argv[1]);

         if(c == 0)
         {
             cout << "(with err from normal distribution): " << endl;
             b = ng.get_b();
         }
         else if(c == 1)
         {
             cout << "(without errors):" << endl;
             b = ng.get_B();
         }
         else if(c == 2)
         {
             cout << "(with err from normal distribution and outliers):" << endl;
             b = ng.get_o();
         }
         else if(c == 3)
         {
             cout << "(with outliers err only): " << endl;
             b = ng.get_O();  // b with outliers err only
         }
         else
         {
             cout << "Uknown type of right side" << endl;
             return -1;
         }
    }
    else
    {
        cout << "(with err from normal distribution): " << endl;
        b = ng.get_o();
    }

    cout << trans(*b) << endl;

    cout << "Dim :" << A->rows() << " x " << A->cols() << endl;
    cout << "Number of outliers: " << ng.get_outliers_count() << endl;
    cout << "Outliers indexes: ";
    ng.print_outliers_indexes();

    cout << "\n#####################################################\n ";
    cout << "                  Adjustmens results                   \n" ;
    cout << "#####################################################\n\n";

    Lsq lsq(A,b);
    lsq.solve();
    lsq.print_solution();

    // Discrete LMS
    DiscreteLms dLms(A,b);

    dLms.solve();

    return 0;
}
