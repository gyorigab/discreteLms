#include <iostream>
#include <fstream>
#include <matvec.h>
#include <discretelms.h>

using namespace std;

typedef GNU_gama::Index Index;
typedef GNU_gama::Mat<double> Mat;
typedef GNU_gama::Vec<double> Vec;

int main(int argc, char *argv[])
{

    ifstream AFileInput, bFileInput;

    AFileInput.open("../data/A");
    bFileInput.open("../data/b");

    Mat A;
    Vec b;

    AFileInput >> A;
    bFileInput >> b;

    DiscreteLms dLms(&A,&b);

    dLms.solve();

    return 0;
}
