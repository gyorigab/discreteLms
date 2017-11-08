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

    // b2 - dve odlahle merania 43,42 (indexovane od 0)
    // b8 - 8 odlahlych merani 43..43-8

    // spravny vysledok x : 0.00318222, 0.00249310

    //AFileInput.open("../discreteLms/data/A");
    //bFileInput.open("../discreteLms/data/b2");

    if(argc > 2)
    {
        AFileInput.open(argv[1]);
        bFileInput.open(argv[2]);
    }

    Mat A;
    Vec b;

    AFileInput >> A;
    bFileInput >> b;

    DiscreteLms dLms(&A,&b);

    dLms.solve();

    return 0;
}
