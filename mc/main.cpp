#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    ifstream AFileInput;

    if(argc > 1)
    {
        AFileInput.open(argv[1]);
    }

    //AFileInput.open("../data_lsq/illc1033.mtx");

    int m,n,nz;

    AFileInput >> m >> n >> nz;

    vector< vector<double> > v;


    for(int i = 0; i<m; i++)
    {
        vector<double> tmp;
        for(int j = 0; j<n; j++)
        {
            tmp.push_back(0);
        }
        v.push_back(tmp);
    }

    int k,l;
    double val;

    while(AFileInput >> k >> l >> val )
    {
        v[k-1][l-1] = val;
    }

    cout << m << " " << n << endl;

    for(int i = 0; i<v.size(); i++)
    {
        for(int j = 0; j<v[i].size(); j++)
        {
            cout << v[i][j] << " ";
        }
        cout << "\n";
    }

    return 0;
}
