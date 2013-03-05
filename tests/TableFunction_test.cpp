#include <iostream>
#include "TableFunction.h"

using namespace std;

int main()
{
    TableFunction lambertw("tables/lambertw_function.dat");
    lambertw.interpolation_model = 5;
    cout << lambertw.map(28.3) << endl;
}