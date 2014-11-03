#include <iostream>
#include <fstream>
#include <vector>
#include "RandomVariableNDArray.h"

using namespace std;

int main()
{
  ifstream in_f("gg.dat");
  double A[261][261];
  for (int i=0; i<261; i++)
    for (int j=0; j<261; j++)
      in_f >> A[i][j];
  vector<long> shape(2);
  shape[0] = 261; shape[1] = 261;
  long indices[2];
  RandomVariableNDArray randNd(A, shape);
  for (long l=0; l<1000000; l++)
  {
    randNd.sampleAccToInvCDF(indices);
    cout << indices[0] << "   " << indices[1] << endl;
  }

}
