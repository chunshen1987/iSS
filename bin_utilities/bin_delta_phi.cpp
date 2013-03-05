#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "Table.h"


using namespace std;

void cos2_col7(vector<double> *data)
{
  double pT = (*data)[5];
  double phi = (*data)[6];
  double 
  (*data)[6] = cos(2*(*data)[6]);
}


int main()
{
  ifstream in_f1("results/samples_211.dat");
  ofstream out_f1("results/samples_211_bin_eta.dat");
  ifstream in_f2("results/samples_321.dat");
  ofstream out_f2("results/samples_321_bin_eta.dat");
  ifstream in_f3("results/samples_2212.dat");
  ofstream out_f3("results/samples_2212_bin_eta.dat");
  void (*cos2_col7_func)(vector<double>*) = cos2_col7;

  Table bins("eta_bins_gauss.dat");
  get_bin_average_and_count(in_f1, out_f1, bins.getColumn(1), 4, cos2_col7_func);
  get_bin_average_and_count(in_f2, out_f2, bins.getColumn(1), 4, cos2_col7_func);
  get_bin_average_and_count(in_f3, out_f3, bins.getColumn(1), 4, cos2_col7_func);
  //vector <double> bins(2);
  //bins[0] = 0; bins[1] = 5.0;
  //get_bin_average_and_count(in_f, out_f, &bins, 5, cos2_col7_func);
}
