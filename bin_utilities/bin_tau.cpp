#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "Table.h"
#include "ParameterReader.h"

using namespace std;

int phi_col_idx;
int calculate_vn_to_order;
int total_number_of_columns;

void operation(vector<double> *data)
// Call-back function. So far only calculate flows.
{
  for (int order=1; order<=calculate_vn_to_order; order++)
      (*data)[total_number_of_columns+order-1] = cos(order*(*data)[phi_col_idx]);
}

void write_format()
// Write a format file corresponding to the new data calculated by operation function.
{
  ofstream vn_format("results/bin_tau_format.dat");
  for (int order=1; order<=calculate_vn_to_order; order++)
  {
    vn_format << "v" << order << " = " << total_number_of_columns+order << endl;
  }

  vn_format << "total_N = " << total_number_of_columns+calculate_vn_to_order+1 << endl;
  vn_format << "dN_dtau = " << total_number_of_columns+calculate_vn_to_order+2 << endl;
  vn_format.close();
}

int main()
{
  // initialization
  char in_filename_buffer[300], out_filename_buffer[300];
  void (*operation_ptr)(vector<double>*) = operation;
  Table bins("tables/tau_bins.dat");
  ParameterReader iSS_para;
  iSS_para.readFromFile("parameters.dat");
  ParameterReader format;
  format.readFromFile("results/samples_format.dat");
  Table chosen_particles("EOS/chosen_particles.dat");

  // get parameters
  total_number_of_columns = format.getVal("total_number_of_columns");
  phi_col_idx = format.getVal("phi");
  int tau_col_idx = format.getVal("tau");
  calculate_vn_to_order = iSS_para.getVal("calculate_vn_to_order");

  // loop over particles
  for (int i=1; i<=chosen_particles.getNumberOfRows(); i++)
  {
    cout << "----- For particle with MC-index " << (int) chosen_particles.get(1,i) << " -----" << endl;
    // prepare files
    sprintf(in_filename_buffer, "results/samples_%d.dat", (int) chosen_particles.get(1,i));
    sprintf(out_filename_buffer, "results/samples_%d_bin_tau.dat", (int) chosen_particles.get(1,i));
    ifstream in_f(in_filename_buffer);
    ofstream out_f(out_filename_buffer);
    get_bin_average_and_count(in_f, out_f, bins.getColumn(1), tau_col_idx, operation, total_number_of_columns+calculate_vn_to_order);
  }

  // write out a format file
  write_format();

  //vector <double> bins(2);
  //bins[0] = 0; bins[1] = 5.0;
  //get_bin_average_and_count(in_f, out_f, &bins, 5, cos2_col7_func);
}
