#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "ParameterReader.h"
#include "Table.h"


using namespace std;

int phi_col_idx; // indexed starting from 0
int calculate_vn_to_order;
int total_number_of_columns;

void operation(vector<double> *data)
// Call-back function. So far only calculate flows.
{
    int write_idx = 1;
    for (int order=1; order<=calculate_vn_to_order; order++)
    {
        (*data)[total_number_of_columns-1+write_idx] = cos(order*(*data)[phi_col_idx]); // -1: indices start with 0
        write_idx ++;
        (*data)[total_number_of_columns-1+write_idx] = sin(order*(*data)[phi_col_idx]); // -1: indices start with 0
        write_idx ++;
    }
}

void write_format()
// Write a format file corresponding to the new data calculated by operation function.
{
  ofstream vn_format("results/bin_pT_format.dat");
  int write_idx = 1;
  for (int order=1; order<=calculate_vn_to_order; order++)
  {
    vn_format << "v" << order << "_real = " << total_number_of_columns+write_idx << endl;
    write_idx ++;
    vn_format << "v" << order << "_imag = " << total_number_of_columns+write_idx << endl;
    write_idx ++;
  }
  vn_format << "total_N = " << total_number_of_columns+write_idx << endl;
  write_idx ++;
  vn_format << "dN_dpt = " << total_number_of_columns+write_idx << endl;
  write_idx ++;
  vn_format.close();
}

int main()
{
  // initialization
  char in_filename_buffer[300], out_filename_buffer[300];
  void (*operation_ptr)(vector<double>*) = operation;
  Table bins("tables/pT_bins.dat");
  ParameterReader iSS_para;
  iSS_para.readFromFile("parameters.dat");
  ParameterReader format;
  format.readFromFile("results/samples_format.dat");
  Table chosen_particles("EOS/chosen_particles.dat");

  // get parameters
  int pT_col_idx = format.getVal("pt");
  total_number_of_columns = format.getVal("total_number_of_columns");
  phi_col_idx = format.getVal("phi") - 1; // indexed starting from 0
  calculate_vn_to_order = iSS_para.getVal("calculate_vn_to_order");

  // loop over particles
  for (int i=1; i<=chosen_particles.getNumberOfRows(); i++)
  {
    cout << "----- For particle with MC-index " << (int) chosen_particles.get(1,i) << " -----" << endl;
    // prepare files
    sprintf(in_filename_buffer, "results/samples_%d.dat", (int) chosen_particles.get(1,i));
    sprintf(out_filename_buffer, "results/samples_%d_bin_pt.dat", (int) chosen_particles.get(1,i));
    ifstream in_f(in_filename_buffer);
    ofstream out_f(out_filename_buffer);
    get_bin_average_and_count(in_f, out_f, bins.getColumn(1), pT_col_idx, operation, total_number_of_columns+calculate_vn_to_order*2);
  }

  // write out a format file
  write_format();

  //vector <double> bins(2);
  //bins[0] = 0; bins[1] = 5.0;
  //get_bin_average_and_count(in_f, out_f, &bins, 5, cos2_col7_func);

}
