// Ver 2.1
// Note the version used in iSS is slightly slower than the one
// used in iS, do not mix them unless necessary.

#ifndef EMISSIONFUNCTION_H
#define EMISSIONFUNCTION_H

#include<string>
#include<vector>

#include "Table.h"
#include "NBD.h"
#include "ParameterReader.h"
#include "main.h"

using namespace std;

class EmissionFunctionArray
{
private:
  Table *pT_tab, *phi_tab, *y_minus_eta_tab;
  int pT_tab_length, phi_tab_length, y_minus_eta_tab_length, y_minus_eta_min_index; // the y_minus_eta_min_index holds the index to the smallest positive y-eta_s value
  bool positive_y_minus_eta_table_only; // true if y_minus_eta_tab has only positive part
  ParameterReader *paraRdr; // used to pass-in parameters
  
  Table *dN_pTdpTdphidy; // dN / (pt dpt dphi dy)
  Table *dN_pTdpTdphidy_max; // store the largest element when summing over xt and eta to get dN / (pt dpt dphi dy); used during sampling.
  double **dN_dxtdetady; // dN / (d^2x_t deta dy) (correction: and dtau)
  double **dN_dxtdetady_pT_max; // dN / (d^2x_t deta dy) is the (weighted) sum of all elements of the [pT_tab_length][phi_tab_length]-sized dN/dall matrix; this array records the maximum value of this dN/dall matrix, for each eta and FO-cell.
  double **dN_dxtdy_4all; // dN / (dxt dy) for all particles
  
  int number_of_chosen_particles;
  int *chosen_particles_01_table; // used in spectra and flow calculations; it has length Nparticle, 0 means miss, 1 means include
  int grouping_particles; // 0/1: Particles with similar mass and chemical potentials will be sampled using the same dN/(dxt deta dy) matrix
  double mc_grouping_tolerance; // Usable only when grouping_particles is 1. If two particles have mass and chemical potentials close within this relative tolerance, they are considered to be identical and will be sampled successively without regenerating the dN / (dxt deta dy) matrix for efficiency.
  int *chosen_particles_sampling_table; // store particle index; the sampling process follows the order specified by this table
  // list for information for all particles
  int Nparticles;
  particle_info* particles;
  // list for information for all fluid cells
  long FO_length;
  FO_surf* FOsurf_ptr;
  int last_particle_idx; // store the last particle index being used by calculate_dNArrays function
  double **trig_phi_table, **hypertrig_y_minus_eta_table;
  inline long determine_number_to_sample(double dN, int model=1, double para1=0, double para2=0, double para3=0, double para4=0, double para5=0);
  NBD nbd; // NBD random sample generator
  bool particles_are_the_same(int, int);
  long *sorted_FZ;
public:
  EmissionFunctionArray(Table* chosen_particle, Table* pt_tab_in, Table* phi_tab_in, Table* eta_tab_in, particle_info* particles_in, int Nparticles, FO_surf* FOsurf_ptr_in, long FO_length_in, ParameterReader* paraRdr_in);
  ~EmissionFunctionArray();

  void calculate_dNArrays(int);
  void calculate_dN_dxtdetady(int);
  void calculate_dN_pTdpTdphidy(int);
  void write_dN_pTdpTdphidy_toFile();
  string dN_pTdpTdphidy_filename; // where to save
  void write_dN_dxtdetady_toFile();
  string dN_dxtdetady_filename;

  void calculate_flows(int to_order, string flow_diff_filename, string flow_inte_filename);
  string flow_differential_filename_old, flow_integrated_filename_old;
  string flow_differential_filename, flow_integrated_filename;

  void calculate_dN_pTdpTdphidy_and_flows_4all_old_output(int perform_sampling=0);
  void calculate_dN_pTdpTdphidy_and_flows_4all(int perform_sampling=0);
  void calculate_dN_dxtdetady_and_sample_4all();
  void shell(); // it all starts here...
  void combine_samples_to_OSCAR();
  string OSCAR_header_filename, OSCAR_output_filename;

  // Sample files
  string samples_filename, samples_control_filename, samples_format_filename; // where samples, its control informations, and its "format file" are stored; the first two can contain a "%d" string to generate multiple files
  // First sampling method
  void sample_using_dN_dxtdetady_smooth_pT_phi();
  void calculate_dN_dtau_using_dN_dxtdetady(double tau0=0, double dtau=0.5, double tau_max=17);
  string dN_dtau_filename;
  void calculate_dN_dphi_using_dN_pTdpTdphidy();
  string dN_dphi_filename;
  void calculate_dN_deta_using_dN_dxtdetady();
  string dN_deta_filename;
  void calculate_dN_dxt_using_dN_dxtdetady();
  string dN_dxt_filename;

  // Second sampling method
  void calculate_dN_dxtdy_4all_particles(); // the following variables need to be set first in order to call this function
  double calculate_total_FZ_energy_flux();
  void sample_using_dN_dxtdy_4all_particles_conventional(); // to be used after calculate_dN_dxtdy_4all_particles

  // Third sampling method
  void sample_using_dN_pTdpTdphidy();
  Table pT_tab4Sampling, phi_tab4Sampling;
  int pT_tab4Sampling_length, phi_tab4Sampling_length;
  double** trig_phi_tab4Sampling;
};

#endif

/*----------------------------------------------------------------------
Change log: See changelog.txt.
----------------------------------------------------------------------*/
