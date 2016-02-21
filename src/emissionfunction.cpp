// Ver 2.1
// Note that all calculations are done at a given particle rapidity y; and all
// "y_minus_eta_s" appearences in the code are y-y_minus_eta_s.

#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_sf_expint.h>
#include<gsl/gsl_sf_lambert.h>

#include "main.h"
#include "readindata.h"
#include "emissionfunction.h"
#include "RandomVariable1DArray.h"
#include "RandomVariable2DArray.h"
#include "RandomVariableNDArray.h"
#include "NBD.h"
#include "Poisson.h"
#include "ParameterReader.h"
#include "arsenal.h"
#include "Stopwatch.h"

#define AMOUNT_OF_OUTPUT 0 // smaller value means less outputs
#define NUMBER_OF_LINES_TO_WRITE   100000 // string buffer for sample files

using namespace std;

// Class EmissionFunctionArray ------------------------------------------
//***************************************************************************
EmissionFunctionArray::EmissionFunctionArray(
    Table* chosen_particles_in, Table* pt_tab_in, Table* phi_tab_in, 
    Table* y_minus_eta_tab_in, particle_info* particles_in, int Nparticles_in, 
    FO_surf* FOsurf_ptr_in, long FO_length_in, ParameterReader* paraRdr_in)
{
  // get info
  pT_tab = pt_tab_in;
  pT_tab_length = pT_tab->getNumberOfRows();
  phi_tab = phi_tab_in; 
  phi_tab_length = phi_tab->getNumberOfRows();

  y_minus_eta_tab = y_minus_eta_tab_in; 
  y_minus_eta_tab_length = y_minus_eta_tab->getNumberOfRows();

  if (y_minus_eta_tab->get(1,1)>-1e-15)
      positive_y_minus_eta_table_only = true;
  else 
      positive_y_minus_eta_table_only = false;

  particles = particles_in;
  Nparticles = Nparticles_in;

  FOsurf_ptr = FOsurf_ptr_in;
  FO_length = FO_length_in;

  paraRdr = paraRdr_in;

  hydro_mode = paraRdr->getVal("hydro_mode");

  F0_IS_NOT_SMALL = paraRdr->getVal("f0_is_not_small");
  USE_OSCAR_FORMAT = paraRdr->getVal("use_OSCAR_format");
  INCLUDE_DELTAF = paraRdr->getVal("include_deltaf_shear");
  INCLUDE_BULK_DELTAF = paraRdr->getVal("include_deltaf_bulk");
  bulk_deltaf_kind = paraRdr->getVal("bulk_deltaf_kind");

  INCLUDE_DIFFUSION_DELTAF = paraRdr->getVal("include_deltaf_diffusion");

  turn_on_rhob = paraRdr->getVal("turn_on_rhob");

  // allocate internal buffer
  dN_pTdpTdphidy = new Table(pT_tab_length, phi_tab_length);
  dN_pTdpTdphidy_max = new Table(pT_tab_length, phi_tab_length);
  dN_pTdpTdphidy_filename = "results/dN_pTdpTdphidy.dat";
  
  dN_dxtdetady = new double*[y_minus_eta_tab_length];
  for (int k=0; k<y_minus_eta_tab_length; k++) 
      dN_dxtdetady[k] = new double[FO_length];

  dN_dxtdetady_pT_max = new double*[y_minus_eta_tab_length];
  for (int k=0; k<y_minus_eta_tab_length; k++)
      dN_dxtdetady_pT_max[k] = new double[FO_length];

  dN_dxtdetady_filename = "results/dN_dxdetady.dat";

  // deal with chosen_particle_xxx tables
  number_of_chosen_particles = chosen_particles_in->getNumberOfRows();
  // first, for spectra and flow calculations
  chosen_particles_01_table = new int[Nparticles];
  for (int n=0; n<Nparticles; n++) chosen_particles_01_table[n]=0;
  for (int m=0; m<number_of_chosen_particles; m++)
  {
    int monval = chosen_particles_in->get(1,m+1);
    for (int n=0; n<Nparticles; n++)
    {
      if (particles[n].monval==monval)
      {
        chosen_particles_01_table[n]=1;
        break;
      }
    }
  }
  // next, for sampling processes
  chosen_particles_sampling_table = new int[number_of_chosen_particles];
  unidentifiedPid_table = new int [number_of_chosen_particles];
  // first copy the chosen_particles table, but now using indecies 
  // instead of monval
  int current_idx = 0;
  int temp_idx = 0;
  for (int m=0; m<number_of_chosen_particles; m++)
  {
    int monval = chosen_particles_in->get(1,m+1);
    for(int n=0; n<Nparticles; n++)
    {
      if (particles[n].monval==monval)
      {
        chosen_particles_sampling_table[current_idx] = n;
        current_idx ++;
        break;
      }
      else if (n == Nparticles - 1)
      {
        unidentifiedPid_table[temp_idx] = monval;
        temp_idx++;
      }
    }
  }
  // check whether all chosen particles are in the particle list
  if(number_of_chosen_particles != current_idx)
  {
     cout << "Warning: not all chosen particles are in the pdg particle list!" 
          << endl;
     cout << "There are " << number_of_chosen_particles - current_idx 
          << " particles can not be found in the pdg particle list!" << endl;
     cout << "Their monte carlo numbers are:" << endl;
     for(int i = 0; i < number_of_chosen_particles - current_idx ; i++)
        cout << unidentifiedPid_table[i] << endl;

     number_of_chosen_particles = current_idx;
  }
  // next re-order them so that particles with similar mass are adjacent
  grouping_particles = paraRdr->getVal("grouping_particles");
  // sort particles according to their mass; bubble-sorting
  if (grouping_particles) 
  {
    for (int m=0; m<number_of_chosen_particles; m++)
      for (int n=0; n<number_of_chosen_particles-m-1; n++)
        if (particles[chosen_particles_sampling_table[n]].mass 
            > particles[chosen_particles_sampling_table[n+1]].mass)
        {
          // swap them
          int particle_idx = chosen_particles_sampling_table[n+1];
          chosen_particles_sampling_table[n+1] = (
                          chosen_particles_sampling_table[n]);
          chosen_particles_sampling_table[n] = particle_idx;
        }
  }

  // for flow calculation
  flow_differential_filename_old = "results/v2data.dat";
  flow_integrated_filename_old = "results/v2data-inte.dat";
  flow_differential_filename = "results/thermal_%d_vndata.dat";
  flow_integrated_filename = "results/thermal_%d_integrated_vndata.dat";
  last_particle_idx = -1;

  // pre-calculate variables
  //cout << "Caching trig(phi) tables... ";
  trig_phi_table = new double*[phi_tab_length];
  for (int j=0; j<phi_tab_length; j++)
  {
    trig_phi_table[j] = new double[2]; // 2: 0,1-> cos,sin
    double phi = phi_tab->get(1,j+1);
    trig_phi_table[j][0] = cos(phi);
    trig_phi_table[j][1] = sin(phi);
  }
  //cout << "done" << endl;

  //cout << "Caching hyperbolictrig(y-y_minus_eta_s) tables... ";
  hypertrig_y_minus_eta_table = new double*[y_minus_eta_tab_length];
  double y_minus_eta_smallest = 100; y_minus_eta_min_index = 0;
  for (int k=0; k<y_minus_eta_tab_length; k++)
  {
    hypertrig_y_minus_eta_table[k] = new double[2]; // 2: 0,1-> cosh,sinh

    // relative to particle_y
    double y_minus_eta_s = y_minus_eta_tab->get(1,k+1);
    // "y_minus_eta_s" here is actually y-y_minus_eta_s
    hypertrig_y_minus_eta_table[k][0] = cosh(-y_minus_eta_s); 
    hypertrig_y_minus_eta_table[k][1] = sinh(-y_minus_eta_s);

    if (y_minus_eta_s>=0.0 && y_minus_eta_s<y_minus_eta_smallest)
    {
        y_minus_eta_smallest = y_minus_eta_s;
        y_minus_eta_min_index = k;
    }
  }
  //cout << "done" << endl;

  samples_filename = "results/samples_%d.dat";
  samples_control_filename = "results/samples_control_%d.dat";
  samples_format_filename = "results/samples_format.dat";

  dN_dtau_filename = "results/dN_dtau_%d.dat";
  dN_dphi_filename = "results/dN_dphi_%d.dat";
  dN_deta_filename = "results/dN_deta_%d.dat";
  dN_dxt_filename = "results/dN_dxt_%d.dat";
  dN_dx_filename = "results/dN_dx_%d.dat";

  OSCAR_header_filename = "OSCAR_header.txt";
  OSCAR_output_filename = "OSCAR.DAT";

  dN_dxtdy_4all = new double*[FO_length];
  for (long l=0; l<FO_length; l++)
      dN_dxtdy_4all[l] = new double[number_of_chosen_particles];
  //sorted_FZ = new long[FO_length];

  // for interpolation for the third way of sampling
  // generate new set of pT and phi table to be interpolated onto
  pT_tab4Sampling.loadTableFromFile("tables/pT_table_for_sampling.dat");
  pT_tab4Sampling_length = pT_tab4Sampling.getNumberOfRows();
  phi_tab4Sampling.loadTableFromFile("tables/phi_table_for_sampling.dat");
  phi_tab4Sampling_length = phi_tab4Sampling.getNumberOfRows();
  // extend pT_tab and phi_tab in order to extract index info for given 
  // pT or phi
  for (int i=0; i<pT_tab_length; i++)
      pT_tab->set(3, i+1, i+1);
  for (int j=0; j<phi_tab_length; j++)
      phi_tab->set(3, j+1, j+1);
  // get the index info of pT's and phi's specified in pT_tab4Sampling 
  // or phi_tab4Sampling in terms of pT_tab or phi_tab; indecies starts with 1.
  for (int i=0; i<pT_tab4Sampling_length; i++)
  {
      double pT = pT_tab4Sampling.get(1, i+1);
      pT_tab4Sampling.set(3, i+1, pT_tab->interp(1,3,pT,2,true));
  }
  for (int j=0; j<phi_tab4Sampling_length; j++)
  {
      double phi = phi_tab4Sampling.get(1, j+1);
      phi_tab4Sampling.set(3, j+1, phi_tab->interp(1,3,phi,2,true));
  }
  // create trig caches
  trig_phi_tab4Sampling = new double*[phi_tab4Sampling_length];
  for (int j=0; j<phi_tab4Sampling_length; j++)
  {
    trig_phi_tab4Sampling[j] = new double[2]; // 2: 0,1-> cos,sin
    double phi = phi_tab4Sampling.get(1,j+1);
    trig_phi_tab4Sampling[j][0] = cos(phi);
    trig_phi_tab4Sampling[j][1] = sin(phi);
  }

  gsl_rng_env_setup();
  gsl_type_random_number = gsl_rng_default;
  gsl_random_r = gsl_rng_alloc(gsl_type_random_number);

  // arrays for bulk delta f coefficients
  if(INCLUDE_BULK_DELTAF == 1 && bulk_deltaf_kind == 0)
  {
      bulkdf_coeff = (
            new Table ("tables/BulkDf_Coefficients_Hadrons_s95p-v0-PCE.dat"));
  }

  // load table for diffusion delta f coeffient
  if(INCLUDE_DIFFUSION_DELTAF == 1)
  {
      load_deltaf_qmu_coeff_table("tables/Coefficients_RTA_diffusion.dat");
  }

  // create arrays for special functions who are needed to compute 
  // particle yields
  initialize_special_function_arrays();
  
}
//***************************************************************************



//***************************************************************************
EmissionFunctionArray::~EmissionFunctionArray()
{
  // Total number of "new" in constructor should equal the total number of 
  // "delete" in the destructor!
  delete dN_pTdpTdphidy;
  delete dN_pTdpTdphidy_max;

  for (int k=0; k<y_minus_eta_tab_length; k++) 
      delete[] dN_dxtdetady[k];
  delete[] dN_dxtdetady;

  for (int k=0; k<y_minus_eta_tab_length; k++)
      delete[] dN_dxtdetady_pT_max[k];
  delete[] dN_dxtdetady_pT_max;


  delete[] chosen_particles_01_table;
  delete[] chosen_particles_sampling_table;
  delete[] unidentifiedPid_table;

  for (int j=0; j<phi_tab_length; j++)
      delete[] trig_phi_table[j];
  delete[] trig_phi_table;

  for (int k=0; k<y_minus_eta_tab_length; k++)
      delete[] hypertrig_y_minus_eta_table[k];
  delete[] hypertrig_y_minus_eta_table;

  for (long l=0; l<FO_length; l++)
      delete[] dN_dxtdy_4all[l];
  delete[] dN_dxtdy_4all;

  //delete[] sorted_FZ;

  for (int j=0; j<phi_tab4Sampling_length; j++)
      delete[] trig_phi_tab4Sampling[j];
  delete[] trig_phi_tab4Sampling;

  if(INCLUDE_BULK_DELTAF == 1 && bulk_deltaf_kind == 0)
      delete bulkdf_coeff;

  if(INCLUDE_DIFFUSION_DELTAF == 1)
  {
      for(int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
      {
          delete [] deltaf_qmu_coeff_tb[i];
      }
      delete [] deltaf_qmu_coeff_tb;
  }

  gsl_rng_free(gsl_random_r);

  // clean arrays for special functions
  for(int i = 0; i < sf_tb_length; i++)
  {
      delete [] sf_bessel_Kn[i];
      if(INCLUDE_DIFFUSION_DELTAF == 1)
          delete [] sf_expint_En[i];
  }
  delete [] sf_bessel_Kn;
  if(INCLUDE_DIFFUSION_DELTAF == 1)
      delete [] sf_expint_En;
  delete [] lambert_W;
}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_dN_dxtdetady(int particle_idx)
// Calculate dN_dxdetady array.
// Note that this function does not calculate dN_pTdpTdphidy array which is
// needed to calcualte flows.
{
  last_particle_idx = particle_idx;
  Stopwatch sw;
  sw.tic();

  int use_pos_dN_only = paraRdr->getVal("use_pos_dN_only");

  particle_info* particle;
  particle = &particles[particle_idx];

  double mass = particle->mass;
  double sign = particle->sign;
  double degen = particle->gspin;
  double baryon = particle->baryon;

  double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;
  
  double *bulkvisCoefficients;
  if(INCLUDE_BULK_DELTAF == 1)
  {
      if(bulk_deltaf_kind == 0)
          bulkvisCoefficients = new double [3];
      else
          bulkvisCoefficients = new double [2];
  }

  FO_surf* surf = &FOsurf_ptr[0];

  // initialize to 0
  for (int k=0; k<y_minus_eta_tab_length; k++)
  for (long l=0; l<FO_length; l++)
  {
    dN_dxtdetady[k][l] = 0.0;
    dN_dxtdetady_pT_max[k][l] = 0.0;
  }

  // create local cache
  double pT_tab_double[pT_tab_length], pT_tab_weight[pT_tab_length];
  double mT_tab_double[pT_tab_length];
  for (int i=0; i<pT_tab_length; i++)
  {
    pT_tab_double[i] = pT_tab->get(1, i+1);
    pT_tab_weight[i] = pT_tab->get(2, i+1);
    mT_tab_double[i] = sqrt(mass*mass + pT_tab_double[i]*pT_tab_double[i]);
  }
  double phi_tab_weight[phi_tab_length];
  for (int j=0; j<phi_tab_length; j++)
  {
    phi_tab_weight[j] = phi_tab->get(2, j+1);
  }

  //---------------------------
  // THE main summation loop
  //---------------------------
  //cout << "------------------------------------- " << endl;
  //cout << "Performing the main summation loop... " << endl;
  double progress_total = FO_length;
  if (AMOUNT_OF_OUTPUT>0) print_progressbar(-1);
  for (long l=0; l<FO_length; l++)
  {
      surf = &FOsurf_ptr[l];

      double Tdec = surf->Tdec;
      double Pdec = surf->Pdec;
      double Edec = surf->Edec;

      double tau = surf->tau;

      double gammaT = surf->u0;
      double ux = surf->u1;
      double uy = surf->u2;

      double mu =  surf->particle_mu[particle_idx];

      double da0 = surf->da0;
      double da1 = surf->da1;
      double da2 = surf->da2;
      double pi00 = surf->pi00;
      double pi01 = surf->pi01;
      double pi02 = surf->pi02;
      double pi11 = surf->pi11;
      double pi12 = surf->pi12;
      double pi22 = surf->pi22;
      double pi33 = surf->pi33;
      double deltaf_prefactor = 0.0;
      if(INCLUDE_DELTAF == 1)
          deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec));
      
      // bulk delta f
      double bulkPi = 0.0;
      if(INCLUDE_BULK_DELTAF == 1)
      {
          if(bulk_deltaf_kind == 0)
              bulkPi = surf->bulkPi;
          else
              bulkPi = surf->bulkPi/hbarC; // unit in fm^-4 
          getbulkvisCoefficients(Tdec, bulkvisCoefficients);
      }

      // diffusion delta f
      double qmu0 = 0.0;
      double qmu1 = 0.0;
      double qmu2 = 0.0;
      double qmu3 = 0.0;
      double deltaf_qmu_coeff = 1.0;
      double prefactor_qmu = 0.0;
      if(INCLUDE_DIFFUSION_DELTAF == 1)
      {
          qmu0 = surf->qmu0;
          qmu1 = surf->qmu1;
          qmu2 = surf->qmu2;
          qmu3 = surf->qmu3;
          
          double mu_B = surf->muB;
          deltaf_qmu_coeff = get_deltaf_qmu_coeff(Tdec, mu_B);
          
          double rho_B = surf->Bn;
          prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
      }

      for (int k=0; k<y_minus_eta_tab_length; k++)
      {
          // use local variables to speed up
          double dN_dxtdetady_tmp = 0.0;
          double dN_dxtdetady_pT_max_tmp = 0.0;

          for (int i=0; i<pT_tab_length; i++)
          {
              double pT = pT_tab_double[i];
              double pT_weight = pT_tab_weight[i];
              double mT = mT_tab_double[i];

              double pt = mT*hypertrig_y_minus_eta_table[k][0];
              double pz = mT*hypertrig_y_minus_eta_table[k][1];

              for (int j=0; j<phi_tab_length; j++)
              {
                  double phi_weight = phi_tab_weight[j];

                  double px = pT*trig_phi_table[j][0];
                  double py = pT*trig_phi_table[j][1];

                  double pT_phi_inte_weight = pT*pT_weight*phi_weight;


                  double pdotu = pt*gammaT - px*ux - py*uy;
                  double expon = (pdotu - mu) / Tdec;
                  // thermal equilibrium distributions
                  double f0 = 1./(exp(expon)+sign);

                  double pdsigma = pt*da0 + px*da1 + py*da2;

                  //viscous corrections
                  double delta_f_shear = 0.0;
                  if(INCLUDE_DELTAF)
                  {
                      double Wfactor = (
                          pt*pt*pi00 - 2.0*pt*px*pi01 - 2.0*pt*py*pi02 
                          + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 
                          + pz*pz*pi33);
                      delta_f_shear = (
                          (1 - F0_IS_NOT_SMALL*sign*f0)*Wfactor
                          *deltaf_prefactor);
                  }

                  double delta_f_bulk = 0.0;
                  if (INCLUDE_BULK_DELTAF== 1)
                  {
                      if(bulk_deltaf_kind == 0)
                      {
                          delta_f_bulk = (
                              -(1. - F0_IS_NOT_SMALL*sign*f0)*bulkPi
                              *(bulkvisCoefficients[0]*mass*mass 
                                + bulkvisCoefficients[1]*pdotu 
                                + bulkvisCoefficients[2]*pdotu*pdotu));
                      }
                      else if (bulk_deltaf_kind == 1)
                      {

                          double E_over_T = pdotu/Tdec;
                          double mass_over_T = mass/Tdec;
                          delta_f_bulk = (
                              -1.0*(1.-sign*f0)/E_over_T*bulkvisCoefficients[0]
                              *(mass_over_T*mass_over_T/3. 
                                - bulkvisCoefficients[1]*E_over_T*E_over_T)
                              *bulkPi);
                      }
                      else if (bulk_deltaf_kind == 2)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = (
                              -1.*(1.-sign*f0)*(-bulkvisCoefficients[0] 
                                  + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                      }
                      else if (bulk_deltaf_kind == 3)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = (
                              -1.0*(1.-sign*f0)/sqrt(E_over_T)
                              *(-bulkvisCoefficients[0] 
                                + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                      }
                      else if (bulk_deltaf_kind == 4)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = (
                              -1.0*(1.-sign*f0)*(bulkvisCoefficients[0] 
                                  - bulkvisCoefficients[1]/E_over_T)*bulkPi);
                      }
                  }
                  
                  // delta f for diffusion
                  double delta_f_qmu = 0.0;
                  if(INCLUDE_DIFFUSION_DELTAF == 1)
                  {
                      double qmufactor = pt*qmu0 - px*qmu1 - py*qmu2 - pz*qmu3;
                      delta_f_qmu = ((1. - sign*f0)
                                     *(prefactor_qmu - baryon/pdotu)
                                     *qmufactor/deltaf_qmu_coeff);
                  }

                  double result;
                  // restrict the size of delta f to be smaller than f_0
                  double ratio_max = 1.0;
                  double deltaf_size = fabs(delta_f_shear + delta_f_bulk 
                                            + delta_f_qmu);
                  double resize_factor = min(1., 
                                             ratio_max/(deltaf_size + 1e-10));
                  result = (prefactor*degen*f0*pdsigma*tau
                         *(1. + (delta_f_shear + delta_f_bulk + delta_f_qmu)
                                *resize_factor));

                  if (use_pos_dN_only && result < 0.)
                      continue;

                  dN_dxtdetady_tmp += result*pT_phi_inte_weight;

                  if (dN_dxtdetady_pT_max_tmp<result*pT)
                      dN_dxtdetady_pT_max_tmp=result*pT;
              } // j
          } // i
          dN_dxtdetady[k][l] = dN_dxtdetady_tmp;
          dN_dxtdetady_pT_max[k][l] = dN_dxtdetady_pT_max_tmp;

      }
      if (AMOUNT_OF_OUTPUT>0)
          print_progressbar((double)(l)/progress_total);
  }
  if (AMOUNT_OF_OUTPUT>0)
      print_progressbar(1);
  //cout << endl << "------------------------------------- " << endl;

  if(INCLUDE_BULK_DELTAF == 1)
      delete [] bulkvisCoefficients;

  sw.toc();
  cout << endl 
       << " -- Calculate_dN_dxtdetady finished in " 
       << sw.takeTime() << " seconds." << endl;
}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_dN_pTdpTdphidy(int particle_idx)
// Calculate only dN_pTdpTdphidy array; tuned for best efficiency. This
// function is ~3 times faster than the generaic version calculate_dNArrays.
// Note that the execution fo this function does NOT provide the dN_dxtdetady
// array which is needed for sampling.
{
  last_particle_idx = particle_idx;
  Stopwatch sw;
  sw.tic();

  int use_pos_dN_only = paraRdr->getVal("use_pos_dN_only");

  particle_info* particle;
  particle = &particles[particle_idx];

  double mass = particle->mass;
  double sign = particle->sign;
  double degen = particle->gspin;
  double baryon = particle->baryon;

  double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

  FO_surf* surf = &FOsurf_ptr[0];
  
  double *bulkvisCoefficients;
  if(INCLUDE_BULK_DELTAF == 1)
  {
      if(bulk_deltaf_kind == 0)
          bulkvisCoefficients = new double [3];
      else
          bulkvisCoefficients = new double [2];
  }

  // for intermedia results
  //cout << "initializing intermedia variables... ";
  double dN_pTdpTdphidy_tab[pT_tab_length][phi_tab_length];
  for (int i=0; i<pT_tab_length; i++)
      for (int j=0; j<phi_tab_length; j++)
          dN_pTdpTdphidy_tab[i][j] = 0.0;
  double dN_pTdpTdphidy_max_tab[pT_tab_length][phi_tab_length];
  for (int i=0; i<pT_tab_length; i++)
      for (int j=0; j<phi_tab_length; j++)
          dN_pTdpTdphidy_max_tab[i][j] = 0.0;


  // create local cache
  double delta_y_minus_eta_tab[y_minus_eta_tab_length];
  for (int k=0; k<y_minus_eta_tab_length; k++)
      delta_y_minus_eta_tab[k] = y_minus_eta_tab->get(2,k+1);

  //---------------------------
  // THE main summation loop
  //---------------------------
  //cout << "------------------------------------- " << endl;
  //cout << "Performing the main summation loop... " << endl;
  double progress_total = pT_tab_length*phi_tab_length;
  if (AMOUNT_OF_OUTPUT>0) print_progressbar(-1);
  for (int i=0; i<pT_tab_length; i++)
  {
      double pT = pT_tab->get(1,i+1);
      double mT = sqrt(mass*mass + pT*pT);

      for (int j=0; j<phi_tab_length; j++)
      {
          double px = pT*trig_phi_table[j][0];
          double py = pT*trig_phi_table[j][1];

          double dN_pTdpTdphidy_tmp = 0.0;
          double dN_pTdpTdphidy_max_tmp = 0.0;

          for (long l=0; l<FO_length; l++)
          {
              surf = &FOsurf_ptr[l];

              double Tdec = surf->Tdec;
              double Pdec = surf->Pdec;
              double Edec = surf->Edec;

              double tau = surf->tau;
              double gammaT = surf->u0;
              double ux = surf->u1;
              double uy = surf->u2;

              double mu = surf->particle_mu[particle_idx];

              double da0 = surf->da0;
              double da1 = surf->da1;
              double da2 = surf->da2;
              double pi00 = surf->pi00;
              double pi01 = surf->pi01;
              double pi02 = surf->pi02;
              double pi11 = surf->pi11;
              double pi12 = surf->pi12;
              double pi22 = surf->pi22;
              double pi33 = surf->pi33;
              double bulkPi = 0.0;
              double deltaf_prefactor = 0.0;
              if(INCLUDE_DELTAF)
                  deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec));
              if(INCLUDE_BULK_DELTAF == 1)
              {
                  if(bulk_deltaf_kind == 0)
                      bulkPi = surf->bulkPi;
                  else 
                      bulkPi = surf->bulkPi/hbarC;   // unit in fm^-4 
                  getbulkvisCoefficients(Tdec, bulkvisCoefficients);
              }

              // diffusion delta f
              double qmu0 = 0.0;
              double qmu1 = 0.0;
              double qmu2 = 0.0;
              double qmu3 = 0.0;
              double deltaf_qmu_coeff = 1.0;
              double prefactor_qmu = 0.0;
              if(INCLUDE_DIFFUSION_DELTAF == 1)
              {
                  qmu0 = surf->qmu0;
                  qmu1 = surf->qmu1;
                  qmu2 = surf->qmu2;
                  qmu3 = surf->qmu3;
                  
                  double mu_B = surf->muB;
                  deltaf_qmu_coeff = get_deltaf_qmu_coeff(Tdec, mu_B);
                  
                  double rho_B = surf->Bn;
                  prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
              }

              for (int k=0; k<y_minus_eta_tab_length; k++)
              {
                  double delta_eta = delta_y_minus_eta_tab[k];

                  double pt = mT*hypertrig_y_minus_eta_table[k][0];
                  double pz = mT*hypertrig_y_minus_eta_table[k][1];

                  double pdotu = pt*gammaT - px*ux - py*uy;
                  double expon = (pdotu - mu) / Tdec;
                  double f0 = 1./(exp(expon)+sign);

                  double pdsigma = pt*da0 + px*da1 + py*da2;

                  //viscous corrections
                  double delta_f_shear = 0.0;
                  if(INCLUDE_DELTAF)
                  {
                      double Wfactor = (
                          pt*pt*pi00 - 2.0*pt*px*pi01 - 2.0*pt*py*pi02 
                          + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 
                          + pz*pz*pi33);
                      delta_f_shear = (
                          (1 - F0_IS_NOT_SMALL*sign*f0)
                          *Wfactor*deltaf_prefactor);
                  }

                  double delta_f_bulk = 0.0;
                  if (INCLUDE_BULK_DELTAF == 1)
                  {
                      if(bulk_deltaf_kind == 0)
                          delta_f_bulk = (
                              -(1. - F0_IS_NOT_SMALL*sign*f0)*bulkPi
                              *(bulkvisCoefficients[0]*mass*mass 
                                + bulkvisCoefficients[1]*pdotu 
                                + bulkvisCoefficients[2]*pdotu*pdotu));
                      else if (bulk_deltaf_kind == 1)
                      {

                          double E_over_T = pdotu/Tdec;
                          double mass_over_T = mass/Tdec;
                          delta_f_bulk = (
                              -1.0*(1.-sign*f0)/E_over_T*bulkvisCoefficients[0]
                              *(mass_over_T*mass_over_T/3. 
                                - bulkvisCoefficients[1]*E_over_T*E_over_T)
                              *bulkPi);
                      }
                      else if (bulk_deltaf_kind == 2)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = (
                              -1.*(1.-sign*f0)
                              *(-bulkvisCoefficients[0] 
                                + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                      }
                      else if (bulk_deltaf_kind == 3)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = (
                              -1.0*(1.-sign*f0)/sqrt(E_over_T)
                              *(-bulkvisCoefficients[0] 
                                + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                      }
                      else if (bulk_deltaf_kind == 4)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = (
                              -1.0*(1.-sign*f0)
                              *(bulkvisCoefficients[0] 
                                - bulkvisCoefficients[1]/E_over_T)*bulkPi);
                      }
                  }

                  double delta_f_qmu = 0.0;
                  if(INCLUDE_DIFFUSION_DELTAF == 1)
                  {
                      double qmufactor = pt*qmu0 - px*qmu1 - py*qmu2 - pz*qmu3;
                      delta_f_qmu = ((1. - sign*f0)
                                     *(prefactor_qmu - baryon/pdotu)
                                     *qmufactor/deltaf_qmu_coeff);
                  }
                  
                  double result;
                  // restrict the size of delta f to be smaller than f_0
                  double ratio_max = 1.0;
                  double deltaf_size = fabs(delta_f_shear + delta_f_bulk
                                            + delta_f_qmu);
                  double resize_factor = min(1., 
                                             ratio_max/(deltaf_size + 1e-10));
                  result = (prefactor*degen*f0*pdsigma*tau
                         *(1. + (delta_f_shear + delta_f_bulk + delta_f_qmu)
                                 *resize_factor));

                  if (use_pos_dN_only && result < 0.)
                      continue;

                  dN_pTdpTdphidy_tmp += result*delta_eta;

                  if (dN_pTdpTdphidy_max_tmp<result)
                      dN_pTdpTdphidy_max_tmp=result;

              } // k
          } // l
          dN_pTdpTdphidy_tab[i][j] = dN_pTdpTdphidy_tmp;
          dN_pTdpTdphidy_max_tab[i][j] = dN_pTdpTdphidy_max_tmp;
          if (AMOUNT_OF_OUTPUT>0)
              print_progressbar((i*phi_tab_length+j)/progress_total);
      }
  }
  if (AMOUNT_OF_OUTPUT>0) print_progressbar(1);
  //cout << endl << "------------------------------------- " << endl;

  //cout << "Copy out results... ";
  for (int i=0; i<pT_tab_length; i++)
      for (int j=0; j<phi_tab_length; j++)
      {
          dN_pTdpTdphidy->set(i+1,j+1,dN_pTdpTdphidy_tab[i][j]);
          dN_pTdpTdphidy_max->set(i+1,j+1,dN_pTdpTdphidy_max_tab[i][j]);
      }
  //cout << "done." << endl;

  if(INCLUDE_BULK_DELTAF == 1)
      delete [] bulkvisCoefficients;

  sw.toc();
  cout << endl 
       << " -- Calculate_dN_pTdpTdphidy finished in " 
       << sw.takeTime() << " seconds." << endl;
}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::write_dN_pTdpTdphidy_toFile()
// Append the dN_pTdpTdphidy results to file.
{
  ofstream of1(dN_pTdpTdphidy_filename.c_str(), ios_base::app);
  dN_pTdpTdphidy->printTable(of1);
  of1.close();
}
//***************************************************************************



//***************************************************************************
void EmissionFunctionArray::write_dN_dxtdetady_toFile()
// Append the dN_pTdpTdphidy results to file.
{
  ofstream of1(dN_dxtdetady_filename.c_str(), ios_base::app);
  for (int k=0; k<y_minus_eta_tab_length; k++)
  {
    for (long l=0; l<FO_length; l++)
    {
         of1 << scientific << setprecision(12) << dN_dxtdetady[k][l] << "   ";
    }
    of1 << endl;
  }
  of1.close();
}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_flows(
    int to_order, string flow_differential_filename_in, 
    string flow_integrated_filename_in)
// Calculate flow from order from_order to to_order and store them to files.
{
  /*
  cout << endl
       <<"*************************************************"
       << endl
       << "Function calculate_flows started... " << endl;*/
  Stopwatch sw;
  sw.tic();

  int from_order = 1;

  int number_of_flows = to_order-from_order+1;

  // line format: pT, mT, dN/(pT dpT), 
  // flow_1_real, flow_1_imag, flow_1_norm, ...
  Table vn_diff(3+number_of_flows*3, pT_tab_length); 
  // line format: order# (starting from 0), numerator_real, numerator_imag, 
  // flow_real, flow_imag, flow_norm
  Table vn_inte(6, to_order+1); 

  double mass = particles[last_particle_idx].mass;

  //---------------------
  // differential flow
  //---------------------
  //cout << "Calculating differential flows... ";

  double normalization[pT_tab_length]; // normalization factor
  for (int i=0; i<pT_tab_length; i++) normalization[i] = 0.0;

  // diff_flow numerators; 2: 0,1->real,imag
  double vn[pT_tab_length][number_of_flows][2]; 
  for (int i=0; i<pT_tab_length; i++)
  for (int t=0; t<number_of_flows; t++)
    {vn[i][t][0]=0; vn[i][t][1]=0;}

  for (int i=0; i<pT_tab_length; i++)
  //for (int i=0; i<1; i++) // for debugging
  {
    double pT = pT_tab->get(1,i+1);
    double mT = sqrt(mass*mass + pT*pT);

    // phi integration
    for(int j=0; j<phi_tab_length; j++)
    //for(int j=0; j<1; j++) // for debugging
    {
      double phi = phi_tab->get(1,j+1), phi_weight = phi_tab->get(2,j+1);
      double dN = dN_pTdpTdphidy->get(i+1,j+1);

      normalization[i] += dN*phi_weight;
      for (int order=from_order; order<=to_order; order++)
      {
        vn[i][order-from_order][0] += dN*phi_weight*cos(order*phi);
        vn[i][order-from_order][1] += dN*phi_weight*sin(order*phi);
      }
    }

    normalization[i] = normalization[i] + 1e-30;
    // store values
    vn_diff.set(1, i+1, pT);
    vn_diff.set(2, i+1, mT-mass);
    vn_diff.set(3, i+1, normalization[i]/(2.0*M_PI)); 
    // 2*pi: azimuthal angle averaged
    for (int t=0; t<number_of_flows; t++)
    {
      vn_diff.set(4+t*3, i+1, vn[i][t][0]/normalization[i]);
      vn_diff.set(5+t*3, i+1, vn[i][t][1]/normalization[i]);
      vn_diff.set(6+t*3, i+1, 
                  sqrt(vn[i][t][0]*vn[i][t][0] + vn[i][t][1]*vn[i][t][1])
                  /normalization[i]);
    }

  }
  //cout << "done." << endl;


  //---------------------
  // integrated flow
  //---------------------
  //cout << "Calculating integrated flows... ";

  double normalizationi = 0;

  // integrated_flow numerators; 2: 0,1->real,imag
  double vni[number_of_flows][2]; 
  for (int t=0; t<number_of_flows; t++)
  {
      vni[t][0]=0; 
      vni[t][1]=0;
  }

  for (int i=0; i<pT_tab_length; i++)
  //for (int i=0; i<1; i++) // for debugging
  {
    double pT = pT_tab->get(1,i+1), pT_weight = pT_tab->get(2,i+1);

    normalizationi += normalization[i]*pT*pT_weight;

    for (int order=from_order; order<=to_order; order++)
    {
      vni[order-from_order][0] += vn[i][order-from_order][0]*pT*pT_weight;
      vni[order-from_order][1] += vn[i][order-from_order][1]*pT*pT_weight;
    }

  }

  // store values
  // To mimic:
  // WRITE(30,941) " ", N, " ",X(N)," ",Y(N),
  //   &  " ",X(N)/X(0)," ",Y(N)/X(0)," ",
  //   &  sqrt(X(N)*X(N)+Y(N)*Y(N))/X(0)
  vn_inte.set(1, 1, 0);
  vn_inte.set(2, 1, normalizationi);
  vn_inte.set(3, 1, 0);
  vn_inte.set(4, 1, 1);
  vn_inte.set(5, 1, 0);
  vn_inte.set(6, 1, 1);

  for (int t=0; t<number_of_flows; t++)
  {
    vn_inte.set(1, t+2, from_order+t);
    vn_inte.set(2, t+2, vni[from_order+t-1][0]);
    vn_inte.set(3, t+2, vni[from_order+t-1][1]);
    vn_inte.set(4, t+2, vni[from_order+t-1][0]/normalizationi);
    vn_inte.set(5, t+2, vni[from_order+t-1][1]/normalizationi);
    vn_inte.set(6, t+2, 
                sqrt(vni[from_order+t-1][0]*vni[from_order+t-1][0] 
                     + vni[from_order+t-1][1]*vni[from_order+t-1][1])
                /normalizationi);
  }
  //cout << "done." << endl;

  // save to files
  //cout << "Writing to files... ";
  ofstream of1(flow_differential_filename_in.c_str(), ios_base::app);
  vn_diff.printTable(of1);
  of1.close();

  ofstream of2(flow_integrated_filename_in.c_str(), ios_base::app);
  vn_inte.printTable(of2);
  of2.close();
  //cout << "done." << endl;

  sw.toc();
  //cout << "calculate_flows finishes " << sw.takeTime() << " seconds." << endl;

}
//***************************************************************************



//***************************************************************************
void EmissionFunctionArray::calculate_dN_pTdpTdphidy_and_flows_4all_old_output(
                                                          int perform_sampling)
// Calculate dN_pTdpTdphidy and flows for all particles given in chosen_particle 
// array.
// This version output dN / (ptdpt dphi dy) matrices and flows in a single file
// as the Azspectra7 did.
// If perform_sampling is 1 then sample_using_dN_pTdpTdphidy will be called.
{

    cout << endl
        << "*****************************************************************"
        << endl
        << "Function calculate_dN_pTdpTdphidy_and_flows_4all(old) started... " 
        << endl;
    Stopwatch sw;
    sw.tic();

    remove(dN_pTdpTdphidy_filename.c_str());
    remove(flow_differential_filename_old.c_str());
    remove(flow_integrated_filename_old.c_str());


    // read in parameters
    int calculate_dN_dphi = paraRdr->getVal("calculate_dN_dphi");
    int to_order = paraRdr->getVal("calculate_vn_to_order");

    // loop over particles
    particle_info* particle = NULL;
    for (int n=0; n<Nparticles; n++)
    {
        particle = &particles[n];
        cout << "Index: " << n << ", Name: " << particle->name 
             << ", Monte-carlo index: " << particle->monval;

        // first, calculate dN_pTdpTdphidy arrays:
        if (chosen_particles_01_table[n]==0)
        {
            cout << " ...skipped." << endl;
            dN_pTdpTdphidy->setAll(0.0);
            last_particle_idx = n; // fake a "calculation"
        }
        else
        {
            // calculate dN_*** arrays
            cout << endl;
            if (n>0 && particles_are_the_same(n,n-1))
            {
                // no need to calculate dN_pTdpTdphidy again
                cout << " -- Using previously calculated dN_pTdpTdphidy... " 
                     << endl;
                last_particle_idx = n; // fake a "calculation"
            }
            else
            {
                cout << " -- Calculating dN_pTdpTdphidy... " << endl;
                calculate_dN_pTdpTdphidy(n);
            }

            // perform sampling
            if (perform_sampling) sample_using_dN_pTdpTdphidy();

            // calculate dN/dphi
            if (calculate_dN_dphi) calculate_dN_dphi_using_dN_pTdpTdphidy();

        }


        // first write dN/(pt dpt dphi dy) array
        write_dN_pTdpTdphidy_toFile();

        // next flows:

        ofstream of1(flow_differential_filename_old.c_str(), ios_base::app);
        of1 << "# Output for particle: " << particle->name << endl;
        of1 << "#                 " << particle->monval << endl;
        of1.close();

        ofstream of2(flow_integrated_filename_old.c_str(), ios_base::app);
        of2 << "# For: " << particle->name << endl;
        of2.close();
        calculate_flows(to_order, flow_differential_filename_old, 
                        flow_integrated_filename_old);

    } // n: particle loop

    sw.toc();
    cout << " -- Calculate_dN_pTdpTdphidy_and_flows_4all finishes " 
         << sw.takeTime() << " seconds." << endl;
}
//***************************************************************************






//***************************************************************************
void EmissionFunctionArray::calculate_dN_pTdpTdphidy_and_flows_4all(
                                                          int perform_sampling)
// Calculate dN_pTdpTdphidy and flows 
// for all particles given in chosen_particle array.
// This version output dN / (ptdpt dphi dy) matrices in one file, 
// but flows in separated ones.
{

    cout << endl
        << "****************************************************************"
        << endl
        << "Function calculate_dN_pTdpTdphidy_and_flows_4all started... " 
        << endl;
    Stopwatch sw;
    sw.tic();

    // read in parameters
    int calculate_dN_dphi = paraRdr->getVal("calculate_dN_dphi");
    int to_order = paraRdr->getVal("calculate_vn_to_order");


    // prepare a huge array to store calculated dN_pTdpTdphidy
    Table* dNs[Nparticles];
    for (int n=0; n<Nparticles; n++) dNs[n]=NULL;

    // loop over chosen particles
    particle_info* particle = NULL;
    for (int m=0; m<number_of_chosen_particles; m++)
    {
        int particle_idx = chosen_particles_sampling_table[m];
        particle = &particles[particle_idx];
        int monval = particle->monval;
        cout << "Index: " << m << ", Name: " << particle->name 
             << ", Monte-carlo index: " << monval << endl;
        // Calculate dN / (ptdpt dphi dy)
        if (m > 0. 
            && particles_are_the_same(particle_idx, 
                                      chosen_particles_sampling_table[m-1]))
        {
           cout << " -- Using dN_pTdpTdphidy from previous calculation... " 
                << endl;
           last_particle_idx = particle_idx; // fake a calculation
        }
        else
        {
            cout << " -- Calculating dN_pTdpTdphidy... " << endl;
            calculate_dN_pTdpTdphidy(particle_idx);
        }

        // perform sampling
        if (perform_sampling)
            sample_using_dN_pTdpTdphidy();

        // Store calculated table
        dNs[particle_idx] = new Table(*dN_pTdpTdphidy);

        // Calcualte dN / dphi
        if (calculate_dN_dphi)
            calculate_dN_dphi_using_dN_pTdpTdphidy();

        char buffer_diff[500], buffer_inte[500];
        sprintf(buffer_diff, flow_differential_filename.c_str(), monval);
        remove(buffer_diff);
        sprintf(buffer_inte, flow_integrated_filename.c_str(), monval);
        remove(buffer_inte);
        calculate_flows(to_order, buffer_diff, buffer_inte);
    }

    // write out dN / (ptdpt dphi dy) matrices
    remove(dN_pTdpTdphidy_filename.c_str());
    ofstream of(dN_pTdpTdphidy_filename.c_str(), ios_base::app);
    Table zero(dN_pTdpTdphidy->getNumberOfCols(), 
               dN_pTdpTdphidy->getNumberOfRows(), 0);
    for (int n=0; n<Nparticles; n++)
    {
        if (dNs[n]==NULL)
        {
            zero.printTable(of);
        }
        else
        {
            dNs[n]->printTable(of);
            delete dNs[n];
        }
    }
    of.close();

    sw.toc();
    cout << " -- Calculate_dN_pTdpTdphidy_and_flows_4all finishes " 
         << sw.takeTime() << " seconds." << endl;

}
//***************************************************************************




//***************************************************************************
void EmissionFunctionArray::sample_using_dN_dxtdetady_smooth_pT_phi()
// This version give smooth distribution in pT and phi.
// Sample according to the dN_dxtdetady and dN_dxtdetady_pT_max array. These
// arrays are asssumed to be already cacluated.
// The particle index is assumed to be already stored in the
// "last_particle_idx" internal variable.
// Emission function below "zero" will be treated as "zero".
// Each line in the sample file contains the following data:
// FZ_cell_idx tau, xpt, ypt, y_minus_eta_s, pT, phi, da0, da1, da2, vx, vy
// Each line in the sample control file records the number of particles
// sampled in one event.
// The format file is a file that can be read by ParameterReader class
// which records the the purpose of each column in the sample file. For
// example, if one line of it reads as "foo=1" means that the variable
// "foo" is the recorded in the 1st column.
// Parameters:
// -- number_of_repeated_sampling: number of successive sampling; 
//    in each sampling the number of particles is determined
//    by dN_dy and the model parameter.
// -- pT_to: sample up to what pT
// -- model parameter:
//    1): The fractional part of dN_dy is used as a probability to determine
//       whether we have 0 or 1 more particle.
{
    Stopwatch sw;
    sw.tic();

    double pT_to = paraRdr->getVal("sample_pT_up_to");
    if (pT_to < 0.)
        pT_to = pT_tab->getLast(1); // use table to determine pT range
    
    // If dN/(dx_t deta dy) is evaluated to be smaller than this value, 
    // then it is replaced by this value.
    double zero = paraRdr->getVal("minimum_emission_function_val"); 
    long number_of_repeated_sampling = (
                    paraRdr->getVal("number_of_repeated_sampling"));

    particle_info* particle = &particles[last_particle_idx];

    double mass = particle->mass;
    double sign = particle->sign;
    double degen = particle->gspin;
    double baryon = particle->baryon;

    double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

    FO_surf* surf = &FOsurf_ptr[0];
  
    double *bulkvisCoefficients;
    if(INCLUDE_BULK_DELTAF == 1)
    {
        if(bulk_deltaf_kind == 0)
            bulkvisCoefficients = new double [3];
        else
            bulkvisCoefficients = new double [2];
    }

    // create local cache
    double delta_y_minus_eta_tab[y_minus_eta_tab_length];
    for (int k=0; k<y_minus_eta_tab_length; k++)
        delta_y_minus_eta_tab[k] = y_minus_eta_tab->get(2,k+1);

    // prepare the inverse CDF
    // put back eta_weight in sampling function
    double **dN_dxtdetady_with_weight;
    dN_dxtdetady_with_weight = new double*[y_minus_eta_tab_length];
    for (int k=0; k<y_minus_eta_tab_length; k++)
        dN_dxtdetady_with_weight[k] = new double[FO_length];

    for (int k=0; k<y_minus_eta_tab_length; k++)
        for (long l=0; l<FO_length; l++)
            dN_dxtdetady_with_weight[k][l] = (
                            dN_dxtdetady[k][l]*delta_y_minus_eta_tab[k]);

    RandomVariable2DArray rand2D(dN_dxtdetady_with_weight, FO_length, 
                                 y_minus_eta_tab_length, zero);

    // first calcualte total number of particles
    double dN_dy = 0.0;
    for (long l=0; l<FO_length; l++)
        for (int k=0; k<y_minus_eta_tab_length; k++)
            dN_dy += dN_dxtdetady_with_weight[k][l];

    // recycle
    for (int k=0; k<y_minus_eta_tab_length; k++)
        delete[] dN_dxtdetady_with_weight[k];
    delete[] dN_dxtdetady_with_weight;

    // prepare for outputs
    // the control file records how many particles are there in each sampling
    char samples_control_filename_buffer[300];
    sprintf(samples_control_filename_buffer, samples_control_filename.c_str(), 
            particle->monval);
    remove(samples_control_filename_buffer);
    ofstream of_control(samples_control_filename_buffer);

    // the sample file contains the actual samples
    char samples_filename_buffer[300];
    sprintf(samples_filename_buffer, samples_filename.c_str(), 
            particle->monval);
    remove(samples_filename_buffer);
    ofstream of_sample(samples_filename_buffer);

    // buffers are used to speed up the output process
    char line_buffer[500]; // only used in text mode
    stringstream sample_str_buffer; // to speed up outputing process
    stringstream control_str_buffer;

    int sampling_model = paraRdr->getVal("dN_dy_sampling_model");
    double sampling_para1 = paraRdr->getVal("dN_dy_sampling_para1");
    double sampling_para2 = paraRdr->getVal("dN_dy_sampling_para2");
    double sampling_para3 = paraRdr->getVal("dN_dy_sampling_para3");
    double sampling_para4 = paraRdr->getVal("dN_dy_sampling_para4");
    double sampling_para5 = paraRdr->getVal("dN_dy_sampling_para5");

    // get y range for sampling
    double y_LB = paraRdr->getVal("y_LB");
    double y_RB = paraRdr->getVal("y_RB");

    cout << " dN_dy=" << dN_dy << ", ";
    double dN = (y_RB-y_LB)*dN_dy;
    cout << "dN=" << dN << "..." << endl;
    if (AMOUNT_OF_OUTPUT>0) print_progressbar(-1);
    long sample_writing_signal = 0;
    long control_writing_signal = 0;
    long maximum_impatience = 5000; // used in pt-phi sampling
    for (long sampling_idx=1; sampling_idx<=number_of_repeated_sampling; 
         sampling_idx++)
    {
        long number_to_sample = determine_number_to_sample(
            dN, sampling_model, sampling_para1, sampling_para2, 
            sampling_para3, sampling_para4, sampling_para5);

        // write to control file
        sprintf(line_buffer, "%lu\n", number_to_sample);
        control_str_buffer << line_buffer;
        control_writing_signal++;
        if (control_writing_signal==NUMBER_OF_LINES_TO_WRITE)
        {
            of_control << control_str_buffer.str();
            control_str_buffer.str("");
            control_writing_signal=0;
        }

        long y_minus_eta_s_idx, FO_idx;
        long sample_idx = 1;
        while (sample_idx <= number_to_sample)
        {
            // first, sample eta and freeze-out cell index
            rand2D.sampleAccToInvCDF(&y_minus_eta_s_idx, &FO_idx);
            surf = &FOsurf_ptr[FO_idx];

            // Table starts with 1
            double y_minus_eta_s = y_minus_eta_tab->get(1,y_minus_eta_s_idx+1); 

            double Tdec = surf->Tdec;
            double inv_Tdec = 1.0/Tdec;
            double Pdec = surf->Pdec;
            double Edec = surf->Edec;

            double tau = surf->tau;

            double gammaT = surf->u0;
            double ux = surf->u1;
            double uy = surf->u2;

            double mu = surf->particle_mu[last_particle_idx];

            double da0 = surf->da0;
            double da1 = surf->da1;
            double da2 = surf->da2;
            double pi00 = surf->pi00;
            double pi01 = surf->pi01;
            double pi02 = surf->pi02;
            double pi11 = surf->pi11;
            double pi12 = surf->pi12;
            double pi22 = surf->pi22;
            double pi33 = surf->pi33;
            double bulkPi = 0.0;
            double deltaf_prefactor = 0.0;
            if(INCLUDE_DELTAF)
                deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec));

            if(INCLUDE_BULK_DELTAF == 1)
            {
                if(bulk_deltaf_kind == 0)
                    bulkPi = surf->bulkPi;
                else
                    bulkPi = surf->bulkPi/hbarC;   // unit in fm^-4 
                getbulkvisCoefficients(Tdec, bulkvisCoefficients);
            }

            // diffusion delta f
            double qmu0 = 0.0;
            double qmu1 = 0.0;
            double qmu2 = 0.0;
            double qmu3 = 0.0;
            double deltaf_qmu_coeff = 1.0;
            double prefactor_qmu = 0.0;
            if(INCLUDE_DIFFUSION_DELTAF == 1)
            {
                qmu0 = surf->qmu0;
                qmu1 = surf->qmu1;
                qmu2 = surf->qmu2;
                qmu3 = surf->qmu3;
                
                double mu_B = surf->muB;
                deltaf_qmu_coeff = get_deltaf_qmu_coeff(Tdec, mu_B);
                
                double rho_B = surf->Bn;
                prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
            }

            // next sample pt and phi
            double pT, mT, phi, px, py;// will-be sampled values
            long tries = 1;
            while (tries<maximum_impatience)
            {

                // refer to calculate_dNArrays function to see how the rate 
                // is calculated
                // Basically it is "just Cooper-Frye"
                pT = drand(0, pT_to); // sample according to pT dpT
                mT = sqrt(mass*mass + pT*pT);

                double pt = (
                    mT*hypertrig_y_minus_eta_table[y_minus_eta_s_idx][0]);
                double pz = (
                    mT*hypertrig_y_minus_eta_table[y_minus_eta_s_idx][1]);

                phi = drand(0, 2*M_PI);
                px = pT*cos(phi);
                py = pT*sin(phi);

                double pdotu = pt*gammaT - px*ux - py*uy;
                double expon = (pdotu - mu) * inv_Tdec;
                double f0 = 1./(exp(expon)+sign);

                double pdsigma = pt*da0 + px*da1 + py*da2;
                  
                double delta_f_shear = 0.0;
                if(INCLUDE_DELTAF)
                {
                    double Wfactor = (
                        pt*pt*pi00 - 2.0*pt*px*pi01 - 2.0*pt*py*pi02 
                        + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 
                        + pz*pz*pi33);
                    delta_f_shear = (
                        (1. - F0_IS_NOT_SMALL*sign*f0)
                        *Wfactor*deltaf_prefactor);
                }
                double delta_f_bulk = 0.0;
                if (INCLUDE_BULK_DELTAF == 1)
                {
                    if(bulk_deltaf_kind == 0)
                        delta_f_bulk = (
                            -(1. - F0_IS_NOT_SMALL*sign*f0)*bulkPi
                            *(bulkvisCoefficients[0]*mass*mass 
                              + bulkvisCoefficients[1]*pdotu 
                              + bulkvisCoefficients[2]*pdotu*pdotu));
                    else if (bulk_deltaf_kind == 1)
                    {

                        double E_over_T = pdotu/Tdec;
                        double mass_over_T = mass/Tdec;
                        delta_f_bulk = (
                            -1.0*(1.-sign*f0)/E_over_T*bulkvisCoefficients[0]
                            *(mass_over_T*mass_over_T/3. 
                              - bulkvisCoefficients[1]*E_over_T*E_over_T)
                            *bulkPi);
                    }
                    else if (bulk_deltaf_kind == 2)
                    {
                        double E_over_T = pdotu/Tdec;
                        delta_f_bulk = (
                            -1.*(1.-sign*f0)
                            *(-bulkvisCoefficients[0] 
                              + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                    }
                    else if (bulk_deltaf_kind == 3)
                    {
                        double E_over_T = pdotu/Tdec;
                        delta_f_bulk = (
                            -1.0*(1.-sign*f0)/sqrt(E_over_T)
                            *(-bulkvisCoefficients[0] 
                              + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                    }
                    else if (bulk_deltaf_kind == 4)
                    {
                        double E_over_T = pdotu/Tdec;
                        delta_f_bulk = (
                            -1.0*(1.-sign*f0)
                            *(bulkvisCoefficients[0] 
                              - bulkvisCoefficients[1]/E_over_T)*bulkPi);
                    }
                }

                // delta f for diffusion
                double delta_f_qmu = 0.0;
                if(INCLUDE_DIFFUSION_DELTAF == 1)
                {
                    double qmufactor = pt*qmu0 - px*qmu1 - py*qmu2 - pz*qmu3;
                    delta_f_qmu = ((1. - sign*f0)
                                   *(prefactor_qmu - baryon/pdotu)
                                   *qmufactor/deltaf_qmu_coeff);
                }
            
                double result;
                // restrict the size of delta f to be smaller than f_0
                double ratio_max = 1.0;
                double deltaf_size = fabs(delta_f_shear + delta_f_bulk
                                          + delta_f_qmu);
                double resize_factor = min(1., 
                                           ratio_max/(deltaf_size + 1e-10));
                result = (prefactor*degen*f0*pdsigma*tau
                       *(1. + (delta_f_shear + delta_f_bulk + delta_f_qmu)
                               *resize_factor));

                // Note that the factor 1.0 used here assumes that the maximum 
                // on the discrete lattice is the same as the maximum of the 
                // actual function, which is of course only an approximation
                double test = (result*pT
                               /dN_dxtdetady_pT_max[y_minus_eta_s_idx][FO_idx]
                               /1.0);

                // for debugging; 1.0: the maximum on the discrete lattice 
                // may not be the maximum for the actual continuous function
                if (AMOUNT_OF_OUTPUT>1 && test > 1.)
                    cout << "WTH?!" << endl; 

                if (drand48() < test) 
                    break;  // accept sample! 

                tries ++;
            }
            if (tries == maximum_impatience && AMOUNT_OF_OUTPUT>5)
                cout << "EmissionFunctionArray::"
                     << "sample_using_dN_dxtdetady_smooth_pT_phi warning: "
                     << "maximum_impatience reached." << endl;

            if (tries == maximum_impatience)
                continue; // resample
            else 
                sample_idx++; // write-out sample

            if (positive_y_minus_eta_table_only)
                y_minus_eta_s = irand(0,1)==0 ? y_minus_eta_s : -y_minus_eta_s;

            double y = y_LB + drand48()*(y_RB-y_LB);
            double eta_s = y - y_minus_eta_s;
            double p_z = mT*sinh(y);
            double E = mT*cosh(y);
            double z = surf->tau*sinh(eta_s);
            double t = surf->tau*cosh(eta_s);
            double vx = ux/gammaT;
            double vy = uy/gammaT;

            // write to sample file
            if (!USE_OSCAR_FORMAT)
            {
                sprintf(line_buffer, 
                        "%lu  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
                        FO_idx, surf->tau, surf->xpt, surf->ypt, y_minus_eta_s, 
                        pT, phi, surf->da0, surf->da1, surf->da2, vx, vy, y, 
                        eta_s, E, p_z, t, z);
            }
            // To be combined to OSCAR
            if (USE_OSCAR_FORMAT)
            {
                sprintf(line_buffer, 
                        "%24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e\n", 
                        px, py, p_z, E, mass, surf->xpt, surf->ypt, z, t);
            }
            sample_str_buffer << line_buffer;
            sample_writing_signal++;
            if (sample_writing_signal==NUMBER_OF_LINES_TO_WRITE)
            {
              of_sample << sample_str_buffer.str();
              sample_str_buffer.str("");
              sample_writing_signal=0;
            }
        }
        if (AMOUNT_OF_OUTPUT>0)
            print_progressbar((double)(sampling_idx)
                                      /number_of_repeated_sampling);
    }
    // flushing buffers
    if (control_writing_signal!=0)
    {
      of_control << control_str_buffer.str();
      control_str_buffer.str("");
      control_writing_signal = 0;
    }
    of_control.close();
    if (sample_writing_signal!=0)
    {
      of_sample << sample_str_buffer.str();
      sample_str_buffer.str("");
      sample_writing_signal = 0;
    }
    of_sample.close();
    if (AMOUNT_OF_OUTPUT>0) print_progressbar(1);

    ofstream of_sample_format(samples_format_filename.c_str());
    // Be careful that ParameterReader class convert all strings to lower case
    // so do NOT use variable names that differ only by cases!
    if (!USE_OSCAR_FORMAT)
    {
        of_sample_format
        << "Total_number_of_columns = " << 18 << endl
        << "FZ_cell_idx = " << 1 << endl
        << "tau = " << 2 << endl
        << "FZ_x = " << 3 << endl
        << "FZ_y = " << 4 << endl
        << "y_minus_eta_s = " << 5 << endl
        << "pT = " << 6 << endl
        << "phi = " << 7 << endl
        << "surf_da0 = " << 8 << endl
        << "surf_da1 = " << 9 << endl
        << "surf_da2 = " << 10 << endl
        << "surf_vx = " << 11 << endl
        << "surf_vy = " << 12 << endl
        << "y = " << 13 << endl
        << "eta_s = " << 14 << endl
        << "E = " << 15 << endl
        << "p_z = " << 16 << endl
        << "t = " << 17 << endl
        << "z = " << 18 << endl;
    }
    if (USE_OSCAR_FORMAT)
    {
        of_sample_format
        << "Total_number_of_columns = " << 9 << endl
        << "t = " << 9 << endl
        << "FZ_x = " << 6 << endl
        << "FZ_y = " << 7 << endl
        << "z = " << 8 << endl
        << "E = " << 4 << endl
        << "px = " << 1 << endl
        << "py = " << 2 << endl
        << "p_z = " << 3 << endl
        << "mass =" << 5 << endl;
    }
    of_sample_format.close();

    if(INCLUDE_BULK_DELTAF == 1)
        delete [] bulkvisCoefficients;

    sw.toc();
    cout << endl 
         << "Sampling finished in " << sw.takeTime() << " seconds." << endl;

}
//***************************************************************************



//***************************************************************************
void EmissionFunctionArray::sample_using_dN_pTdpTdphidy()
// Sample according to the dN_pTdpTdphidy and dN_pTdpTdphidy_max array. These
// arrays are asssumed to be already cacluated.
// The particle index is assumed to be already stored in the
// "last_particle_idx" internal variable.
// Emission function below "zero" will be treated as "zero".
// Each line in the sample file contains the following data:
// FZ_cell_idx tau, xpt, ypt, y_minus_eta_s, pT, phi, da0, da1, da2, vx, vy
// Each line in the sample control file records the number of particles
// sampled in one event.
// The format file is a file that can be read by ParameterReader class
// which records the the purpose of each column in the sample file. For
// example, if one line of it reads as "foo=1" means that the variable
// "foo" is the recorded in the 1st column.
// Parameters:
// -- number_of_repeated_sampling: number of successive sampling; in each 
//    sampling the number of particles is determined by dN_dy and the 
//    model parameter.
// -- pT_to: sample up to what pT
// -- model parameter:
//    1): The fractional part of dN_dy is used as a probability to determine
//       whether we have 0 or 1 more particle.
{
    Stopwatch sw;
    sw.tic();

    double pT_to = paraRdr->getVal("sample_pT_up_to");
    if (pT_to<0)
        pT_to = pT_tab->getLast(1); // use table to determine pT range
    // If dN/(dx_t deta dy) is evaluated to be smaller than this value, 
    // then it is replaced by this value.
    double zero = paraRdr->getVal("minimum_emission_function_val"); 
    long number_of_repeated_sampling = paraRdr->getVal(
                                                "number_of_repeated_sampling");

    particle_info* particle = &particles[last_particle_idx];

    double mass = particle->mass;
    double sign = particle->sign;
    double degen = particle->gspin;
    double baryon = particle->baryon;

    double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

    FO_surf* surf = &FOsurf_ptr[0];

    double *bulkvisCoefficients;
    if(INCLUDE_BULK_DELTAF == 1)
    {
        if(bulk_deltaf_kind == 0)
            bulkvisCoefficients = new double [3];
        else 
            bulkvisCoefficients = new double [2];
    }

    //------------------------------------------------------------------------
    // Prepare the inverse CDF
    // Step 1) interp dN and dN_max matrices from pT_tab and phi_tab to 
    // pT_tab4Sampling and phi_tab4Sampling allocation
    double **dN_pTdpTdphidy_with_weight_4Sampling = (
                                         new double* [pT_tab4Sampling_length]);
    for (int i=0; i<pT_tab4Sampling_length; i++)
        dN_pTdpTdphidy_with_weight_4Sampling[i] = (
                                         new double [phi_tab4Sampling_length]);

    double **dN_pTdpTdphidy_max_4Sampling = (
                                         new double* [pT_tab4Sampling_length]);
    for (int i=0; i<pT_tab4Sampling_length; i++)
        dN_pTdpTdphidy_max_4Sampling[i] = new double [phi_tab4Sampling_length];

    // interpolation
    for (int i=0; i<pT_tab4Sampling_length; i++)
    {
        for (int j=0; j<phi_tab4Sampling_length; j++)
        {
            double pT = pT_tab4Sampling.get(1, i+1); // get pT
            double pT_weight = pT_tab4Sampling.get(2, i+1); // get pT weight
            double pT_index = pT_tab4Sampling.get(3, i+1); // get pT index
            double phi_weight = phi_tab4Sampling.get(2, j+1); // get phi weight
            double phi_index = phi_tab4Sampling.get(3, j+1); // get phi index
            dN_pTdpTdphidy_with_weight_4Sampling[i][j] = (
                dN_pTdpTdphidy->interp2(pT_index, phi_index, 4)
                *pT*pT_weight*phi_weight); 
            // interp2: parameter 2 -> allow extrapolation
            dN_pTdpTdphidy_max_4Sampling[i][j] = (
                dN_pTdpTdphidy_max->interp2(pT_index, phi_index, 4));
        }
    }
        
    // create random variable using inverse CDF
    RandomVariable2DArray rand2D(dN_pTdpTdphidy_with_weight_4Sampling, 
                                 phi_tab4Sampling_length, 
                                 pT_tab4Sampling_length, zero);

    //-----------------------------------------------------------------------
    // start sampling
    // first calcualte total number of particles
    double dN_dy = 0.0;
    for (long i=0; i<pT_tab4Sampling_length; i++)
        for (int j=0; j<phi_tab4Sampling_length; j++)
            dN_dy += dN_pTdpTdphidy_with_weight_4Sampling[i][j];

    // dN_dy *= 100000;

    // create caches
    double mT_tab[pT_tab4Sampling_length]; // mT table
    for (long i=0; i<pT_tab4Sampling_length; i++)
        mT_tab[i] = sqrt(mass*mass + pT_tab4Sampling.get(1,i+1)
                                     *pT_tab4Sampling.get(1,i+1));

    // prepare for outputs
    // the control file records how many particles are there in each sampling
    char samples_control_filename_buffer[300];
    sprintf(samples_control_filename_buffer, samples_control_filename.c_str(), 
            particle->monval);
    remove(samples_control_filename_buffer);
    ofstream of_control(samples_control_filename_buffer);

    // the sample file contains the actual samples
    char samples_filename_buffer[300];
    sprintf(samples_filename_buffer, samples_filename.c_str(), 
            particle->monval);
    remove(samples_filename_buffer);
    ofstream of_sample(samples_filename_buffer);

    // buffers are used to speed up the output process
    char line_buffer[500]; // only used in text mode
    stringstream sample_str_buffer; // to speed up outputing process
    stringstream control_str_buffer;

    int sampling_model = paraRdr->getVal("dN_dy_sampling_model");
    double sampling_para1 = paraRdr->getVal("dN_dy_sampling_para1");
    double sampling_para2 = paraRdr->getVal("dN_dy_sampling_para2");
    double sampling_para3 = paraRdr->getVal("dN_dy_sampling_para3");
    double sampling_para4 = paraRdr->getVal("dN_dy_sampling_para4");
    double sampling_para5 = paraRdr->getVal("dN_dy_sampling_para5");

    // get y range for sampling
    double y_LB = paraRdr->getVal("y_LB");
    double y_RB = paraRdr->getVal("y_RB");

    cout << " dN_dy=" << dN_dy << ", ";
    double dN = (y_RB-y_LB)*dN_dy;
    cout << "dN=" << dN << "..." << endl;
    if (AMOUNT_OF_OUTPUT>0) print_progressbar(-1);
    long sample_writing_signal = 0;
    long control_writing_signal = 0;
    long maximum_impatience = 5000; // used in pt-phi sampling

    long total_tries = 0, total_violation = 0;
    //double largest_violation = -1.0;
    
    for (long sampling_idx=1; sampling_idx<=number_of_repeated_sampling; 
         sampling_idx++)
    {
        long number_to_sample = determine_number_to_sample(
            dN, sampling_model, sampling_para1, sampling_para2, sampling_para3, 
            sampling_para4, sampling_para5);

        // write to control file
        sprintf(line_buffer, "%lu\n", number_to_sample);
        control_str_buffer << line_buffer;
        control_writing_signal++;
        if (control_writing_signal==NUMBER_OF_LINES_TO_WRITE)
        {
            of_control << control_str_buffer.str();
            control_str_buffer.str("");
            control_writing_signal=0;
        }

        long pT_idx, phi_idx;
        long FO_idx, y_minus_eta_s_idx;

        long sample_idx = 1;
        while (sample_idx <= number_to_sample)
        {
            // first, sample eta and freeze-out cell index
            rand2D.sampleAccToInvCDF(&pT_idx, &phi_idx);
            double pT = pT_tab4Sampling.get(1,pT_idx+1);
            double phi = phi_tab4Sampling.get(1,phi_idx+1);

            double mT = mT_tab[pT_idx];
            double px = pT*trig_phi_tab4Sampling[phi_idx][0];
            double py = pT*trig_phi_tab4Sampling[phi_idx][1];

            // use local variable to substitute global ones
            double dN_max_sampling = (
                            dN_pTdpTdphidy_max_4Sampling[pT_idx][phi_idx]);

            long tries = 1;
            while (tries<maximum_impatience)
            {
                bool found_sample = false;

                // get a surface index
                FO_idx = floor(drand48()*(FO_length-1e-30));

                surf = &FOsurf_ptr[FO_idx];

                double Tdec = surf->Tdec;
                double inv_Tdec = 1.0/Tdec;
                double Pdec = surf->Pdec;
                double Edec = surf->Edec;

                double tau = surf->tau;

                double gammaT = surf->u0;
                double ux = surf->u1;
                double uy = surf->u2;

                double mu = surf->particle_mu[last_particle_idx];

                double da0 = surf->da0;
                double da1 = surf->da1;
                double da2 = surf->da2;
                double pi00 = surf->pi00;
                double pi01 = surf->pi01;
                double pi02 = surf->pi02;
                double pi11 = surf->pi11;
                double pi12 = surf->pi12;
                double pi22 = surf->pi22;
                double pi33 = surf->pi33;
                double bulkPi = 0.0;
                double deltaf_prefactor = 0.0;
                if(INCLUDE_DELTAF)
                    deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec));

                if(INCLUDE_BULK_DELTAF == 1)
                {
                    if(bulk_deltaf_kind == 0)
                        bulkPi = surf->bulkPi;
                    else 
                        bulkPi = surf->bulkPi/hbarC;   // unit in fm^-4 
                    getbulkvisCoefficients(Tdec, bulkvisCoefficients);
                }

                // diffusion delta f
                double qmu0 = 0.0;
                double qmu1 = 0.0;
                double qmu2 = 0.0;
                double qmu3 = 0.0;
                double deltaf_qmu_coeff = 1.0;
                double prefactor_qmu = 0.0;
                if(INCLUDE_DIFFUSION_DELTAF == 1)
                {
                    qmu0 = surf->qmu0;
                    qmu1 = surf->qmu1;
                    qmu2 = surf->qmu2;
                    qmu3 = surf->qmu3;
                    
                    double mu_B = surf->muB;
                    deltaf_qmu_coeff = get_deltaf_qmu_coeff(Tdec, mu_B);
                    
                    double rho_B = surf->Bn;
                    prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
                }

                y_minus_eta_s_idx = y_minus_eta_min_index;
                double pt_max = (
                        mT*hypertrig_y_minus_eta_table[y_minus_eta_s_idx][0]);
                double pz_max = (
                        mT*hypertrig_y_minus_eta_table[y_minus_eta_s_idx][1]);

                double pdotu_max = pt_max*gammaT - px*ux - py*uy;
                double expon_max = (pdotu_max - mu) * inv_Tdec;
                double f0_max = 1./(exp(expon_max)+sign);

                double pdsigma_max = pt_max*da0 + px*da1 + py*da2;

                //viscous corrections
                double delta_f_shear = 0.0;
                if(INCLUDE_DELTAF)
                {
                    double Wfactor_max = (
                        pt_max*pt_max*pi00 - 2.0*pt_max*px*pi01 
                        - 2.0*pt_max*py*pi02 + px*px*pi11 
                        + 2.0*px*py*pi12 + py*py*pi22 
                        + pz_max*pz_max*pi33);
                    delta_f_shear = (
                        (1. - F0_IS_NOT_SMALL*sign*f0_max)
                        *Wfactor_max*deltaf_prefactor);
                }

                double delta_f_bulk = 0.0;
                if (INCLUDE_BULK_DELTAF == 1)
                {
                    if(bulk_deltaf_kind == 0)
                        delta_f_bulk = (
                            -(1. - F0_IS_NOT_SMALL*sign*f0_max)*bulkPi
                            *(bulkvisCoefficients[0]*mass*mass 
                              + bulkvisCoefficients[1]*pdotu_max 
                              + bulkvisCoefficients[2]*pdotu_max*pdotu_max));
                    else if (bulk_deltaf_kind == 1)
                    {

                        double E_over_T = pdotu_max/Tdec;
                        double mass_over_T = mass/Tdec;
                        delta_f_bulk = (
                            -1.0*(1.-sign*f0_max)/E_over_T
                            *bulkvisCoefficients[0]
                            *(mass_over_T*mass_over_T/3. 
                              - bulkvisCoefficients[1]*E_over_T*E_over_T)
                            *bulkPi);
                    }
                    else if (bulk_deltaf_kind == 2)
                    {
                        double E_over_T = pdotu_max/Tdec;
                        delta_f_bulk = (
                            -1.*(1.-sign*f0_max)
                            *(-bulkvisCoefficients[0] 
                              + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                    }
                    else if (bulk_deltaf_kind == 3)
                    {
                        double E_over_T = pdotu_max/Tdec;
                        delta_f_bulk = (
                            -1.0*(1.-sign*f0_max)/sqrt(E_over_T)
                            *(-bulkvisCoefficients[0] 
                              + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                    }
                    else if (bulk_deltaf_kind == 4)
                    {
                        double E_over_T = pdotu_max/Tdec;
                        delta_f_bulk = (
                            -1.0*(1.-sign*f0_max)
                            *(bulkvisCoefficients[0] 
                              - bulkvisCoefficients[1]/E_over_T)*bulkPi);
                    }
                }

                // delta f for diffusion
                double delta_f_qmu = 0.0;
                if(INCLUDE_DIFFUSION_DELTAF == 1)
                {
                    double qmufactor = (
                        pt_max*qmu0 - px*qmu1 - py*qmu2 - pz_max*qmu3);
                    delta_f_qmu = ((1. - sign*f0_max)
                                   *(prefactor_qmu - baryon/pdotu_max)
                                   *qmufactor/deltaf_qmu_coeff);
                }
                
                double results_max;
                delta_f_shear = max(delta_f_shear, 0.0);
                delta_f_bulk = max(delta_f_bulk, 0.0);
                delta_f_qmu = max(delta_f_qmu, 0.0);
                results_max = (prefactor*degen*f0_max
                         *(1. + delta_f_shear + delta_f_bulk + delta_f_qmu)
                         *pdsigma_max*tau);

                // for debugging; 1.0: the maximum on the discrete lattice 
                // may not be the maximum for the actual continuous function
                if (results_max/dN_max_sampling > (1.0+1e-6) )
                {
                    total_violation++;
                } 
                
                if (drand48() > results_max/dN_max_sampling)
                {
                    tries++; 
                    continue;
                } // discard this freeze-out cell
                
                for (int sub_try=0; sub_try<y_minus_eta_tab_length/2; 
                     sub_try++)
                {
                    
                    // get a y-eta_s index
                    y_minus_eta_s_idx = (floor(drand48()
                                         *(y_minus_eta_tab_length-1e-30))); 
                    double pt = (
                        mT*hypertrig_y_minus_eta_table[y_minus_eta_s_idx][0]);
                    double pz = (
                        mT*hypertrig_y_minus_eta_table[y_minus_eta_s_idx][1]);

                    double pdotu = (pt*gammaT - px*ux - py*uy);
                    double expon = (pdotu - mu) * inv_Tdec;
                    double f0 = 1./(exp(expon)+sign);

                    double pdsigma = pt*da0 + px*da1 + py*da2;

                    //viscous corrections
                    double delta_f_shear = 0.0;
                    if(INCLUDE_DELTAF)
                    {
                        double Wfactor = (
                            pt*pt*pi00 - 2.0*pt*px*pi01 - 2.0*pt*py*pi02 
                            + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 
                            + pz*pz*pi33);
                        delta_f_shear = (
                            (1. - F0_IS_NOT_SMALL*sign*f0)
                            *Wfactor*deltaf_prefactor);
                    }
                    double delta_f_bulk = 0.0;
                    if (INCLUDE_BULK_DELTAF == 1)
                    {
                        if(bulk_deltaf_kind == 0)
                            delta_f_bulk = (
                                -(1. - F0_IS_NOT_SMALL*sign*f0)*bulkPi
                                *(bulkvisCoefficients[0]*mass*mass 
                                  + bulkvisCoefficients[1]*pdotu 
                                  + bulkvisCoefficients[2]*pdotu*pdotu));
                        else if (bulk_deltaf_kind == 1)
                        {

                            double E_over_T = pdotu/Tdec;
                            double mass_over_T = mass/Tdec;
                            delta_f_bulk = (
                                -1.0*(1.-sign*f0)/E_over_T
                                *bulkvisCoefficients[0]
                                *(mass_over_T*mass_over_T/3. 
                                  - bulkvisCoefficients[1]*E_over_T*E_over_T)
                                *bulkPi);
                        }
                        else if (bulk_deltaf_kind == 2)
                        {
                            double E_over_T = pdotu/Tdec;
                            delta_f_bulk = (
                                -1.*(1.-sign*f0)
                                *(-bulkvisCoefficients[0] 
                                  + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                        }
                        else if (bulk_deltaf_kind == 3)
                        {
                            double E_over_T = pdotu/Tdec;
                            delta_f_bulk = (
                                -1.0*(1.-sign*f0)/sqrt(E_over_T)
                                *(-bulkvisCoefficients[0] 
                                  + bulkvisCoefficients[1]*E_over_T)*bulkPi);
                        }
                        else if (bulk_deltaf_kind == 4)
                        {
                            double E_over_T = pdotu/Tdec;
                            delta_f_bulk = (
                                -1.0*(1.-sign*f0)
                                *(bulkvisCoefficients[0] 
                                  - bulkvisCoefficients[1]/E_over_T)*bulkPi);
                        }
                    }

                    // delta f for diffusion
                    double delta_f_qmu = 0.0;
                    if(INCLUDE_DIFFUSION_DELTAF == 1)
                    {
                        double qmufactor = (
                                        pt*qmu0 - px*qmu1 - py*qmu2 - pz*qmu3);
                        delta_f_qmu = ((1. - sign*f0)
                                       *(prefactor_qmu - baryon/pdotu)
                                       *qmufactor/deltaf_qmu_coeff);
                    }

                    double result;
                    // restrict the size of delta f to be smaller than f_0
                    double ratio_max = 1.0;
                    double deltaf_size = fabs(delta_f_shear + delta_f_bulk
                                              + delta_f_qmu);
                    double resize_factor = min(1., 
                                              ratio_max/(deltaf_size + 1e-10));
                    result = (prefactor*degen*f0*pdsigma*tau*(1. 
                              + (delta_f_shear + delta_f_bulk + delta_f_qmu)
                                *resize_factor));

                    if ( result/results_max > (1.0+1e-6))
                    {
                        total_violation++;
                    } 
                    
                    if (drand48()<result/results_max)
                    {
                        found_sample=true; 
                        break;
                    } // accept sample! 

                    tries ++;

                }

                if (found_sample) break;
            }

            total_tries += tries;

            if (tries >= maximum_impatience && AMOUNT_OF_OUTPUT > 5)
                cout << "EmissionFunctionArray::"
                     << "sample_using_dN_pTdpTdphidy warning: "
                     << "maximum_impatience reached." << endl;

            if (tries >= maximum_impatience)
                continue; // resample
            else
                sample_idx++; // write-out sample

            double y_minus_eta_s = y_minus_eta_tab->get(1,y_minus_eta_s_idx+1);
            if (positive_y_minus_eta_table_only)
                y_minus_eta_s = irand(0,1)==0 ? y_minus_eta_s : -y_minus_eta_s;

            double y = y_LB + drand48()*(y_RB-y_LB);
            double eta_s = y - y_minus_eta_s;
            double p_z = mT*sinh(y);
            double E = mT*cosh(y);
            double z = surf->tau*sinh(eta_s);
            double t = surf->tau*cosh(eta_s);

            // write to sample file
            if (!USE_OSCAR_FORMAT)
            {
                sprintf(line_buffer, 
                        "%lu  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
                        FO_idx, surf->tau, surf->xpt, surf->ypt, y_minus_eta_s, 
                        pT, phi, surf->da0, surf->da1, surf->da2, 
                        surf->u1/surf->u0, surf->u2/surf->u0, y, eta_s, E, 
                        p_z, t, z);
            }
            // To be combined to OSCAR
            if (USE_OSCAR_FORMAT)
            {
                sprintf(line_buffer, 
                        "%24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e\n", 
                        px, py, p_z, E, mass, surf->xpt, surf->ypt, z, t);
            }
            sample_str_buffer << line_buffer;
            sample_writing_signal++;
            if (sample_writing_signal==NUMBER_OF_LINES_TO_WRITE)
            {
              of_sample << sample_str_buffer.str();
              sample_str_buffer.str("");
              sample_writing_signal=0;
            }

        }
        if (AMOUNT_OF_OUTPUT > 0)
            print_progressbar((double)(sampling_idx)
                                      /number_of_repeated_sampling);
    }

    // flushing buffers
    if (control_writing_signal!=0)
    {
      of_control << control_str_buffer.str();
      control_str_buffer.str("");
      control_writing_signal = 0;
    }
    of_control.close();
    if (sample_writing_signal!=0)
    {
      of_sample << sample_str_buffer.str();
      sample_str_buffer.str("");
      sample_writing_signal = 0;
    }
    of_sample.close();
    if (AMOUNT_OF_OUTPUT>0) print_progressbar(1);

    cout << endl << "Average number of tries: " 
         << total_tries/number_of_repeated_sampling/dN << endl;
    cout << endl << "Average number of violations: " 
         << float(total_violation)/number_of_repeated_sampling/dN << endl;
    
    ofstream of_sample_format(samples_format_filename.c_str());
    // Be careful that ParameterReader class convert all strings to lower case
    // so do NOT use variable names that differ only by cases!
    if (!USE_OSCAR_FORMAT)
    {
        of_sample_format
        << "Total_number_of_columns = " << 18 << endl
        << "FZ_cell_idx = " << 1 << endl
        << "tau = " << 2 << endl
        << "FZ_x = " << 3 << endl
        << "FZ_y = " << 4 << endl
        << "y_minus_eta_s = " << 5 << endl
        << "pT = " << 6 << endl
        << "phi = " << 7 << endl
        << "surf_da0 = " << 8 << endl
        << "surf_da1 = " << 9 << endl
        << "surf_da2 = " << 10 << endl
        << "surf_vx = " << 11 << endl
        << "surf_vy = " << 12 << endl
        << "y = " << 13 << endl
        << "eta_s = " << 14 << endl
        << "E = " << 15 << endl
        << "p_z = " << 16 << endl
        << "t = " << 17 << endl
        << "z = " << 18 << endl;
    }
    if (USE_OSCAR_FORMAT)
    {
        of_sample_format
        << "Total_number_of_columns = " << 9 << endl
        << "t = " << 9 << endl
        << "FZ_x = " << 6 << endl
        << "FZ_y = " << 7 << endl
        << "z = " << 8 << endl
        << "E = " << 4 << endl
        << "px = " << 1 << endl
        << "py = " << 2 << endl
        << "p_z = " << 3 << endl
        << "mass =" << 5 << endl;
    }
    of_sample_format.close();

    // delete the intermedia variables (very important)
    for (int i=0; i<pT_tab4Sampling_length; i++)
        delete[] dN_pTdpTdphidy_with_weight_4Sampling[i];
    delete[] dN_pTdpTdphidy_with_weight_4Sampling;

    for (int i=0; i<pT_tab4Sampling_length; i++)
        delete[] dN_pTdpTdphidy_max_4Sampling[i];
    delete[] dN_pTdpTdphidy_max_4Sampling;

    if(INCLUDE_BULK_DELTAF == 1)
        delete [] bulkvisCoefficients;

    sw.toc();
    cout << endl 
         << "Sampling finished in " << sw.takeTime() << " seconds." << endl;
}
//***************************************************************************


//***************************************************************************
inline long EmissionFunctionArray::determine_number_to_sample(
    double dN_dy_in, int model, double para1, double para2, 
    double para3, double para4, double para5)
// From a non-integer averaged particles number dN, return an actual interger
// particle number that can be used in sampling.
{
    if (dN_dy_in<0)
    {
        cout << "EmissionFunctionArray::"
             << "determine_number_to_sample error: "
             << "dN_dy should be positive but receives " << dN_dy_in << endl;
        exit(-1);
    }
    double dN_dy = dN_dy_in;
    long number_to_sample;
    long dN_dy_int = long(dN_dy); // integer part
    double dN_dy_fraction = dN_dy - dN_dy_int; // fractional part
    double p,k; // parameters in NBD

    // determine actual number of particles
    switch (model) // explained in parameters.dat
    {
        case 1: // with possibly 1 more particle
            number_to_sample = dN_dy_int;
            // lucky! 1 more particle...
            if (drand48() < dN_dy_fraction)
                number_to_sample++; 
            break;
        case 10: // use NBD
            k = para1*dN_dy_fraction;
            p = 1.0/(1.0+para1);
            if (k<1e-15)
                number_to_sample = dN_dy_int;
            else 
                number_to_sample = (
                    dN_dy_int + gsl_ran_negative_binomial(gsl_random_r, p, k));
            break;
        case 20:
            k = para1*dN_dy;
            p = 1.0/(1.0+para1);
            if (k<1e-15) 
                number_to_sample = dN_dy_int;
            else 
                number_to_sample = (
                                gsl_ran_negative_binomial(gsl_random_r, p, k));
            break;
        case 30:  // use Poisson distribution
            if (dN_dy < 1e-15)
                number_to_sample = 0;
            else
                //number_to_sample = poissonDistribution.rand(dN_dy);
                number_to_sample = gsl_ran_poisson(gsl_random_r, dN_dy);
            break;
        default:
            cout << "EmissionFunctionArray::"
                 << "determine_number_to_sample error: "
                 << "the model specified (" << model << ") "
                 << "to determine the number of particles is not found" 
                 << endl;
            exit(-1);
    }

    return number_to_sample;
}
//***************************************************************************



//***************************************************************************
void EmissionFunctionArray::calculate_dN_dxtdetady_and_sample_4all()
// Calculate dN_dxtdetady then sample.
{

    cout << endl
        << "***********************************************************"
        << endl
        << "Function calculate_dN_dxtdetady_and_sample_4all started... " 
        << endl;
    Stopwatch sw;
    sw.tic();

    // read in parameters
    int calculate_dN_dtau = paraRdr->getVal("calculate_dN_dtau");
    int calculate_dN_deta = paraRdr->getVal("calculate_dN_deta");
    int calculate_dN_dxt = paraRdr->getVal("calculate_dN_dxt");
    int calculate_dN_dx = paraRdr->getVal("calculate_dN_dx");

    // loop over chosen particles
    particle_info* particle = NULL;
    for (int m=0; m<number_of_chosen_particles; m++)
    {
        int particle_idx = chosen_particles_sampling_table[m];
        particle = &particles[particle_idx];
        cout << "Index: " << m << ", Name: " << particle->name 
             << ", Monte-carlo index: " << particle->monval << endl;
        if (m > 0. && 
            particles_are_the_same(particle_idx, 
                                   chosen_particles_sampling_table[m-1]))
        {
           cout << " -- Using dN_dxdetady from previous calculation... " 
                << endl;
           last_particle_idx = particle_idx; // fake a calculation
        }
        else
        {
            cout << " -- Calculating dN_dxdetady... " << endl;
            calculate_dN_dxtdetady(particle_idx);
            write_dN_dxtdetady_toFile();
        }

        if (calculate_dN_dtau)
        {
            double tau0 = paraRdr->getVal("bin_tau0");
            double dtau = paraRdr->getVal("bin_dtau");
            double tau_max = paraRdr->getVal("bin_tau_max");
            calculate_dN_dtau_using_dN_dxtdetady(tau0, dtau, tau_max);
        }
        
        if (calculate_dN_dx)
        {
            double x_min = paraRdr->getVal("bin_x_min");
            double dx = paraRdr->getVal("bin_dx");
            double x_max = paraRdr->getVal("bin_x_max");
            calculate_dN_dx_using_dN_dxtdetady(x_min, x_max, dx);
        }

        if (calculate_dN_deta)
        {
            calculate_dN_deta_using_dN_dxtdetady();
        }

        if (calculate_dN_dxt)
        {
            calculate_dN_dxt_using_dN_dxtdetady();
        }

        cout << " -- Sampling using dN_dxtdetady: ";
        sample_using_dN_dxtdetady_smooth_pT_phi();

    }

    sw.toc();
    cout << " -- calculate_dN_dxtdetady_and_sample_4all finishes " 
         << sw.takeTime() << " seconds." << endl;
}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_dN_dtau_using_dN_dxtdetady(
                                      double tau0, double dtau, double tau_max)
// Calculate dN/dtau. Should be called after calculate_dN_dxtdetady. 
// The emission function will be binned into bins tau0:dtau:tau_max. 
// The result will be written directly to file.
// Each line of the output file has the format:
// tau-at-the-center-of-the-bin mean-tau dN_dtau count
{
    // create local cache
    double delta_y_minus_eta_tab[y_minus_eta_tab_length];
    for (int k=0; k<y_minus_eta_tab_length; k++)
        delta_y_minus_eta_tab[k] = y_minus_eta_tab->get(2,k+1);

    // use a buffer to store summed data
    long number_of_bins = (tau_max-tau0)/dtau;
    vector<double> sum_dN(number_of_bins, 0);
    vector<double> sum_tau(number_of_bins, 0);
    vector<double> count(number_of_bins,0);

    // construct bins
    vector<double> bins;
    for (double tau=tau0; tau<tau_max; tau+=dtau)
        bins.push_back(tau);

    // summing to bins
    for (int k=0; k<y_minus_eta_tab_length; k++)
    for (long l=0; l<FO_length; l++)
    {
        double tau = FOsurf_ptr[l].tau;
        long idx = binarySearch(&bins, tau, true);
        if (idx==-1)
            continue; // skip those not falling in any bins
        sum_dN[idx] += dN_dxtdetady[k][l]*delta_y_minus_eta_tab[k];
        sum_tau[idx] += tau;
        count[idx] ++;

    }

    // average them and output
    particle_info* particle = &particles[last_particle_idx];
    char dN_dtau_filename_buffer[300];
    sprintf(dN_dtau_filename_buffer, dN_dtau_filename.c_str(), 
            particle->monval);
    ofstream of(dN_dtau_filename_buffer);
    for (long idx=0; idx<number_of_bins; idx++)
        formatedPrint(of, 4, tau0+(idx+0.5)*dtau, 
                      sum_tau[idx]/(count[idx]+1e-30), sum_dN[idx]/dtau, 
                      count[idx]);
    of.close();

}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_dN_dx_using_dN_dxtdetady(
                                         double x_min, double x_max, double dx)
// Calculate dN/dx. Should be called after calculate_dN_dxtdetady. The emission
// function will be binned into bins x_min:x_max:dx. The result will be written
// directly to file.
// Each line of the output file has the format:
// x-at-the-center-of-the-bin mean-x dN_dx count
{
    // create local cache
    double delta_y_minus_eta_tab[y_minus_eta_tab_length];
    for (int k=0; k<y_minus_eta_tab_length; k++)
        delta_y_minus_eta_tab[k] = y_minus_eta_tab->get(2,k+1);

    // use a buffer to store summed data
    long number_of_bins = (x_max - x_min)/dx;
    vector<double> sum_dN1(number_of_bins, 0);   // dN/dydx1 @ |x2| < 0.5 fm
    vector<double> sum_dN2(number_of_bins, 0);   // dN/dydx2 @ |x1| < 0.5 fm
    vector<double> sum_x1(number_of_bins, 0);
    vector<double> sum_x2(number_of_bins, 0);
    vector<double> count1(number_of_bins,0);
    vector<double> count2(number_of_bins,0);

    // construct bins
    vector<double> bins;
    for (double x = x_min; x < x_max; x+=dx)
        bins.push_back(x);

    // summing to bins
    for (int k=0; k<y_minus_eta_tab_length; k++)
    for (long l=0; l<FO_length; l++)
    {
        double x_local = FOsurf_ptr[l].xpt;
        double y_local = FOsurf_ptr[l].ypt;
        if(fabs(y_local) < 0.5)
        {
            long idx = binarySearch(&bins, x_local, true);
            if (idx==-1)
                continue; // skip those not falling in any bins
            sum_dN1[idx] += dN_dxtdetady[k][l]*delta_y_minus_eta_tab[k];
            sum_x1[idx] += x_local;
            count1[idx] ++;
        }

        if(fabs(x_local) < 0.5)
        {
            long idx = binarySearch(&bins, y_local, true);
            if (idx==-1)
                continue; // skip those not falling in any bins
            sum_dN2[idx] += dN_dxtdetady[k][l]*delta_y_minus_eta_tab[k];
            sum_x2[idx] += y_local;
            count2[idx] ++;
        }
    }

    // average them and output
    particle_info* particle = &particles[last_particle_idx];
    char dN_dx_filename_buffer[300];
    sprintf(dN_dx_filename_buffer, dN_dx_filename.c_str(), particle->monval);
    ofstream of(dN_dx_filename_buffer);
    for (long idx=0; idx<number_of_bins; idx++) 
        formatedPrint(of, 7, x_min+(idx+0.5)*dx, 
                      sum_x1[idx]/(count1[idx]+1e-30), sum_dN1[idx]/dx, 
                      count1[idx], sum_x2[idx]/(count2[idx]+1e-30), 
                      sum_dN2[idx]/dx, count2[idx]);
    of.close();

}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_dN_dphi_using_dN_pTdpTdphidy()
// Calculate dN/dphi. Using the dN_pTdpTdphidy array.
// The output has three columns: phi dN/dphi phi_weight
{

    vector<double> dN_dphi(phi_tab_length,0);
    for (int i=0; i<pT_tab_length; i++)
    {
        double pT = pT_tab->get(1,i+1);
        double pT_weight = pT_tab->get(2,i+1);
        for (int j=0; j<phi_tab_length; j++)
        {
            dN_dphi[j] += dN_pTdpTdphidy->get(i+1, j+1)*pT*pT_weight;
        }

    }

    // output
    particle_info* particle = &particles[last_particle_idx];
    char dN_dphi_filename_buffer[300];
    sprintf(dN_dphi_filename_buffer, dN_dphi_filename.c_str(),
            particle->monval);
    ofstream of(dN_dphi_filename_buffer);
    for (int j=0; j<phi_tab_length; j++)
        formatedPrint(of, 3, phi_tab->get(1,j+1), dN_dphi[j], 
                      phi_tab->get(2,j+1));
    of.close();

}
//***************************************************************************




//***************************************************************************
void EmissionFunctionArray::calculate_dN_deta_using_dN_dxtdetady()
// Calculate dN/deta. Using the dN_dxtdetady array.
// The output has three columns: eta dN/deta eta_weight
{

    vector<double> dN_deta(y_minus_eta_tab_length,0);
    for (int k=0; k<y_minus_eta_tab_length; k++)
    for (long l=0; l<FO_length; l++)
    {
        dN_deta[k] += dN_dxtdetady[k][l];
    }

    // output
    particle_info* particle = &particles[last_particle_idx];
    char dN_deta_filename_buffer[300];
    sprintf(dN_deta_filename_buffer, dN_deta_filename.c_str(),
            particle->monval);
    ofstream of(dN_deta_filename_buffer);
    for (int k=0; k<y_minus_eta_tab_length; k++) 
        formatedPrint(of, 3, y_minus_eta_tab->get(1,k+1), dN_deta[k], 
                      dN_deta[k]*y_minus_eta_tab->get(2,k+1));
    of.close();

}


//***************************************************************************
void EmissionFunctionArray::calculate_dN_dxt_using_dN_dxtdetady()
{

    vector<double> dN_dxt(FO_length,0);
    for (long l=0; l<FO_length; l++)
    for (int k=0; k<y_minus_eta_tab_length; k++)
    {
        dN_dxt[l] += dN_dxtdetady[k][l]*y_minus_eta_tab->get(2, k+1);
    }

    // output
    particle_info* particle = &particles[last_particle_idx];
    char dN_dxt_filename_buffer[300];
    sprintf(dN_dxt_filename_buffer, dN_dxt_filename.c_str(), particle->monval);
    ofstream of(dN_dxt_filename_buffer);
    for (long l=0; l<FO_length; l++)
        formatedPrint(of, 2, (double)(l), dN_dxt[l]);
    of.close();
}



//***************************************************************************
bool EmissionFunctionArray::particles_are_the_same(int idx1, int idx2)
{
    if (particles[idx1].sign!=particles[idx2].sign)
        return false;
    if (particles[idx1].gspin!=particles[idx2].gspin)
        return false;
    double tolerance = paraRdr->getVal("grouping_tolerance");
    if (abs((particles[idx1].mass-particles[idx2].mass)
            /(particles[idx2].mass+1e-30)) > tolerance)
        return false;
    for (long l=0; l<FO_length; l++)
    {
        double chem1 = FOsurf_ptr[l].particle_mu[idx1];
        double chem2 = FOsurf_ptr[l].particle_mu[idx2];
        if (abs((chem1-chem2)/(chem2+1e-30)) > tolerance)
        {
          return false;
        }

    }
    return true;
}



//***************************************************************************
void EmissionFunctionArray::shell()
{
    // read in parameters
    int calculate_vn = paraRdr->getVal("calculate_vn");
    int historic_format = paraRdr->getVal("use_historic_flow_output_format");
    int MC_sampling = paraRdr->getVal("MC_sampling");

    int perform_sampling_during_calculation = 0;
    if (MC_sampling==3)
    {
        calculate_vn=1;
        perform_sampling_during_calculation = 1;
    }

    if (calculate_vn)
    {
        if (historic_format)
            calculate_dN_pTdpTdphidy_and_flows_4all_old_output(
                            perform_sampling_during_calculation);
        else
            calculate_dN_pTdpTdphidy_and_flows_4all(
                            perform_sampling_during_calculation);
    }

    if (MC_sampling==1)
    {
        calculate_dN_dxtdetady_and_sample_4all();
        if (USE_OSCAR_FORMAT)
            combine_samples_to_OSCAR();
    }
    else if (MC_sampling==2)
    {
        calculate_dN_dxtdy_4all_particles();
        sample_using_dN_dxtdy_4all_particles_conventional();
        if (USE_OSCAR_FORMAT)
            combine_samples_to_OSCAR();
    }
    else if (MC_sampling==3)
    {
        if (USE_OSCAR_FORMAT)
            combine_samples_to_OSCAR();
    }
}



//***************************************************************************
void EmissionFunctionArray::combine_samples_to_OSCAR()
{
    Stopwatch sw;
    sw.tic();
    cout << " -- Now combine sample files to OSCAR file..." << endl;

    char line_buffer[500];

    // open file for output
    remove(OSCAR_output_filename.c_str());
    ofstream oscar(OSCAR_output_filename.c_str());

    // write header first
    ifstream header(OSCAR_header_filename.c_str());
    if (!header.is_open())
    {
        cout << endl 
            << "combine_samples_to_OSCAR error: OSCAR header file " 
            << OSCAR_header_filename.c_str() << " not found." << endl;
        exit(-1);
    }
    while (true)
    {
        header.getline(line_buffer, 500);
        if (!header.eof()) oscar << line_buffer << endl;
        else break;
    }
    header.close();

    // open control and sample files
    vector<ifstream*> controls(number_of_chosen_particles); // control files
    vector<ifstream*> samples(number_of_chosen_particles); // sample files
    for (int m=0; m<number_of_chosen_particles; m++)
    {
        char filename_buffer[300];
        int monval = particles[chosen_particles_sampling_table[m]].monval;
        // control files first
        sprintf(filename_buffer, samples_control_filename.c_str(), monval);
        controls[m] = new ifstream ;
        controls[m]->open(filename_buffer);
        if (!controls[m]->is_open())
        {
            cout << endl 
                 << "combine_samples_to_OSCAR error: control file " 
                 << filename_buffer << " not found." << endl;
            exit(-1);
        }
        sprintf(filename_buffer, samples_filename.c_str(), monval);
        samples[m] = new ifstream ;
        samples[m]->open(filename_buffer);
        if (!samples[m]->is_open())
        {
            cout << endl 
                 << "combine_samples_to_OSCAR error: sample file " 
                 << filename_buffer << " not found." << endl;
            exit(-1);
        }
    }

    // big loop for generating OSCAR
    if (AMOUNT_OF_OUTPUT>0) print_progressbar(-1);
    long number_of_repeated_sampling = paraRdr->getVal(
                    "number_of_repeated_sampling");
    for (long sample_idx=1; sample_idx<=number_of_repeated_sampling; 
         sample_idx++)
    {
        // read-in number of particles for each species
        int number_of_particles[number_of_chosen_particles];
        long total_number_of_particles = 0;
        for (int m=0; m<number_of_chosen_particles; m++)
        {
            (*controls[m]) >> number_of_particles[m];
            total_number_of_particles += number_of_particles[m];
        }
        // sub-header for each event
        oscar << setw(10) << sample_idx << "  " 
              << setw(10) << total_number_of_particles << "  " 
              << setw(8) << 0.0 << "  " << setw(8) << 0.0 << endl;

        // now copy each line from samples file to OSCAR file
        long ipart = 1;
        for (int m=0; m<number_of_chosen_particles; m++)
        {
            int monval = particles[chosen_particles_sampling_table[m]].monval;
            for (long ii=1; ii<=number_of_particles[m]; ii++)
            {
                oscar << setw(10) << ipart << "  " 
                      << setw(10) << monval << "  ";
                samples[m]->getline(line_buffer, 500);
                oscar << line_buffer << endl;
                ipart ++;
            }
        }
        if (AMOUNT_OF_OUTPUT>0)
            print_progressbar((double)sample_idx/number_of_repeated_sampling);
    }

    sw.toc();
    cout << endl
         << " -- combine_samples_to_OSCAR samples finishes " 
         << sw.takeTime() << " seconds."
         << endl;

}


//--------------------------------------------------------------------------
//
// The following functions are added starting from Ver 2.0
//
//--------------------------------------------------------------------------


//***************************************************************************
void EmissionFunctionArray::calculate_dN_dxtdy_4all_particles()
/*
 The p_integral_table is a table stores results for the integral
   int( p^2 / ( exp(sqrt(p^2+m^2)-mu) + 1), p=0..inf)
 and m_integral_table is a table stores results for the integral
   int( p^2 / ( exp(sqrt(p^2+m^2)-mu) - 1), p=0..inf)
 for different m and mu values.
 The m values changes with row indices and it starts with
 integral_table_m0 with step integral_table_dm.
 The m-mu values changes with column indices and it starts with
 integral_table_m_minus_mu0 with step integral_table_dm_minus_mu.
 Note that here m and mu represent in fact m/T and mu/T so they are
 unitless.
*/
{
    Stopwatch sw;
    sw.tic();
  
    double *bulkvisCoefficients;
    if(INCLUDE_BULK_DELTAF == 1)
    {
        if(bulk_deltaf_kind == 0)
            bulkvisCoefficients = new double [3];
        else
            bulkvisCoefficients = new double [2];
    }

    // sort freeze-out temperature for furture use
    //for (long l=0; l<FO_length; l++)
    //    sorted_FZ[l] = l; // natural order
    //// now bubble sort temperature; smaller temperature goes first
    //long sort_to = FO_length-1;
    //while (sort_to>0)
    //{
    //    long last_operation = 0;
    //    for (long l2=0; l2<sort_to; l2++)
    //    {
    //        if (FOsurf_ptr[sorted_FZ[l2]].Tdec 
    //                        > FOsurf_ptr[sorted_FZ[l2+1]].Tdec)
    //        {
    //            // swap
    //            long tmp = sorted_FZ[l2];
    //            sorted_FZ[l2] = sorted_FZ[l2+1];
    //            sorted_FZ[l2+1] = tmp;
    //            last_operation = l2;
    //        }
    //    }
    //    sort_to = last_operation;
    //}

    // read parameters
    double tolerance = paraRdr->getVal("grouping_tolerance");

    // now loop over all freeze-out cells and particles
    FO_surf *surf; particle_info* particle;
    double unit_factor = 1.0/pow(hbarC,3);  // unit: convert to unitless

    double integral_laststep[number_of_chosen_particles];

    if (AMOUNT_OF_OUTPUT > 0)
        print_progressbar(-1);

    // loop over all the fluid cells
    for (long l = 0; l < FO_length; l++)
    {
        //surf = &FOsurf_ptr[sorted_FZ[l]];
        surf = &FOsurf_ptr[l];
        double temp = surf->Tdec;
        double tau = surf->tau;

        double gammaT = surf->u0;
        double ux = surf->u1;
        double uy = surf->u2;
        double uz = surf->u3;   // uz = tau*u^\eta

        double da0 = surf->da0;
        double da1 = surf->da1;
        double da2 = surf->da2;
        double da3 = surf->da3;
        
        double dsigma_dot_u = tau*(da0*gammaT + ux*da1 + uy*da2 + uz*da3/tau);

        // bulk delta f contribution
        double bulkPi = 0.0;
        if(INCLUDE_BULK_DELTAF == 1)
        {
            if(bulk_deltaf_kind == 0)
                bulkPi = surf->bulkPi;
            else
                bulkPi = surf->bulkPi/hbarC; // unit in fm^-4 
            getbulkvisCoefficients(temp, bulkvisCoefficients);
        }

        // diffusion delta f
        double dsigma_dot_q = 0.0;
        double deltaf_qmu_coeff = 1.0;
        double prefactor_qmu = 0.0;
        if(INCLUDE_DIFFUSION_DELTAF == 1)
        {
            double qmu0 = surf->qmu0;
            double qmu1 = surf->qmu1;
            double qmu2 = surf->qmu2;
            double qmu3 = surf->qmu3;
            dsigma_dot_q = tau*(da0*qmu0 + da1*qmu1 + da2*qmu2 + da3*qmu3/tau);
            
            double mu_B = surf->muB;
            deltaf_qmu_coeff = get_deltaf_qmu_coeff(temp, mu_B);
            
            double rho_B = surf->Bn;
            double Edec = surf->Edec;
            double Pdec = surf->Pdec;
            prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
        }

        // calculate dN / (dxt dy) for all particles
        int last_particle_sign=0;
        int last_particle_degen=0;
        int last_particle_baryon = 2;
        double last_particle_mass=-1;
        double last_particle_mu=-1;
        double total_N = 0;
        for (long n = 0; n < number_of_chosen_particles; n++)
        {
            long real_particle_idx = chosen_particles_sampling_table[n];
            particle = &particles[real_particle_idx];

            int sign = particle->sign;
            int degen = particle->gspin;
            double mass = particle->mass;
            int baryon = particle->baryon;
            double mu = surf->particle_mu[real_particle_idx];

            double prefactor = degen/(2.*M_PI*M_PI);

            if ( n > 0 && last_particle_sign == sign 
                 && last_particle_degen == degen 
                 && abs((last_particle_mass - mass)
                         /(last_particle_mass+1e-30)) < tolerance 
                 && abs((last_particle_mu - mu)
                         /(last_particle_mu+1e-30)) < tolerance
            )
            {
                // skip calculation for the current particle
                integral_laststep[n] = total_N;
            }
            else
            {
                // calculate dN / (dxt dy)
                last_particle_sign = sign;
                last_particle_degen = degen;
                last_particle_mass = mass;
                last_particle_mu = mu;
                last_particle_baryon = baryon;

                double* results_ptr = new double [5];
                calculate_dN_analytic(particle, mu, temp, results_ptr);

                double N_eq = (
                        unit_factor*prefactor*dsigma_dot_u*results_ptr[0]);

                double deltaN_bulk = 0.0;
                if(INCLUDE_BULK_DELTAF == 1)
                {
                    deltaN_bulk = (unit_factor*prefactor*dsigma_dot_u
                                   *(- bulkPi*bulkvisCoefficients[0])
                                   *(- bulkvisCoefficients[1]*results_ptr[1]
                                     + results_ptr[2]));
                }

                double deltaN_qmu = 0.0;
                if(INCLUDE_DIFFUSION_DELTAF == 1)
                {
                    deltaN_qmu = (unit_factor*prefactor
                                  *dsigma_dot_q/deltaf_qmu_coeff
                                  *(- prefactor_qmu*results_ptr[3] 
                                    - baryon*results_ptr[4]));
                }

                total_N = N_eq + deltaN_bulk + deltaN_qmu;

                integral_laststep[n] = total_N;

                delete [] results_ptr;
            }
            dN_dxtdy_4all[l][n] = integral_laststep[n];
        }

        if (AMOUNT_OF_OUTPUT > 0)
            print_progressbar((double)(l)/FO_length);
    }
    if (AMOUNT_OF_OUTPUT > 0)
            print_progressbar(1);

    if ((int)(paraRdr->getVal("output_dN_dxtdy_4all")))
    {
        Table to_write(dN_dxtdy_4all, FO_length, number_of_chosen_particles);
        char dN_dxt_filename_buffer[300];
        // 0 means "all"
        sprintf(dN_dxt_filename_buffer, dN_dxt_filename.c_str(), 0); 
        ofstream of(dN_dxt_filename_buffer);
        to_write.printTable(of);
        of.close();
    }

    if(INCLUDE_BULK_DELTAF == 1)
        delete [] bulkvisCoefficients;

    sw.toc();
    cout << endl << " -- Calculate_dN_dxtdy_4all_particles finished in " 
         << sw.takeTime() << " seconds." << endl;

}


//***************************************************************************
void EmissionFunctionArray::calculate_dN_analytic(
      particle_info* particle, double mu, double Temperature, double* results)
/* calculate particle yield using analytic formula as a series sum of Bessel 
   functions. The particle yield emitted from a given fluid cell is 
   \dsigma_mu u^\mu times the return value from this function
*/
{
   double N_eq = 0.0;                  // equilibrium contribution
   double deltaN_bulk_term1 = 0.0;     // contribution from bulk delta f
   double deltaN_bulk_term2 = 0.0;     // contribution from bulk delta f
   double deltaN_qmu_term1 = 0.0;      // contribution from baryon diffusion 
   double deltaN_qmu_term2 = 0.0;      // contribution from baryon diffusion 

   int truncate_order = 10;            // truncation order in taylor expansion

   int sign = particle->sign;
   double mass = particle->mass;
   double beta = 1./Temperature;

   double lambda = exp(beta*mu);  // fugacity factor

   // compute the sum in the series
   for(int n = 1; n <= truncate_order; n++)
   {
      double arg = n*mass*beta;  // argument inside bessel functions

      double theta = pow(-sign, n-1);
      double fugacity = pow(lambda, n);
      //double K_2 = gsl_sf_bessel_Kn(2, arg);
      double K_2 = get_special_function_K2(arg);

      N_eq += theta/n*fugacity*K_2;

      if(INCLUDE_BULK_DELTAF == 1 && bulk_deltaf_kind == 1)
      {
          //double K_1 = gsl_sf_bessel_K1(arg);
          double K_1 = get_special_function_K1(arg);
          deltaN_bulk_term1 += theta*fugacity*(mass*beta*K_1 + 3./n*K_2);
          deltaN_bulk_term2 += theta*fugacity*K_1;
      }

      if(INCLUDE_DIFFUSION_DELTAF == 1)
      {
          deltaN_qmu_term1 += theta/n*fugacity*K_2;

          double *sf_expint_En_ptr = new double [sf_expint_truncate_order-1];
          get_special_function_En(arg, sf_expint_En_ptr);

          double mbeta = mass*beta;
          double I_1_n = 0.0;
          double I_1_1 = exp(-arg)/arg*(2./(arg*arg) + 2./arg - 1./2.);
          //double I_1_2 = 3./8.*gsl_sf_expint_E2(arg);
          double I_1_2 = 3./8.*sf_expint_En_ptr[0];
          I_1_n = I_1_1 + I_1_2;

          double double_factorial = 1.;  // record (2k-5)!! 
          double factorial = 2.;         // record k! start with 2!
          double factor_2_to_k_power = 4.;      // record 2^k start with 2^2
          for(int k = 3; k <= sf_expint_truncate_order; k++)
          {
              double_factorial *= (2*k - 5);
              factorial *= k;
              factor_2_to_k_power *= 2;
              //double I_1_k = (
              //    3.*double_factorial/factor_2_to_k_power/factorial
              //    *gsl_sf_expint_En(2*k-2, arg));
              double I_1_k = (
                  3.*double_factorial/factor_2_to_k_power/factorial
                  *sf_expint_En_ptr[k-2]);
              I_1_n += I_1_k;
          }
          I_1_n = -(mbeta*mbeta*mbeta)*I_1_n;
          deltaN_qmu_term2 += n*theta*fugacity*I_1_n;

          delete [] sf_expint_En_ptr;
      }
   }

   // equilibrium contribution
   double prefactor_Neq = mass*mass*Temperature;
   N_eq = prefactor_Neq*N_eq;

   // contribution from bulk viscosity
   if(INCLUDE_BULK_DELTAF == 1 && bulk_deltaf_kind == 1)
   {
       deltaN_bulk_term1 = mass*mass/beta*deltaN_bulk_term1;
       deltaN_bulk_term2 = mass*mass*mass/3.*deltaN_bulk_term2;
   }
   else
   {
       deltaN_bulk_term1 = 0.0;
       deltaN_bulk_term2 = 0.0;
   }

   // contribution from baryon diffusion
   if(INCLUDE_DIFFUSION_DELTAF == 1)
   {
       deltaN_qmu_term1 = mass*mass/(beta*beta)*deltaN_qmu_term1;
       deltaN_qmu_term2 = 1./(3.*beta*beta*beta)*deltaN_qmu_term2;
   }

   // final results
   results[0] = N_eq;
   results[1] = deltaN_bulk_term1;
   results[2] = deltaN_bulk_term2;
   results[3] = deltaN_qmu_term1;
   results[4] = deltaN_qmu_term2;
}

//***************************************************************************
double EmissionFunctionArray::calculate_total_FZ_energy_flux()
/*
 Return the total energy flux on the freeze-out surface at eta_s=0.
*/
{
    Stopwatch sw;
    sw.tic();
    cout << " Function calculate_total_FZ_energy_flux started..." << endl;

    FO_surf *surf;
    double total_energy = 0;
    for (long l=0; l<FO_length; l++)
    {
        surf = &FOsurf_ptr[l];

        double p = surf->Pdec;
        double e = surf->Edec;

        double u0 = surf->u0;
        double u1 = surf->u1;
        double u2 = surf->u2;

        double da0 = surf->da0;
        double da1 = surf->da1;
        double da2 = surf->da2;

        double pi00 = surf->pi00;
        double pi01 = surf->pi01;
        double pi02 = surf->pi02;

        total_energy += ((e*u0*u0 + pi00)*da0 
                + ((e+p)*u0*u1 + pi01)*da1 + ((e+p)*u0*u2 + pi02)*da2);
    }

    return total_energy;
}




//***************************************************************************
void EmissionFunctionArray::sample_using_dN_dxtdy_4all_particles_conventional()
/*
 Sample using the already calculated dN_dxtdy_4all array.
 Different model can be used in the sampling.
*/
{
    Stopwatch sw_total;
    sw_total.tic();
    cout << " Function sample_using_dN_dxtdy_4all_particles started..." 
         << endl;

    // reusable local variables
    FO_surf *surf; particle_info* particle;
    double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;
  
    double *bulkvisCoefficients;
    if(INCLUDE_BULK_DELTAF == 1)
    {
        if(bulk_deltaf_kind == 0)
            bulkvisCoefficients = new double [3];
        else
            bulkvisCoefficients = new double [2];
    }

    // load pre-calculated table
    // (x,y) that y*exp(-y) = x; y<=1
    TableFunction z_exp_m_z("tables/z_exp_m_z.dat"); 
    z_exp_m_z.interpolation_model = 5;

    // control variables
    int sampling_model = paraRdr->getVal("dN_dy_sampling_model");
    double sampling_para1 = paraRdr->getVal("dN_dy_sampling_para1");
    double sampling_para2 = paraRdr->getVal("dN_dy_sampling_para2");
    double sampling_para3 = paraRdr->getVal("dN_dy_sampling_para3");
    double sampling_para4 = paraRdr->getVal("dN_dy_sampling_para4");
    double sampling_para5 = paraRdr->getVal("dN_dy_sampling_para5");

    long number_of_repeated_sampling = paraRdr->getVal(
                                                "number_of_repeated_sampling");
    double pT_to = paraRdr->getVal("sample_pT_up_to");
    if (pT_to<0)
        pT_to = pT_tab->getLast(1); // use table to determine pT range
    double y_minus_eta_s_range = paraRdr->getVal("sample_y_minus_eta_s_range");

    if (sampling_model<=100)
    {
        // use "conventional" sampling
        // reusable variables
        // dN_dxtdy for 1 particle
        vector<double> dN_dxtdy_single_particle(FO_length, 0); 
        cout << endl << "Sampling using dN/dy with "
            << "sample_using_dN_dxtdy_4all_particles function." 
            << endl << endl;
        for (long n=0; n<number_of_chosen_particles; n++)
        {
            long real_particle_idx = chosen_particles_sampling_table[n];
            particle = &particles[real_particle_idx];
            double mass = particle->mass;
            int sign = particle->sign;
            int degen = particle->gspin;
            int baryon = particle->baryon;
            cout << "Index: " << n << ", Name: " << particle->name 
                 << ", Monte-carlo index: " << particle->monval << endl;

            // prepare for outputs
            // the control file records how many particles are there 
            // in each sampling
            char samples_control_filename_buffer[300];
            sprintf(samples_control_filename_buffer, 
                    samples_control_filename.c_str(), particle->monval);
            remove(samples_control_filename_buffer);
            ofstream of_control(samples_control_filename_buffer);

            // the sample file contains the actual samples
            char samples_filename_buffer[300];
            sprintf(samples_filename_buffer, samples_filename.c_str(), 
                    particle->monval);
            remove(samples_filename_buffer);
            ofstream of_sample(samples_filename_buffer);

            // buffers are used to speed up the output process
            char line_buffer[500]; // only used in text mode
            stringstream sample_str_buffer; // to speed up outputing process
            stringstream control_str_buffer;

            // get y range for sampling
            double y_LB = paraRdr->getVal("y_LB");
            double y_RB = paraRdr->getVal("y_RB");

            // prepare the inverse CDF
            for (long l=0; l<FO_length; l++)
                dN_dxtdy_single_particle[l] = dN_dxtdy_4all[l][n];
            RandomVariable1DArray rand1D(&dN_dxtdy_single_particle, 0);

            // first get total number of particles
            double dN_dy = rand1D.return_sum();

            cout << " -- Sampling using dN_dy=" << dN_dy << ", ";

            double dN;
            if(hydro_mode != 2)
                dN = (y_RB-y_LB)*dN_dy;
            else  // for (3+1)-d case, dN_dy is total N (summing over all etas)
                dN = dN_dy;
            cout << "dN=" << dN << "..." << endl;

            Stopwatch sw;
            sw.tic();

            // The following variables are for "dynamic maximum" treatment
            // The maximum used in pdf sampling is "maximum_guess", this value 
            // is guaranteed (proven) to be larger than the emission function 
            // but may over estimate the actual maximum 
            // (hard to find analytically).
            // The code first use the guessed maximum (maximum_guess) for 
            // adjust_maximum_after number of sampling-tries; 
            // during which process the largest ratio between the emission 
            // function and maximum_guess is stored in maximum_ratio variable.
            // Next, if use_dynamic_maximum is set to 1, then the code adjust 
            // maximum_guess to maximum_guess*maximum_ratio*adjust_maximum_to 
            // for the rest of sampling. Note that you want to set 
            // adjust_maximum_to to be slightly larger than 1 to avoid errors.
            int use_dynamic_maximum = paraRdr->getVal("use_dynamic_maximum");
            long adjust_maximum_after = paraRdr->getVal("adjust_maximum_after");
            double adjust_maximum_to = paraRdr->getVal("adjust_maximum_to");
            long number_of_tries = 0, number_of_success = 0;
            double maximum_ratio = 0;
            bool has_adjusted_maximum = false;
            double actual_adjusted_maximum_factor = 1.0;

            long sample_writing_signal = 0;
            long control_writing_signal = 0;
            long maximum_impatience = 5000; // used in etas-pt-phi sampling

            if (AMOUNT_OF_OUTPUT>0)
                print_progressbar(-1);

            for (long repeated_sampling_idx = 1; 
                 repeated_sampling_idx <= number_of_repeated_sampling; 
                 repeated_sampling_idx++)
            {
                long number_to_sample = determine_number_to_sample(
                    dN, sampling_model, sampling_para1, sampling_para2, 
                    sampling_para3, sampling_para4, sampling_para5);

                // write to control file
                sprintf(line_buffer, "%lu\n", number_to_sample);
                control_str_buffer << line_buffer;
                control_writing_signal++;
                if (control_writing_signal == NUMBER_OF_LINES_TO_WRITE)
                {
                    of_control << control_str_buffer.str();
                    if (use_dynamic_maximum && !has_adjusted_maximum)
                    {
                        if (number_of_tries > adjust_maximum_after)
                        {
                            // adjust maximum
                            actual_adjusted_maximum_factor = (
                                            maximum_ratio*adjust_maximum_to);
                            has_adjusted_maximum = true;
                        }
                    }
                    control_str_buffer.str("");
                    control_writing_signal=0;
                }

                long sample_idx = 1;
                while (sample_idx <= number_to_sample)
                {
                    // first, sample eta and freeze-out cell index
                    //long FO_idx =  sorted_FZ[rand1D.rand()];
                    long FO_idx =  rand1D.rand();

                    surf = &FOsurf_ptr[FO_idx];

                    double Tdec = surf->Tdec;
                    double inv_Tdec = 1.0/Tdec;
                    double Pdec = surf->Pdec;
                    double Edec = surf->Edec;

                    double tau = surf->tau;
                    double u0 = surf->u0;
                    double u1 = surf->u1;
                    double u2 = surf->u2;
                    double u3 = surf->u3;

                    double mu = surf->particle_mu[real_particle_idx];

                    double da0 = surf->da0;
                    double da1 = surf->da1;
                    double da2 = surf->da2;
                    double da3 = surf->da3;

                    double eta_s;
                    if(hydro_mode == 2)
                        eta_s = surf->eta;
                    else
                        eta_s = 0.0;

                    double pi00 = surf->pi00;
                    double pi01 = surf->pi01;
                    double pi02 = surf->pi02;
                    double pi03 = surf->pi03;
                    double pi11 = surf->pi11;
                    double pi12 = surf->pi12;
                    double pi13 = surf->pi13;
                    double pi22 = surf->pi22;
                    double pi23 = surf->pi23;
                    double pi33 = surf->pi33;

                    double deltaf_prefactor = 0.0;
                    if(INCLUDE_DELTAF)
                        deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec));
                    
                    double bulkPi = 0.0;
                    if(INCLUDE_BULK_DELTAF == 1)
                    {
                        if(bulk_deltaf_kind == 0)
                            bulkPi = surf->bulkPi;
                        else
                            bulkPi = surf->bulkPi/hbarC;   // unit in fm^-4 
                        getbulkvisCoefficients(Tdec, bulkvisCoefficients);
                    }

                    // diffusion delta f
                    double qmu0 = 0.0;
                    double qmu1 = 0.0;
                    double qmu2 = 0.0;
                    double qmu3 = 0.0;
                    double deltaf_qmu_coeff = 1.0;
                    double prefactor_qmu = 0.0;
                    if(INCLUDE_DIFFUSION_DELTAF == 1)
                    {
                        qmu0 = surf->qmu0;
                        qmu1 = surf->qmu1;
                        qmu2 = surf->qmu2;
                        qmu3 = surf->qmu3;
                        
                        double mu_B = surf->muB;
                        deltaf_qmu_coeff = get_deltaf_qmu_coeff(Tdec, mu_B);
                        
                        double rho_B = surf->Bn;
                        prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
                    }

                    // calculate maximum value for p*dsigma f, 
                    // used in PDF accept/reject sampling
                    double u_dot_dsigma = tau*(
                                    u0*da0 + u1*da1 + u2*da2 + u3*da3/tau);
                    double dsigma_sq = tau*tau*(
                              da0*da0 - da1*da1 - da2*da2 - da3*da3/tau/tau);
                    double dsigmaT = sqrt(fabs(dsigma_sq 
                                               - u_dot_dsigma*u_dot_dsigma));
                    //double dsigmaT = sqrt(da0*da0*(u1*u1+u2*u2) 
                    //                      + da1*da1*(1+u1*u1) 
                    //                      + da2*da2*(1+u2*u2) 
                    //                      + 2*da0*da1*u0*u1 + 2*da0*da2*u0*u2 
                    //                      + 2*da1*da2*u1*u2);
                    
                    double dsigma_all = abs(u_dot_dsigma) + dsigmaT;

                    double f0_mass = 1./(exp((mass - mu)*inv_Tdec) + sign);

                    // ideal first, p*dsigma pT f < dsgima_all*E^2*f0
                    // Here pT is involved because instead of sampling pT^2 
                    // uniformly by using d(pT^2) uniformly, 
                    // we sample pT*d(pT) where pT is uniformly to avoid 
                    // taking sqrt, which saves time
                    // solve (1 -+ f0) = A/(beta E) 
                    // (A= power of E, for ideal case A=2), which gives
                    // (beta*E-A)*exp(beta*E-A) = +- A*exp(beta*mu-A) --- (*)
                    int A;
                    A = 1;
                    // ideal part for the guess of the maximum
                    double guess_ideal = 0; 
                    if (sign==1) // fermion
                    {
                        // choose upper sign; (*) always has a solution
                        // which gives the maximum
                        double Emax = Tdec*(
                                        get_special_function_lambertW(
                                                A*exp(inv_Tdec*mu - A)) + A);
                        if (Emax < mass)
                            Emax = mass; // maximum in [mass, inf]
                        guess_ideal = (
                            pow(Emax, A)/(exp((Emax - mu)*inv_Tdec) + sign));
                    }
                    else // boson
                    {
                        // choose lower sign; (*) has a solution only 
                        // when A*exp(beta*mu-A)<=1/e
                        double rhs = A*exp(inv_Tdec*mu - A);
                        if (rhs>0.3678794) // 1/e = 0.367879441171442
                        {
                            // no solution; maximum is attained at E=mass
                            guess_ideal = pow(mass, A)*f0_mass;
                        }
                        else
                        {
                            // has a solution; maximum is attained either at 
                            // the solution E=Emax or at E=mass. Note that 
                            // there are two Emax solutions to (*) and we want 
                            // the larger one.
                            double Emax = Tdec*(A - z_exp_m_z.map(rhs));
                            if (Emax < mass)
                            {
                                guess_ideal = pow(mass, A)*f0_mass;
                            }
                            else
                            {
                                double guess_ideal1 = (pow(Emax, A)
                                          /(exp((Emax - mu)*inv_Tdec) + sign));
                                double guess_ideal2 = pow(mass, A)*f0_mass;
                                guess_ideal = (guess_ideal1>guess_ideal2 ? 
                                                guess_ideal1 : guess_ideal2);
                            }
                        }
                    }

                    // next viscous part
                    // p*dsigma pT f < dsgima_all*tmp_factor*sqrt(3)
                    //                 *E^3*f0*trace_Pi2/(2*T^2*(e+p))
                    double trace_Pi2 = (
                        pi00*pi00 + pi11*pi11 + pi22*pi22 + pi33*pi33 
                        - 2*pi01*pi01 - 2*pi02*pi02 -2.*pi03*pi03
                        + 2*pi12*pi12 + 2.*pi13*pi13 + 2.*pi23*pi23);
                    double tmp_factor = 1;
                    if (F0_IS_NOT_SMALL && sign==-1)
                    {
                        // (1+f0) <= 2*f0
                        tmp_factor = 2.0;
                    }
                    else
                    {
                        // (1-f0) or (1-0*f0) <= 1*f0
                        tmp_factor = 1.0;
                    }
                    // viscous case, solve (1 -+ f0) = A/(beta E) for A=3. 
                    // Ref. ideal case
                    A = 3;
                    // ideal part for the guess of the maximum
                    double guess_viscous = 0; 
                    if (sign == 1) // fermion
                    {
                        // choose upper sign; (*) always has a solution 
                        // which gives the maximum
                        double Emax = (Tdec*(
                            get_special_function_lambertW(
                                    A*exp(inv_Tdec*mu - A)) + A));
                        if (Emax < mass)
                            Emax = mass;
                        guess_viscous = (pow(Emax, A)
                                    /(exp((Emax - mu)*inv_Tdec) + sign));
                    }
                    else // boson
                    {
                        double rhs = A*exp(inv_Tdec*mu - A);
                        if (rhs>0.3678794) // 1/e = 0.367879441171442
                        {
                            guess_viscous = pow(mass, A)*f0_mass;
                        }
                        else
                        {
                            double Emax = Tdec*(A - z_exp_m_z.map(rhs));
                            if (Emax < mass)
                            {
                                guess_viscous = pow(mass, A)*f0_mass;
                            }
                            else
                            {
                                double guess_viscous1 = (pow(Emax, A)
                                          /(exp((Emax - mu)*inv_Tdec) + sign));
                                double guess_viscous2 = pow(mass, A)*f0_mass;
                                guess_viscous = (
                                    guess_viscous1>guess_viscous2 ? 
                                    guess_viscous1 : guess_viscous2);
                            }
                        }
                    }
                    guess_viscous *= (tmp_factor/2.0*inv_Tdec*inv_Tdec
                                      *sqrt(trace_Pi2)/(Edec+Pdec));

                    // bulk delta f
                    double guess_bulk = 0.0;
                    if(INCLUDE_BULK_DELTAF == 1)
                    {
                        guess_bulk = (fabs(bulkPi*bulkvisCoefficients[0])
                                      *mass*mass*inv_Tdec/3.
                                      *f0_mass*(1. - sign*f0_mass));
                    }

                    // baryon diffusion delta f
                    double guess_qmu = 0.0;
                    if(INCLUDE_DIFFUSION_DELTAF == 1)
                    {
                        double qmu_sq = (qmu0*qmu0 - qmu1*qmu1 
                                        - qmu2*qmu2 - qmu3*qmu3);
                        double qmu_mag_over_kappa_hat = (
                                        sqrt(fabs(qmu_sq))/deltaf_qmu_coeff);
                        // term 1
                        A = 2;
                        // ideal part for the guess of the maximum
                        double guess_G2max = 0; 
                        if (sign == 1) // fermion
                        {
                            // choose upper sign; (*) always has a solution 
                            // which gives the maximum
                            double Emax = (Tdec*(
                                get_special_function_lambertW(
                                                A*exp(inv_Tdec*mu - A)) + A));
                            if (Emax < mass)
                                Emax = mass;
                            guess_G2max = (pow(Emax, A)
                                        /(exp((Emax - mu)*inv_Tdec) + sign));
                        }
                        else // boson
                        {
                            double rhs = A*exp(inv_Tdec*mu - A);
                            if (rhs>0.3678794) // 1/e = 0.367879441171442
                            {
                                guess_G2max = pow(mass, A)*f0_mass;
                            }
                            else
                            {
                                double Emax = Tdec*(A - z_exp_m_z.map(rhs));
                                if (Emax < mass)
                                {
                                    guess_G2max = pow(mass, A)*f0_mass;
                                }
                                else
                                {
                                    double guess_G2max_temp1 = (
                                        pow(Emax, A)/(exp((Emax - mu)*inv_Tdec) 
                                                      + sign));
                                    double guess_G2max_temp2 = (
                                                    pow(mass, A)*f0_mass);
                                    guess_G2max = (
                                        guess_G2max_temp1>guess_G2max_temp2 ? 
                                        guess_G2max_temp1 : guess_G2max_temp2);
                                }
                            }
                        }

                        double guess_G1max = guess_ideal;

                        guess_qmu = prefactor_qmu*guess_G2max;
                        if(baryon > 0)
                            guess_qmu += baryon*guess_G1max;

                        guess_qmu *= tmp_factor*qmu_mag_over_kappa_hat;
                    }

                    
                    // combine
                    double maximum_guess = (prefactor*degen*dsigma_all
                                            *(guess_ideal + guess_viscous 
                                              + guess_bulk + guess_qmu));

                    // next sample pt and phi
                    // will-be sampled values
                    double pT, mT, phi, px, py; 
                    double rapidity_y, y_minus_eta_s;
                    long tries = 1;
                    double accept_prob_max = 0.0;
                    while (tries < maximum_impatience)
                    {
                        // refer to calculate_dNArrays function to see how 
                        // the rate is calculated
                        // Basically it is "just Cooper-Frye"

                        // sample according to pT dpT
                        pT = sqrt(drand(0, pT_to*pT_to)); 
                        mT = sqrt(mass*mass + pT*pT);

                        phi = drand(0, 2*M_PI);
                        px = pT*cos(phi);
                        py = pT*sin(phi);

                        y_minus_eta_s = drand(-y_minus_eta_s_range, 
                                              y_minus_eta_s_range);

                        double p0 = mT*cosh(y_minus_eta_s);  // p0 = p^tau
                        double p3 = mT*sinh(y_minus_eta_s);  // p3 = tau p^eta

                        double pdotu = p0*u0 - px*u1 - py*u2 - p3*u3;
                        double expon = (pdotu - mu)*inv_Tdec;
                        double f0 = 1./(exp(expon)+sign);

                        double pdsigma = p0*da0 + px*da1 + py*da2 + p3*da3/tau;

                        double Wfactor = (
                            p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 
                            - 2.0*p0*p3*pi03
                            + px*px*pi11 + 2.0*px*py*pi12 + 2.0*px*p3*pi13
                            + py*py*pi22 + 2.0*py*p3*pi23
                            + p3*p3*pi33);

                        double delta_f_shear = (
                            (1. - F0_IS_NOT_SMALL*sign*f0)*Wfactor
                            *deltaf_prefactor);

                        double delta_f_bulk = 0.0;
                        if (INCLUDE_BULK_DELTAF== 1)
                        {
                            if(bulk_deltaf_kind == 0)
                            {
                                delta_f_bulk = (
                                    -(1. - F0_IS_NOT_SMALL*sign*f0)*bulkPi
                                    *(bulkvisCoefficients[0]*mass*mass 
                                      + bulkvisCoefficients[1]*pdotu 
                                      + bulkvisCoefficients[2]*pdotu*pdotu));
                            }
                            else if (bulk_deltaf_kind == 1)
                            {

                                double E_over_T = pdotu/Tdec;
                                double mass_over_T = mass/Tdec;
                                delta_f_bulk = (
                                    -1.0*(1.-sign*f0)/E_over_T
                                    *bulkvisCoefficients[0]
                                    *(mass_over_T*mass_over_T/3. 
                                      - bulkvisCoefficients[1]
                                        *E_over_T*E_over_T)
                                    *bulkPi);
                            }
                            else if (bulk_deltaf_kind == 2)
                            {
                                double E_over_T = pdotu/Tdec;
                                delta_f_bulk = (
                                    -1.*(1.-sign*f0)*(-bulkvisCoefficients[0] 
                                        + bulkvisCoefficients[1]*E_over_T)
                                    *bulkPi);
                            }
                            else if (bulk_deltaf_kind == 3)
                            {
                                double E_over_T = pdotu/Tdec;
                                delta_f_bulk = (
                                    -1.0*(1.-sign*f0)/sqrt(E_over_T)
                                    *(-bulkvisCoefficients[0] 
                                      + bulkvisCoefficients[1]*E_over_T)
                                    *bulkPi);
                            }
                            else if (bulk_deltaf_kind == 4)
                            {
                                double E_over_T = pdotu/Tdec;
                                delta_f_bulk = (
                                    -1.0*(1.-sign*f0)*(bulkvisCoefficients[0] 
                                        - bulkvisCoefficients[1]/E_over_T)
                                    *bulkPi);
                            }
                        }

                        // delta f for diffusion
                        double delta_f_qmu = 0.0;
                        if(INCLUDE_DIFFUSION_DELTAF == 1)
                        {
                            double qmufactor = (
                                    p0*qmu0 - px*qmu1 - py*qmu2 - p3*qmu3);
                            delta_f_qmu = ((1. - sign*f0)
                                           *(prefactor_qmu - baryon/pdotu)
                                           *qmufactor/deltaf_qmu_coeff);
                        }

                        double result;
                        // restrict the size of delta f to be smaller than f_0
                        double ratio_max = 1.0 - 1e-6;
                        double deltaf_size = fabs(delta_f_shear + delta_f_bulk 
                                                  + delta_f_qmu);
                        double resize_factor = (
                            min(1., ratio_max/(deltaf_size + 1e-10)));
                        result = (prefactor*degen*f0*pdsigma*tau
                                  *(1. + (delta_f_shear + delta_f_bulk 
                                          + delta_f_qmu)*resize_factor));

                        double accept_prob = (
                            result/(actual_adjusted_maximum_factor
                                    *maximum_guess));
                        if (AMOUNT_OF_OUTPUT > 5)
                        {
                            if(accept_prob > accept_prob_max)
                                accept_prob_max = accept_prob;
                        }

                        // to track success rate
                        number_of_tries ++;
                        if (accept_prob > maximum_ratio)
                                maximum_ratio = accept_prob;

                        // for debugging
                        if (accept_prob > 1 && AMOUNT_OF_OUTPUT > 1)
                            cout << "EmissionFunctionArray::"
                                 << "sample_using_dN_dxtdy_4all_particles "
                                 << "warning: emission function is bigger "
                                 << "than 1: " << accept_prob << endl;

                        if (drand48() < accept_prob)
                            break; // accept sample!

                        tries++;
                    }
                    if (tries == maximum_impatience && AMOUNT_OF_OUTPUT > 5)
                    {
                        cout << "EmissionFunctionArray::"
                             << "sample_using_dN_dxtdy_4all_particles warning:"
                             << "maximum_impatience reached." << endl;
                        cout << "accept_prob_max = " << accept_prob_max << ", "
                             << "maximum_guess = " << maximum_guess << endl;
                        cout << "guess_ideal = " << guess_ideal << ", "
                             << "guess_shear = " << guess_viscous << ", "
                             << "guess_bulk = " << guess_bulk << endl;
                    }

                    if (tries == maximum_impatience)
                        continue; // resample
                    else 
                        sample_idx++; // write-out sample

                    number_of_success++; // to track success rate

                    if(hydro_mode != 2)
                    {
                        rapidity_y = y_LB + drand48()*(y_RB-y_LB);
                        eta_s = rapidity_y - y_minus_eta_s;
                    }
                    else
                    {
                        rapidity_y = y_minus_eta_s + eta_s;
                    }

                    double p_z = mT*sinh(rapidity_y);
                    double E = mT*cosh(rapidity_y);
                    double z = surf->tau*sinh(eta_s);
                    double t = surf->tau*cosh(eta_s);

                    // write to sample file
                    if (!USE_OSCAR_FORMAT)
                    {
                        sprintf(line_buffer, 
                                "%lu  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
                                FO_idx, surf->tau, surf->xpt, surf->ypt, 
                                y_minus_eta_s, pT, phi, surf->da0, surf->da1, 
                                surf->da2, surf->u1/surf->u0, 
                                surf->u2/surf->u0, rapidity_y, eta_s, E, p_z, 
                                t, z);
                    }
                    // To be combined to OSCAR
                    if (USE_OSCAR_FORMAT)
                    {
                        sprintf(line_buffer, 
                                "%24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e\n", 
                                px, py, p_z, E, mass, surf->xpt, 
                                surf->ypt, z, t);
                    }
                    sample_str_buffer << line_buffer;
                    sample_writing_signal++;
                    if (sample_writing_signal==NUMBER_OF_LINES_TO_WRITE)
                    {
                        of_sample << sample_str_buffer.str();
                        sample_str_buffer.str("");
                        sample_writing_signal=0;
                    }

                }

                // check whether to adjust maximum
                if (use_dynamic_maximum && !has_adjusted_maximum)
                {
                    if (number_of_tries > adjust_maximum_after)
                    {
                        // adjust maximum
                        actual_adjusted_maximum_factor = (
                                        maximum_ratio*adjust_maximum_to);
                        has_adjusted_maximum = true;
                    }
                }

                if (AMOUNT_OF_OUTPUT>0)
                    print_progressbar((double)(repeated_sampling_idx)
                                              /number_of_repeated_sampling);
            }
            // flushing buffers
            if (control_writing_signal!=0)
            {
                of_control << control_str_buffer.str();
                control_str_buffer.str("");
                control_writing_signal = 0;
            }
            of_control.close();
            if (sample_writing_signal!=0)
            {
                of_sample << sample_str_buffer.str();
                sample_str_buffer.str("");
                sample_writing_signal = 0;
            }
            of_sample.close();
            if (AMOUNT_OF_OUTPUT>0) print_progressbar(1);

            if (AMOUNT_OF_OUTPUT>3)
            {
                cout << endl << " -- -- Number of tries: " << number_of_tries 
                     << ", number of success: " << number_of_success << endl
                     << " -- -- Success rate: " 
                     << (double) number_of_success/number_of_tries << ", " 
                     << "maximum accept rate: " << maximum_ratio << endl;
            }

            sw.toc();
            cout << endl << "Sampling finished in " 
                 << sw.takeTime() << " seconds." << endl;

        } // n; particle loop

        ofstream of_sample_format(samples_format_filename.c_str());
        // Be careful that ParameterReader class convert all strings to lower 
        // case so do NOT use variable names that differ only by cases!
        if (!USE_OSCAR_FORMAT)
        {
            of_sample_format
            << "Total_number_of_columns = " << 18 << endl
            << "FZ_cell_idx = " << 1 << endl
            << "tau = " << 2 << endl
            << "FZ_x = " << 3 << endl
            << "FZ_y = " << 4 << endl
            << "y_minus_eta_s = " << 5 << endl
            << "pT = " << 6 << endl
            << "phi = " << 7 << endl
            << "surf_da0 = " << 8 << endl
            << "surf_da1 = " << 9 << endl
            << "surf_da2 = " << 10 << endl
            << "surf_vx = " << 11 << endl
            << "surf_vy = " << 12 << endl
            << "y = " << 13 << endl
            << "eta_s = " << 14 << endl
            << "E = " << 15 << endl
            << "p_z = " << 16 << endl
            << "t = " << 17 << endl
            << "z = " << 18 << endl;
        }
        if (USE_OSCAR_FORMAT)
        {
            of_sample_format
            << "Total_number_of_columns = " << 9 << endl
            << "t = " << 9 << endl
            << "FZ_x = " << 6 << endl
            << "FZ_y = " << 7 << endl
            << "z = " << 8 << endl
            << "E = " << 4 << endl
            << "px = " << 1 << endl
            << "py = " << 2 << endl
            << "p_z = " << 3 << endl
            << "mass =" << 5 << endl;
        }
        of_sample_format.close();

    }
    else
    {
        cout << "EmissionFunctionArray::"
             << "sample_using_dN_dxtdy_4all_particles error: sampling model " 
             << sampling_model << " is not supported." << endl;
        exit(-1);
    }
    
    if(INCLUDE_BULK_DELTAF == 1)
        delete [] bulkvisCoefficients;

    sw_total.toc();
    cout << endl 
         << "sample_using_dN_dxtdy_4all_particles finished in " 
         << sw_total.takeTime() << " seconds." << endl;
}

void EmissionFunctionArray::getbulkvisCoefficients(
                double Tdec, double* bulkvisCoefficients)
{
   double Tdec_fm = Tdec/hbarC;  // [1/fm]
   double Tdec_fm_power[11];    // cache the polynomial power of Tdec_fm
   Tdec_fm_power[1] = Tdec_fm;
   for(int ipower = 2; ipower < 11; ipower++)
       Tdec_fm_power[ipower] = Tdec_fm_power[ipower-1]*Tdec_fm;
   if(bulk_deltaf_kind == 0)       // 14 moment expansion
   {
        // load from file
        
        //B0 [fm^3/GeV^3]
        bulkvisCoefficients[0] = (
                        bulkdf_coeff->interp(1, 2, Tdec_fm, 5)/pow(hbarC, 3));

        // D0 [fm^3/GeV^2]
        bulkvisCoefficients[1] = (
                        bulkdf_coeff->interp(1, 3, Tdec_fm, 5)/pow(hbarC, 2));

        // E0 [fm^3/GeV^3]
        bulkvisCoefficients[2] = (
                        bulkdf_coeff->interp(1, 4, Tdec_fm, 5)/pow(hbarC, 3));

        // parameterization for mu = 0
        // B0[fm^3/GeV^3]
        //bulkvisCoefficients[0] = (
        //              exp(-15.04512474*Tdec_fm + 11.76194266)/pow(hbarC, 3)); 
        // D0 [fm^3/GeV^2]
        //bulkvisCoefficients[1] = (
        //              exp( -12.45699277*Tdec_fm + 11.4949293)/hbarC/hbarC);  
        // E0 [fm^3/GeV^3]
        //bulkvisCoefficients[2] = (
        //            -exp(-14.45087586*Tdec_fm + 11.62716548)/pow(hbarC, 3));
   }
   else if(bulk_deltaf_kind == 1)  // relaxation type
   {
       // parameterization from JF
       // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV, 
       // do not trust it beyond
       bulkvisCoefficients[0] = (  642096.624265727 
                                 - 8163329.49562861*Tdec_fm_power[1] 
                                 + 47162768.4292073*Tdec_fm_power[2] 
                                 - 162590040.002683*Tdec_fm_power[3] 
                                 + 369637951.096896*Tdec_fm_power[4] 
                                 - 578181331.809836*Tdec_fm_power[5] 
                                 + 629434830.225675*Tdec_fm_power[6] 
                                 - 470493661.096657*Tdec_fm_power[7] 
                                 + 230936465.421*Tdec_fm_power[8] 
                                 - 67175218.4629078*Tdec_fm_power[9] 
                                 + 8789472.32652964*Tdec_fm_power[10]);

       bulkvisCoefficients[1] = (  1.18171174036192 
                                 - 17.6740645873717*Tdec_fm_power[1]
                                 + 136.298469057177*Tdec_fm_power[2] 
                                 - 635.999435106846*Tdec_fm_power[3] 
                                 + 1918.77100633321*Tdec_fm_power[4] 
                                 - 3836.32258307711*Tdec_fm_power[5] 
                                 + 5136.35746882372*Tdec_fm_power[6] 
                                 - 4566.22991441914*Tdec_fm_power[7] 
                                 + 2593.45375240886*Tdec_fm_power[8] 
                                 - 853.908199724349*Tdec_fm_power[9]
                                 + 124.260460450113*Tdec_fm_power[10]);
   }
   else if (bulk_deltaf_kind == 2)
   {
       // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV
       // do not trust it beyond
       bulkvisCoefficients[0] = (  21091365.1182649 
                                 - 290482229.281782*Tdec_fm_power[1] 
                                 + 1800423055.01882*Tdec_fm_power[2] 
                                 - 6608608560.99887*Tdec_fm_power[3] 
                                 + 15900800422.7138*Tdec_fm_power[4] 
                                 - 26194517161.8205*Tdec_fm_power[5] 
                                 + 29912485360.2916*Tdec_fm_power[6] 
                                 - 23375101221.2855*Tdec_fm_power[7] 
                                 + 11960898238.0134*Tdec_fm_power[8] 
                                 - 3618358144.18576*Tdec_fm_power[9] 
                                 + 491369134.205902*Tdec_fm_power[10]);

       bulkvisCoefficients[1] = (  4007863.29316896 
                                 - 55199395.3534188*Tdec_fm_power[1] 
                                 + 342115196.396492*Tdec_fm_power[2]
                                 - 1255681487.77798*Tdec_fm_power[3] 
                                 + 3021026280.08401*Tdec_fm_power[4] 
                                 - 4976331606.85766*Tdec_fm_power[5] 
                                 + 5682163732.74188*Tdec_fm_power[6] 
                                 - 4439937810.57449*Tdec_fm_power[7] 
                                 + 2271692965.05568*Tdec_fm_power[8] 
                                 - 687164038.128814*Tdec_fm_power[9] 
                                 + 93308348.3137008*Tdec_fm_power[10]);
   }
   else if (bulk_deltaf_kind == 3)
   {
       bulkvisCoefficients[0] = (  160421664.93603
                                 - 2212807124.97991*Tdec_fm_power[1] 
                                 + 13707913981.1425*Tdec_fm_power[2]
                                 - 50204536518.1767*Tdec_fm_power[3] 
                                 + 120354649094.362*Tdec_fm_power[4]
                                 - 197298426823.223*Tdec_fm_power[5] 
                                 + 223953760788.288*Tdec_fm_power[6]
                                 - 173790947240.829*Tdec_fm_power[7] 
                                 + 88231322888.0423*Tdec_fm_power[8] 
                                 - 26461154892.6963*Tdec_fm_power[9] 
                                 + 3559805050.19592*Tdec_fm_power[10]);
       bulkvisCoefficients[1] = (  33369186.2536556 
                                 - 460293490.420478*Tdec_fm_power[1] 
                                 + 2851449676.09981*Tdec_fm_power[2] 
                                 - 10443297927.601*Tdec_fm_power[3] 
                                 + 25035517099.7809*Tdec_fm_power[4] 
                                 - 41040777943.4963*Tdec_fm_power[5] 
                                 + 46585225878.8723*Tdec_fm_power[6] 
                                 - 36150531001.3718*Tdec_fm_power[7] 
                                 + 18353035766.9323*Tdec_fm_power[8] 
                                 - 5504165325.05431*Tdec_fm_power[9] 
                                 + 740468257.784873*Tdec_fm_power[10]);
   }
   else if (bulk_deltaf_kind == 4)
   {
       bulkvisCoefficients[0] = (  1167272041.90731 
                                 - 16378866444.6842*Tdec_fm_power[1] 
                                 + 103037615761.617*Tdec_fm_power[2] 
                                 - 382670727905.111*Tdec_fm_power[3] 
                                 + 929111866739.436*Tdec_fm_power[4] 
                                 - 1540948583116.54*Tdec_fm_power[5] 
                                 + 1767975890298.1*Tdec_fm_power[6] 
                                 - 1385606389545*Tdec_fm_power[7] 
                                 + 709922576963.213*Tdec_fm_power[8] 
                                 - 214726945096.326*Tdec_fm_power[9] 
                                 + 29116298091.9219*Tdec_fm_power[10]);
       bulkvisCoefficients[1] = (  5103633637.7213 
                                 - 71612903872.8163*Tdec_fm_power[1] 
                                 + 450509014334.964*Tdec_fm_power[2] 
                                 - 1673143669281.46*Tdec_fm_power[3] 
                                 + 4062340452589.89*Tdec_fm_power[4] 
                                 - 6737468792456.4*Tdec_fm_power[5] 
                                 + 7730102407679.65*Tdec_fm_power[6] 
                                 - 6058276038129.83*Tdec_fm_power[7] 
                                 + 3103990764357.81*Tdec_fm_power[8] 
                                 - 938850005883.612*Tdec_fm_power[9] 
                                 + 127305171097.249*Tdec_fm_power[10]);
   }
   return;
}

void EmissionFunctionArray::load_deltaf_qmu_coeff_table(string filename)
{
    ifstream table(filename.c_str());                                         
    deltaf_qmu_coeff_table_length_T = 150;                                    
    deltaf_qmu_coeff_table_length_mu = 100;                                   
    delta_qmu_coeff_table_T0 = 0.05;                                          
    delta_qmu_coeff_table_mu0 = 0.0;                                          
    delta_qmu_coeff_table_dT = 0.001;                                         
    delta_qmu_coeff_table_dmu = 0.007892;                                     
    
    deltaf_qmu_coeff_tb = new double* [deltaf_qmu_coeff_table_length_T];      
    for(int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)                  
        deltaf_qmu_coeff_tb[i] = new double [deltaf_qmu_coeff_table_length_mu];

    // load 2D table
    double dummy;                                                             
    for(int j = 0; j < deltaf_qmu_coeff_table_length_mu; j++)
        for(int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)               
            table >> dummy >> dummy >> deltaf_qmu_coeff_tb[i][j];               
    table.close();                                                            
}                                                                             

double EmissionFunctionArray::get_deltaf_qmu_coeff(double T, double muB)
{
    int idx_T = (int)((T - delta_qmu_coeff_table_T0)/delta_qmu_coeff_table_dT);
    int idx_mu = (int)((muB - delta_qmu_coeff_table_mu0)
                       /delta_qmu_coeff_table_dmu);
    double x_fraction = ((T - delta_qmu_coeff_table_T0)
                         /delta_qmu_coeff_table_dT - idx_T);
    double y_fraction = ((muB - delta_qmu_coeff_table_mu0)
                         /delta_qmu_coeff_table_dmu - idx_mu); 
    
    //avoid overflow
    if(idx_mu > deltaf_qmu_coeff_table_length_mu - 2)
        return(1e30);
    if(idx_T > deltaf_qmu_coeff_table_length_T - 2)
        return(1e30);
    
    // avoid underflow
    if(idx_mu < 0)
        return(1e30);
    if(idx_T < 0)
        return(1e30);
    
    double f1 = deltaf_qmu_coeff_tb[idx_T][idx_mu];
    double f2 = deltaf_qmu_coeff_tb[idx_T][idx_mu+1];
    double f3 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu+1];
    double f4 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu];
    double coeff = (  f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction 
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));
    return(coeff*hbarC);
}

void EmissionFunctionArray::initialize_special_function_arrays()
{
    cout << "Initializing special function arrays ... ";
    sf_expint_truncate_order = 10;
    sf_x_min = 0.5;
    sf_x_max = 400;
    sf_dx = 0.05;
    sf_tb_length = (int)((sf_x_max - sf_x_min)/sf_dx) + 1;
    sf_bessel_Kn = new double* [sf_tb_length];
    if(INCLUDE_DIFFUSION_DELTAF == 1)
        sf_expint_En = new double* [sf_tb_length];
    for(int i = 0; i < sf_tb_length; i++)
    {
        sf_bessel_Kn[i] = new double [2];

        double sf_x = sf_x_min + i*sf_dx;

        if(INCLUDE_BULK_DELTAF == 1)
            sf_bessel_Kn[i][0] = gsl_sf_bessel_K1(sf_x);     // store K_1
        else
            sf_bessel_Kn[i][0] = 0.0;

        sf_bessel_Kn[i][1] = gsl_sf_bessel_Kn(2, sf_x);  // store K_2

        if(INCLUDE_DIFFUSION_DELTAF == 1)
        {
            sf_expint_En[i] = new double [sf_expint_truncate_order-1];
            sf_expint_En[i][0] = gsl_sf_expint_E2(sf_x);     // store E_2
            for(int k = 1; k < sf_expint_truncate_order-1; k++)
            {
                sf_expint_En[i][k] = gsl_sf_expint_En(2*k+2, sf_x);
            }
        }
    }

    // pre-tabulate the lambert W function
    lambert_x_min = 0;
    lambert_x_max = 200.0;
    lambert_dx = 0.005;
    lambert_tb_length = (int)((lambert_x_max - lambert_x_min)/lambert_dx) + 1;
    lambert_W = new double [lambert_tb_length];
    //ofstream check("lambertw_function.dat");
    for(int i = 0; i < lambert_tb_length; i++)
    {
        double lambert_x_local = lambert_x_min + i*lambert_dx;
        lambert_W[i] = gsl_sf_lambert_W0(lambert_x_local);
        //check << scientific << setw(18) << setprecision(8)
        //      << lambert_x_local << "   " << lambert_W[i] << endl;
    }
    //check.close();
    cout << "done!" << endl;
}

double EmissionFunctionArray::get_special_function_K2(double arg)
{
    double results;
    if(arg < sf_x_min || arg > sf_x_max-sf_dx)
    {
        if(AMOUNT_OF_OUTPUT > 5)
        {
            cout << "EmissionFunctionArray::get_special_function_K2: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results = gsl_sf_bessel_Kn(2, arg);
    }
    else
    {
        int idx = (int)((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        results = ((1. - fraction)*sf_bessel_Kn[idx][1] 
                    + fraction*sf_bessel_Kn[idx+1][1]);
    }
    return(results);
}

double EmissionFunctionArray::get_special_function_K1(double arg)
{
    double results;
    if(arg < sf_x_min || arg > sf_x_max-sf_dx)
    {
        if(AMOUNT_OF_OUTPUT > 5)
        {
            cout << "EmissionFunctionArray::get_special_function_K1: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results = gsl_sf_bessel_K1(arg);
    }
    else
    {
        int idx = (int)((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        results = ((1. - fraction)*sf_bessel_Kn[idx][0] 
                    + fraction*sf_bessel_Kn[idx+1][0]);
    }
    return(results);
}

void EmissionFunctionArray::get_special_function_En(double arg, 
                                                    double* results)
{
    if(arg < sf_x_min || arg > sf_x_max-sf_dx)
    {
        if(AMOUNT_OF_OUTPUT > 5)
        {
            cout << "EmissionFunctionArray::get_special_function_En: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results[0] = gsl_sf_expint_E2(arg);
        for(int i = 1; i < sf_expint_truncate_order-1; i++)
        {
            results[i] = gsl_sf_expint_En(2*i+2, arg);
        }
    }
    else
    {
        int idx = (int)((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        for(int i = 0; i < sf_expint_truncate_order-1; i++)
        {
            results[i] = ((1. - fraction)*sf_expint_En[idx][i] 
                          + fraction*sf_expint_En[idx+1][i]);
        }
    }
}

double EmissionFunctionArray::get_special_function_lambertW(double arg)
{
    double results;
    if(arg < lambert_x_min || arg > lambert_x_max-lambert_dx)
    {
        if(AMOUNT_OF_OUTPUT > 5)
        {
            cout << "EmissionFunctionArray::get_special_function_lambertW: "
                 << "out of the table bound!" << endl;
            cout << "lambert_x_min = " << lambert_x_min 
                 << ", lambert_x_max = " << lambert_x_max
                 << ", arg = " << arg << endl;
        }
        results = gsl_sf_lambert_W0(arg);
    }
    else
    {
        int idx = (int)((arg - lambert_x_min)/lambert_dx);
        double fraction = (arg - lambert_x_min - idx*lambert_dx)/lambert_dx;
        results = (1. - fraction)*lambert_W[idx] + fraction*lambert_W[idx+1];
    }
    return(results);
}

