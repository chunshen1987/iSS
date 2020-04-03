// Copyright @ 2012 Zhi Qiu and Chun Shen
// Ver 2.1
// Note the version used in iSS is slightly slower than the one
// used in iS, do not mix them unless necessary.

#ifndef SRC_EMISSIONFUNCTION_H_
#define SRC_EMISSIONFUNCTION_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
#include <vector>
#include <array>
#include <memory>

#include "Table.h"
#include "TableFunction.h"
#include "ParameterReader.h"
#include "particle_decay.h"
#include "data_struct.h"
#include "data_struct.h"
#include "Random.h"
#include "pretty_ostream.h"

class EmissionFunctionArray {
 private:
    int hydro_mode;   // switch for (2+1)-d or (3+1)-d hypersurface
    int flag_PCE_;
    int flag_restrict_deltaf;
    double deltaf_max_ratio;

    pretty_ostream messager;

    const std::string path_;
    const std::string table_path_;
    const AfterburnerType afterburner_type_;

    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

    int MC_sampling;
    int number_of_repeated_sampling;
    int local_charge_conservation;

    Table *pT_tab, *phi_tab, *y_minus_eta_tab;
    int pT_tab_length, phi_tab_length;

    // the y_minus_eta_min_index holds the index 
    // to the smallest positive y-eta_s value
    int y_minus_eta_tab_length, y_minus_eta_min_index;
    // true if y_minus_eta_tab has only positive part
    bool positive_y_minus_eta_table_only;

    ParameterReader *paraRdr; // used to pass-in parameters
    int USE_OSCAR_FORMAT;
    int USE_GZIP_FORMAT;
    int INCLUDE_DELTAF, INCLUDE_BULK_DELTAF, INCLUDE_DIFFUSION_DELTAF;
    int bulk_deltaf_kind;

    int turn_on_rhob;

    Table *dN_pTdpTdphidy;  // dN / (pt dpt dphi dy)
    // store the largest element when summing over xt and eta
    // to get dN / (pt dpt dphi dy); used during sampling.
    Table *dN_pTdpTdphidy_max;

    double **dN_dxtdetady;  // dN / (d^2x_t deta dy) (correction: and dtau)
    // dN / (d^2x_t deta dy) is the (weighted) sum of all elements of
    // the [pT_tab_length][phi_tab_length]-sized dN/dall matrix;
    // this array records the maximum value of this dN/dall matrix,
    // for each eta and FO-cell.
    double **dN_dxtdetady_pT_max;
    double **dN_dxtdy_4all;  // dN / (dxt dy) for all particles

    // dN/(dxt dy) for one particle species
    std::vector<double> dN_dxtdy_for_one_particle_species;

    int number_of_chosen_particles;
    // used in spectra and flow calculations; 
    // it has length Nparticle, 0 means miss, 1 means include
    std::vector<int> chosen_particles_01_table;

    // 0/1: Particles with similar mass and chemical potentials 
    // will be sampled using the same dN/(dxt deta dy) matrix
    int grouping_particles;

    // Usable only when grouping_particles is 1.
    // If two particles have mass and chemical potentials
    // close within this relative tolerance, they are considered to be
    // identical and will be sampled successively without regenerating
    // the dN / (dxt deta dy) matrix for efficiency.
    double mc_grouping_tolerance;

    // store particle index; 
    // the sampling process follows the order specified by this table
    std::vector<int> chosen_particles_sampling_table;

    // store chosen particle monte carlo number that are not found in pdg.dat
    std::vector<int> unidentifiedPid_table;

    // list for information for all particles
    int Nparticles;
    std::vector<particle_info> particles;

    // list for information for all fluid cells
    long FO_length;
    const std::vector<FO_surf> &FOsurf_ptr;

    // store the last particle index being used by calculate_dNArrays function
    int last_particle_idx;

    double **trig_phi_table, **hypertrig_y_minus_eta_table;
    inline long determine_number_to_sample(
                    double dN, int model=1, double para1=0);

    const gsl_rng_type *gsl_type_random_number;
    gsl_rng *gsl_random_r;

    bool particles_are_the_same(int, int);
    //long *sorted_FZ;

    //array for bulk delta f coefficients
    Table *bulkdf_coeff;

    // table parameter for diffusion deltaf coefficient
    int deltaf_qmu_coeff_table_length_T;
    int deltaf_qmu_coeff_table_length_mu;
    double delta_qmu_coeff_table_T0, delta_qmu_coeff_table_mu0;
    double delta_qmu_coeff_table_dT, delta_qmu_coeff_table_dmu;
    double **deltaf_qmu_coeff_tb;       // kappa_B [1/fm^3]

    // arrays to speed up computing particle yield
    double sf_dx, sf_x_min, sf_x_max;
    int sf_tb_length;
    int sf_expint_truncate_order;
    double** sf_bessel_Kn;
    double** sf_expint_En;

    double lambert_x_min, lambert_x_max, lambert_dx;
    std::vector<double> lambert_W;

    int flag_output_samples_into_files;
    int flag_store_samples_in_memory;

    //! particle decay
    int flag_perform_decays;
    particle_decay *decayer_ptr;

 public:
    EmissionFunctionArray(std::shared_ptr<RandomUtil::Random> ran_gen,
                          Table* chosen_particle, Table* pt_tab_in,
                          Table* phi_tab_in, Table* eta_tab_in,
                          std::vector<particle_info> particles_in,
                          const std::vector<FO_surf> &FOsurf_ptr_in,
                          int flag_PCE_in, ParameterReader* paraRdr_in,
                          std::string path_in, std::string table_path,
                          AfterburnerType afterburner_type);
    ~EmissionFunctionArray();

    std::vector< std::vector<iSS_Hadron>* >* Hadron_list;

    void initialize_special_function_arrays();
    double get_special_function_K1(double arg);
    double get_special_function_K2(double arg);
    void get_special_function_En(double arg, std::vector<double> &results);
    double get_special_function_lambertW(double arg);

    void calculate_dNArrays(int);
    void calculate_dN_dxtdetady(int);
    void calculate_dN_pTdpTdphidy(int);
    void write_dN_pTdpTdphidy_toFile();
    std::string dN_pTdpTdphidy_filename;  // where to save
    void write_dN_dxtdetady_toFile();
    std::string dN_dxtdetady_filename;

    void calculate_flows(int to_order, std::string flow_diff_filename,
                         std::string flow_inte_filename);
    std::string flow_differential_filename_old, flow_integrated_filename_old;
    std::string flow_differential_filename, flow_integrated_filename;

    void calculate_dN_pTdpTdphidy_and_flows_4all_old_output(
                                                    int perform_sampling = 0);
    void calculate_dN_pTdpTdphidy_and_flows_4all(int perform_sampling = 0);
    void calculate_dN_dxtdetady_and_sample_4all();

    void shell();  // it all starts here...

    void combine_samples_to_OSCAR();
    std::string OSCAR_header_filename, OSCAR_output_filename;
    void combine_samples_to_gzip_file();

    // Sample files
    // where samples, its control informations, and its "format file" 
    // are stored; 
    // the first two can contain a "%d" string to generate multiple files
    std::string samples_format_filename;

    // First sampling method
    void sample_using_dN_dxtdetady_smooth_pT_phi();
    void calculate_dN_dtau_using_dN_dxtdetady(
                    double tau0 = 0, double dtau = 0.5, double tau_max = 17);
    void calculate_dN_dphi_using_dN_pTdpTdphidy();
    void calculate_dN_deta_using_dN_dxtdetady();
    void calculate_dN_dxt_using_dN_dxtdetady();
    void calculate_dN_dx_using_dN_dxtdetady(
                    double x_min, double x_max, double dx);

    // Second sampling method
    void calculate_dN_analytic(const particle_info* particle, double mu,
                               double Temperature,
                               std::array<double, 5> &results);

    // the following variables need to be set first in order to
    // call this function
    void calculate_dN_dxtdy_4all_particles();

    void calculate_dN_dxtdy_for_one_particle_species(const int particle_idx);

    double calculate_total_FZ_energy_flux();
    // to be used after calculate_dN_dxtdy_4all_particles
    void sample_using_dN_dxtdy_4all_particles_conventional();

    // Third sampling method
    void sample_using_dN_pTdpTdphidy();
    Table pT_tab4Sampling, phi_tab4Sampling;
    int pT_tab4Sampling_length, phi_tab4Sampling_length;
    double** trig_phi_tab4Sampling;

    void getbulkvisCoefficients(double Tdec,
                                std::array<double, 3> &bulkvisCoefficients);
    void load_deltaf_qmu_coeff_table(std::string filename);
    double get_deltaf_qmu_coeff(double T, double muB);

    void check_samples_in_memory();
    int get_number_of_sampled_events() {return(Hadron_list->size());};
    int get_number_of_particles(int iev) {
        return((*Hadron_list)[iev]->size());
    };

    iSS_Hadron get_hadron(int iev, int ipart) {
        return((*(*Hadron_list)[iev])[ipart]);
    };

    std::vector<iSS_Hadron>* get_hadron_list_iev(const int iev) {
        return(Hadron_list->at(iev));
    }

    void perform_resonance_feed_down(
                std::vector< std::vector<iSS_Hadron>* >* input_particle_list);
    int compute_number_of_sampling_needed(int number_of_particles_needed);
    double estimate_ideal_maximum(
        int sign, double mass, double Tdec, double mu, double f0_mass,
        TableFunction &z_exp_m_z);
    double estimate_shear_viscous_maximum(
        int sign, double mass, double Tdec, double mu, double f0_mass,
        TableFunction &z_exp_m_z, double pi_size);
    double estimate_diffusion_maximum(
        int sign, int baryon, double mass, double Tdec, double mu,
        double f0_mass, TableFunction &z_exp_m_z,
        double prefactor_qmu, double guess_ideal, double q_size);
    double get_deltaf_bulk(
        double mass, double pdotu, double bulkPi, double Tdec, int sign,
        double f0, const std::array<double, 3> bulkvisCoefficients);
    int sample_momemtum_from_a_fluid_cell(
        const double mass, const double degen, const int sign,
        const int baryon, const int strange, const int charge,
        const double pT_to, const double y_minus_eta_s_range,
        const double maximum_guess, const FO_surf *surf,
        const std::array<double, 3> bulkvisCoefficients,
        const double deltaf_qmu_coeff,
        double &pT, double &phi, double &y_minus_eta_s);
    double estimate_maximum(
        const FO_surf *surf, const int real_particle_idx, const double mass,
        const double sign, const double degen,
        const int baryon, const int strange, const int charge,
        TableFunction &z_exp_m_z,
        const std::array<double, 3> bulkvisCoefficients,
        const double deltaf_qmu_coeff);
    std::string add_one_sampled_particle(
        const int repeated_sampling_idx, 
        const unsigned long FO_idx, const FO_surf *surf,
        const int particle_monval, const double mass,
        const double pT, const double phi,
        const double y_minus_eta_s, const double eta_s);
};

#endif  // SRC_EMISSIONFUNCTION_H_

/*----------------------------------------------------------------------
Change log: See changelog.txt.
----------------------------------------------------------------------*/
