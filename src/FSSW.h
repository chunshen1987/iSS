// Copyright @ 2020 Chun Shen
// The Fastest Spectra Sampler in the West

#ifndef SRC_FSSW_H_
#define SRC_FSSW_H_

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
#include "Random.h"
#include "pretty_ostream.h"
#include "MomentumSamplerShell.h"

class FSSW {
 private:
    int hydro_mode;   // switch for (2+1)-d or (3+1)-d hypersurface
    int flag_PCE_;

    pretty_ostream messager_;

    const std::string path_;
    const std::string table_path_;
    const AfterburnerType afterburner_type_;

    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;
    std::shared_ptr<MomentumSamplerShell> momentum_sampler_ptr_;

    int number_of_repeated_sampling;
    int local_charge_conservation;

    ParameterReader *paraRdr; // used to pass-in parameters
    int USE_OSCAR_FORMAT;
    int USE_GZIP_FORMAT;
    int USE_BINARY_FORMAT;
    int INCLUDE_DELTAF, INCLUDE_BULK_DELTAF, INCLUDE_DIFFUSION_DELTAF;
    int bulk_deltaf_kind;

    // dN/(dxt dy) for one particle species
    std::vector<double> dN_dxtdy_for_one_particle_species;

    int number_of_chosen_particles;

    // store particle index;
    // the sampling process follows the order specified by this table
    std::vector<int> chosen_particles_sampling_table;

    // list for information for all particles
    std::vector<particle_info> particles;

    // list for information for all fluid cells
    long FO_length;
    const std::vector<FO_surf_LRF> &FOsurf_ptr;

    inline long determine_number_to_sample(
                    double dN, int model=1, double para1=0);

    const gsl_rng_type *gsl_type_random_number;
    gsl_rng *gsl_random_r;

    bool particles_are_the_same(int, int);

    //array for bulk delta f coefficients
    Table *bulkdf_coeff;

    // table parameter for diffusion deltaf coefficient
    int deltaf_qmu_coeff_table_length_T;
    int deltaf_qmu_coeff_table_length_mu;
    double delta_qmu_coeff_table_T0, delta_qmu_coeff_table_mu0;
    double delta_qmu_coeff_table_dT, delta_qmu_coeff_table_dmu;
    double **deltaf_qmu_coeff_tb;       // kappa_B [1/fm^3]

    // table parameter for bulk deltaf coefficient 14 moments
    int deltaf_bulk_coeff_14mom_table_length_T_;
    int deltaf_bulk_coeff_14mom_table_length_mu_;
    double deltaf_bulk_coeff_14mom_table_T0_;
    double deltaf_bulk_coeff_14mom_table_mu0_;
    double deltaf_bulk_coeff_14mom_table_dT_;
    double deltaf_bulk_coeff_14mom_table_dmu_;
    double **deltaf_bulk_coeff_14mom_c0_tb_;
    double **deltaf_bulk_coeff_14mom_c1_tb_;
    double **deltaf_bulk_coeff_14mom_c2_tb_;


    // arrays to speed up computing particle yield
    double sf_dx, sf_x_min, sf_x_max;
    int sf_tb_length;
    int sf_expint_truncate_order;
    double** sf_bessel_Kn;
    double** sf_expint_En;

    int flag_output_samples_into_files;
    int flag_store_samples_in_memory;

    //! particle decay
    int flag_perform_decays;
    particle_decay *decayer_ptr;

 public:
    FSSW(std::shared_ptr<RandomUtil::Random> ran_gen,
           Table* chosen_particles_in,
           std::vector<particle_info> particles_in,
           const std::vector<FO_surf_LRF> &FOsurf_ptr_in, int flag_PCE,
           ParameterReader* paraRdr_in, string path, string table_path,
           AfterburnerType afterburner_type);
    ~FSSW();

    std::vector< std::vector<iSS_Hadron>* >* Hadron_list;

    void initialize_special_function_arrays();
    double get_special_function_K1(double arg);
    double get_special_function_K2(double arg);
    double get_special_function_K3(double arg);
    void get_special_function_En(double arg, std::vector<double> &results);

    void shell();  // it all starts here...

    void combine_samples_to_OSCAR();
    std::string OSCAR_header_filename, OSCAR_output_filename;
    void combine_samples_to_gzip_file();
    void combine_samples_to_binary_file();

    void calculate_dN_analytic(const particle_info* particle, double mu,
                               double Temperature,
                               std::array<double, 6> &results);

    // the following variables need to be set first in order to
    // call this function
    void calculate_dN_dxtdy_for_one_particle_species(const int particle_idx);
    void sample_using_dN_dxtdy_4all_particles_conventional();

    void getbulkvisCoefficients(const double Tdec,
                                std::array<double, 3> &bulkvisCoefficients);
    void getbulkvisCoefficients(const double Tdec, const double mu_B,
                                std::array<double, 3> &bulkvisCoefficients);

    void load_bulk_deltaf_14mom_table(string filepath);
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
    double get_deltaf_bulk(
        const double mass, const double pdotu, const double bulkPi,
        const double Tdec, const int sign, const int baryon,
        const double f0, const std::array<double, 3> bulkvisCoefficients);
    int sample_momemtum_from_a_fluid_cell(
        const double mass, const int sign,
        const int baryon, const int strange, const int charge,
        const FO_surf_LRF *surf,
        const std::array<double, 3> bulkvisCoefficients,
        const double deltaf_qmu_coeff,
        double &pT, double &phi, double &y_minus_eta_s);
    void add_one_sampled_particle(
        const int repeated_sampling_idx, const FO_surf_LRF *surf,
        const int particle_monval, const double mass,
        const double pT, const double phi,
        const double y_minus_eta_s, const double eta_s);
    void boost_vector_back_to_lab_frame(iSS_data::Vec4 &p_LRF,
                                        iSS_data::Vec4 &p_lab,
                                        iSS_data::Vec4 &umu) const;
};

#endif  // SRC_FSSW_H_
