#ifndef ISS_H
#define ISS_H

#include <string>
#include <vector>
#include <memory>

#include "data_struct.h"
#include "ParameterReader.h"
#include "readindata.h"
#include "emissionfunction.h"
#include "Random.h"
#include "pretty_ostream.h"
#include "FSSW.h"

class iSS {
 private:
    const std::string path_;
    const std::string table_path_;
    const std::string particle_table_path_;
    const std::string surface_filename_;

    std::vector<FO_surf> FOsurf_array_;
    std::vector<FO_surf_LRF> FOsurf_LRF_array_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;

    int Nparticle;
    int flag_PCE_;
    AfterburnerType afterburner_type_;

    long randomSeed_;

    std::vector<particle_info> particle;

    std::unique_ptr<EmissionFunctionArray> efa_;
    std::unique_ptr<FSSW> spectra_sampler_;

    pretty_ostream messager;

 public:
    iSS(std::string path, std::string table_path="iSS_tables",
        std::string particle_table_path="iSS_tables",
        std::string inputfile="iSS_parameters.dat",
        std::string surface_filename="surface.dat");
    ~iSS();

    ParameterReader *paraRdr_ptr;

    void set_random_seed();
    void set_random_seed(int randomSeed_in);

    void perform_checks();
    void construct_Tmunu();

    int shell();
    int read_in_FO_surface();
    int read_in_FO_surface(std::vector<FO_surf> &surf_vec);
    int generate_samples();

    int get_number_of_sampled_events() {
        if (paraRdr_ptr->getVal("MC_sampling") == 4) {
            return(spectra_sampler_->get_number_of_sampled_events());
        } else {
            return(efa_->get_number_of_sampled_events());
        }
    };

    int get_number_of_particles(int iev) {
        if (paraRdr_ptr->getVal("MC_sampling") == 4) {
            return(spectra_sampler_->get_number_of_particles(iev));
        } else {
            return(efa_->get_number_of_particles(iev));
        }
    };

    iSS_Hadron get_hadron(int iev, int ipart) {
        if (paraRdr_ptr->getVal("MC_sampling") == 4) {
            return(spectra_sampler_->get_hadron(iev, ipart));
        } else {
            return(efa_->get_hadron(iev, ipart));
        }
    };

    std::vector<iSS_Hadron>* get_hadron_list_iev(const int iev) {
        if (paraRdr_ptr->getVal("MC_sampling") == 4) {
            return(spectra_sampler_->get_hadron_list_iev(iev));
        } else {
            return(efa_->get_hadron_list_iev(iev));
        }
    }

    void clear();

    // this function transform all the variables to local rest frame of the
    // fluid cell and trasform them to the t-z coordinate
    void transform_to_local_rest_frame(
                std::vector<FO_surf> &FOsurf_ptr,
                std::vector<FO_surf_LRF> &FOsurf_LRF_ptr);
};


#endif  // ISS_H
