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

class iSS {
 private:
    const std::string path_;
    const std::string table_path_;
    const std::string particle_table_path_;

    std::vector<FO_surf> FOsurf_ptr;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

    int Nparticle;
    int flag_PCE_;
    AfterburnerType afterburner_type_;

    long randomSeed;

    std::vector<particle_info> particle;

    std::unique_ptr<EmissionFunctionArray> efa;

    pretty_ostream messager;

 public:
    iSS(std::string path, std::string table_path="iSS_tables",
        std::string particle_table_path="iSS_tables",
        std::string inputfile="iSS_parameters.dat");
    ~iSS();

    ParameterReader *paraRdr_ptr;

    void set_random_seed();
    void set_random_seed(int randomSeed_in);

    int shell();
    int read_in_FO_surface();
    int generate_samples();

    int get_number_of_sampled_events() {
        return(efa->get_number_of_sampled_events());
    };

    int get_number_of_particles(int iev) {
        return(efa->get_number_of_particles(iev));
    };

    iSS_Hadron get_hadron(int iev, int ipart) {
        return(efa->get_hadron(iev, ipart));
    };

    std::vector<iSS_Hadron>* get_hadron_list_iev(const int iev) {
        return(efa->get_hadron_list_iev(iev));
    }

    void clear();
};


#endif  // ISS_H
