#ifndef ISS_H
#define ISS_H

#include <string>

#include "./ParameterReader.h"
#include "./readindata.h"
#include "./emissionfunction.h"

using namespace std;

class iSS {
 private:
    string path;
    
    int FO_length;
    FO_surf* FOsurf_ptr;

    int Nparticle;
    int flag_PCE;

    long randomSeed;

    particle_info *particle;
    
    EmissionFunctionArray *efa;

 public:
    iSS();
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
};


#endif  // ISS_H
