#ifndef ISS_H
#define ISS_H

#include "./ParameterReader.h"
#include "./readindata.h"
#include "./emissionfunction.h"

using namespace std;

class iSS {
 private:
    string path;
    
    ParameterReader *paraRdr_ptr;
    int FO_length;
    FO_surf* FOsurf_ptr;

    int Nparticle;
    int flag_PCE;

    particle_info *particle;

 public:
    iSS(ParameterReader *paraRdr_in);
    ~iSS();

    int shell();
    int read_in_FO_surface();
    int generate_samples();
};


#endif  // ISS_H
