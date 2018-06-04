// Copyright @ 2012 Chun Shen and Zhi Qiu
#ifndef SRC_READINDATA_H_
#define SRC_READINDATA_H_

#include <fstream>
#include <string>

#include "data_struct.h"
#include "ParameterReader.h"

using namespace std;

class read_FOdata {
 private:
    ParameterReader* paraRdr;
    string path;
    int mode;

    // flag to determine whether the EoS is partial chemical equilibrium or not
    int flag_PCE;
    int turn_on_bulk;       // switch to read in bulk viscous pressure
    int turn_on_rhob;       // switch to read in net baryon density
    int turn_on_diff;       // switch to read in diffusion current

    int n_eta_skip;
    int IEOS_music;

 public:
    read_FOdata(ParameterReader* paraRdr_in, string path);
    ~read_FOdata();
    int get_number_of_freezeout_cells();
    int get_flag_PCE() {return(flag_PCE);}
    void read_in_freeze_out_data(int length, FO_surf* surf_ptr);
    int read_in_chemical_potentials(string path, int FO_length,
                                    FO_surf* surf_ptr,
                                    particle_info* particle_ptr);
    void read_decdat(int length, FO_surf* surf_ptr);
    void read_surfdat(int length, FO_surf* surf_ptr);
    void read_FOsurfdat_VISH2p1(int length, FO_surf* surf_ptr);
    void read_FOsurfdat_MUSIC(int length, FO_surf* surf_ptr);
    void read_FOsurfdat_MUSIC_boost_invariant(int length,
                                              FO_surf* surf_ptr);
    void read_FOsurfdat_hydro_analysis_boost_invariant(int length,
                                                       FO_surf* surf_ptr);
    void read_decdat_mu(int FO_length, int N_stable, double** particle_mu);
    void read_chemical_potentials_music(int FO_length, FO_surf* FOsurf_ptr,
                                        int N_stable, double** particle_mu);
    int read_resonances_list(particle_info* particle);
    void calculate_particle_mu_PCE(int Nparticle, FO_surf* FOsurf_ptr,
                                   int FO_length, particle_info* particle,
                                   double** particle_mu);
    void regulate_Wmunu(double* u, double** Wmunu, double** Wmunu_regulated);
};

#endif  // SRC_READINDATA_H_

/*
ver 1.2 10-01-2012
fixed potential bug in adding the anti-baryon decay list. 
Previously, the decay particle lists for anti-baryons are taken to be the same
as their correspond baryons. Only for the purpose of calculating particles'
chemical potential, this will be correct if baryon and anti-baryon have the
same chemical potentials. Now, I added some more comparisons and produced
the correct decay lists for anti-baryons. With correct Monte-Carlo numbers
for both baryons and mesons.
(major modifications are in read_resonance() function.)

Ver 1.1 04-12-2012
fixed calculation of the chemical potential for unstable particles in
calculate_particle_mu() function.

Ver 1.0 04-11-2012
Change the structure of read in chemical potential. For chemical potential
table less than FO_length, copy the last chemical potential values to the
rest fluid cell.
*/
