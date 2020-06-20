// Copyright @ 2012 Chun Shen and Zhi Qiu
#ifndef SRC_READINDATA_H_
#define SRC_READINDATA_H_

#include <fstream>
#include <string>
#include <vector>

#include "data_struct.h"
#include "ParameterReader.h"
#include "pretty_ostream.h"

class read_FOdata {
 private:
    ParameterReader* paraRdr;

    const std::string path_;
    const std::string table_path_;
    const std::string particle_table_path_;

    int mode;
    bool surface_in_binary;
    bool quantum_statistics_;

    // flag to determine whether the EoS is partial chemical equilibrium or not
    int flag_PCE_;
    int turn_on_bulk;       // switch to read in bulk viscous pressure
    int turn_on_rhob;       // switch to read in net baryon density
    int turn_on_diff;       // switch to read in diffusion current

    int n_eta_skip;
    int iEOS_MUSIC_;

    AfterburnerType afterburner_type_;

    pretty_ostream messager;

 public:
    read_FOdata(ParameterReader* paraRdr_in, std::string path,
                std::string table_path, std::string particle_table_path);
    ~read_FOdata();
    int get_number_of_freezeout_cells();
    int get_IEOS_music() const {return(iEOS_MUSIC_);}
    AfterburnerType get_afterburner_type() const {return(afterburner_type_);}
    int get_number_of_lines_of_binary_surface_file(string filename);
    int get_flag_PCE() {return(flag_PCE_);}
    void read_in_freeze_out_data(std::vector<FO_surf> &surf_ptr);
    void read_in_chemical_potentials(std::vector<FO_surf> &surf_ptr,
                                     std::vector<particle_info> &particle_ptr);
    void read_decdat(std::vector<FO_surf> &surf_ptr);
    void read_surfdat(std::vector<FO_surf> &surf_ptr);
    void read_FOsurfdat_VISH2p1(std::vector<FO_surf> &surf_ptr);
    void read_FOsurfdat_MUSIC(std::vector<FO_surf> &surf_ptr);
    void read_FOsurfdat_MUSIC_boost_invariant(std::vector<FO_surf> &surf_ptr);
    void read_FOsurfdat_hydro_analysis_boost_invariant(
                                        std::vector<FO_surf> &surf_ptr);
    void read_decdat_mu(int FO_length, int N_stable, double** particle_mu);
    void read_chemical_potentials_music(int FO_length,
                                        std::vector<FO_surf> &FOsurf_ptr,
                                        int N_stable, double** particle_mu);
    int read_resonances_list(std::vector<particle_info> &particle);
    void calculate_particle_mu_PCE(int Nparticle,
                                   std::vector<FO_surf> &FOsurf_ptr,
                                   int FO_length,
                                   std::vector<particle_info> &particle,
                                   double** particle_mu);
    void regulate_surface_cells(std::vector<FO_surf> &surf_ptr);
    void regulate_Wmunu(double u[4], double Wmunu[4][4],
                        double Wmunu_regulated[4][4]);
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
