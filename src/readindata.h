//ver 1.1
#ifndef READINDATA_H
#define READINDATA_H

#include "main.h"
#include "ParameterReader.h"
#include<fstream>

using namespace std;

typedef struct
{
  int monval;     // Monte Carlo number according PDG
  string name;
  double mass;
  double width;
  int gspin;      // spin degeneracy
  int baryon;
  int strange;
  int charm;
  int bottom;
  int gisospin;   // isospin degeneracy
  int charge;
  int decays;     // amount of decays listed for this resonance
  int stable;     // defines whether this particle is considered as stable
  int decays_Npart[Maxdecaychannel];
  double decays_branchratio[Maxdecaychannel];
  int decays_part[Maxdecaychannel][Maxdecaypart];
  int sign;                   // Bose-Einstein or Dirac-Fermi statistics
}particle_info;

typedef struct
{
   double tau, xpt, ypt, eta;
   double da0, da1, da2, da3;
   double u0, u1, u2, u3;
   double Edec, Tdec, Pdec;
   double Bn, muB, muS;
   double pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
   double bulkPi;
   double qmu0, qmu1, qmu2, qmu3;
   double particle_mu[Maxparticle];
}FO_surf;

class read_FOdata
{
    private:
        ParameterReader* paraRdr;
        string path;
        int mode;

        int turn_on_bulk;       // switch to read in bulk viscous pressure
        int turn_on_rhob;       // switch to read in net baryon density
        int turn_on_diff;       // switch to read in diffusion current

        int n_eta_skip;
        int IEOS_music;

    public:
        read_FOdata(ParameterReader* paraRdr_in, string path);
        ~read_FOdata();
        
        int get_number_of_freezeout_cells();
        void read_in_freeze_out_data(int length, FO_surf* surf_ptr);
        int read_in_chemical_potentials(
            string path, int FO_length, FO_surf* surf_ptr, 
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
        void calculate_particle_mu(int Nparticle, 
                                   FO_surf* FOsurf_ptr, int FO_length, 
                                   particle_info* particle);
};

#endif


/*
ver 1.2 10-01-2012
fixed potential bug in adding the anti-baryon decay list. Previously, the decay particle lists for anti-baryons are taken to be the same as their correspond baryons. Only for the purpose of calculating particles' chemical potential, this will be correct if baryon and anti-baryon have the same chemical potentials. Now, I added some more comparisons and produced the correct decay lists for anti-baryons. With correct Monte-Carlo numbers for both baryons and mesons. (major modifications are in read_resonance() function.)

Ver 1.1 04-12-2012
fixed calculation of the chemical potential for unstable particles in calculate_particle_mu() function.

Ver 1.0 04-11-2012
Change the structure of read in chemical potential. For chemical potential table less than FO_length, copy the last chemical potential values to the rest fluid cell.

*/
