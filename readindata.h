//ver 1.1
#ifndef READINDATA_H
#define READINDATA_H

#include "main.h"
#include<fstream>
using namespace std;

typedef struct
{
   double tau, xpt, ypt;
   double da0, da1, da2;
   double vx, vy;
   double Edec, Tdec, Pdec;
   double Bn, muB, muS;
   double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
   double particle_mu[Maxparticle];
}FO_surf;

int get_filelength(string filepath);
void read_decdat(string path, int length, FO_surf* surf_ptr);
void read_surfdat(string path, int length, FO_surf* surf_ptr);
void read_decdat_mu(string path, int FO_length, int N_stable, double** particle_mu);
int read_resonance(particle_info* particle);
void calculate_particle_mu(int Nparticle, FO_surf* FOsurf_ptr, int FO_length, particle_info* particle, double** particle_mu);

#endif


/*
ver 1.2 10-01-2012
fixed potential bug in adding the anti-baryon decay list. Previously, the decay particle lists for anti-baryons are taken to be the same as their correspond baryons. Only for the purpose of calculating particles' chemical potential, this will be correct if baryon and anti-baryon have the same chemical potentials. Now, I added some more comparisons and produced the correct decay lists for anti-baryons. With correct Monte-Carlo numbers for both baryons and mesons. (major modifications are in read_resonance() function.)

Ver 1.1 04-12-2012
fixed calculation of the chemical potential for unstable particles in calculate_particle_mu() function.

Ver 1.0 04-11-2012
Change the structure of read in chemical potential. For chemical potential table less than FO_length, copy the last chemical potential values to the rest fluid cell.

*/
