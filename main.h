#ifndef MAIN_H
#define MAIN_H

#include<string>
using namespace std;

const double hbarC=0.197327053;  //GeV*fm

const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

typedef struct
{
  int monval;			// Montecarlo number according PDG
  string name;
  double mass;
  double width;
  int gspin;			// spin degeneracy
  int baryon;
  int strange;
  int charm;
  int bottom;
  int gisospin;			// isospin degeneracy
  int charge;
  int decays;			// amount of decays listed for this resonance
  int stable;			// defines whether this particle is considered as stable
  int decays_Npart[Maxdecaychannel];
  double decays_branchratio[Maxdecaychannel];
  int decays_part[Maxdecaychannel][Maxdecaypart];
  int sign;       //Bose-Einstein or Dirac-Fermi statistics
}particle_info;

#endif
