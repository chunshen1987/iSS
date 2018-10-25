#ifndef DATA_STRUCT_H_
#define DATA_STRUCT_H_

#include <string>
#include <array>
#include <vector>

const double hbarC=0.197327053;  //GeV*fm

const int Maxparticle=500;            //size of array for storage of the particles
const int Maxdecaychannel=20;
const int Maxdecaypart=5;
const std::string table_path="iSS_tables";

typedef std::array<double, 4> Vec4;
typedef std::array<double, 10> ViscousVec;

typedef struct {
    int monval;     // Monte Carlo number according PDG
    std::string name;
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
} particle_info;

typedef struct {
    double tau, xpt, ypt, eta;
    double da0, da1, da2, da3;
    double u0, u1, u2, u3;
    double Edec, Tdec, Pdec;
    double Bn, muB, muS, muC;
    double pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
    double bulkPi;
    double qmu0, qmu1, qmu2, qmu3;
    double *particle_mu_PCE;
} FO_surf;

typedef struct {
    int decay_Npart;
    double branching_ratio;
    int decay_part[5];
} decay_channel_info;

typedef struct {
    int monval;     // Monte Carlo number according PDG
    std::string name;
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
    std::vector<decay_channel_info*> decay_channels;
    int sign;                   // Bose-Einstein or Dirac-Fermi statistics
} particle_decay_info;

struct iSS_Hadron {
     int pid;
     double mass;
     double E, px, py, pz;
     double t, x, y, z;
};

#endif  // DATA_STRUCT_H_
