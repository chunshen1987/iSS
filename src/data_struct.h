#ifndef DATA_STRUCT_H_
#define DATA_STRUCT_H_

#include <string>
#include <array>
#include <vector>

namespace iSS_data {
    const double hbarC=0.197327053;  //GeV*fm

    typedef std::array<double, 4> Vec4;
    typedef std::array<double, 10> ViscousVec;

    const int NUMBER_OF_LINES_TO_WRITE = 100000;
    const int AMOUNT_OF_OUTPUT = 0;
}


enum class AfterburnerType {
    UrQMD,
    SMASH,
    JAM,
    PDG_Decay,
};


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
} particle_info;


typedef struct {
    float tau, xpt, ypt, eta;
    float da0, da1, da2, da3;
    float u0, u1, u2, u3;
    float Edec, Tdec, Pdec;
    float Bn, muB, muS, muC;
    float pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
    float bulkPi;
    float qmu0, qmu1, qmu2, qmu3;
    std::vector<float> particle_mu_PCE;
} FO_surf;


struct iSS_Hadron {
     int pid;
     double mass;
     double E, px, py, pz;
     double t, x, y, z;
};

#endif  // DATA_STRUCT_H_
