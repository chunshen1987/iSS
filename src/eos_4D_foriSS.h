// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_4D_H_
#define SRC_EOS_4D_H_

#include <cmath>
#include <array>
#include <vector>

#include "pretty_ostream.h"

class EOS_4D {
 private:
    // variables for header infos.
    float mubtilde0, muqtilde0, mustilde0, Ttilde0;
    float dmubtilde, dmuqtilde, dmustilde, dTtilde;
    int N_mub, N_muq, N_mus, N_T;

    // max values for the tilde variables
    float T_tilde_max;
    float mub_tilde_max;
    float muq_tilde_max;
    float mus_tilde_max;
    //double Ttilde, mubtilde, muqtilde, mustilde;

    // useful constants
    //const double alphaNf = (8/45.0 + 7/60.0*3.0)*M_PI*M_PI;
    //const double OneoveralphaNf = 1./alphaNf;

    const double hbarc3 = 0.19733*0.19733*0.19733;

    // 1D tables.
    std::vector<float> pressure_vec;
    std::vector<float> temp_vec;
    std::vector<float> mub_vec;
    std::vector<float> muq_vec;
    std::vector<float> mus_vec;

    std::vector<std::vector<float>> dfCoeffs_;

    // method to read/mainupalate header info and data
    void read_eos_binary(std::string filepath,
                         std::vector<float> &out, int header_size=12);
    void read_dfCoeffs_binary(std::string filepath, const int dfType);

    void get_eos_max_values();

    void read_header_binary(std::string filepath, int header_size=12);

    // Shift in the index corresponds to the header size.
    int index(int i_T, int i_mub, int i_muq, int i_mus) const;

    void FourDLInterp(const std::vector<float> &data,
                      const std::array<float, 4> &TildeVar,
                      std::array<float, 5> &ResArr,
                      bool compute_derivatives=false) const;
    void get_tilde_variables(double e, double nb, double nq, double ns,
                             std::array<float, 4> &TildeVar) const;
    pretty_ostream messenger;

 public:
    EOS_4D();
    ~EOS_4D();

    void initialize_eos();
    void initialize_dfCoeffs(const int dfType=0);

    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muB        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muS        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muQ        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_pressure   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;

    void getThermalVariables(const double epsilon, const double rhob,
                             const double rhoq, const double rhos,
                             std::vector<float> &thermalVec) const;
    void getDeltafCoeffs(const double epsilon, const double rhob,
                         const double rhoq, const double rhos,
                         std::vector<double> &deltafVec) const;
};

#endif  // SRC_EOS_4D_H_
