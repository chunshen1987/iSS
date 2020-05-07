// Copyright @ 2020 Chun Shen

#ifndef SRC_SPIN_POLARIZATION_H_
#define SRC_SPIN_POLARIZATION_H_

#include<vector>

#include "data_struct.h"
#include "Table.h"
#include "ParameterReader.h"

class SpinPolarization {
 private:
    const std::string path_;
    const std::string table_path_;
    const std::vector<FO_surf> &FOsurf_ptr_;
    const std::vector<particle_info> &particle_info_;
    const ParameterReader &paraRdr_;

    const int NpT  = 15;
    const int Nphi = 48;
    const int Ny   = 30;

    std::vector<double> pT_arr, phi_arr, y_arr;
    std::vector<double> cos_phi_arr, sin_phi_arr;
    double ***St_pTdpTdphidy_;
    double ***Sx_pTdpTdphidy_;
    double ***Sy_pTdpTdphidy_;
    double ***Sz_pTdpTdphidy_;

 public:
    SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                     const std::vector<particle_info> &particles,
                     std::string path, std::string table_path,
                     ParameterReader &paraRdr);
    ~SpinPolarization();

    void compute_spin_polarization();
    void compute_spin_polarization_for_a_given_p(
        const particle_info &POI_info, const iSS_data::Vec4 &pmu,
        iSS_data::Vec4 &Smu);
};

#endif  // SRC_SPIN_POLARIZATION_H_
