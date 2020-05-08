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

    const int NpT_  = 15;
    const int Nphi_ = 48;
    const int Ny_   = 30;

    std::vector<double> pT_arr_, phi_arr_, y_arr_;
    std::vector<double> cos_phi_arr_, sin_phi_arr_;

    std::vector<iSS_data::Vec4> Smu_pT_;
    std::vector<iSS_data::Vec4> Smu_phi_;
    std::vector<iSS_data::Vec4> Smu_y_;

    std::vector< std::vector<iSS_data::Vec4> > Smu_pTdpTdphi_;

    double ***dN_pTdpTdphidy_;
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
        iSS_data::Vec4 &Smu, double &dN);
    void compute_integrated_spin_polarizations();
    void output_integrated_spin_polarizations(int POI_monval);
};

#endif  // SRC_SPIN_POLARIZATION_H_
