// Copyright @ 2020 Chun Shen

#ifndef SRC_SPIN_POLARIZATION_H_
#define SRC_SPIN_POLARIZATION_H_

#include<vector>
#include<cstring>

#include "data_struct.h"
#include "Table.h"
#include "ParameterReader.h"

class SpinPolarization {
 private:
    const std::string path_;
    const std::string table_path_;
    const std::vector<FO_surf> &FOsurf_ptr_;
    const std::vector<particle_info> &particle_info_;
    std::vector<std::string> vorticity_typenames_;
    std::vector<std::string> thermal_shear_typenames_;
    const ParameterReader &paraRdr_;

    const int NpT_  = 30;
    const int Nphi_ = 48;
    const int Ny_   = 51;

    std::vector<double> pT_arr_, phi_arr_, y_arr_;
    std::vector<double> cos_phi_arr_, sin_phi_arr_;

    double dN_;
    std::vector<double> dN_pT_;
    std::vector<double> dN_phi_;
    std::vector<double> dN_y_;

    iSS_data::Vec4 Smu_;
    std::vector<iSS_data::Vec4> Smu_pT_;
    std::vector<iSS_data::Vec4> Smu_phi_;
    std::vector<iSS_data::Vec4> Smu_y_;

    iSS_data::Vec4 SmuLRF_;
    std::vector<iSS_data::Vec4> SmuLRF_pT_;
    std::vector<iSS_data::Vec4> SmuLRF_phi_;
    std::vector<iSS_data::Vec4> SmuLRF_y_;

    std::vector< std::vector<double> > dN_pTdpTdphi_;
    std::vector< std::vector<iSS_data::Vec4> > Smu_pTdpTdphi_;
    std::vector< std::vector<iSS_data::Vec4> > SmuLRF_pTdpTdphi_;

    std::vector<double> Rspin_y_;                     // Rspin(y)
    std::vector< std::vector<double> > dN_dpTdy_;
    std::vector< std::vector<double> > Rspin_pTy_;    // Rspin(y, pT)

    double ***dN_pTdpTdphidy_;
    double ****Smu_pTdpTdphidy_;
    double ****SmuLRF_pTdpTdphidy_;

 public:
    SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                     const std::vector<particle_info> &particles,
                     std::string path, std::string table_path,
                     ParameterReader &paraRdr);
    ~SpinPolarization();

    void compute_spin_polarization_shell();
    void compute_spin_polarization(
            const int POI_monval, const int irap_type, const int ivor_type,
            const int ithShearType, const int Flag_MuIP, const int Flag_SIP);
    void compute_spin_polarization_for_a_given_p(
        const particle_info &POI_info, const iSS_data::Vec4 &pmu,
        const int ivor_type, const int ithShearType,
        const int Flag_MuIP, const int Flag_SIP,
        iSS_data::Vec4 &Smu, iSS_data::Vec4 &SmuLRF, double &dN);
    void compute_integrated_spin_polarizations(double ****SmuMat,
                                               double ****SmuLRFMat);
    void output_integrated_spin_polarizations(
            const int POI_monval, const std::string rap_typename,
            const std::string vorticity_typename,
            const std::string thermal_shear_typenames,
            const int Flag_MuIP, const int Flag_SIP);
};

#endif  // SRC_SPIN_POLARIZATION_H_
