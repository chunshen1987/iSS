// Copyright @ 2020 Chun Shen

#include "spin_polarization.h"
#include "arsenal.h"
#include "data_struct.h"
#include<iostream>
#include<sstream>
#include<cstdlib>
#include<fstream>
#include<iomanip>

#ifndef _OPENMP
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

using std::cout;
using std::endl;
using std::scientific;
using std::setprecision;
using std::setw;
using iSS_data::hbarC;

SpinPolarization::SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                                   const std::vector<particle_info> &particles,
                                   std::string path, std::string table_path,
                                   ParameterReader &paraRdr) :
        path_(path), table_path_(table_path),
        FOsurf_ptr_(FOsurf_ptr), particle_info_(particles),
        paraRdr_(paraRdr) {

    double dpT = 3.0/NpT_;
    pT_arr_.resize(NpT_);
    for (int i = 0; i < NpT_; i++)
        pT_arr_[i] = (i + 0.5)*dpT;
    double dphi = 2.*M_PI/Nphi_;
    phi_arr_.resize(Nphi_);
    cos_phi_arr_.resize(Nphi_);
    sin_phi_arr_.resize(Nphi_);
    for (int i = 0; i < Nphi_; i++) {
        phi_arr_[i] = 0.0 + i*dphi;
        cos_phi_arr_[i] = cos(phi_arr_[i]);
        sin_phi_arr_[i] = sin(phi_arr_[i]);
    }
    double y_size = 10.0;
    double dy = y_size/(Ny_ - 1);
    y_arr_.resize(Ny_);
    for (int i = 0; i < Ny_; i++)
        y_arr_[i] = - y_size/2. + i*dy;

    dN_pT_.resize(NpT_, 0.);
    dN_phi_.resize(Nphi_, 0.);
    dN_y_.resize(Ny_, 0.);
    Smu_pT_.resize(NpT_, {0., 0., 0., 0.});
    Smu_phi_.resize(Nphi_, {0., 0., 0., 0.});
    Smu_y_.resize(Ny_, {0., 0., 0., 0.});
    SmuLRF_pT_.resize(NpT_, {0., 0., 0., 0.});
    SmuLRF_phi_.resize(Nphi_, {0., 0., 0., 0.});
    SmuLRF_y_.resize(Ny_, {0., 0., 0., 0.});

    Rspin_y_.resize(Ny_, 0.);

    for (int i = 0; i < NpT_; i++) {
        std::vector<iSS_data::Vec4> Smu(Nphi_, {0., 0., 0., 0.});
        Smu_pTdpTdphi_.push_back(Smu);
        std::vector<iSS_data::Vec4> SmuLRF(Nphi_, {0., 0., 0., 0.});
        SmuLRF_pTdpTdphi_.push_back(SmuLRF);

        std::vector<double> dNtmp(Nphi_, 0.);
        dN_pTdpTdphi_.push_back(dNtmp);
    }
    for (int i = 0; i < Ny_; i++) {
        std::vector<double> Rtmp(NpT_, 0.);
        Rspin_pTy_.push_back(Rtmp);
        std::vector<double> dNdpTtmp(NpT_, 0.);
        dN_dpTdy_.push_back(dNdpTtmp);
    }

    dN_pTdpTdphidy_ = create_a_3D_Matrix(Ny_, NpT_, Nphi_, 0.);
    Smu_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuLRF_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuSIP1_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuSIP1LRF_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuSIP2_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuSIP2LRF_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuFull_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);
    SmuFullLRF_pTdpTdphidy_ = create_a_4D_Matrix(Ny_, NpT_, Nphi_, 4, 0.);

    vorticity_typenames_.push_back("KineticSP");
    vorticity_typenames_.push_back("Kinetic");
    vorticity_typenames_.push_back("Thermal");
    vorticity_typenames_.push_back("Temperature");
    thermal_shear_typenames_.push_back("BBPP");
    thermal_shear_typenames_.push_back("LY");
    rapidity_typenames_.push_back("rapidity");
    rapidity_typenames_.push_back("pseudorapidity");
}


SpinPolarization::~SpinPolarization() {
    delete_a_3D_Matrix(dN_pTdpTdphidy_, Ny_, NpT_);
    delete_a_4D_Matrix(Smu_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuLRF_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuSIP1_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuSIP1LRF_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuSIP2_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuSIP2LRF_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuFull_pTdpTdphidy_, Ny_, NpT_, Nphi_);
    delete_a_4D_Matrix(SmuFullLRF_pTdpTdphidy_, Ny_, NpT_, Nphi_);
}


void SpinPolarization::compute_spin_polarization_shell() {
    #pragma omp parallel for
    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            cout << "OpenMP: using " << omp_get_num_threads()
                 << " threads." << endl;
        }
    }

    // by default we compute Polarization in rapidity and pseudorapidity
    int PolRapType = static_cast<int>(paraRdr_.getVal("polarizationRapType",
                                                      2));
    int ivor_type = 2;  // thermal voriticity

    std::array<int, 2> POI_list = {3122, -3122};   // Lambda and Anti-Lambda
    for (const auto &POI_monval: POI_list) {
        if (PolRapType < 2) {
            compute_spin_polarization(POI_monval, PolRapType, ivor_type);
        } else {
            for (int irap = 0; irap < PolRapType; irap++) {
                compute_spin_polarization(POI_monval, irap, ivor_type);
            }
        }
    }
}


void SpinPolarization::compute_spin_polarization(
        const int POI_monval, const int irap_type, const int ivor_type) {
    // first clean up previous results
    set_val_in_3D_Matrix(dN_pTdpTdphidy_, Ny_, NpT_, Nphi_, 0.0);
    set_val_in_4D_Matrix(Smu_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuLRF_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuSIP1_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuSIP1LRF_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuSIP2_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuSIP2LRF_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuFull_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);
    set_val_in_4D_Matrix(SmuFullLRF_pTdpTdphidy_, Ny_, NpT_, Nphi_, 4, 0.0);

    // find the information about the particle of interest
    particle_info POI_info;
    for (const auto &particle_i: particle_info_) {
        if (particle_i.monval == POI_monval) {
            POI_info = particle_i;
            break;
        }
    }

    cout << "Computing spin polarization for " << POI_info.name
         << ", Monte-carlo index: " << POI_info.monval << endl;
    if (irap_type == 0) {
        cout << "Rapidity type : rapidity" << endl;
    } else {
        cout << "Rapidity type : pseudo-rapidity" << endl;
    }
    cout << "spin polarization tensor type : "
         << vorticity_typenames_[ivor_type] << endl;

    const double mass  = POI_info.mass;
    const double prefactor = 1./(2.*mass);

    for (int iy = 0; iy < Ny_; iy++) {
        cout << "progress: " << iy << "/" << Ny_ << endl;
        double y = y_arr_[iy];
        double cosh_y = cosh(y);
        double sinh_y = sinh(y);
        for (int ipT = 0; ipT < NpT_; ipT++) {
            double p_perp = pT_arr_[ipT];
            double p0, pz;
            if (irap_type == 0) {
                // rapidity
                double m_perp = sqrt(mass*mass + p_perp*p_perp);
                p0 = m_perp*cosh_y;
                pz = m_perp*sinh_y;
            } else {
                // pseudo-rapidity
                p0 = sqrt(p_perp*cosh_y*p_perp*cosh_y + mass*mass);
                pz = p_perp*sinh_y;
            }
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                double px = p_perp*cos_phi_arr_[iphi];
                double py = p_perp*sin_phi_arr_[iphi];

                iSS_data::Vec4 pmu = {static_cast<float>(p0),
                                      static_cast<float>(px),
                                      static_cast<float>(py),
                                      static_cast<float>(pz)};
                iSS_data::Vec4 Smu = {0., 0., 0., 0.};
                iSS_data::Vec4 Smu_LRF = {0., 0., 0., 0.};
                iSS_data::Vec4 SmuSIP1 = {0., 0., 0., 0.};
                iSS_data::Vec4 SmuSIP1_LRF = {0., 0., 0., 0.};
                iSS_data::Vec4 SmuSIP2 = {0., 0., 0., 0.};
                iSS_data::Vec4 SmuSIP2_LRF = {0., 0., 0., 0.};
                iSS_data::Vec4 SmuFull = {0., 0., 0., 0.};
                iSS_data::Vec4 SmuFull_LRF = {0., 0., 0., 0.};
                double dN = 0.;
                compute_spin_polarization_for_a_given_p(
                    POI_info, pmu, ivor_type, dN,
                    Smu, Smu_LRF, SmuSIP1, SmuSIP1_LRF, SmuSIP2, SmuSIP2_LRF,
                    SmuFull, SmuFull_LRF);

                dN_pTdpTdphidy_[iy][ipT][iphi] = dN;
                for (int mu = 0; mu < 4; mu++) {
                    Smu_pTdpTdphidy_[iy][ipT][iphi][mu] = prefactor*Smu[mu];
                    SmuLRF_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*Smu_LRF[mu]);
                    SmuSIP1_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*SmuSIP1[mu]);
                    SmuSIP1LRF_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*SmuSIP1_LRF[mu]);
                    SmuSIP2_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*SmuSIP2[mu]);
                    SmuSIP2LRF_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*SmuSIP2_LRF[mu]);
                    SmuFull_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*SmuFull[mu]);
                    SmuFullLRF_pTdpTdphidy_[iy][ipT][iphi][mu] = (
                                                prefactor*SmuFull_LRF[mu]);
                }
            }
        }
    }

    // vorticity
    compute_integrated_spin_polarizations(Smu_pTdpTdphidy_,
                                          SmuLRF_pTdpTdphidy_);
    output_integrated_spin_polarizations(
        POI_monval, rapidity_typenames_[irap_type],
        vorticity_typenames_[ivor_type], thermal_shear_typenames_[0], 0, 0);
    // vorticity + SIP(BBPP)
    compute_integrated_spin_polarizations(SmuSIP1_pTdpTdphidy_,
                                          SmuSIP1LRF_pTdpTdphidy_);
    output_integrated_spin_polarizations(
        POI_monval, rapidity_typenames_[irap_type],
        vorticity_typenames_[ivor_type], thermal_shear_typenames_[0], 0, 1);
    // vorticity + SIP(LY)
    compute_integrated_spin_polarizations(SmuSIP2_pTdpTdphidy_,
                                          SmuSIP2LRF_pTdpTdphidy_);
    output_integrated_spin_polarizations(
        POI_monval, rapidity_typenames_[irap_type],
        vorticity_typenames_[ivor_type], thermal_shear_typenames_[1], 0, 1);
    // vorticity + SIP(LY) + muIP
    compute_integrated_spin_polarizations(SmuFull_pTdpTdphidy_,
                                          SmuFullLRF_pTdpTdphidy_);
    output_integrated_spin_polarizations(
        POI_monval, rapidity_typenames_[irap_type],
        vorticity_typenames_[ivor_type], thermal_shear_typenames_[1], 1, 1);
}


void SpinPolarization::compute_integrated_spin_polarizations(
                            double ****SmuMat, double ****SmuLRFMat) {
    cout << "computing the integrated spin polarization ... " << endl;

    const double dpT = pT_arr_[1] - pT_arr_[0];
    const double dphi = phi_arr_[1] - phi_arr_[0];
    const double dy = y_arr_[1] - y_arr_[0];
    // clean up previous results
    Smu_ = {0.};
    SmuLRF_ = {0.};
    dN_ = 0.;
    for (int ipT = 0; ipT < NpT_; ipT++) {
        dN_pT_[ipT] = 0.;
        Smu_pT_[ipT] = {0.};
        SmuLRF_pT_[ipT] = {0.};
    }

    for (int iphi = 0; iphi < Nphi_; iphi++) {
        dN_phi_[iphi] = 0.;
        Smu_phi_[iphi] = {0.};
        SmuLRF_phi_[iphi] = {0.};
    }

    for (int iy = 0; iy < Ny_; iy++) {
        dN_y_[iy] = 0.;
        Smu_y_[iy] = {0.};
        SmuLRF_y_[iy] = {0.};
        Rspin_y_[iy] = 0.;
        for (int ipT = 0; ipT < NpT_; ipT++) {
            dN_dpTdy_[iy][ipT] = 0.;
            Rspin_pTy_[iy][ipT] = 0.;
        }
    }

    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iphi = 0; iphi < Nphi_; iphi++) {
            dN_pTdpTdphi_[ipT][iphi] = 0.0;
            for (int i = 0; i < 4; i++) {
                Smu_pTdpTdphi_[ipT][iphi][i] = 0.0;
                SmuLRF_pTdpTdphi_[ipT][iphi][i] = 0.0;
            }
        }
    }

    // compute S^mu
    double dN = 0.;
    for (int iy = 0; iy < Ny_; iy++) {
        if (std::abs(y_arr_[iy]) > 1.) continue;
        for (int ipT = 0; ipT < NpT_; ipT++) {
            if (pT_arr_[ipT] < 0.5) continue;
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                dN_ += dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                dN  += dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                for (int mu = 0; mu < 4; mu++) {
                    Smu_[mu] += SmuMat[iy][ipT][iphi][mu]*pT_arr_[ipT];
                    SmuLRF_[mu] += SmuLRFMat[iy][ipT][iphi][mu]*pT_arr_[ipT];
                }
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        Smu_[i] /= dN;
        SmuLRF_[i] /= dN;
    }
    dN_ *= dpT*dphi*dy;

    // compute S^mu(pT)
    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iy = 0; iy < Ny_; iy++) {
            if (std::abs(y_arr_[iy]) > 1.) continue;
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                dN_pT_[ipT] += dN_pTdpTdphidy_[iy][ipT][iphi];
                for (int mu = 0; mu < 4; mu++) {
                    Smu_pT_[ipT][mu] += SmuMat[iy][ipT][iphi][mu];
                    SmuLRF_pT_[ipT][mu] += SmuLRFMat[iy][ipT][iphi][mu];
                }
            }
        }
        for (int i = 0; i < 4; i++) {
            Smu_pT_[ipT][i] /= dN_pT_[ipT];
            SmuLRF_pT_[ipT][i] /= dN_pT_[ipT];
        }
        dN_pT_[ipT] *= dphi*dy;
    }

    // compute S^mu(phi)
    for (int iphi = 0; iphi < Nphi_; iphi++) {
        double dN_phi = 0.;
        for (int iy = 0; iy < Ny_; iy++) {
            if (std::abs(y_arr_[iy]) > 1.) continue;
            for (int ipT = 0; ipT < NpT_; ipT++) {
                if (pT_arr_[ipT] < 0.5) continue;
                dN_phi_[iphi] += dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                dN_phi        += dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                for (int mu = 0; mu < 4; mu++) {
                    Smu_phi_[iphi][mu] += (
                            SmuMat[iy][ipT][iphi][mu]*pT_arr_[ipT]);
                    SmuLRF_phi_[iphi][mu] += (
                            SmuLRFMat[iy][ipT][iphi][mu]*pT_arr_[ipT]);
                }
            }
        }
        for (int i = 0; i < 4; i++) {
            Smu_phi_[iphi][i] /= dN_phi;
            SmuLRF_phi_[iphi][i] /= dN_phi;
        }
        dN_phi_[iphi] *= dpT*dy;
    }

    // compute S^mu(y)
    for (int iy = 0; iy < Ny_; iy++) {
        double dN_y = 0.;
        for (int ipT = 0; ipT < NpT_; ipT++) {
            if (pT_arr_[ipT] < 0.5) continue;
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                dN_y_[iy] += dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                dN_y      += dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                for (int mu = 0; mu < 4; mu++) {
                    Smu_y_[iy][mu] += SmuMat[iy][ipT][iphi][mu]*pT_arr_[ipT];
                    SmuLRF_y_[iy][mu] += (
                            SmuLRFMat[iy][ipT][iphi][mu]*pT_arr_[ipT]);
                }
            }
        }
        for (int i = 0; i < 4; i++) {
            Smu_y_[iy][i] /= dN_y;
            SmuLRF_y_[iy][i] /= dN_y;
        }
        dN_y_[iy] *= dpT*dphi;
    }

    // compute S^mu(pT, phi)
    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iphi = 0; iphi < Nphi_; iphi++) {
            for (int iy = 0; iy < Ny_; iy++) {
                if (std::abs(y_arr_[iy]) > 1.) continue;
                dN_pTdpTdphi_[ipT][iphi] += dN_pTdpTdphidy_[iy][ipT][iphi];
                for (int mu = 0; mu < 4; mu++) {
                    Smu_pTdpTdphi_[ipT][iphi][mu] += SmuMat[iy][ipT][iphi][mu];
                    SmuLRF_pTdpTdphi_[ipT][iphi][mu] += (
                                            SmuLRFMat[iy][ipT][iphi][mu]);
                }
            }
            for (int i = 0; i < 4; i++) {
                Smu_pTdpTdphi_[ipT][iphi][i] /= dN_pTdpTdphi_[ipT][iphi];
                SmuLRF_pTdpTdphi_[ipT][iphi][i] /= dN_pTdpTdphi_[ipT][iphi];
            }
            dN_pTdpTdphi_[ipT][iphi] *= dy;
        }
    }

    // compute Rspin(pT, y) and Rspin(y)
    for (int iy = 0; iy < Ny_; iy++) {
        double dN_y = 0.;
        for (int ipT = 0; ipT < NpT_; ipT++) {
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                double dNtmp = dN_pTdpTdphidy_[iy][ipT][iphi]*pT_arr_[ipT];
                double pxhat = cos(phi_arr_[iphi]);
                double pyhat = sin(phi_arr_[iphi]);
                double Rspin = (- SmuMat[iy][ipT][iphi][1]*pyhat
                                + SmuMat[iy][ipT][iphi][2]*pxhat);
                dN_dpTdy_[iy][ipT]  += dNtmp;
                Rspin_pTy_[iy][ipT] += Rspin*pT_arr_[ipT];

                if (pT_arr_[ipT] > 0.5) {
                    dN_y += dNtmp;
                    Rspin_y_[iy] += Rspin*pT_arr_[ipT];
                }
            }
            Rspin_pTy_[iy][ipT] /= dN_dpTdy_[iy][ipT];
            dN_dpTdy_[iy][ipT] *= dphi;
        }
        Rspin_y_[iy] /= dN_y;
    }
}


void SpinPolarization::output_integrated_spin_polarizations(
        const int POI_monval, const std::string rap_typename,
        const std::string vorticity_typename,
        const std::string thermal_shear_typenames,
        const int Flag_MuIP, const int Flag_SIP) {
    cout << "output spin polarization results to files ... " << endl;
    std::ofstream of;
    std::stringstream fileTypeName;

    fileTypeName << vorticity_typename << "_" << rap_typename
                 << "_" << POI_monval;
    if (Flag_MuIP == 1)
        fileTypeName << "_" << "wMuIP";
    if (Flag_SIP == 1) 
        fileTypeName << "_wSIP_" << thermal_shear_typenames;

    std::stringstream Smu_filename;
    Smu_filename << path_ << "/Smu_" << fileTypeName.str() << ".dat";
    remove(Smu_filename.str().c_str());
    of.open(Smu_filename.str().c_str(), std::ios::out);
    of << "# dN S^t  S^x  S^y  S^z S^t_LRF S^x_LRF S^y_LRF S^z_LRF" << endl;
    of << scientific << setw(10) << setprecision(6) << dN_ << "  ";
    for (int i = 0; i < 4; i++) {
        of << scientific << setw(10) << setprecision(6)
           << Smu_[i] << "  ";
    }
    for (int i = 0; i < 4; i++) {
        of << scientific << setw(10) << setprecision(6)
           << SmuLRF_[i] << "  ";
    }
    of << endl;
    of.close();

    std::stringstream SmupT_filename;
    SmupT_filename << path_ << "/Smu_pT_" << fileTypeName.str() << ".dat";
    remove(SmupT_filename.str().c_str());
    of.open(SmupT_filename.str().c_str(), std::ios::out);
    of << "# pT[GeV]  dN/(pTdpT)[GeV^-2]  S^t(pT)  S^x(pT)  S^y(pT)  S^z(pT)  "
       << "S^t_LRF(pT)  S^x_LRF(pT)  S^y_LRF(pT)  S^z_LRF(pT)" << endl;
    for (int ipT = 0; ipT < NpT_; ipT++) {
        of << scientific << setw(10) << setprecision(6)
           << pT_arr_[ipT] << "  " << dN_pT_[ipT] << "  ";
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << Smu_pT_[ipT][i] << "  ";
        }
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << SmuLRF_pT_[ipT][i] << "  ";
        }
        of << endl;
    }
    of.close();

    std::stringstream Smuphi_filename;
    Smuphi_filename << path_ << "/Smu_phi_" << fileTypeName.str() << ".dat";
    remove(Smuphi_filename.str().c_str());
    of.open(Smuphi_filename.str().c_str(), std::ios::out);
    of << "# phi  dN/dphi  S^t(phi)  S^x(phi)  S^y(phi)  S^z(phi)  "
       << "S^t_LRF(phi)  S^x_LRF(phi)  S^y_LRF(phi)  S^z_LRF(phi)" << endl;
    for (int iphi = 0; iphi < Nphi_; iphi++) {
        of << scientific << setw(10) << setprecision(6)
           << phi_arr_[iphi] << "  " << dN_phi_[iphi] << "  ";
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << Smu_phi_[iphi][i] << "  ";
        }
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << SmuLRF_phi_[iphi][i] << "  ";
        }
        of << endl;
    }
    of.close();

    std::stringstream Smuy_filename;
    Smuy_filename << path_ << "/Smu_y_" << fileTypeName.str() << ".dat";
    remove(Smuy_filename.str().c_str());
    of.open(Smuy_filename.str().c_str(), std::ios::out);
    of << "# y  dN/dy  S^t(y)  S^x(y)  S^y(y)  S^z(y)  "
       << "S^t_LRF(y)  S^x_LRF(y)  S^y_LRF(y)  S^z_LRF(y)" << endl;
    for (int iy = 0; iy < Ny_; iy++) {
        of << scientific << setw(10) << setprecision(6)
           << y_arr_[iy] << "  " << dN_y_[iy] << "  ";
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << Smu_y_[iy][i] << "  ";
        }
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << SmuLRF_y_[iy][i] << "  ";
        }
        of << endl;
    }
    of.close();

    std::stringstream Smu_dpTdphi_filename;
    Smu_dpTdphi_filename << path_ << "/Smu_dpTdphi_" << fileTypeName.str()
                         << ".dat";
    remove(Smu_dpTdphi_filename.str().c_str());
    of.open(Smu_dpTdphi_filename.str().c_str(), std::ios::out);
    of << "# pT[GeV]  phi  dN/(pTdpTdphi)[GeV^-2}  "
       << "S^t(pT, phi)  S^x(pT, phi)  S^y(pT, phi)  S^z(pT, phi)  "
       << "S^t_LRF(pT, phi)  S^x_LRF(pT, phi)  S^y_LRF(pT, phi)  "
       << "S^z_LRF(pT, phi)" << endl;
    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iphi = 0; iphi < Nphi_; iphi++) {
            of << scientific << setw(10) << setprecision(6)
               << pT_arr_[ipT] << "  " << phi_arr_[iphi] << "  "
               << dN_pTdpTdphi_[ipT][iphi] << "  ";
            for (int i = 0; i < 4; i++) {
                of << scientific << setw(10) << setprecision(6)
                   << Smu_pTdpTdphi_[ipT][iphi][i] << "  ";
            }
            for (int i = 0; i < 4; i++) {
                of << scientific << setw(10) << setprecision(6)
                   << SmuLRF_pTdpTdphi_[ipT][iphi][i] << "  ";
            }
            of << endl;
        }
    }
    of.close();

    std::stringstream Rspin_pTy_filename;
    Rspin_pTy_filename << path_ << "/Rspin_pTy_" << fileTypeName.str()
                       << ".dat";
    remove(Rspin_pTy_filename.str().c_str());
    of.open(Rspin_pTy_filename.str().c_str(), std::ios::out);
    of << "# y  pT[GeV]  dN/(dpTdy) [GeV^-1]  Rspin" << endl;
    for (int iy = 0; iy < Ny_; iy++) {
        for (int ipT = 0; ipT < NpT_; ipT++) {
            of << scientific << setw(10) << setprecision(6)
               << y_arr_[iy] << "  " << pT_arr_[ipT] << "  "
               << dN_dpTdy_[iy][ipT] << "  "
               << Rspin_pTy_[iy][ipT] << endl;
        }
    }
    of.close();

    std::stringstream Rspin_y_filename;
    Rspin_y_filename << path_ << "/Rspin_y_" << fileTypeName.str() << ".dat";
    remove(Rspin_y_filename.str().c_str());
    of.open(Rspin_y_filename.str().c_str(), std::ios::out);
    of << "# y  Rspin" << endl;
    for (int iy = 0; iy < Ny_; iy++) {
        of << scientific << setw(10) << setprecision(6)
           << y_arr_[iy] << "  " << Rspin_y_[iy] << endl;
    }
    of.close();

    //std::stringstream Smu_dpTdphidy_filename;
    //Smu_dpTdphidy_filename << path_ << "/Smu_dpTdphidy_" << fileTypeName.str()
    //                       << ".dat";
    //remove(Smu_dpTdphidy_filename.str().c_str());
    //of.open(Smu_dpTdphidy_filename.str().c_str(), std::ios::out);
    //of << "# y  pT[GeV]  phi  dN/(pTdpTdphidy)[GeV^-2]  S^t  S^x  S^y  S^z  "
    //   << "S^t_LRF  S^x_LRF  S^y_LRF  S^z_LRF" << endl;
    //for (int iy = 0; iy < Ny_; iy++) {
    //    for (int ipT = 0; ipT < NpT_; ipT++) {
    //        for (int iphi = 0; iphi < Nphi_; iphi++) {
    //            of << scientific << setw(10) << setprecision(6)
    //               << y_arr_[iy] << "  " << pT_arr_[ipT] << "  "
    //               << phi_arr_[iphi] << "  "
    //               << dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << St_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << Sx_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << Sy_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << Sz_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << StLRF_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << SxLRF_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << SyLRF_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << SzLRF_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
    //               << endl;
    //        }
    //    }
    //}
    //of.close();
}


void SpinPolarization::compute_spin_polarization_for_a_given_p(
        const particle_info &POI_info, const iSS_data::Vec4 &pmu,
        const int ivor_type, double &dN,
        iSS_data::Vec4 &Smu, iSS_data::Vec4 &SmuLRF,
        iSS_data::Vec4 &SmuSIP1, iSS_data::Vec4 &SmuSIP1LRF,
        iSS_data::Vec4 &SmuSIP2, iSS_data::Vec4 &SmuSIP2LRF,
        iSS_data::Vec4 &SmuFull, iSS_data::Vec4 &SmuFullLRF) {
    const double hbarC3 = hbarC*hbarC*hbarC;
    double Smu_tmp[4] = {0., 0., 0., 0.};
    double SmuSIP1_tmp[4] = {0., 0., 0., 0.};
    double SmuSIP2_tmp[4] = {0., 0., 0., 0.};
    double SmuFull_tmp[4] = {0., 0., 0., 0.};
    #pragma omp parallel for reduction(+: Smu_tmp[:4], SmuSIP1_tmp[:4], SmuSIP2_tmp[:4], SmuFull_tmp[:4], dN)
    for (unsigned int i = 0; i < FOsurf_ptr_.size(); i++) {
        const FO_surf &surf = FOsurf_ptr_[i];
        const float tau = surf.tau;
        const float ut = surf.u0*surf.cosh_eta + surf.u3*surf.sinh_eta;
        const float ux = surf.u1;
        const float uy = surf.u2;
        const float uz = surf.u0*surf.sinh_eta + surf.u3*surf.cosh_eta;
        const float ptau = pmu[0]*surf.cosh_eta - pmu[3]*surf.sinh_eta;
        const float tau_peta = pmu[3]*surf.cosh_eta - pmu[0]*surf.sinh_eta;   // tau*p^eta
        const double mu = (  POI_info.baryon*surf.muB
                           + POI_info.strange*surf.muS
                           + POI_info.charge*surf.muC);
        const double pdotu = pmu[0]*ut - pmu[1]*ux - pmu[2]*uy - pmu[3]*uz;
        const double pdsigma = (tau*(  ptau*surf.da0 + pmu[1]*surf.da1
                                     + pmu[2]*surf.da2 + tau_peta*surf.da3/tau)
                                /hbarC3);    // [GeV^-2]
        const double expon = (pdotu - mu)/surf.Tdec;
        const double f0 = 1./(exp(expon) + POI_info.sign);
        const double prefactor = pdsigma*f0*(1. - f0);

        dN += pdsigma*f0;

        // thermal vorticity
        const float omega_tx = surf.vorticity_arr[6*ivor_type + 0];
        const float omega_ty = surf.vorticity_arr[6*ivor_type + 1];
        const float omega_tz = surf.vorticity_arr[6*ivor_type + 2];
        const float omega_xy = surf.vorticity_arr[6*ivor_type + 3];
        const float omega_xz = surf.vorticity_arr[6*ivor_type + 4];
        const float omega_yz = surf.vorticity_arr[6*ivor_type + 5];

        float sigma[10];
        for (int j = 0; j < 10; j++)
            sigma[j] = surf.vorticity_arr[24 + j];

        float Smu_th[4] = {0., 0., 0., 0.};
        Smu_th[0] = prefactor*(
            (omega_yz*pmu[1] - omega_xz*pmu[2] + omega_xy*pmu[3]));
        Smu_th[1] = prefactor*(
            -(- omega_yz*pmu[0] - omega_ty*pmu[3] + omega_tz*pmu[2]));
        Smu_th[2] = prefactor*(
            -(omega_tx*pmu[3] - omega_tz*pmu[1] + omega_xz*pmu[0]));
        Smu_th[3] = prefactor*(
            -(- omega_tx*pmu[2] + omega_ty*pmu[1] - omega_xy*pmu[0]));

        for (int imu = 0; imu < 4; imu++) {
            Smu_tmp[imu] += Smu_th[imu];
        }

        // thermal shear tensor
        float Smu_SIP1[4] = {0., 0., 0., 0.};
        float pdotsigma[4];
        const double prefactorSIP = prefactor/pdotu;
        // Becattini
        pdotsigma[0] = (  pmu[0]*sigma[0] - pmu[1]*sigma[1]
                        - pmu[2]*sigma[2] - pmu[3]*sigma[3]);
        pdotsigma[1] = (  pmu[0]*sigma[1] - pmu[1]*sigma[4]
                        - pmu[2]*sigma[5] - pmu[3]*sigma[6]);
        pdotsigma[2] = (  pmu[0]*sigma[2] - pmu[1]*sigma[5]
                        - pmu[2]*sigma[7] - pmu[3]*sigma[8]);
        pdotsigma[3] = (  pmu[0]*sigma[3] - pmu[1]*sigma[6]
                        - pmu[2]*sigma[8] - pmu[3]*sigma[9]);
        Smu_SIP1[0] = 0.;
        Smu_SIP1[1] = prefactorSIP*(- pmu[2]*pdotsigma[3]
                                    + pmu[3]*pdotsigma[2]);
        Smu_SIP1[2] = prefactorSIP*(  pmu[1]*pdotsigma[3]
                                    - pmu[3]*pdotsigma[1]);
        Smu_SIP1[3] = prefactorSIP*(- pmu[1]*pdotsigma[2]
                                    + pmu[2]*pdotsigma[1]);

        for (int imu = 0; imu < 4; imu++) {
            SmuSIP1_tmp[imu] += Smu_th[imu] + Smu_SIP1[imu];
        }

        // Liu-Yin
        float Smu_SIP2[4] = {0., 0., 0., 0.};
        double p_perpMu[4] = {
            pmu[0] - pdotu*ut,
            pmu[1] - pdotu*ux,
            pmu[2] - pdotu*uy,
            pmu[3] - pdotu*uz,
        };
        pdotsigma[0] = (  p_perpMu[0]*sigma[0] - p_perpMu[1]*sigma[1]
                        - p_perpMu[2]*sigma[2] - p_perpMu[3]*sigma[3]);
        pdotsigma[1] = (  p_perpMu[0]*sigma[1] - p_perpMu[1]*sigma[4]
                        - p_perpMu[2]*sigma[5] - p_perpMu[3]*sigma[6]);
        pdotsigma[2] = (  p_perpMu[0]*sigma[2] - p_perpMu[1]*sigma[5]
                        - p_perpMu[2]*sigma[7] - p_perpMu[3]*sigma[8]);
        pdotsigma[3] = (  p_perpMu[0]*sigma[3] - p_perpMu[1]*sigma[6]
                        - p_perpMu[2]*sigma[8] - p_perpMu[3]*sigma[9]);

        Smu_SIP2[0] = - prefactorSIP*(
              (ux*pmu[2] - uy*pmu[1])*pdotsigma[3]
            + (uy*pmu[3] - uz*pmu[2])*pdotsigma[1]
            - (ux*pmu[3] - uz*pmu[1])*pdotsigma[2]
        );
        Smu_SIP2[1] = prefactorSIP*(
            - (uy*pmu[3] - uz*pmu[2])*pdotsigma[0]
            - (ut*pmu[2] - uy*pmu[0])*pdotsigma[3]
            + (ut*pmu[3] - uz*pmu[0])*pdotsigma[2]
        );
        Smu_SIP2[2] = prefactorSIP*(
              (ut*pmu[1] - ux*pmu[0])*pdotsigma[3]
            - (ut*pmu[3] - uz*pmu[0])*pdotsigma[1]
            + (ux*pmu[3] - uz*pmu[1])*pdotsigma[0]
        );
        Smu_SIP2[3] = prefactorSIP*(
            - (ut*pmu[1] - ux*pmu[0])*pdotsigma[2]
            + (ut*pmu[2] - uy*pmu[0])*pdotsigma[1]
            - (ux*pmu[2] - uy*pmu[1])*pdotsigma[0]
        );

        for (int imu = 0; imu < 4; imu++) {
            SmuSIP2_tmp[imu] += Smu_th[imu] + Smu_SIP2[imu];
        }

        // muB induced Polarization
        float Smu_muIP[4] = {0., 0., 0., 0.};
        float DmuB_over_T[4];
        for (int j = 0; j < 4; j++)
            DmuB_over_T[j] = surf.vorticity_arr[34 + j];
        const double prefactorMuIP = prefactor*hbarC*POI_info.baryon/pdotu;
        Smu_muIP[0] = prefactorMuIP*(
              (ux*pmu[2] - uy*pmu[1])*DmuB_over_T[3]
            + (uy*pmu[3] - uz*pmu[2])*DmuB_over_T[1]
            - (ux*pmu[3] - uz*pmu[1])*DmuB_over_T[2]
        );
        Smu_muIP[1] = - prefactorMuIP*(
            - (uy*pmu[3] - uz*pmu[2])*DmuB_over_T[0]
            - (ut*pmu[2] - uy*pmu[0])*DmuB_over_T[3]
            + (ut*pmu[3] - uz*pmu[0])*DmuB_over_T[2]
        );
        Smu_muIP[2] = - prefactorMuIP*(
              (ut*pmu[1] - ux*pmu[0])*DmuB_over_T[3]
            - (ut*pmu[3] - uz*pmu[0])*DmuB_over_T[1]
            + (ux*pmu[3] - uz*pmu[1])*DmuB_over_T[0]
        );
        Smu_muIP[3] = - prefactorMuIP*(
            - (ut*pmu[1] - ux*pmu[0])*DmuB_over_T[2]
            + (ut*pmu[2] - uy*pmu[0])*DmuB_over_T[1]
            - (ux*pmu[2] - uy*pmu[1])*DmuB_over_T[0]
        );

        for (int imu = 0; imu < 4; imu++) {
            SmuFull_tmp[imu] += Smu_th[imu] + Smu_SIP2[imu] + Smu_muIP[imu];
        }
    }

    for (int i = 0; i < 4; i++) {
        Smu[i] = Smu_tmp[i];
        SmuSIP1[i] = SmuSIP1_tmp[i];
        SmuSIP2[i] = SmuSIP2_tmp[i];
        SmuFull[i] = SmuFull_tmp[i];
    }

    // transform to the particle LRF
    transformSmuToLRF(POI_info.mass, pmu, Smu_tmp, SmuLRF);
    transformSmuToLRF(POI_info.mass, pmu, SmuSIP1_tmp, SmuSIP1LRF);
    transformSmuToLRF(POI_info.mass, pmu, SmuSIP2_tmp, SmuSIP2LRF);
    transformSmuToLRF(POI_info.mass, pmu, SmuFull_tmp, SmuFullLRF);
}


void SpinPolarization::transformSmuToLRF(
        const double mass, const iSS_data::Vec4 &pmu,
        double Smu[], iSS_data::Vec4 &SmuLRF) const {
    // transform to the particle LRF
    double p_dot_S = 0.;
    for (int i = 1; i < 4; i++) {
        p_dot_S += pmu[i]*Smu[i];
    }
    SmuLRF[0] = pmu[0]/mass*Smu[0] - p_dot_S/mass;
    for (int i = 1; i < 4; i++) {
        SmuLRF[i] = Smu[i] - p_dot_S*pmu[i]/(pmu[0]*(pmu[0] + mass));
    }
}

