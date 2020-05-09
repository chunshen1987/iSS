// Copyright @ 2020 Chun Shen

#include "spin_polarization.h"
#include "arsenal.h"
#include "data_struct.h"
#include<iostream>
#include<sstream>
#include<cstdlib>
#include<fstream>
#include<iomanip>

using std::cout;
using std::endl;
using std::scientific;
using std::setprecision;
using std::setw;


SpinPolarization::SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                                   const std::vector<particle_info> &particles,
                                   std::string path, std::string table_path,
                                   ParameterReader &paraRdr) :
        path_(path), table_path_(table_path),
        FOsurf_ptr_(FOsurf_ptr), particle_info_(particles),
        paraRdr_(paraRdr) {

    double dpT = 3.0/(NpT_ - 1);
    pT_arr_.resize(NpT_);
    for (int i = 0; i < NpT_; i++)
        pT_arr_[i] = 0.0 + i*dpT;
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

    Smu_pT_.resize(NpT_, {0., 0., 0., 0.});
    Smu_phi_.resize(Nphi_, {0., 0., 0., 0.});
    Smu_y_.resize(Ny_, {0., 0., 0., 0.});

    for (int i = 0; i < NpT_; i++) {
        std::vector<iSS_data::Vec4> Smu(Nphi_, {0., 0., 0., 0.});
        Smu_pTdpTdphi_.push_back(Smu);
    }

    dN_pTdpTdphidy_ = create_a_3D_Matrix(Ny_, NpT_, Nphi_, 0.);
    St_pTdpTdphidy_ = create_a_3D_Matrix(Ny_, NpT_, Nphi_, 0.);
    Sx_pTdpTdphidy_ = create_a_3D_Matrix(Ny_, NpT_, Nphi_, 0.);
    Sy_pTdpTdphidy_ = create_a_3D_Matrix(Ny_, NpT_, Nphi_, 0.);
    Sz_pTdpTdphidy_ = create_a_3D_Matrix(Ny_, NpT_, Nphi_, 0.);

    vorticity_typenames_.push_back("KineticSP");
    vorticity_typenames_.push_back("Kinetic");
    vorticity_typenames_.push_back("Thermal");
    vorticity_typenames_.push_back("Temperature");
}


SpinPolarization::~SpinPolarization() {
    delete_a_3D_Matrix(dN_pTdpTdphidy_, Ny_, NpT_);
    delete_a_3D_Matrix(St_pTdpTdphidy_, Ny_, NpT_);
    delete_a_3D_Matrix(Sx_pTdpTdphidy_, Ny_, NpT_);
    delete_a_3D_Matrix(Sy_pTdpTdphidy_, Ny_, NpT_);
    delete_a_3D_Matrix(Sz_pTdpTdphidy_, Ny_, NpT_);
}


void SpinPolarization::compute_spin_polarization_shell() {
    const int POI_monval = 3122;  // Lambda
    for (unsigned int itype = 0; itype < vorticity_typenames_.size(); itype++) {
        compute_spin_polarization(POI_monval, itype);
        output_integrated_spin_polarizations(POI_monval,
                                             vorticity_typenames_[itype]);
    }
}


void SpinPolarization::compute_spin_polarization(const int POI_monval,
                                                 const int itype) {
    // first clean up previous results
    set_val_in_3D_Matrix(dN_pTdpTdphidy_, Ny_, NpT_, Nphi_, 0.0);
    set_val_in_3D_Matrix(St_pTdpTdphidy_, Ny_, NpT_, Nphi_, 0.0);
    set_val_in_3D_Matrix(Sx_pTdpTdphidy_, Ny_, NpT_, Nphi_, 0.0);
    set_val_in_3D_Matrix(Sy_pTdpTdphidy_, Ny_, NpT_, Nphi_, 0.0);
    set_val_in_3D_Matrix(Sz_pTdpTdphidy_, Ny_, NpT_, Nphi_, 0.0);

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
    cout << "spin polarization tensor type : " << vorticity_typenames_[itype]
         << endl;

    const double mass  = POI_info.mass;
    const double prefactor = -1./(8.*mass);

    for (int iy = 0; iy < Ny_; iy++) {
        cout << "progress: " << iy << "/" << Ny_ << endl;
        double y = y_arr_[iy];
        double cosh_y = cosh(y);
        double sinh_y = sinh(y);
        for (int ipT = 0; ipT < NpT_; ipT++) {
            double p_perp = pT_arr_[ipT];
            double m_perp = sqrt(mass*mass + p_perp*p_perp);
            double pt = m_perp*cosh_y;
            double pz = m_perp*sinh_y;
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                double px = p_perp*cos_phi_arr_[iphi];
                double py = p_perp*sin_phi_arr_[iphi];

                iSS_data::Vec4 pmu = {pt, px, py, pz};
                iSS_data::Vec4 Smu = {0., 0., 0., 0.};
                double dN = 0.;
                compute_spin_polarization_for_a_given_p(POI_info, pmu, itype,
                                                        Smu, dN);
                dN_pTdpTdphidy_[iy][ipT][iphi] = dN;
                St_pTdpTdphidy_[iy][ipT][iphi] = prefactor*Smu[0];
                Sx_pTdpTdphidy_[iy][ipT][iphi] = prefactor*Smu[1];
                Sy_pTdpTdphidy_[iy][ipT][iphi] = prefactor*Smu[2];
                Sz_pTdpTdphidy_[iy][ipT][iphi] = prefactor*Smu[3];
            }
        }
    }
    compute_integrated_spin_polarizations();
}


void SpinPolarization::compute_integrated_spin_polarizations() {
    cout << "computing the integrated spin polarization ... " << endl;

    // clean up previous results
    Smu_ = {0.};
    for (int ipT = 0; ipT < NpT_; ipT++)
        Smu_pT_[ipT] = {0.};
    for (int iphi = 0; iphi < Nphi_; iphi++)
        Smu_phi_[iphi] = {0.};
    for (int iy = 0; iy < Ny_; iy++)
        Smu_y_[iy] = {0.};
    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iphi = 0; iphi < Nphi_; iphi++) {
            for (int i = 0; i < 4; i++)
                Smu_pTdpTdphi_[ipT][iphi][i] = 0.0;
        }
    }

    // compute S^mu(pT)
    double dN = 0.;
    for (int iy = 0; iy < Ny_; iy++) {
        if (std::abs(y_arr_[iy]) > 1.) continue;
        for (int ipT = 0; ipT < NpT_; ipT++) {
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                dN      += dN_pTdpTdphidy_[iy][ipT][iphi];
                Smu_[0] += St_pTdpTdphidy_[iy][ipT][iphi];
                Smu_[1] += Sx_pTdpTdphidy_[iy][ipT][iphi];
                Smu_[2] += Sy_pTdpTdphidy_[iy][ipT][iphi];
                Smu_[3] += Sz_pTdpTdphidy_[iy][ipT][iphi];
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        Smu_[i] /= dN;
    }

    // compute S^mu(pT)
    for (int ipT = 0; ipT < NpT_; ipT++) {
        double dN_pT = 0.;
        for (int iy = 0; iy < Ny_; iy++) {
            if (std::abs(y_arr_[iy]) > 1.) continue;
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                dN_pT           += dN_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pT_[ipT][0] += St_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pT_[ipT][1] += Sx_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pT_[ipT][2] += Sy_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pT_[ipT][3] += Sz_pTdpTdphidy_[iy][ipT][iphi];
            }
        }
        for (int i = 0; i < 4; i++) {
            Smu_pT_[ipT][i] /= dN_pT;
        }
    }

    // compute S^mu(phi)
    for (int iphi = 0; iphi < Nphi_; iphi++) {
        double dN_phi = 0.;
        for (int iy = 0; iy < Ny_; iy++) {
            if (std::abs(y_arr_[iy]) > 1.) continue;
            for (int ipT = 0; ipT < NpT_; ipT++) {
                dN_phi            += dN_pTdpTdphidy_[iy][ipT][iphi];
                Smu_phi_[iphi][0] += St_pTdpTdphidy_[iy][ipT][iphi];
                Smu_phi_[iphi][1] += Sx_pTdpTdphidy_[iy][ipT][iphi];
                Smu_phi_[iphi][2] += Sy_pTdpTdphidy_[iy][ipT][iphi];
                Smu_phi_[iphi][3] += Sz_pTdpTdphidy_[iy][ipT][iphi];
            }
        }
        cout << "dN_phi = " << dN_phi << endl;
        for (int i = 0; i < 4; i++) {
            Smu_phi_[iphi][i] /= dN_phi;
        }
    }

    // compute S^mu(y)
    for (int iy = 0; iy < Ny_; iy++) {
        double dN_y = 0.;
        for (int ipT = 0; ipT < NpT_; ipT++) {
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                dN_y          += dN_pTdpTdphidy_[iy][ipT][iphi];
                Smu_y_[iy][0] += St_pTdpTdphidy_[iy][ipT][iphi];
                Smu_y_[iy][1] += Sx_pTdpTdphidy_[iy][ipT][iphi];
                Smu_y_[iy][2] += Sy_pTdpTdphidy_[iy][ipT][iphi];
                Smu_y_[iy][3] += Sz_pTdpTdphidy_[iy][ipT][iphi];
            }
        }
        for (int i = 0; i < 4; i++) {
            Smu_y_[iy][i] /= dN_y;
        }
    }

    // compute S^mu(pT, phi)
    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iphi = 0; iphi < Nphi_; iphi++) {
            double dN_pTdpTdphi = 0.;
            for (int iy = 0; iy < Ny_; iy++) {
                if (std::abs(y_arr_[iy]) > 1.) continue;
                dN_pTdpTdphi                 += dN_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pTdpTdphi_[ipT][iphi][0] += St_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pTdpTdphi_[ipT][iphi][1] += Sx_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pTdpTdphi_[ipT][iphi][2] += Sy_pTdpTdphidy_[iy][ipT][iphi];
                Smu_pTdpTdphi_[ipT][iphi][3] += Sz_pTdpTdphidy_[iy][ipT][iphi];
            }
            for (int i = 0; i < 4; i++) {
                Smu_pTdpTdphi_[ipT][iphi][i] /= dN_pTdpTdphi;
            }
        }
    }
}


void SpinPolarization::output_integrated_spin_polarizations(
        const int POI_monval, const std::string vorticity_typename) {
    cout << "output spin polarization results to files ... " << endl;
    std::ofstream of;

    std::stringstream Smu_filename;
    Smu_filename << path_ << "/Smu_" << vorticity_typename << "_"
                 << POI_monval << ".dat";
    remove(Smu_filename.str().c_str());
    of.open(Smu_filename.str().c_str(), std::ios::out);
    of << "# S^t(pT)  S^x(pT)  S^y(pT)  S^z(pT)" << endl;
    for (int i = 0; i < 4; i++)
        of << scientific << setw(10) << setprecision(6) << Smu_[i] << "  ";
    of << endl;
    of.close();

    std::stringstream SmupT_filename;
    SmupT_filename << path_ << "/Smu_pT_" << vorticity_typename << "_"
                   << POI_monval << ".dat";
    remove(SmupT_filename.str().c_str());
    of.open(SmupT_filename.str().c_str(), std::ios::out);
    of << "# pT[GeV]  S^t(pT)  S^x(pT)  S^y(pT)  S^z(pT)" << endl;
    for (int ipT = 0; ipT < NpT_; ipT++) {
        of << scientific << setw(10) << setprecision(6)
           << pT_arr_[ipT] << "  ";
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << Smu_pT_[ipT][i] << "  ";
        }
        of << endl;
    }
    of.close();

    std::stringstream Smuphi_filename;
    Smuphi_filename << path_ << "/Smu_phi_" << vorticity_typename << "_"
                    << POI_monval << ".dat";
    remove(Smuphi_filename.str().c_str());
    of.open(Smuphi_filename.str().c_str(), std::ios::out);
    of << "# phi  S^t(phi)  S^x(phi)  S^y(phi)  S^z(phi)" << endl;
    for (int iphi = 0; iphi < Nphi_; iphi++) {
        of << scientific << setw(10) << setprecision(6)
           << phi_arr_[iphi] << "  ";
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << Smu_phi_[iphi][i] << "  ";
        }
        of << endl;
    }
    of.close();

    std::stringstream Smuy_filename;
    Smuy_filename << path_ << "/Smu_y_" << vorticity_typename << "_"
                  << POI_monval << ".dat";
    remove(Smuy_filename.str().c_str());
    of.open(Smuy_filename.str().c_str(), std::ios::out);
    of << "# y  S^t(y)  S^x(y)  S^y(y)  S^z(y)" << endl;
    for (int iy = 0; iy < Ny_; iy++) {
        of << scientific << setw(10) << setprecision(6)
           << y_arr_[iy] << "  ";
        for (int i = 0; i < 4; i++) {
            of << scientific << setw(10) << setprecision(6)
               << Smu_y_[iy][i] << "  ";
        }
        of << endl;
    }
    of.close();

    std::stringstream Smu_dpTdphi_filename;
    Smu_dpTdphi_filename << path_ << "/Smu_dpTdphi_" << vorticity_typename
                         << "_" << POI_monval << ".dat";
    remove(Smu_dpTdphi_filename.str().c_str());
    of.open(Smu_dpTdphi_filename.str().c_str(), std::ios::out);
    of << "# pT[GeV]  phi  S^t(pT, phi)  S^x(pT, phi)  S^y(pT, phi)  "
       << "S^z(pT, phi)" << endl;
    for (int ipT = 0; ipT < NpT_; ipT++) {
        for (int iphi = 0; iphi < Nphi_; iphi++) {
            of << scientific << setw(10) << setprecision(6)
               << pT_arr_[ipT] << "  " << phi_arr_[iphi] << "  ";
            for (int i = 0; i < 4; i++) {
                of << scientific << setw(10) << setprecision(6)
                   << Smu_pTdpTdphi_[ipT][iphi][i] << "  ";
            }
            of << endl;
        }
    }
    of.close();

    std::stringstream Smu_dpTdphidy_filename;
    Smu_dpTdphidy_filename << path_ << "/Smu_dpTdphidy_" << vorticity_typename
                           << "_" << POI_monval << ".dat";
    remove(Smu_dpTdphidy_filename.str().c_str());
    of.open(Smu_dpTdphidy_filename.str().c_str(), std::ios::out);
    of << "# y  pT[GeV]  phi  S^t  S^x  S^y  S^z" << endl;
    for (int iy = 0; iy < Ny_; iy++) {
        for (int ipT = 0; ipT < NpT_; ipT++) {
            for (int iphi = 0; iphi < Nphi_; iphi++) {
                of << scientific << setw(10) << setprecision(6)
                   << y_arr_[iy] << "  " << pT_arr_[ipT] << "  "
                   << phi_arr_[iphi] << "  "
                   << St_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
                   << Sx_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
                   << Sy_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
                   << Sz_pTdpTdphidy_[iy][ipT][iphi]/dN_pTdpTdphidy_[iy][ipT][iphi] << "  "
                   << endl;
            }
        }
    }
    of.close();
}


void SpinPolarization::compute_spin_polarization_for_a_given_p(
        const particle_info &POI_info, const iSS_data::Vec4 &pmu,
        const int itype, iSS_data::Vec4 &Smu, double &dN) {
    double Smu_tmp[4] = {0., 0., 0., 0.};
    #pragma omp parallel for reduction(+: Smu_tmp[:4], dN)
    for (unsigned int i = 0; i < FOsurf_ptr_.size(); i++) {
        const FO_surf &surf = FOsurf_ptr_[i];
        const float tau = surf.tau;
        const float ptau = pmu[0]*surf.cosh_eta - pmu[3]*surf.sinh_eta;
        const float tau_peta = pmu[3]*surf.cosh_eta - pmu[0]*surf.sinh_eta;   // tau*p^eta
        const double mu = (  POI_info.baryon*surf.muB
                           + POI_info.strange*surf.muS
                           + POI_info.charge*surf.muC);
        const double pdotu = (  ptau*surf.u0 - pmu[1]*surf.u1
                              - pmu[2]*surf.u2 - tau_peta*surf.u3);
        const double pdsigma = tau*(  ptau*surf.da0 + pmu[1]*surf.da1
                                    + pmu[2]*surf.da2 + tau_peta*surf.da3/tau);
        const double expon = (pdotu - mu)/surf.Tdec;
        const double f0 = 1./(exp(expon) + POI_info.sign);

        const float omega_tx = surf.vorticity_arr[6*itype + 0];
        const float omega_ty = surf.vorticity_arr[6*itype + 1];
        const float omega_tz = surf.vorticity_arr[6*itype + 2];
        const float omega_xy = surf.vorticity_arr[6*itype + 3];
        const float omega_xz = surf.vorticity_arr[6*itype + 4];
        const float omega_yz = surf.vorticity_arr[6*itype + 5];

        const double prefactor = pdsigma*f0*(1. - f0)*2.;
        dN     += pdsigma*f0;
        Smu_tmp[0] += prefactor*(- omega_yz*pmu[1] + omega_xz*pmu[2]
                                 - omega_xy*pmu[3]);
        Smu_tmp[1] += prefactor*(- omega_ty*pmu[3] + omega_tz*pmu[2]
                                 - omega_yz*pmu[0]);
        Smu_tmp[2] += prefactor*(- omega_tz*pmu[1] + omega_tx*pmu[3]
                                 - omega_xz*pmu[0]);
        Smu_tmp[3] += prefactor*(omega_ty*pmu[1] - omega_tx*pmu[2]
                                 - omega_xy*pmu[0]);
    }

    // transform to the particle LRF
    iSS_data::Vec4 Smu_LRF;
    const double mass = POI_info.mass;
    double p_dot_S = 0.;
    for (int i = 1; i < 4; i++)
        p_dot_S += pmu[i]*Smu_tmp[i];
    Smu_LRF[0] = pmu[0]/mass*Smu_tmp[0] - p_dot_S/mass;
    for (int i = 1; i < 4; i++)
        Smu_LRF[i] = Smu_tmp[i] - p_dot_S*pmu[i]/(pmu[0]*(pmu[0] + mass));

    Smu = Smu_LRF;
}
