// Copyright @ 2020 Chun Shen

#include "spin_polarization.h"
#include "arsenal.h"
#include "data_struct.h"
#include<iostream>

using std::cout;
using std::endl;
using iSS_data::hbarC;


SpinPolarization::SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                                   const std::vector<particle_info> &particles,
                                   std::string path, std::string table_path,
                                   ParameterReader &paraRdr) :
        path_(path), table_path_(table_path),
        FOsurf_ptr_(FOsurf_ptr), particle_info_(particles),
        paraRdr_(paraRdr) {

    double dpT = 3.0/(NpT - 1);
    pT_arr.resize(NpT);
    for (int i = 0; i < NpT; i++)
        pT_arr[i] = 0.0 + i*dpT;
    double dphi = 2.*M_PI/Nphi;
    phi_arr.resize(Nphi);
    cos_phi_arr.resize(Nphi);
    sin_phi_arr.resize(Nphi);
    for (int i = 0; i < Nphi; i++) {
        phi_arr[i] = 0.0 + i*dphi;
        cos_phi_arr[i] = cos(phi_arr[i]);
        sin_phi_arr[i] = sin(phi_arr[i]);
    }
    double y_size = 10.0;
    double dy = -y_size/2./(Ny - 1);
    y_arr.resize(Ny);
    for (int i = 0; i < Ny; i++)
        y_arr[i] = - y_size/2. + i*dy;

    St_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
    Sx_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
    Sy_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
    Sz_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
}


SpinPolarization::~SpinPolarization() {
    delete_a_3D_Matrix(St_pTdpTdphidy_, NpT, Nphi);
    delete_a_3D_Matrix(Sx_pTdpTdphidy_, NpT, Nphi);
    delete_a_3D_Matrix(Sy_pTdpTdphidy_, NpT, Nphi);
    delete_a_3D_Matrix(Sz_pTdpTdphidy_, NpT, Nphi);
}


void SpinPolarization::compute_spin_polarization() {
    int POI_monval = 3122;  // Lambda

    set_val_in_3D_Matrix(St_pTdpTdphidy_, NpT, Nphi, Ny, 0.0);
    set_val_in_3D_Matrix(Sx_pTdpTdphidy_, NpT, Nphi, Ny, 0.0);
    set_val_in_3D_Matrix(Sy_pTdpTdphidy_, NpT, Nphi, Ny, 0.0);
    set_val_in_3D_Matrix(Sz_pTdpTdphidy_, NpT, Nphi, Ny, 0.0);

    particle_info POI_info;
    // find the information about the particle of interest
    for (const auto &particle_i: particle_info_) {
        if (particle_i.monval == POI_monval) {
            POI_info = particle_i;
            break;
        }
    }

    cout << "Computing spin polarization for " << POI_info.name
         << ", Monte-carlo index: " << POI_info.monval << endl;

    const double mass  = POI_info.mass;

    double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

    for (int iy = 0; iy < Ny; iy++) {
        double y = y_arr[iy];
        double cosh_y = cosh(y);
        double sinh_y = sinh(y);
        for (int ipT = 0; ipT < NpT; ipT++) {
            double pT = pT_arr[ipT];
            double mT = sqrt(mass*mass + pT*pT);
            double pt = mT*cosh_y;
            double pz = mT*sinh_y;
            for (int iphi = 0; iphi < Nphi; iphi++) {
                double px = pT*cos_phi_arr[iphi];
                double py = pT*sin_phi_arr[iphi];

                iSS_data::Vec4 pmu = {pt, px, py, pz};
                iSS_data::Vec4 Smu = {0., 0., 0., 0.};
                compute_spin_polarization_for_a_given_p(POI_info, pmu, Smu);
                St_pTdpTdphidy_[ipT][iphi][iy] = Smu[0];
                Sx_pTdpTdphidy_[ipT][iphi][iy] = Smu[1];
                Sy_pTdpTdphidy_[ipT][iphi][iy] = Smu[2];
                Sz_pTdpTdphidy_[ipT][iphi][iy] = Smu[3];
            }
        }
    }
}


void SpinPolarization::compute_spin_polarization_for_a_given_p(
        const particle_info &POI_info, const iSS_data::Vec4 &pmu,
        iSS_data::Vec4 &Smu) {
    double yield = 0.;
    for (unsigned int i = 0; i < FOsurf_ptr_.size(); i++) {
        const FO_surf &surf = FOsurf_ptr_[i];
        const float tau = surf.tau;
        const float ptau = pmu[0]*surf.cosh_eta - pmu[3]*surf.sinh_eta;
        const float tau_peta = pmu[3]*surf.cosh_eta - pmu[0]*surf.sinh_eta;   // tau*p^eta
        double mu = (  POI_info.baryon*surf.muB + POI_info.strange*surf.muS
                     + POI_info.charge*surf.muC);
        const double pdotu = (  ptau*surf.u0 - pmu[1]*surf.u1
                              - pmu[2]*surf.u2 - tau_peta*surf.u3);
        const double pdsigma = tau*(  ptau*surf.da0 + pmu[1]*surf.da1
                                    + pmu[2]*surf.da2 + tau_peta*surf.da3/tau);
        const double expon = (pdotu - mu)/surf.Tdec;
        const double f0 = 1./(exp(expon) + POI_info.sign);

        const float omega_tx = surf.vorticity_arr[0];
        const float omega_ty = surf.vorticity_arr[1];
        const float omega_tz = surf.vorticity_arr[2];
        const float omega_xy = surf.vorticity_arr[3];
        const float omega_xz = surf.vorticity_arr[4];
        const float omega_yz = surf.vorticity_arr[5];

        const double prefactor = pdsigma*f0*(1. - f0)*2.;
        yield  += pdsigma*f0;
        Smu[0] += prefactor*(- omega_yz*pmu[1] + omega_xz*pmu[2]
                             - omega_xy*pmu[3]);
        Smu[1] += prefactor*(- omega_ty*pmu[3] + omega_tz*pmu[2]
                             - omega_yz*pmu[0]);
        Smu[2] += prefactor*(- omega_tz*pmu[1] + omega_tx*pmu[3]
                             - omega_xz*pmu[0]);
        Smu[3] += prefactor*(omega_ty*pmu[1] - omega_tx*pmu[2]
                             - omega_xy*pmu[0]);
    }
    const double denorm = 1./(8.*POI_info.mass*yield);
    for (int i = 0; i < 4; i++)
        Smu[i] = - Smu[i]*denorm;
}
