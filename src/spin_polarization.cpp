// Copyright @ 2020 Chun Shen

#include "spin_polarization.h"
#include "arsenal.h"
#include<iostream>


using std::cout;
using std::endl;

SpinPolarization::SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                                   const std::vector<particle_info> &particles,
                                   std::string path, std::string table_path) :
        path_(path), table_path_(table_path),
        FOsurf_ptr_(FOsurf_ptr), particle_info_(particles) {

    Sx_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
    Sy_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
    Sz_pTdpTdphidy_ = create_a_3D_Matrix(NpT, Nphi, Ny, 0.);
}


SpinPolarization::~SpinPolarization() {
    delete_a_3D_Matrix(Sx_pTdpTdphidy_, NpT, Nphi);
    delete_a_3D_Matrix(Sy_pTdpTdphidy_, NpT, Nphi);
    delete_a_3D_Matrix(Sz_pTdpTdphidy_, NpT, Nphi);
}

void SpinPolarization::compute_spin_polarization() {
    int POI_monval = 3122;  // Lambda

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

}

