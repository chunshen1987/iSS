// Copyright @ 2020 Chun Shen

#include "spin_polarization.h"
#include "arsenal.h"

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

}

