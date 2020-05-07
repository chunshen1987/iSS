// Copyright @ 2020 Chun Shen

#ifndef SRC_SPIN_POLARIZATION_H_
#define SRC_SPIN_POLARIZATION_H_

#include<vector>

#include "data_struct.h"
#include "Table.h"

class SpinPolarization {
 private:
    const std::string path_;
    const std::string table_path_;
    const std::vector<FO_surf> &FOsurf_ptr_;
    const std::vector<particle_info> particle_info_;

    const int NpT  = 15;
    const int Nphi = 48;
    const int Ny   = 30;

    double ***Sx_pTdpTdphidy_;
    double ***Sy_pTdpTdphidy_;
    double ***Sz_pTdpTdphidy_;

 public:
    SpinPolarization(const std::vector<FO_surf> &FOsurf_ptr,
                     const std::vector<particle_info> &particles,
                     std::string path, std::string table_path);
    ~SpinPolarization();

    void compute_spin_polarization();

};

#endif  // SRC_SPIN_POLARIZATION_H_
