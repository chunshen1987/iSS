// Copyright @ Chun Shen 2020

#include "MomentumSamplerBase.h"

#include <stdlib.h>

#include <cmath>
#include <iomanip>
#include <iostream>

MomentumSamplerBase::MomentumSamplerBase() : tol_(1e-8), eps_(1e-16) {}

void MomentumSamplerBase::update_cache(const double m_minus_mu_tilde) {
    m_minus_mu_tilde_ = m_minus_mu_tilde;
    CDF_0_mtilde_ = CDF_0(m_minus_mu_tilde);
    CDF_1_mtilde_ = CDF_1(m_minus_mu_tilde);
    CDF_2_mtilde_ = CDF_2(m_minus_mu_tilde);
}

double MomentumSamplerBase::Sample_a_momentum(double m, double T, double mu) {
    T = std::max(eps_, T);
    const double m_tilde = m / T;
    const double mu_tilde = mu / T;
    if (std::abs((m_tilde - mu_tilde) - m_minus_mu_tilde_) > tol_)
        update_cache(m_tilde - mu_tilde);
    double m_term =
        (CDF_2_mtilde_ + 2. * mu_tilde * CDF_1_mtilde_
         + (mu_tilde * mu_tilde - m_tilde * m_tilde / 2.) * CDF_0_mtilde_);

    const int idx_max = tb_len_ - 1;
    double CDF_max =
        (CDF_2_[idx_max] + 2. * mu_tilde * CDF_1_[idx_max]
         + (mu_tilde * mu_tilde - m_tilde * m_tilde / 2.) * CDF_0_[idx_max]
         - m_term);

    const int idx_min = static_cast<int>(
        (m_tilde - mu_tilde - Etilde_[0]) / (Etilde_[1] - Etilde_[0]));
    if (idx_min < 0 || idx_min >= idx_max) {
        std::cout << "[Error: MomentumSampler] out of range, "
                  << "Etilde_min = " << Etilde_[0]
                  << ", Etilde_max = " << Etilde_[idx_max] << ", m = " << m
                  << " GeV, T = " << T << " GeV, mu = " << mu
                  << " GeV. mtilde - mutilde = " << m_tilde - mu_tilde
                  << std::endl;
        exit(1);
    }

    double p_sample = -1.;
    double E_sample = -1.;
    double accept_ratio = 1.;
    do {
        double r = ran_gen_ptr_->rand_uniform() * CDF_max;
        double Etilde_sample =
            inverse_CDF(r, idx_min, CDF_max, m_term, m_tilde, mu_tilde);
        // convert back to particle momentum
        E_sample = T * Etilde_sample + mu;
        p_sample = sqrt(E_sample * E_sample - m * m);
        accept_ratio =
            (p_sample / E_sample) / (1. - m * m / (2. * E_sample * E_sample));
    } while (ran_gen_ptr_->rand_uniform() > accept_ratio);
    return (p_sample);
}

double MomentumSamplerBase::inverse_CDF(
    const double r, const int idx_0, const double CDF_max, const double m_term,
    const double m_tilde, const double mu_tilde) const {
    int idx_min = idx_0;
    int idx_max = tb_len_ - 1;
    double r_min =
        (CDF_2_[idx_min] + 2. * mu_tilde * CDF_1_[idx_min]
         + (mu_tilde * mu_tilde - m_tilde * m_tilde / 2.) * CDF_0_[idx_min]
         - m_term);
    double r_max = CDF_max;
    while ((idx_max - idx_min) > 1) {
        int idx_mid = static_cast<int>((idx_max + idx_min) / 2 + 0.1);
        double r_mid =
            (CDF_2_[idx_mid] + 2. * mu_tilde * CDF_1_[idx_mid]
             + (mu_tilde * mu_tilde - m_tilde * m_tilde / 2.) * CDF_0_[idx_mid]
             - m_term);
        if (r < r_mid) {
            idx_max = idx_mid;
            r_max = r_mid;
        } else {
            idx_min = idx_mid;
            r_min = r_mid;
        }
    }
    double E0 = Etilde_[idx_min];
    if (E0 < m_tilde - mu_tilde) {
        E0 = m_tilde - mu_tilde;
        r_min = 0.;
    }
    double res =
        (E0
         + (Etilde_[idx_max] - E0) / std::max(eps_, (r_max - r_min))
               * (r - r_min));
    return (res);
}
