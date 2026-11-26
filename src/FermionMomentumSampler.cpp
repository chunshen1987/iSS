// Copyright @ Chun Shen 2020

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "FermionMomentumSampler.h"

FermionMomentumSampler::FermionMomentumSampler(
                std::shared_ptr<RandomUtil::Random> ran_gen,
                double m0, int trunc_order) {
    m0_ = m0;
    trunc_order_ = trunc_order;
    ran_gen_ptr_ = ran_gen;
    update_cache(1.0);
    double E_tilde_min = m0_;
    double E_tilde_max = E_tilde_min + 40.;
    double dE_tilde = 0.02;
    int npoints = (E_tilde_max - E_tilde_min)/dE_tilde + 1;
    Etilde_.resize(npoints);
    CDF_0_.resize(npoints);
    CDF_1_.resize(npoints);
    CDF_2_.resize(npoints);
    tb_len_ = npoints;
    for (int i = 0; i < npoints; i++) {
        const double E_tilde_local = E_tilde_min + i*dE_tilde;
        Etilde_[i] = E_tilde_local;
        CDF_0_[i] = CDF_0(E_tilde_local);
        CDF_1_[i] = CDF_1(E_tilde_local);
        CDF_2_[i] = CDF_2(E_tilde_local);
    }
}


double FermionMomentumSampler::CDF_0(const double Etilde) const {
    return(CDF_0(Etilde, m0_));
}


double FermionMomentumSampler::CDF_0(const double Etilde,
                                     const double mtilde) const {
    double res = 0.;
    if (trunc_order_ > 5) {
        res = - exp(mtilde)*log((1. + exp(-Etilde))/(1. + exp(-mtilde)));
    } else {
        int sign = 1;
        for (int n = 0; n < trunc_order_; n++) {
            int n1 = n + 1;
            res += (sign/n1)*exp(-mtilde*n)*(1. - exp((mtilde - Etilde)*n1));
            sign *= -1;
        }
    }
    return(res);
}


double FermionMomentumSampler::CDF_1(const double Etilde) const {
    return(CDF_1(Etilde, m0_));
}


double FermionMomentumSampler::CDF_1(const double Etilde,
                                     const double mtilde) const {
    double res = 0.;
    double sign = 1.;
    for (int n = 0; n < trunc_order_; n++) {
        int n1 = n + 1;
        res += (sign/(n1*n1)*exp(-mtilde*n)*(
            (mtilde*n1 + 1) - exp((mtilde - Etilde)*n1)*(Etilde*n1 + 1)));
        sign *= -1.;
    }
    return(res);
}


double FermionMomentumSampler::CDF_2(const double Etilde) const {
    return(CDF_2(Etilde, m0_));
}


double FermionMomentumSampler::CDF_2(const double Etilde,
                                     const double mtilde) const {
    double res = 0.;
    for (int n = 0; n < trunc_order_; n++) {
        int n1 = n + 1;
        res += (1./(n1*n1*n1)*exp(-mtilde*n)*(
              (mtilde*n1*(mtilde*n1 + 2) + 2)
            - exp((mtilde - Etilde)*n1)*(Etilde*n1*(Etilde*n1 + 2) + 2))
        );
    }
    return(res);
}
