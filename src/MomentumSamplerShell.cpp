// Copyright @ Chun Shen 2020

#include <memory>
#include "MomentumSamplerShell.h"

MomentumSamplerShell::MomentumSamplerShell(
                    std::shared_ptr<RandomUtil::Random> ran_gen_ptr) : 
            mtilde_1_(30.), mtilde_2_(50.) {
    fermion_sampler_0_ = std::shared_ptr<FermionMomentumSampler>(
            new FermionMomentumSampler(ran_gen_ptr, 0., 10));
    fermion_sampler_1_ = std::shared_ptr<FermionMomentumSampler>(
            new FermionMomentumSampler(ran_gen_ptr, mtilde_1_., 2));
    fermion_sampler_2_ = std::shared_ptr<FermionMomentumSampler>(
            new FermionMomentumSampler(ran_gen_ptr, mtilde_2_, 1));
    boson_sampler_0_ = std::shared_ptr<BosonMomentumSampler>(
            new BosonMomentumSampler(ran_gen_ptr, 0.05, 10));
    boson_sampler_1_ = std::shared_ptr<BosonMomentumSampler>(
            new BosonMomentumSampler(ran_gen_ptr, mtilde_1_., 2));
    boson_sampler_2_ = std::shared_ptr<BosonMomentumSampler>(
            new BosonMomentumSampler(ran_gen_ptr, mtilde_2_, 1));
}


double MomentumSamplerShell::Sample_a_momentum(double m, double T,
                                               double mu, int sign) const {
    const double m0tilde = m/T - mu/T;
    double p_sample = 0.;
    if (sign == -1) {
        // Bosons
        if (m0tilde < mtilde_1_) {
            p_sample = boson_sampler_0_.Sample_a_momentum(m, T, mu);
        } else if (m0tilde < mtilde_2_) {
            p_sample = boson_sampler_1_.Sample_a_momentum(m, T, mu);
        } else {
            p_sample = boson_sampler_2_.Sample_a_momentum(m, T, mu);
        }
    } else {
        // Fermions
        if (m0tilde < mtilde_1_) {
            p_sample = fermion_sampler_0_.Sample_a_momentum(m, T, mu);
        } else if (m0tilde < mtilde_2_) {
            p_sample = fermion_sampler_1_.Sample_a_momentum(m, T, mu);
        } else {
            p_sample = fermion_sampler_2_.Sample_a_momentum(m, T, mu);
        }
    }
    return(p_sample);
}
