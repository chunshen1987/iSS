// Copyright @ Chun Shen 2020

#ifndef SRC_MomentumSamplerShell_
#define SRC_MomentumSamplerShell_

#include <memory>
#include "Random.h"
#include "BosonMomentumSampler.h"
#include "FermionMomentumSampler.h"

class MomentumSamplerShell {
 private:
     const double mtilde_1_;
     const double mtilde_2_;
     std::shared_ptr<FermionMomentumSampler> fermion_sampler_0_;
     std::shared_ptr<FermionMomentumSampler> fermion_sampler_1_;
     std::shared_ptr<FermionMomentumSampler> fermion_sampler_2_;
     std::shared_ptr<BosonMomentumSampler> boson_sampler_0_;
     std::shared_ptr<BosonMomentumSampler> boson_sampler_1_;
     std::shared_ptr<BosonMomentumSampler> boson_sampler_2_;

 public:
    MomentumSamplerShell(std::shared_ptr<RandomUtil::Random> ran_gen_ptr);
    ~MomentumSamplerShell() {};

    double Sample_a_momentum(double m, double T, double mu, int sign) const;

};

#endif  // SRC_MomentumSamplerShell_
