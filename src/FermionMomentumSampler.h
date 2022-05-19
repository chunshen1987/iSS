// Copyright @ Chun Shen 2020

#ifndef SRC_FermionMomentumSampler_
#define SRC_FermionMomentumSampler_

#include "MomentumSamplerBase.h"

class FermionMomentumSampler : public MomentumSamplerBase {
 private:
    double m0_;
    int trunc_order_;

 public:
    FermionMomentumSampler(std::shared_ptr<RandomUtil::Random> ran_gen,
                           double m0, int trunc_order);
    ~FermionMomentumSampler() {};

    double CDF_0(const double Etilde) const;
    double CDF_1(const double Etilde) const;
    double CDF_2(const double Etilde) const;
    double CDF_0(const double Etilde, const double mtilde) const;
    double CDF_1(const double Etilde, const double mtilde) const;
    double CDF_2(const double Etilde, const double mtilde) const;
};

#endif  // SRC_FermionMomentumSampler_
