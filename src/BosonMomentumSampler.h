// Copyright @ Chun Shen 2020

#ifndef SRC_BosonMomentumSampler_
#define SRC_BosonMomentumSampler_

#include "MomentumSamplerBase.h"

class BosonMomentumSampler : public MomentumSamplerBase {
  private:
    double m0_;
    int trunc_order_;

  public:
    BosonMomentumSampler(
        std::shared_ptr<RandomUtil::Random> ran_gen, double m0,
        int trunc_order);
    ~BosonMomentumSampler() {};

    double CDF_0(const double Etilde) const;
    double CDF_1(const double Etilde) const;
    double CDF_2(const double Etilde) const;
    double CDF_0(const double Etilde, const double mtilde) const;
    double CDF_1(const double Etilde, const double mtilde) const;
    double CDF_2(const double Etilde, const double mtilde) const;
};

#endif  // SRC_BosonMomentumSampler_
