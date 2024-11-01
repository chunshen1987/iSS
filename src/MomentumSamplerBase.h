// Copyright @ Chun Shen 2020

#ifndef SRC_MomentumSamplerBase_
#define SRC_MomentumSamplerBase_

#include <memory>
#include <vector>

#include "Random.h"

class MomentumSamplerBase {
  private:
  public:
    const double tol_;
    const double eps_;

    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;

    //! cache the previous value
    double m_minus_mu_tilde_;
    double CDF_0_mtilde_;
    double CDF_1_mtilde_;
    double CDF_2_mtilde_;

    int tb_len_;
    std::vector<double> Etilde_;
    std::vector<double> CDF_0_;
    std::vector<double> CDF_1_;
    std::vector<double> CDF_2_;

    MomentumSamplerBase();
    virtual ~MomentumSamplerBase() {};

    virtual double CDF_0(const double Etilde) const { return (0.0 * Etilde); }
    virtual double CDF_1(const double Etilde) const { return (0.0 * Etilde); }
    virtual double CDF_2(const double Etilde) const { return (0.0 * Etilde); }

    void update_cache(const double m_minus_mu_tilde);

    double Sample_a_momentum(double m, double T, double mu);

    double inverse_CDF(
        const double r, const int idx_0, const double CDF_max,
        const double m_term, const double m_tilde, const double mu_tilde) const;
};

#endif  // SRC_MomentumSamplerBase_
