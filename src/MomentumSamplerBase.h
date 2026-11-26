// Copyright @ Chun Shen 2020

#ifndef SRC_MomentumSamplerBase_
#define SRC_MomentumSamplerBase_

#include "Random.h"
#include <array>
#include <cmath>
#include <vector>
#include <memory>

class MomentumSamplerBase {
 protected:
    template <bool First, bool Second, bool AlternateFirst = false>
    static std::array<double, 2> evaluate_cdf_terms(
            double Etilde, double mtilde, int order) {
        std::array<double, 2> result = {0., 0.};
        double sign = 1.;
        for (int n = 0; n < order; ++n) {
            const int n1 = n + 1;
            const double scale = std::exp(-mtilde*n);
            const double tail = std::exp((mtilde - Etilde)*n1);
            if (First) {
                result[0] += (sign/(n1*n1)*scale*(
                    (mtilde*n1 + 1) - tail*(Etilde*n1 + 1)));
            }
            if (Second) {
                result[1] += (1./(n1*n1*n1)*scale*(
                      (mtilde*n1*(mtilde*n1 + 2) + 2)
                    - tail*(Etilde*n1*(Etilde*n1 + 2) + 2)));
            }
            // The fermion first term alternates. Its second term has unit sign.
            if (AlternateFirst) sign *= -1.;
        }
        return result;
    }

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
    // Each inverse CDF probe uses all three terms for the same energy bin.
    std::vector<std::array<double, 3>> cdfTable_;

    MomentumSamplerBase();
    virtual ~MomentumSamplerBase() {};

    virtual double CDF_0(const double Etilde) const {return(0.0*Etilde);}
    virtual double CDF_1(const double Etilde) const {return(0.0*Etilde);}
    virtual double CDF_2(const double Etilde) const {return(0.0*Etilde);}

    virtual std::array<double, 2> CDF_12(double Etilde) const {
        return {{CDF_1(Etilde), CDF_2(Etilde)}};
    }

    void update_cache(const double m_minus_mu_tilde);

    double Sample_a_momentum(double m, double T, double mu);

    double inverse_CDF(const double r, const int idx_0,
                       const double CDF_max, const double m_term,
                       const double m_tilde, const double mu_tilde) const;

};

#endif  // SRC_MomentumSamplerBase_
