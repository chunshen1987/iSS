#include <cmath>
#include <memory>

#include "MomentumSamplerShell.h"
#include "doctest.h"

namespace {

double integrate_momentum_density(double upper, double mass, double temperature,
                                  double chemical_potential, int sign) {
    // Midpoint quadrature of p^2 f(E) in the local rest frame. The reference
    // avoids the sampler's proposal distribution and inverse CDF table.
    const int intervals = 8192;
    const double dp = upper/intervals;
    double sum = 0.;
    for (int i = 0; i < intervals; ++i) {
        const double p = (i + 0.5)*dp;
        const double energy = std::sqrt(p*p + mass*mass);
        const double boltzmann = std::exp((chemical_potential - energy)
                                          /temperature);
        sum += p*p*boltzmann/(1. + sign*boltzmann);
    }
    return dp*sum;
}

void check_thermal_distribution(double mass, double mu, int sign) {
    const int samples = 32768;
    const double temperature = 0.15;
    const double cuts[] = {0.30, 0.75};
    auto rng = std::make_shared<RandomUtil::Random>(917L);
    MomentumSamplerShell sampler(rng);
    CAPTURE(mass);
    CAPTURE(mu);
    CAPTURE(sign);
    int counts[] = {0, 0};
    for (int i = 0; i < samples; ++i) {
        const double p = sampler.Sample_a_momentum(mass, temperature, mu, sign);
        REQUIRE(std::isfinite(p));
        REQUIRE(p >= 0.);
        counts[0] += (p < cuts[0]);
        counts[1] += (p < cuts[1]);
    }

    // At 8 GeV the occupation is suppressed by more than exp(-50).
    const double norm = integrate_momentum_density(8., mass, temperature, mu, sign);
    for (int i = 0; i < 2; ++i) {
        CAPTURE(cuts[i]);
        const double expected = integrate_momentum_density(
            cuts[i], mass, temperature, mu, sign)/norm;
        const double measured = static_cast<double>(counts[i])/samples;
        CAPTURE(expected);
        CAPTURE(measured);
        const double sigma = std::sqrt(expected*(1. - expected)/samples);
        // A binomial bound leaves room for Monte Carlo fluctuations.
        CHECK(std::abs(measured - expected) <= 6.*sigma);
    }
}

}  // namespace

TEST_CASE("Thermal momentum samples follow Bose and Fermi occupation laws") {
    check_thermal_distribution(0.14, -0.04, -1);
    check_thermal_distribution(0.94, 0.24, 1);
}

TEST_CASE("Light fermion momentum follows the Fermi occupation law") {
    check_thermal_distribution(0.14, -0.04, 1);
}

TEST_CASE("A momentum draw is independent of a nearby cached mass") {
    auto rng = std::make_shared<RandomUtil::Random>(17L);
    MomentumSamplerShell fresh(rng);
    MomentumSamplerShell warmed(rng);
    const double mass = 0.14;
    const double temperature = 0.15;
    const double mu = -0.04;
    for (const int sign : {-1, 1}) {
        CAPTURE(sign);
        RandomUtil::Random warmup_rng(83L);
        warmed.Sample_a_momentum(mass + 1.e-9, temperature, mu, sign,
                                 warmup_rng);
        RandomUtil::Random fresh_rng(71L);
        RandomUtil::Random warmed_rng(71L);
        const double expected = fresh.Sample_a_momentum(
            mass, temperature, mu, sign, fresh_rng);
        const double actual = warmed.Sample_a_momentum(
            mass, temperature, mu, sign, warmed_rng);
        CHECK(actual == expected);
        CHECK(warmed_rng.rand_uint32() == fresh_rng.rand_uint32());
    }
}
