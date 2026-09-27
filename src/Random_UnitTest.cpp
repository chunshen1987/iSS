// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include <array>
#include <vector>

#include "Random.h"

TEST_CASE("Test fixed random seed") {
    RandomUtil::Random ran_test1(1);
    RandomUtil::Random ran_test2(1);
    CHECK(ran_test1.get_seed() == 1);
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() == ran_test2.rand_uniform());
}

TEST_CASE("Test device random seed") {
    RandomUtil::Random ran_test1(-1);
    RandomUtil::Random ran_test2(-1);
    CHECK(ran_test1.get_seed() != ran_test2.get_seed());
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() != ran_test2.rand_uniform());
}

TEST_CASE("Test random number range") {
    RandomUtil::Random ran_test1(-1, 1., 2.);
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() >= 1.);
}

TEST_CASE("Philox4x32-10 known-answer vectors (Random123)") {
    using RandomUtil::philox4x32_10;
    const std::array<uint32_t, 4> r0 = philox4x32_10({{0, 0, 0, 0}}, {{0, 0}});
    CHECK(r0[0] == 0x6627e8d5u); CHECK(r0[1] == 0xe169c58du);
    CHECK(r0[2] == 0xbc57ac4cu); CHECK(r0[3] == 0x9b00dbd8u);
    const uint32_t f = 0xffffffffu;
    const std::array<uint32_t, 4> r1 = philox4x32_10({{f, f, f, f}}, {{f, f}});
    CHECK(r1[0] == 0x408f276du); CHECK(r1[1] == 0x41c83b0eu);
    CHECK(r1[2] == 0xa20bc7c6u); CHECK(r1[3] == 0x6d5451fdu);
    const std::array<uint32_t, 4> r2 = philox4x32_10(
        {{0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u}},
        {{0xa4093822u, 0x299f31d0u}});
    CHECK(r2[0] == 0xd16cfe09u); CHECK(r2[1] == 0x94fdccebu);
    CHECK(r2[2] == 0x5001e420u); CHECK(r2[3] == 0x24126ea1u);
}

TEST_CASE("A keyed stream depends on its address only") {
    RandomUtil::Random ran(7);
    RandomUtil::Random reference(7);
    const double mt0 = ran.rand_uniform();
    CHECK(mt0 == reference.rand_uniform());
    ran.set_stream(1, 2, 3, 4, 5);
    std::vector<double> first;
    for (int i = 0; i < 9; i++) first.push_back(ran.rand_uniform());
    ran.unset_stream();
    // the Mersenne Twister goes on where it was left
    CHECK(ran.rand_uniform() == reference.rand_uniform());
    ran.set_stream(1, 2, 3, 4, 6);          // another address
    const double other = ran.rand_uniform();
    ran.set_stream(1, 2, 3, 4, 5);          // the same address again
    for (int i = 0; i < 9; i++) {
        const double u = ran.rand_uniform();
        CHECK(u == first[i]);
        CHECK(u > 0.);
        CHECK(u < 1.);
    }
    CHECK(other != first[0]);
    RandomUtil::Random other_seed(8);       // the address, not the seed, decides
    other_seed.set_stream(1, 2, 3, 4, 5);
    CHECK(other_seed.rand_uniform() == first[0]);
}

TEST_CASE("Poisson by inversion") {
    const int n = 200000;
    for (double lambda : {0.001, 0.3, 2.5, 17., 480.}) {
        double sum = 0., sum2 = 0.;
        long previous = 0;
        bool monotonic = true;
        for (int i = 0; i < n; i++) {
            const double u = (i + 0.5)/n;
            const long k = RandomUtil::poisson_inverse(u, lambda);
            monotonic = monotonic && (k >= previous);
            previous = k;
            sum += k;
            sum2 += static_cast<double>(k)*k;
        }
        const double mean = sum/n;
        const double var = sum2/n - mean*mean;
        CHECK(monotonic);
        CHECK(mean == doctest::Approx(lambda).epsilon(1e-3));
        CHECK(var == doctest::Approx(lambda).epsilon(2e-2));
        // nearly equal means: almost always the same k for the same u
        int changed = 0;
        for (int i = 0; i < n; i++) {
            const double u = (i + 0.5)/n;
            changed += (RandomUtil::poisson_inverse(u, lambda)
                        != RandomUtil::poisson_inverse(u, lambda*(1. + 1e-6)));
        }
        CHECK(changed < 0.001*n);
    }
    CHECK(RandomUtil::poisson_inverse(0.5, 0.) == 0);
    CHECK(RandomUtil::poisson_inverse(1e-9, 0.01) == 0);
    CHECK(RandomUtil::poisson_inverse(1. - 1e-12, 0.01) >= 1);
}
