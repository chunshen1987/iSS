// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "Random.h"

TEST_CASE("Test fixed random seed") {
    RandomUtil::Random ran_test1(1);
    RandomUtil::Random ran_test2(1);
    CHECK(ran_test1.get_seed() == 1);
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() == ran_test2.rand_uniform());
}

TEST_CASE("Test clock random seed replay") {
    RandomUtil::Random clock_random(-1);
    const int seed = clock_random.get_seed();
    CHECK(seed > 0);
    CHECK(seed <= (1 << 30));
    RandomUtil::Random replay(seed);
    for (int i = 0; i < 10; i++)
        CHECK(clock_random.rand_uniform() == replay.rand_uniform());
}

TEST_CASE("Test random number range") {
    RandomUtil::Random ran_test1(-1, 1., 2.);
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() >= 1.);
}
