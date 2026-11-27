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

TEST_CASE("Event streams are reproducible and use every seed component") {
    RandomUtil::Random reference(11u, 29u, 3u);
    RandomUtil::Random same_event(11u, 29u, 3u);
    RandomUtil::Random next_event(11u, 29u, 4u);
    RandomUtil::Random different_high_key(12u, 29u, 3u);
    RandomUtil::Random different_low_key(11u, 30u, 3u);
    bool event_differs = false;
    bool high_key_differs = false;
    bool low_key_differs = false;
    for (int i = 0; i < 16; ++i) {
        const auto draw = reference.rand_uint32();
        CHECK(same_event.rand_uint32() == draw);
        event_differs |= (next_event.rand_uint32() != draw);
        high_key_differs |= (different_high_key.rand_uint32() != draw);
        low_key_differs |= (different_low_key.rand_uint32() != draw);
    }
    CHECK(event_differs);
    CHECK(high_key_differs);
    CHECK(low_key_differs);
}

TEST_CASE("Test random number range") {
    RandomUtil::Random ran_test1(7L, 1., 2.);
    for (int i = 0; i < 100; ++i) {
        const double draw = ran_test1.rand_uniform();
        CHECK(draw >= 1.);
        CHECK(draw < 2.);
    }
}
