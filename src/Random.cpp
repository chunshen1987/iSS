// Copyright @ Chun Shen 2018

#include "Random.h"

namespace RandomUtil {

Random::Random(long seed_in, double min, double max) :
    rand_uniform_dist(min, max) {
    seed = seed_in;
    if (seed == -1) {
        std::random_device ran_dev;
        seed = ran_dev();
    }
    ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));
}


Random::Random(std::uint32_t key_hi, std::uint32_t key_lo,
               std::uint32_t stream) :
    seed(static_cast<long>(stream)), rand_uniform_dist(0., 1.) {
    std::seed_seq seq{key_hi, key_lo, stream};
    ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seq));
}

}
