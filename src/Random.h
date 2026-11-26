// Copyright @ Chun Shen 2018

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <random>
#include <memory>
#include <cstdint>

namespace RandomUtil {

class Random {
 private:
    long seed;
    std::unique_ptr<std::mt19937> ran_generator;
    std::uniform_real_distribution<double> rand_uniform_dist;

 public:
    Random(long seed_in, double min = 0.0, double max = 1.0);
    //! An independent stream seeded from (key, stream) through
    //! std::seed_seq, e.g. one stream per event for parallel sampling.
    Random(std::uint32_t key_hi, std::uint32_t key_lo, std::uint32_t stream);
    double rand_uniform() {return(rand_uniform_dist(*ran_generator));}
    //! a raw 32-bit draw, used to derive the keys of sub-streams
    std::uint32_t rand_uint32() {
        return(static_cast<std::uint32_t>((*ran_generator)()));
    }
    int get_seed() const {return(seed);}
};

}

#endif  // SRC_RANDOM_H_
