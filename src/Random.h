// Copyright @ Chun Shen 2018

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <array>
#include <cstdint>
#include <random>
#include <memory>

#include "Philox.h"

namespace RandomUtil {

class Random {
 private:
    long seed;
    std::random_device ran_dev;
    std::unique_ptr<std::mt19937> ran_generator;
    std::uniform_real_distribution<double> rand_uniform_dist;

    // keyed stream (set_stream): Philox4x32-10 at (key, counter), the last
    // counter word counting the blocks of four numbers drawn
    bool stream_on_ = false;
    std::array<uint32_t, 2> stream_key_{{0, 0}};
    std::array<uint32_t, 4> stream_ctr_{{0, 0, 0, 0}};
    std::array<uint32_t, 4> stream_buf_{{0, 0, 0, 0}};
    int stream_used_ = 4;

    double stream_uniform() {
        if (stream_used_ == 4) {
            stream_buf_ = philox4x32_10(stream_ctr_, stream_key_);
            stream_ctr_[3]++;
            stream_used_ = 0;
        }
        const double u = philox_to_uniform(stream_buf_[stream_used_++]);
        return rand_uniform_dist.a() + (rand_uniform_dist.b()
                                        - rand_uniform_dist.a())*u;
    }

 public:
    Random(long seed_in, double min = 0.0, double max = 1.0);
    double rand_uniform() {
        if (stream_on_) return(stream_uniform());
        return(rand_uniform_dist(*ran_generator));
    }
    int get_seed() const {return(seed);}

    //! From now on rand_uniform() returns the numbers of the Philox stream
    //! addressed by (key0, key1; c0, c1, c2): the same address gives the same
    //! numbers, whatever was drawn before. The Mersenne Twister is untouched.
    void set_stream(uint32_t key0, uint32_t key1,
                    uint32_t c0, uint32_t c1, uint32_t c2) {
        stream_key_ = {{key0, key1}};
        stream_ctr_ = {{c0, c1, c2, 0}};
        stream_used_ = 4;
        stream_on_ = true;
    }
    //! back to the Mersenne Twister, where it was left
    void unset_stream() {stream_on_ = false;}
};

//! Poisson(lambda) by inversion: the smallest k with P(k; lambda) >= u. It is
//! monotonic in u and continuous in lambda, so one uniform gives the same k
//! for nearly equal lambdas (correlated sampling).
long poisson_inverse(const double u, const double lambda);

}

#endif  // SRC_RANDOM_H_
