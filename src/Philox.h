// Philox4x32-10 counter-based random number generator
// (Salmon, Moraes, Dror, Shaw, "Parallel random numbers: as easy as 1, 2, 3",
// SC11). Same output as Random123's philox4x32 with 10 rounds.
//
// A counter-based generator returns the random numbers for a (key, counter)
// pair directly, so a draw can be addressed by what it is for (a cell, a
// species, a hadron) instead of by its position in one long sequence.

#ifndef SRC_PHILOX_H_
#define SRC_PHILOX_H_

#include <array>
#include <cstdint>

namespace RandomUtil {

inline std::array<uint32_t, 4> philox4x32_10(std::array<uint32_t, 4> ctr,
                                             std::array<uint32_t, 2> key) {
    const uint32_t M0 = 0xD2511F53u, M1 = 0xCD9E8D57u;
    const uint32_t W0 = 0x9E3779B9u, W1 = 0xBB67AE85u;
    for (int round = 0; round < 10; round++) {
        if (round > 0) {
            key[0] += W0;
            key[1] += W1;
        }
        const uint64_t p0 = static_cast<uint64_t>(M0)*ctr[0];
        const uint64_t p1 = static_cast<uint64_t>(M1)*ctr[2];
        const uint32_t hi0 = static_cast<uint32_t>(p0 >> 32);
        const uint32_t lo0 = static_cast<uint32_t>(p0);
        const uint32_t hi1 = static_cast<uint32_t>(p1 >> 32);
        const uint32_t lo1 = static_cast<uint32_t>(p1);
        ctr = {hi1 ^ ctr[1] ^ key[0], lo1, hi0 ^ ctr[3] ^ key[1], lo0};
    }
    return ctr;
}

//! a 32-bit random word -> a double in (0, 1)
inline double philox_to_uniform(uint32_t x) {
    return (static_cast<double>(x) + 0.5)*(1.0/4294967296.0);
}

}  // namespace RandomUtil

#endif  // SRC_PHILOX_H_
