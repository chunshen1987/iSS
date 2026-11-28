#ifndef SRC_URQMD_TRANSPORT_H_
#define SRC_URQMD_TRANSPORT_H_

#include <memory>
#include <vector>

#include "data_struct.h"

class UrqmdTransport {
 public:
    struct Result {
        std::vector<iSS_Hadron> particles;
        int collision_attempts = 0;
    };

    static bool available() noexcept;

    // Times are in fm/c. The duration is divided into equal steps no larger
    // than max_step.
    UrqmdTransport(double duration, double max_step, int seed);
    ~UrqmdTransport();

    // An empty sample returns an empty result and leaves the RNG streams intact.
    Result evolve(const std::vector<iSS_Hadron> &event);

 private:
    struct State;
    std::unique_ptr<State> state_;
};

#endif  // SRC_URQMD_TRANSPORT_H_
