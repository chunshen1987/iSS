#ifndef SRC_URQMD_EVENT_H_
#define SRC_URQMD_EVENT_H_

#include <memory>
#include <vector>

#include "gsl/gsl_rng.h"

#include "data_struct.h"

namespace iSS_data {

// Layout shared with the bind(C) record in urqmd_bridge.f.
struct UrqmdParticle {
    double t, x, y, z, energy, px, py, pz, mass;
    double decay_time, formation_time, cross_section_factor;
    int type, twice_isospin, charge;
};

class UrqmdEventConverter {
 public:
    UrqmdEventConverter();
    std::vector<UrqmdParticle> convert(const std::vector<iSS_Hadron> &event);

 private:
    std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)> rng_;
};

int urqmd_pdg_id(int type, int twice_isospin);

}  // namespace iSS_data

#endif  // SRC_URQMD_EVENT_H_
