#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "urqmd_transport.h"

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    try {
        require(UrqmdTransport::available(), "UrQMD support is unavailable");
        UrqmdTransport transport(15., 5., 4242);

        // At rest, the charged rho should decay into pi+ and pi0.
        const std::vector<iSS_Hadron> rho{{213, .775f, .775f, 0, 0, 0, 1, 0, 0, 0}};
        auto decay = transport.evolve(rho);
        auto &daughters = decay.particles;
        std::sort(daughters.begin(), daughters.end(),
                  [](const iSS_Hadron &a, const iSS_Hadron &b) { return a.pid < b.pid; });
        require(daughters.size() == 2 && daughters[0].pid == 111 && daughters[1].pid == 211,
                "rho decay products were not exported");
        require(std::abs(daughters[0].E + daughters[1].E - rho[0].E) < 1e-6,
                "rho decay did not conserve energy");

        // Approach at low center-of-mass energy with a small impact parameter.
        const float energy = std::hypot(.938f, .35f);
        auto scattered = transport.evolve({
            {2212, .938f, energy, .35f, 0, 0, 1, -2, .15f, 0},
            {2212, .938f, energy, -.35f, 0, 0, 1, 2, -.15f, 0},
        });
        require(scattered.collision_attempts == 1 && scattered.particles.size() == 2,
                "proton pair did not scatter");
        require(std::abs(scattered.particles[0].py) > 1e-5
                || std::abs(scattered.particles[0].pz) > 1e-5,
                "UrQMD did not change the proton momenta");

        require(transport.evolve({}).particles.empty(), "empty sample returned particles");
        auto pion = transport.evolve({{-211, .138f, .138f, 0, 0, 0, 1, 3, 7, 13}});
        require(pion.particles.size() == 1 && pion.particles[0].pid == -211
                && pion.particles[0].t == 15.f && pion.particles[0].x == 3.f,
                "state from an earlier event leaked into the pion event");

        UrqmdTransport replay(15., 5., 4242);
        auto repeated = replay.evolve(rho).particles;
        std::sort(repeated.begin(), repeated.end(),
                  [](const iSS_Hadron &a, const iSS_Hadron &b) { return a.pid < b.pid; });
        require(repeated.size() == daughters.size() && repeated[0].px == daughters[0].px
                && repeated[0].py == daughters[0].py && repeated[0].pz == daughters[0].pz,
                "a new wrapper did not reproduce the first seeded event");
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
