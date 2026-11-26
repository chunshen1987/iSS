// Run from the repository root so the decayer can load iSS_tables.
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "particle_decay.h"

int main() {
    auto random = std::make_shared<RandomUtil::Random>(1);
    particle_decay decayer(random, AfterburnerType::UrQMD);
    bool saw_two_body = false;
    bool saw_three_body = false;
    for (const int pid : {113, 223}) {
        iSS_Hadron mother{};
        mother.pid = pid;
        mother.mass = decayer.get_particle_mass(pid);
        mother.E = mother.mass;
        for (int trial = 0; trial < 20; ++trial) {
            std::vector<iSS_Hadron> daughters;
            decayer.perform_decays(&mother, &daughters);
            saw_two_body |= daughters.size() == 2;
            saw_three_body |= daughters.size() == 3;
            std::array<double, 4> sum{};
            for (const auto &daughter : daughters) {
                sum[0] += daughter.E;
                sum[1] += daughter.px;
                sum[2] += daughter.py;
                sum[3] += daughter.pz;
            }
            const std::array<double, 4> expected = {mother.E, 0., 0., 0.};
            for (std::size_t i = 0; i < sum.size(); ++i) {
                if (!std::isfinite(sum[i]) || std::abs(sum[i] - expected[i]) > 1e-6) {
                    std::cerr << "Invalid decay at rest for " << pid << '\n';
                    return 1;
                }
            }
        }
    }
    if (!saw_two_body || !saw_three_body) return 1;
    std::cout << "Two- and three-body decays at rest conserve four-momentum.\n";
}
