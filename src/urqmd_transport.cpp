#include "urqmd_transport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

#include "urqmd_event.h"

#ifdef ISS_WITH_URQMD
extern "C" {
void iss_urqmd_setup();
int iss_urqmd_run(int count, const iSS_data::UrqmdParticle *initial,
                  int seed, int steps, double step_size, int index,
                  int *collisions);
void iss_urqmd_particles(iSS_data::UrqmdParticle *final);
}
#endif

namespace {

constexpr int max_seed = 2147483398;  // UrQMD's smaller ran2 modulus minus one.

}  // namespace

struct UrqmdTransport::State {
    State(int steps_in, double step_size_in, int seed)
        : steps(steps_in), step_size(step_size_in), next_seed(seed) {}

    iSS_data::UrqmdEventConverter converter;
    int steps;
    double step_size;
    int next_seed;
    int next_event = 1;
};

bool UrqmdTransport::available() noexcept {
#ifdef ISS_WITH_URQMD
    return true;
#else
    return false;
#endif
}

UrqmdTransport::UrqmdTransport(double duration, double max_step, int seed) {
    if (!available()) {
        throw std::runtime_error("in-memory UrQMD requires a build with ISS_ENABLE_URQMD=ON");
    }
    if (!std::isfinite(duration) || duration <= 0.
        || !std::isfinite(max_step) || max_step <= 0.) {
        throw std::invalid_argument("UrQMD duration and step size must be finite and positive");
    }
    const double steps = std::max(1., std::ceil(duration/max_step));
    if (steps >= std::numeric_limits<int>::max() || duration/steps == 0.) {
        throw std::invalid_argument("UrQMD duration and step size give an unsupported step count");
    }
    if (seed < 1 || seed > max_seed) {
        throw std::invalid_argument("UrQMD seed must be in the range 1 to 2147483398");
    }
    state_.reset(new State(static_cast<int>(steps), duration/steps, seed));
}

UrqmdTransport::~UrqmdTransport() = default;

UrqmdTransport::Result UrqmdTransport::evolve(const std::vector<iSS_Hadron> &event) {
    if (event.empty()) return {};
#ifdef ISS_WITH_URQMD
    // UrQMD stores the event, tables, and RNG in process-wide COMMON blocks.
    static std::mutex urqmd_mutex;
    static std::once_flag initialize;
    const std::lock_guard<std::mutex> lock(urqmd_mutex);
    if (state_->next_event == std::numeric_limits<int>::max()) {
        throw std::runtime_error("UrQMD event index exhausted");
    }
    auto particles = state_->converter.convert(event);
    std::call_once(initialize, iss_urqmd_setup);

    Result result;
    const int count = iss_urqmd_run(static_cast<int>(particles.size()),
                                   particles.data(), state_->next_seed,
                                   state_->steps, state_->step_size,
                                   state_->next_event, &result.collision_attempts);
    if (count < 0) throw std::runtime_error("UrQMD rejected the event");
    ++state_->next_event;
    state_->next_seed = state_->next_seed == max_seed ? 1 : state_->next_seed + 1;
    particles.resize(count);
    if (count != 0) iss_urqmd_particles(particles.data());

    result.particles.reserve(count);
    for (const auto &p : particles) {
        const double values[] = {p.mass, p.energy, p.px, p.py, p.pz,
                                 p.t, p.x, p.y, p.z};
        for (double value : values) {
            if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<float>::max()) {
                throw std::runtime_error("UrQMD returned an invalid particle record");
            }
        }
        result.particles.push_back({
            iSS_data::urqmd_pdg_id(p.type, p.twice_isospin),
            static_cast<float>(p.mass), static_cast<float>(p.energy),
            static_cast<float>(p.px), static_cast<float>(p.py),
            static_cast<float>(p.pz), static_cast<float>(p.t),
            static_cast<float>(p.x), static_cast<float>(p.y),
            static_cast<float>(p.z),
        });
    }
    return result;
#else
    throw std::runtime_error("in-memory UrQMD is unavailable");
#endif
}
