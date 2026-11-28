#include "urqmd_output.h"
#include "urqmd_event.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "sample_output.h"

namespace iSS_data {
namespace {

constexpr int max_baryon_type = 55;
constexpr int first_meson_type = 100;
constexpr std::size_t max_particles = 100000;  // nmax in UrQMD's coms.f

struct ParticleType {
    int type;
    int twice_isospin;
};

const std::unordered_map<long long, ParticleType> &particle_types() {
    // PDG mapping from osc2u/src/pdg2ityp.f (Steffen A. Bass).
    static const std::unordered_map<long long, ParticleType> types = {
        // Neutron
        {  2112, {  1, -1}},
        // Proton
        {  2212, {  1,  1}},
        // N*
        { 12112, {  2, -1}}, { 12212, {  2,  1}},
        {  1214, {  3, -1}}, {  2124, {  3,  1}},
        { 22112, {  4, -1}}, { 22212, {  4,  1}},
        { 32112, {  5, -1}}, { 32212, {  5,  1}},
        {  2116, {  6, -1}}, {  2216, {  6,  1}},
        { 12116, {  7, -1}}, { 12216, {  7,  1}},
        { 21214, {  8, -1}}, { 22124, {  8,  1}},
        { 42112, {  9, -1}}, { 42212, {  9,  1}},
        { 31214, { 10, -1}}, { 32124, { 10,  1}},
        { 41214, { 11, -1}}, { 42124, { 11,  1}},
        { 11218, { 12, -1}}, { 12128, { 12,  1}},
        { 51214, { 13, -1}}, { 52124, { 13,  1}},
        { 21218, { 14, -1}}, { 22128, { 14,  1}},
        {  2110, { 15, -1}}, {  2210, { 15,  1}},
        { 12110, { 16, -1}}, { 12210, { 16,  1}},
        // Delta
        {  1114, { 17, -3}}, {  2114, { 17, -1}}, {  2214, { 17,  1}}, {  2224, { 17,  3}},
        { 31114, { 18, -3}}, { 32114, { 18, -1}}, { 32214, { 18,  1}}, { 32224, { 18,  3}},
        {  1112, { 19, -3}}, {  1212, { 19, -1}}, {  2122, { 19,  1}}, {  2222, { 19,  3}},
        { 11114, { 20, -3}}, { 12114, { 20, -1}}, { 12214, { 20,  1}}, { 12224, { 20,  3}},
        { 11112, { 21, -3}}, { 11212, { 21, -1}}, { 12122, { 21,  1}}, { 12222, { 21,  3}},
        {  1116, { 22, -3}}, {  1216, { 22, -1}}, {  2126, { 22,  1}}, {  2226, { 22,  3}},
        { 21112, { 23, -3}}, { 21212, { 23, -1}}, { 22122, { 23,  1}}, { 22222, { 23,  3}},
        { 21114, { 24, -3}}, { 22114, { 24, -1}}, { 22214, { 24,  1}}, { 22224, { 24,  3}},
        { 11116, { 25, -3}}, { 11216, { 25, -1}}, { 12126, { 25,  1}}, { 12226, { 25,  3}},
        {  1118, { 26, -3}}, {  2118, { 26, -1}}, {  2218, { 26,  1}}, {  2228, { 26,  3}},
        // Lambda
        {  3122, { 27,  0}},
        { 13122, { 28,  0}},
        {  3124, { 29,  0}},
        { 23122, { 30,  0}},
        { 33122, { 31,  0}},
        { 13124, { 32,  0}},
        { 43122, { 33,  0}},
        { 53122, { 34,  0}},
        {  3126, { 35,  0}},
        { 13126, { 36,  0}},
        { 23124, { 37,  0}},
        {  3128, { 38,  0}},
        { 23126, { 39,  0}},
        // Sigma
        {  3112, { 40, -2}}, {  3212, { 40,  0}}, {  3222, { 40,  2}},
        {  3114, { 41, -2}}, {  3214, { 41,  0}}, {  3224, { 41,  2}},
        { 13112, { 42, -2}}, { 13212, { 42,  0}}, { 13222, { 42,  2}},
        { 13114, { 43, -2}}, { 13214, { 43,  0}}, { 13224, { 43,  2}},
        { 23112, { 44, -2}}, { 23212, { 44,  0}}, { 23222, { 44,  2}},
        {  3116, { 45, -2}}, {  3216, { 45,  0}}, {  3226, { 45,  2}},
        { 13116, { 46, -2}}, { 13216, { 46,  0}}, { 13226, { 46,  2}},
        { 23114, { 47, -2}}, { 23214, { 47,  0}}, { 23224, { 47,  2}},
        {  3118, { 48, -2}}, {  3218, { 48,  0}}, {  3228, { 48,  2}},
        // Xi
        {  3312, { 49, -1}}, {  3322, { 49,  1}},
        {  3314, { 50, -1}}, {  3324, { 50,  1}},
        { 13312, { 51, -1}}, { 13322, { 51,  1}},
        { 13314, { 52, -1}}, { 13324, { 52,  1}},
        { 23314, { 53, -1}}, { 23324, { 53,  1}},
        {  3316, { 54, -1}}, {  3326, { 54,  1}},
        // Omega
        {  3334, { 55,  0}},
        // gamma
        {    22, {100,  0}},
        // pion
        {  -211, {101, -2}}, {   111, {101,  0}}, {   211, {101,  2}},
        // eta
        {   221, {102,  0}},
        // omega
        {   223, {103,  0}},
        // rho
        {  -213, {104, -2}}, {   113, {104,  0}}, {   213, {104,  2}},
        // f0(980)
        { 10221, {105,  0}},
        // kaon
        {   311, {106, -1}}, {   321, {106,  1}},
        // eta'
        {   331, {107,  0}},
        // k*(892)
        {   313, {108, -1}}, {   323, {108,  1}},
        // phi
        {   333, {109,  0}},
        // k0*(1430)
        { 10311, {110, -1}}, { 10321, {110,  1}},
        // a0(980)
        {-10211, {111, -2}}, { 10111, {111,  0}}, { 10211, {111,  2}},
        // f0(1370)
        { 20221, {112,  0}},
        // k1(1270)
        { 10313, {113, -1}}, { 10323, {113,  1}},
        // a1(1260)
        {-20213, {114, -2}}, { 20113, {114,  0}}, { 20213, {114,  2}},
        // f1(1285)
        { 20223, {115,  0}},
        // f1'(1510)
        { 40223, {116,  0}},
        // k2*(1430)
        {   315, {117, -1}}, {   325, {117,  1}},
        // a2(1329)
        {  -215, {118, -2}}, {   115, {118,  0}}, {   215, {118,  2}},
        // f2(1270)
        {   225, {119,  0}},
        // f2'(1525)
        {   335, {120,  0}},
        // k1(1400)
        { 20313, {121, -1}}, { 20323, {121,  1}},
        // b1
        {-10213, {122, -2}}, { 10113, {122,  0}}, { 10213, {122,  2}},
        // h1
        { 10223, {123,  0}},
        // h1 (1380)
        { 10333, {124,  0}},
        // K* (1410)
        { 30313, {125, -1}}, { 30323, {125,  1}},
        // rho (1450)
        {-40213, {126, -2}}, { 40113, {126,  0}}, { 40213, {126,  2}},
        // omega (1420)
        { 50223, {127,  0}},
        // phi(1680)
        {100333, {128,  0}},
        // k*(1680)
        { 40313, {129, -1}}, { 40323, {129,  1}},
        // rho(1700)
        {-30213, {130, -2}}, { 30113, {130,  0}}, { 30213, {130,  2}},
        // omega(1600)
        { 60223, {131,  0}},
        // phi(1850)
        {   337, {132,  0}},
        // D
        {   421, {133, -1}}, {   411, {133,  1}},
        // D*
        { 10421, {134, -1}}, { 10411, {134,  1}},
        // J
    };
    return types;
}

ParticleType urqmd_type(int pdg) {
    const auto &types = particle_types();
    auto match = types.find(pdg);
    if (match != types.end()) return match->second;
    if (pdg < 0) {
        match = types.find(-static_cast<long long>(pdg));
        if (match != types.end()) {
            return {-match->second.type, -match->second.twice_isospin};
        }
    }
    throw std::invalid_argument("unsupported PDG ID " + std::to_string(pdg));
}

int electric_charge(const ParticleType &particle) {
    const int type = std::abs(particle.type);
    const int sign = particle.type > 0 ? 1 : -1;
    const int baryon = type <= max_baryon_type ? 1 : 0;
    int strange = 0;
    switch (type) {
        case 106: case 108: case 110: case 113: case 117:
        case 121: case 125: case 129:
            strange = -1;
            break;
    }
    if (type >= 27 && type <= 48) strange = 1;
    if (type >= 49 && type <= 54) strange = 2;
    if (type == 55) strange = 3;
    // Use the bundled UrQMD charge convention for charmed mesons.
    if (type == 138 || type == 139) return sign;
    const int charm = (type == 133 || type == 134) ? 1 : 0;
    return (particle.twice_isospin + sign*(baryon - strange + charm))/2;
}

double decay_width(int type) {
    // osc2u's widres and widmes tables store most literals as Fortran REAL
    // before promotion. The final seven widths come from the UrQMD table.
    static constexpr float baryons[] = {
        0.f,
        .200f, .125f, .150f, .150f, .140f, .120f, .100f, .110f,
        .150f, .500f, .550f, .250f, .550f, .550f, .470f,
        .115f, .200f, .180f, .300f, .240f, .280f, .250f, .150f, .250f, .250f,
        0.f, .050f, .016f, .150f, .035f, .060f, .300f, .150f, .080f,
        .095f, .100f, .200f, .200f,
        0.f, .036f, .100f, .060f, .090f, .120f, .120f, .220f, .180f,
        0.f, .009f, .050f, .024f, .060f, .020f, 0.f,
    };
    static constexpr float mesons[] = {
        0.f, 0.f, 0.f, .00843f, .151f, .100f, 0.f, .000201f,
        .050f, .00443f, .287f, .100f, .200f, .090f, .400f, .024f,
        .350f, .100f, .107f, .185f, .076f, .174f, .142f, .360f,
        .080f, .227f, .310f, .174f, .150f, .323f, .235f, .280f,
        .400f, 0.f, .002f, 0.f, 0.f, 0.f, 0.f, 0.f,
    };
    return type <= max_baryon_type ? baryons[type - 1]
                                   : mesons[type - first_meson_type];
}

struct Collision {
    int projectile_mass;
    int projectile_charge;
    int target_mass;
    int target_charge;
    int frame;
    double beam_energy;
    double cm_energy;
    double beam_momentum;
    double beta_nn;
    double beta_lab;
    double beta_projectile;
};

Collision read_collision(const std::string &filename) {
    std::ifstream header(filename);
    if (!header) sample_output_error(filename, "cannot open collision metadata");
    std::string line;
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(header, line)) {
            sample_output_error(filename, "missing OSCAR collision metadata");
        }
    }
    // Skip the two A8 model/version fields and their separators.
    std::istringstream fields(line.size() >= 20 ? line.substr(20) : "");
    fields.imbue(std::locale::classic());
    Collision collision{};
    char open_p, comma_p, close_p, plus, open_t, comma_t, close_t;
    std::string frame, energy;
    int test_particles;
    if (!(fields >> open_p >> collision.projectile_mass >> comma_p
                 >> collision.projectile_charge >> close_p >> plus
                 >> open_t >> collision.target_mass >> comma_t
                 >> collision.target_charge >> close_t >> frame
                 >> energy >> test_particles)
        || open_p != '(' || comma_p != ',' || close_p != ')' || plus != '+'
        || open_t != '(' || comma_t != ',' || close_t != ')'
        || collision.projectile_mass <= 0 || collision.target_mass < 0) {
        sample_output_error(filename, "invalid OSCAR collision metadata");
    }
    std::replace(energy.begin(), energy.end(), 'D', 'E');
    std::replace(energy.begin(), energy.end(), 'd', 'e');
    std::istringstream beam_energy(energy);
    beam_energy.imbue(std::locale::classic());
    if (!(beam_energy >> collision.beam_energy) || !beam_energy.eof()
        || !std::isfinite(collision.beam_energy) || collision.beam_energy < 0.) {
        sample_output_error(filename, "invalid beam energy");
    }

    constexpr double nucleon_mass = 0.938f;
    const double target_mass = collision.target_mass == 0 ? 0. : nucleon_mass;
    const double kinetic_energy = collision.beam_energy;
    collision.cm_energy = std::sqrt(nucleon_mass*nucleon_mass
        + 2.*(kinetic_energy + nucleon_mass)*target_mass + target_mass*target_mass);
    collision.beam_momentum = std::sqrt(kinetic_energy*(kinetic_energy + 2.*nucleon_mass));
    if (!std::isfinite(collision.cm_energy) || !std::isfinite(collision.beam_momentum)) {
        sample_output_error(filename, "beam energy exceeds the supported range");
    }
    const double pcm = std::max(1e-10,
        collision.beam_momentum*target_mass/collision.cm_energy);
    const double beta_target_cm = pcm/std::hypot(pcm, target_mass);
    const double beta_projectile_cm = pcm/std::hypot(pcm, nucleon_mass);
    const double beta_beam = collision.beam_momentum
                             /std::hypot(collision.beam_momentum, nucleon_mass);
    if (frame == "eqsp") {
        collision.beta_lab = beta_target_cm;
        collision.beta_projectile = -beta_projectile_cm;
    } else if (frame == "tar") {
        collision.frame = 1;
        collision.beta_nn = -beta_target_cm;
        collision.beta_projectile = -beta_beam;
    } else if (frame == "pro") {
        collision.frame = 2;
        collision.beta_nn = beta_target_cm;
        collision.beta_lab = beta_beam;
    } else {
        sample_output_error(filename, "unknown OSCAR reference frame");
    }
    return collision;
}

// Serialize the Fortran Ew.d fields used by the file-14 reader.
void write_exponential(std::ostream &out, double value, int width, int digits) {
    char scientific[32];
    std::snprintf(scientific, sizeof(scientific), "%.*E", digits - 1,
                  std::abs(value));
    const char *exponent_start = std::strchr(scientific, 'E');
    const char *fraction = exponent_start - (digits - 1);
    const int exponent = value == 0. ? 0 : std::atoi(exponent_start + 1) + 1;
    char field[32];
    const char *sign = std::signbit(value) ? "-" : "";
    if (std::abs(exponent) < 100) {
        std::snprintf(field, sizeof(field), "%s0.%c%.*sE%+03d",
                      sign, scientific[0], digits - 1, fraction, exponent);
    } else {
        std::snprintf(field, sizeof(field), "%s0.%c%.*s%+04d",
                      sign, scientific[0], digits - 1, fraction, exponent);
    }
    out << std::setw(width) << field;
}

void write_header(std::ostream &out, const Collision &collision,
                  std::size_t event_index, std::size_t count) {
    out << "UQMD   version:       10000   1000  10001  output_file  14\n"
        << "projectile:  (mass, char) " << std::setw(4) << collision.projectile_mass
        << std::setw(4) << collision.projectile_charge
        << "   target:  (mass, char) " << std::setw(4) << collision.target_mass
        << std::setw(4) << collision.target_charge << ' ' << '\n'
        << "transformation betas (NN,lab,pro)   " << std::fixed << std::setprecision(7)
        << std::setw(11) << collision.beta_nn << std::setw(11) << collision.beta_lab
        << std::setw(11) << collision.beta_projectile << '\n'
        << "impact_parameter_real/min/max(fm):    0.00  0.00  0.00"
           "  total_cross_section(mbarn):       0.00\n"
        << "equation_of_state:    0  E_lab(GeV/u):";
    write_exponential(out, collision.beam_energy, 11, 4);
    out << "  sqrt(s)(GeV):";
    write_exponential(out, collision.cm_energy, 11, 4);
    out << "  p_lab(GeV/u):";
    write_exponential(out, collision.beam_momentum, 11, 4);
    out << '\n' << "event# " << std::setw(9) << event_index
        << " random seed:           0 readin   total_time(fm/c):        0"
           " Delta(t)_O(fm/c):        0.000\n";

    // Integration defaults from osc2u.f, with the selected reference frame.
    const int options[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, collision.frame, 0, 2, 1,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    };
    const double parameters[] = {
        1., .52, .5, .3, 0., .37, 0., .093, .35, .25, 0., .5,
        .27, .49, .27, 1., 1.6, .85, 1.55, 0., 0., 1., 0., 0.,
        .9, 50., 1., 1., .4, 1.5, 1.6, 0., 2.5, .1, 3., .275,
        .42, 1.08, .8, .5, 0., .55, 5., .8, .5, 800000., 1., 2.,
    };
    for (int i = 0; i < 45; ++i) {
        if (i % 15 == 0) out << "op";
        out << std::setw(3) << options[i] << "  ";
        if (i % 15 == 14) out << '\n';
    }
    for (int i = 0; i < 48; ++i) {
        if (i % 12 == 0) out << "pa";
        write_exponential(out, parameters[i], 11, 4);
        out << "  ";
        if (i % 12 == 11) out << '\n';
    }
    out << "pvec: r0              rx              ry              rz"
           "              p0              px              py              pz"
           "              m          ityp 2i3 chg lcl#  ncl or \n"
        << std::setw(12) << count << std::setw(12) << 0 << '\n'
        << "       0       0       0       0       0       0       0       0\n";
}

struct EventParticle {
    const iSS_Hadron *sample;
    ParticleType type;
};

bool is_baryon(const EventParticle &particle) {
    return std::abs(particle.type.type) <= max_baryon_type;
}

void validate_sample(const iSS_Hadron &sample, const ParticleType &type) {
    const float values[] = {sample.t, sample.x, sample.y, sample.z, sample.E,
                            sample.px, sample.py, sample.pz, sample.mass};
    for (float value : values) {
        if (!std::isfinite(value)) throw std::invalid_argument("nonfinite particle record");
    }
    if (!(sample.E > 0.) || !(sample.mass >= 0.)
        || (sample.mass == 0. && decay_width(std::abs(type.type)) > 0.)) {
        throw std::invalid_argument("invalid particle energy or mass");
    }
}

UrqmdParticle convert_particle(const EventParticle &particle, gsl_rng *rng) {
    const auto &sample = *particle.sample;
    const double width = decay_width(std::abs(particle.type.type));
    double decay_time = 1e34;
    if (width > 1e-10) {
        // Match osc2u's single-precision ran2 draw and time dilation.
        constexpr float scale = 1.f/2147483563.f;
        const double uniform = std::min(scale*static_cast<float>(gsl_rng_get(rng)),
                                        1.f - 1.2e-7f);
        decay_time = -std::log(1. - uniform)/width
                     *(static_cast<double>(sample.E)/sample.mass)*0.197327f;
    }
    decay_time += sample.t;

    // vni_procev.f starts the cascade at t=0 and writes thad (the original
    // freeze-out time) in the formation-time field.
    double time = sample.t;
    const double energy = sample.E;
    double x = sample.x;
    double y = sample.y;
    double z = sample.z;
    if (time > 1e-5) {
        x -= sample.px/energy*time;
        y -= sample.py/energy*time;
        z -= sample.pz/energy*time;
        time = 0.;
    }
    return {time, x, y, z, energy, sample.px, sample.py, sample.pz,
            sample.mass, decay_time, sample.t, 0., particle.type.type,
            particle.type.twice_isospin, electric_charge(particle.type)};
}

void write_particle(std::ostream &out, const UrqmdParticle &particle) {
    const double values[] = {
        particle.t, particle.x, particle.y, particle.z, particle.energy,
        particle.px, particle.py, particle.pz, particle.mass,
    };
    for (double value : values) write_exponential(out, value, 16, 8);
    out << std::setw(11) << particle.type
        << std::setw(3) << particle.twice_isospin
        << std::setw(3) << particle.charge
        << std::setw(9) << 0 << std::setw(5) << 0 << std::setw(10) << 99;
    write_exponential(out, particle.decay_time, 16, 8);
    write_exponential(out, particle.formation_time, 16, 8);
    write_exponential(out, particle.cross_section_factor, 16, 8);
    out << '\n';
}

}  // namespace

UrqmdEventConverter::UrqmdEventConverter()
    : rng_(gsl_rng_alloc(gsl_rng_ran2), gsl_rng_free) {
    gsl_rng_set(rng_.get(), 123987);
    gsl_rng_get(rng_.get());  // osc2u's SSEED discards the first draw.
}

std::vector<UrqmdParticle> UrqmdEventConverter::convert(
        const std::vector<iSS_Hadron> &event) {
    if (event.size() > max_particles) {
        throw std::invalid_argument("event exceeds UrQMD's 100000 particle limit");
    }
    std::vector<EventParticle> ordered;
    ordered.reserve(event.size());
    for (const auto &sample : event) {
        const auto type = urqmd_type(sample.pid);
        validate_sample(sample, type);
        ordered.push_back({&sample, type});
    }
    std::stable_partition(ordered.begin(), ordered.end(), is_baryon);

    std::vector<UrqmdParticle> particles;
    particles.reserve(event.size());
    for (const auto &particle : ordered) {
        particles.push_back(convert_particle(particle, rng_.get()));
    }
    return particles;
}

int urqmd_pdg_id(int type, int twice_isospin) {
    // Types above 1000 encode the original PYTHIA PDG ID with an offset.
    if (type > 1000) return type - 1000;
    if (type < -1000) return type + 1000;

    // Use the same PDG assignments as input. In particular, the legacy UrQMD
    // pdgid table omits several resonances and has duplicate kaon IDs.
    using Isospins = std::unordered_map<int, int>;
    static const auto pdg_ids = [] {
        std::unordered_map<int, Isospins> ids;
        for (const auto &entry : particle_types()) {
            ids[entry.second.type].emplace(entry.second.twice_isospin, entry.first);
        }
        return ids;
    }();
    const auto multiplet = pdg_ids.find(std::abs(type));
    if (multiplet != pdg_ids.end()) {
        const auto state = multiplet->second.find(type < 0 ? -twice_isospin
                                                             : twice_isospin);
        if (state != multiplet->second.end()) {
            return type < 0 ? -state->second : state->second;
        }
    }
    throw std::runtime_error("unknown UrQMD type " + std::to_string(type)
                             + " with twice isospin " + std::to_string(twice_isospin));
}

void write_urqmd_samples(const std::string &filename,
                         const std::string &oscar_header_filename,
                         const std::vector<std::vector<iSS_Hadron>*> &events,
                         std::size_t first_event) {
    const auto collision = read_collision(oscar_header_filename);
    std::ofstream out(filename);
    if (!out) sample_output_error(filename, "cannot open file");
    out.imbue(std::locale::classic());
    UrqmdEventConverter converter;
    for (std::size_t event_index = 0; event_index < events.size(); ++event_index) {
        const auto &event = *events[event_index];
        if (event.empty()) continue;  // UrQMD initializes its clock from particle 1.
        try {
            const auto particles = converter.convert(event);
            write_header(out, collision, first_event + event_index, particles.size());
            for (const auto &particle : particles) write_particle(out, particle);
        } catch (const std::invalid_argument &error) {
            sample_output_error(filename, error.what());
        }
        if (!out) sample_output_error(filename, "write failed");
    }
    out.close();
    if (!out) sample_output_error(filename, "write failed");
}

}  // namespace iSS_data
