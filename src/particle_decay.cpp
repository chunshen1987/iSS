// Copyright Chun Shen @ 2016

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "data_struct.h"
#include "particle_decay.h"

using namespace std;

particle_decay::particle_decay(shared_ptr<RandomUtil::Random> ran_gen,
                               AfterburnerType afterburner_type,
                               string table_path) :
        table_path_(table_path), afterburner_type_(afterburner_type) {
    ran_gen_ptr = ran_gen;
    resonance_table.clear();
    read_resonances_list();
    prepare_channel_kinematics();
}

particle_decay::~particle_decay() {
    if (n_unreachable_.load() > 0) {
        cerr << "[Warning] particle_decay: " << n_unreachable_.load()
             << " resonances were below all of their decay channels and"
             << " were kept as final particles" << endl;
    }
    if (n_mass_fallback_.load() > 0) {
        cerr << "[Warning] particle_decay: daughter masses of "
             << n_mass_fallback_.load() << " decays were set without sampling"
             << endl;
    }
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        for (int j = 0; j < resonance_table[i]->decays; j++) {
            delete resonance_table[i]->decay_channels[j];
        }
        delete resonance_table[i];
    }
    resonance_table.clear();
}

//! This function reads in resonance decay table
int particle_decay::read_resonances_list() {
    double eps = 1e-15;
    cout << " -- Read in particle resonance decay table..." << endl;
    std::string reso_filename;
    if (afterburner_type_ == AfterburnerType::SMASH) {
        reso_filename = table_path_ + "/pdg-SMASH.dat";
    } else if (afterburner_type_ == AfterburnerType::UrQMD) {
        reso_filename = table_path_ + "/pdg-urqmd_v3.3+.dat";
    } else {
        reso_filename = table_path_ + "/pdg-s95pv1.dat";
    }
    std::ifstream resofile(reso_filename.c_str());
    if (!resofile.good()) {
        cout << "[Error] Can not found pdg file: " << reso_filename << endl;
        exit(1);
    }

    int dummy_int;
    int particle_monval;
    resofile >> particle_monval;
    while (!resofile.eof()) {
        // add one resonance
        particle_info *temp_resonance = new particle_info;
        temp_resonance->monval = particle_monval;
        resofile >> temp_resonance->name;
        resofile >> temp_resonance->mass;
        resofile >> temp_resonance->width;
        resofile >> temp_resonance->gspin;        // spin degeneracy
        resofile >> temp_resonance->baryon;
        resofile >> temp_resonance->strange;
        resofile >> temp_resonance->charm;
        resofile >> temp_resonance->bottom;
        resofile >> temp_resonance->gisospin;     // isospin degeneracy
        resofile >> temp_resonance->charge;
        resofile >> temp_resonance->decays;
        for (int j = 0; j < temp_resonance->decays; j++) {
            // read in its decay channels
            decay_channel_info *temp_decay_channel = new decay_channel_info;
            resofile >> dummy_int;
            resofile >> temp_decay_channel->decay_Npart;
            resofile >> temp_decay_channel->branching_ratio;
            resofile >> temp_decay_channel->decay_part[0];
            resofile >> temp_decay_channel->decay_part[1];
            resofile >> temp_decay_channel->decay_part[2];
            resofile >> temp_decay_channel->decay_part[3];
            resofile >> temp_decay_channel->decay_part[4];
            temp_resonance->decay_channels.push_back(temp_decay_channel);
        }

        // decide whether particle is stable under strong interactions
        if (temp_resonance->decay_channels[0]->decay_Npart == 1) {
            temp_resonance->stable = 1;
        } else {
            temp_resonance->stable = 0;
        }
        resonance_table.push_back(temp_resonance);

        if (temp_resonance->baryon == 1) {
            // add anti-particle entry for baryon
            particle_info *temp_anti_resonance = new particle_info;
            temp_anti_resonance->monval = -temp_resonance->monval;
            ostringstream antiname;
            antiname << "Anti-" << temp_resonance->name;
            temp_anti_resonance->name = antiname.str();
            temp_anti_resonance->mass = temp_resonance->mass;
            temp_anti_resonance->width = temp_resonance->width;
            temp_anti_resonance->gspin = temp_resonance->gspin;
            temp_anti_resonance->baryon = -temp_resonance->baryon;
            temp_anti_resonance->strange = -temp_resonance->strange;
            temp_anti_resonance->charm = -temp_resonance->charm;
            temp_anti_resonance->bottom = -temp_resonance->bottom;
            temp_anti_resonance->gisospin = temp_resonance->gisospin;
            temp_anti_resonance->charge = -temp_resonance->charge;
            temp_anti_resonance->decays = temp_resonance->decays;
            temp_anti_resonance->stable = temp_resonance->stable;
            for (int j = 0; j < temp_resonance->decays; j++) {
                // add anti-particle decay channels
                decay_channel_info *temp_anti_decay_channel = (
                                                    new decay_channel_info);
                temp_anti_decay_channel->decay_Npart =
                        temp_resonance->decay_channels[j]->decay_Npart;
                temp_anti_decay_channel->branching_ratio =
                        temp_resonance->decay_channels[j]->branching_ratio;
                for (int k = 0; k < 5; k++) {
                    int decay_part_monval =
                            temp_resonance->decay_channels[j]->decay_part[k];
                    if (decay_part_monval == 0) {
                        // a null entry
                        temp_anti_decay_channel->decay_part[k] = 0;
                    } else {
                        // find the index for decay particle in the
                        // current resonance table
                        int idx;
                        int reso_table_length = resonance_table.size();
                        for (idx = 0; idx < reso_table_length; idx++) {
                            if (resonance_table[idx]->monval
                                == decay_part_monval) {
                                break;
                            }
                        }
                        double temp_br =
                            temp_resonance->decay_channels[j]->branching_ratio;
                        if (idx == reso_table_length
                            && temp_resonance->stable == 0 && temp_br > eps) {
                            cout << "Error: can not find decay particle index "
                                 << "for anti-baryon!" << endl;
                            cout << "particle monval : " 
                                 << decay_part_monval << endl;
                            exit(1);
                        }
                        if (resonance_table[idx]->baryon == 0
                            && resonance_table[idx]->charge == 0 
                            && resonance_table[idx]->strange == 0) {
                            temp_anti_decay_channel->decay_part[k] = (
                                                        decay_part_monval);
                        } else {
                            temp_anti_decay_channel->decay_part[k] = (
                                                        - decay_part_monval);
                        }
                    }
                }
                temp_anti_resonance->decay_channels.push_back(
                                                    temp_anti_decay_channel);
            }
            resonance_table.push_back(temp_anti_resonance);
        }
        resofile >> particle_monval;
    }
    resofile.close();  // close the pdg file
    int Nparticle = resonance_table.size();
    for (int i = 0; i < Nparticle; i++) {
        // determine the quantum statistics for bosons and fermions
        if (resonance_table[i]->baryon == 0) {
            resonance_table[i]->sign = -1;
        } else {
            resonance_table[i]->sign = 1;
        }
    }
    return(Nparticle);
}


//! This is a test function to check whether the resonance table is read in
//! correctly
void particle_decay::check_resonance_table() {
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        cout << "name: " << resonance_table[i]->name << endl;
        cout << "monval: " << resonance_table[i]->monval << endl;
        cout << "mass: " << resonance_table[i]->mass << endl;
        cout << "stable: " << resonance_table[i]->stable << endl;
        cout << "# of decay channels: " << resonance_table[i]->decays << endl;
        for (int j = 0; j < resonance_table[i]->decays; j++) {
            cout << j << ": # of daughters: "
                 << resonance_table[i]->decay_channels[j]->decay_Npart << endl;
            cout << j << ": branching_ratio: "
                 << resonance_table[i]->decay_channels[j]->branching_ratio
                 << endl;
            cout << j << ": decay part: ";
            for (int k = 0;
                 k < resonance_table[i]->decay_channels[j]->decay_Npart;
                 k++) {
                 cout << resonance_table[i]->decay_channels[j]->decay_part[k]
                      << "  ";
            }
            cout << endl;
        }
    }
}

//! This function returns particle width in GeV
double particle_decay::get_particle_width(iSS_Hadron *part) {
    double width = 0.0;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (part->pid == resonance_table[i]->monval) {
            width = resonance_table[i]->width;
            break;
        }
    }
    return(width);
}

//! This function checks whether the particle is stable.  A resonance whose
//! mass is below all of its decay channels counts as stable too:
//! perform_decays() can only hand it back unchanged, so a cascade that
//! queued it again would never end.  Such particles are counted and
//! reported at the end of the run.
int particle_decay::check_particle_stable(iSS_Hadron *part) {
    const int idx = find_index(part->pid);
    if (idx < 0) return(0);  // perform_decays() stops on it
    if (resonance_table[idx]->stable == 1) return(1);
    if (!has_reachable_channel(idx, part->mass)) {
        n_unreachable_.fetch_add(1, std::memory_order_relaxed);
        return(1);
    }
    return(0);
}


//! This function returns the electric charge of particle
int particle_decay::get_particle_charge(int monval) {
    int charge = 0;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (monval == resonance_table[i]->monval) {
            charge = resonance_table[i]->charge;
            break;
        }
    }
    return(charge);
}

//! This function returns the baryon number of particle
int particle_decay::get_particle_baryon_number(int monval) {
    int baryon = 0;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (monval == resonance_table[i]->monval) {
            baryon = resonance_table[i]->baryon;
            break;
        }
    }
    return(baryon);
}

//! This function returns the strange number of particle
int particle_decay::get_particle_strange_number(int monval) {
    int strange = 0;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (monval == resonance_table[i]->monval) {
            strange = resonance_table[i]->strange;
            break;
        }
    }
    return(strange);
}


namespace {

//! a channel with daughters below their pole masses is reachable only above
//! its threshold by more than this much [GeV]; far above the rounding of
//! the single-precision masses of iSS_Hadron
const double kMassMargin = 1e-5;

//! widths above this [GeV] are strong decays (see is_broad())
const double kBroadWidth = 1e-4;

//! two-body breakup momentum of M -> m1 m2 (0 below threshold)
double breakup_momentum(double M, double m1, double m2) {
    double t = (M*M - (m1 + m2)*(m1 + m2))*(M*M - (m1 - m2)*(m1 - m2));
    return t > 0. ? sqrt(t)/(2.*M) : 0.;
}

//! boost (E, px, py, pz) by the velocity (vx, vy, vz)
void boost(double vx, double vy, double vz,
           double &E, double &px, double &py, double &pz) {
    double v2 = vx*vx + vy*vy + vz*vz;
    if (v2 < 1e-30) return;
    double gamma = 1./sqrt(1. - v2);
    double vp = vx*px + vy*py + vz*pz;
    double coef = (gamma - 1.)*vp/v2 + gamma*E;
    E = gamma*(E + vp);
    px += coef*vx;
    py += coef*vy;
    pz += coef*vz;
}

//! The mass at which to decay the mother: that of its stored four-momentum.
//! iSS_Hadron keeps E and p in float, so for a mother with Lorentz factor
//! gamma, sqrt(E^2 - p^2) differs from its mass by about gamma^2 times the
//! float precision (1e-3 at gamma = 100).  Daughters that share this mass
//! in the rest frame add up to the mother's four-momentum.  Right at a
//! threshold, where it may be below the daughters' masses, the mother's
//! mass is kept.
double rest_frame_mass(const iSS_Hadron &mother, double daughter_masses) {
    const double E = mother.E;
    const double p2 = (static_cast<double>(mother.px)*mother.px
                       + static_cast<double>(mother.py)*mother.py
                       + static_cast<double>(mother.pz)*mother.pz);
    const double m2 = E*E - p2;
    if (m2 > daughter_masses*daughter_masses) return sqrt(m2);
    return mother.mass;
}

//! Boosts (E, px, py, pz) from the rest frame of the mother, of mass M, to
//! the frame where it has its stored four-momentum P.  Written with P/M,
//! not the velocity, so that it keeps its precision at large gamma; a set
//! of daughters with total (M, 0) ends up with total P.
void boost_from_rest_frame(const iSS_Hadron &mother, double M, double &E,
                           double &px, double &py, double &pz) {
    const double PE = mother.E;
    const double Px = mother.px, Py = mother.py, Pz = mother.pz;
    const double Pp = Px*px + Py*py + Pz*pz;
    const double E_lab = (PE*E + Pp)/M;
    const double coef = (E + Pp/(PE + M))/M;
    E = E_lab;
    px += coef*Px;
    py += coef*Py;
    pz += coef*Pz;
}

}  // namespace


int particle_decay::find_index(int monval) const {
    auto it = table_index_.find(monval);
    return it == table_index_.end() ? -1 : it->second;
}


//! Decays through the strong interaction: its mass can be sampled when it
//! is a daughter.  Widths below 0.1 MeV are electromagnetic or weak decays
//! (pi0, eta, Sigma0: eV to keV), kept at their pole masses.
bool particle_decay::is_broad(int idx) const {
    return (resonance_table[idx]->stable == 0
            && resonance_table[idx]->width > kBroadWidth);
}


//! The lowest mass a broad daughter may be given: the lowest at which every
//! one of its own decay channels is still reachable in the sense of
//! channel_reachable(), so that a daughter below its pole keeps the
//! branching ratios of the table and can always decay; never above its
//! pole mass.  The two margins put it strictly inside the reachable range:
//! one is the margin of channel_reachable(), the other leaves room for the
//! daughters' masses in sample_daughter_masses().  Narrow and stable
//! particles: their pole masses.
double particle_decay::compute_lowest_mass(int idx, vector<int> &state) {
    if (state[idx] == 2) return lowest_mass_[idx];
    const particle_info *info = resonance_table[idx];
    const double pole = static_cast<float>(info->mass);
    if (state[idx] == 1) return pole;  // a loop in the decay chains
    if (!is_broad(idx)) {
        lowest_mass_[idx] = pole;
        state[idx] = 2;
        return pole;
    }
    state[idx] = 1;
    double lowest = 0.;
    for (int j = 0; j < info->decays; j++) {
        const decay_channel_info *ch = info->decay_channels[j];
        int n = std::abs(ch->decay_Npart);
        if (n < 2 || ch->branching_ratio <= 0.) continue;
        double sum = 0.;
        for (int k = 0; k < n; k++) {
            int d = find_index(ch->decay_part[k]);
            if (d < 0) continue;  // reported by prepare_channel_kinematics
            sum += compute_lowest_mass(d, state);
        }
        lowest = std::max(lowest, sum);
    }
    lowest = std::min(lowest + 2.*kMassMargin, pole);
    lowest_mass_[idx] = lowest;
    state[idx] = 2;
    return lowest;
}


void particle_decay::prepare_channel_kinematics() {
    const int ntable = resonance_table.size();
    table_index_.clear();
    for (int i = 0; i < ntable; i++) {
        // the first entry wins, as in the linear searches
        table_index_.emplace(resonance_table[i]->monval, i);
    }
    lowest_mass_.assign(resonance_table.size(), 0.);
    vector<int> state(resonance_table.size(), 0);
    for (int i = 0; i < ntable; i++) compute_lowest_mass(i, state);

    channel_kinematics_.assign(resonance_table.size(),
                               vector<ChannelKinematics>());
    for (int i = 0; i < ntable; i++) {
        const particle_info *info = resonance_table[i];
        if (info->stable == 1) continue;
        for (int j = 0; j < info->decays; j++) {
            const decay_channel_info *ch = info->decay_channels[j];
            ChannelKinematics kin;
            kin.threshold_pole = 0.;
            kin.threshold_min = 0.;
            kin.has_broad_daughter = false;
            int n = std::abs(ch->decay_Npart);
            if (n > 5) {
                cout << "[Error] particle_decay: " << info->monval
                     << " has a channel with " << n << " daughters" << endl;
                exit(1);
            }
            for (int k = 0; k < n; k++) {
                int d = find_index(ch->decay_part[k]);
                if (d < 0) {
                    if (ch->branching_ratio > 0. && n > 1) {
                        cout << "[Error] particle_decay: daughter "
                             << ch->decay_part[k] << " of " << info->monval
                             << " is not in the particle table" << endl;
                        exit(1);
                    }
                    continue;
                }
                kin.daughter_index.push_back(d);
                // float, as in iSS_Hadron and the decay routines
                kin.threshold_pole += static_cast<float>(
                                            resonance_table[d]->mass);
                // as the lower limits in sample_daughter_masses()
                kin.threshold_min += (is_broad(d)
                    ? lowest_mass_[d]
                    : static_cast<float>(resonance_table[d]->mass));
                if (is_broad(d)) kin.has_broad_daughter = true;
            }
            channel_kinematics_[i].push_back(kin);
        }
    }

    // A species that cannot decay at its pole mass, or at the lowest mass
    // it can be given as a daughter, stays in the final state
    // (check_particle_stable()): say so once, when the table is read.
    for (int i = 0; i < ntable; i++) {
        const particle_info *info = resonance_table[i];
        if (info->stable == 1) continue;
        const double pole = static_cast<float>(info->mass);
        if (!has_reachable_channel(i, pole)) {
            cout << "[Warning] particle_decay: " << info->monval
                 << " has no decay channel open at its pole mass " << pole
                 << " GeV; it will be kept as a final particle" << endl;
        } else if (!has_reachable_channel(i, lowest_mass_[i])) {
            cout << "[Warning] particle_decay: " << info->monval
                 << " has no decay channel open at " << lowest_mass_[i]
                 << " GeV, the lowest mass it can be given as a daughter"
                 << endl;
        }
    }
}


//! Whether a mother of mass M can decay through a channel: either the
//! channel is open at the daughters' pole masses, strictly, as the decay
//! routines need, or its broad daughters can be given lower masses
//! (sample_daughter_masses()).
bool particle_decay::channel_reachable(const ChannelKinematics &channel,
                                       double M) const {
    if (channel.daughter_index.size() < 2) return false;
    if (M > channel.threshold_pole) return true;
    return (channel.has_broad_daughter
            && M - kMassMargin > channel.threshold_min);
}


bool particle_decay::has_reachable_channel(int idx, double M) const {
    const particle_info *info = resonance_table[idx];
    const vector<ChannelKinematics> &kin = channel_kinematics_[idx];
    for (int i = 0; i < info->decays; i++) {
        if (info->decay_channels[i]->branching_ratio > 0.
                && channel_reachable(kin[i], M)) {
            return true;
        }
    }
    return false;
}


//! Masses for the daughters of a channel that is closed at their pole
//! masses: every broad daughter gets a mass from its (non-relativistic)
//! Breit-Wigner distribution, truncated below at compute_lowest_mass(), and
//! the set is accepted when the decay is open.  This is the Breit-Wigner
//! tail through which the channel is reached; stable and narrow daughters
//! keep their pole masses.
void particle_decay::sample_daughter_masses(double M,
        const ChannelKinematics &channel, vector<double> &masses,
        RandomUtil::Random &ran) {
    const int n = masses.size();
    const double M_top = M - kMassMargin;
    vector<double> lo(n), hi(n);
    double sum_lo = 0.;
    for (int k = 0; k < n; k++) {
        int d = channel.daughter_index[k];
        lo[k] = is_broad(d) ? lowest_mass_[d] : masses[k];
        sum_lo += lo[k];
    }
    for (int k = 0; k < n; k++) hi[k] = M_top - (sum_lo - lo[k]);
    vector<double> trial(masses);
    for (int itry = 0; itry < 100000; itry++) {
        double sum = 0.;
        for (int k = 0; k < n; k++) {
            int d = channel.daughter_index[k];
            if (is_broad(d)) {
                double m0 = resonance_table[d]->mass;
                double half_width = 0.5*resonance_table[d]->width;
                double a = atan2(lo[k] - m0, half_width);
                double b = atan2(hi[k] - m0, half_width);
                double theta = a + (b - a)*ran.rand_uniform();
                trial[k] = std::min(std::max(m0 + half_width*tan(theta),
                                             lo[k]), hi[k]);
            }
            sum += static_cast<float>(trial[k]);
        }
        if (sum < M_top) {
            masses = trial;
            return;
        }
    }
    // the allowed region is far in the tails: share the energy above the
    // lowest masses equally among the unstable daughters
    n_mass_fallback_.fetch_add(1, std::memory_order_relaxed);
    int n_broad = 0;
    for (int k = 0; k < n; k++)
        if (is_broad(channel.daughter_index[k])) n_broad++;
    for (int k = 0; k < n; k++) {
        if (is_broad(channel.daughter_index[k])) {
            masses[k] = lo[k] + 0.5*(M_top - sum_lo)/n_broad;
        } else {
            masses[k] = lo[k];
        }
    }
}


//! Decays a resonance into one of its channels.  The channel is drawn with
//! the probabilities of the table's branching ratios, normalised over the
//! channels the mother can reach: those open at the daughters' pole masses
//! and those that open when its unstable daughters take masses below their
//! poles (sample_daughter_masses()).  A channel that neither reaches is
//! skipped, and its share goes to the others; a mother that reaches none is
//! handed back unchanged, and check_particle_stable() calls it stable, so
//! that the cascades keep it as a final particle.  A negative number of
//! daughters in the table (the UrQMD-derived lists use it to mark channels
//! closed at the pole masses) is read as its absolute value.  Channels with
//! two and three daughters use the original routines and random-number
//! sequence; four and five use perform_n_body_decay().
void particle_decay::perform_decays(
            iSS_Hadron *mother, vector<iSS_Hadron>* daughter_list) {
    perform_decays(mother, *daughter_list, *ran_gen_ptr);
}


void particle_decay::perform_decays(iSS_Hadron *mother,
                                    vector<iSS_Hadron> &daughter_list,
                                    RandomUtil::Random &ran) {
    const int idx = find_index(mother->pid);
    if (idx < 0) {
        cout << "[Error] particle_decay::perform_decays: particle "
             << mother->pid << " is not in the particle table" << endl;
        exit(1);
    }
    const particle_info *mother_decay_info = resonance_table[idx];
    if (mother_decay_info->stable == 1) {
        // the particle is a stable particle
        daughter_list.push_back(*mother);
        return;
    }
    const double M = mother->mass;
    const vector<ChannelKinematics> &kin = channel_kinematics_[idx];
    const int N_decay_channel = mother_decay_info->decays;
    // per-thread buffers: a decay would otherwise allocate three vectors
    static thread_local vector<double> weight;
    static thread_local vector<double> masses;
    static thread_local vector<iSS_Hadron> daughters;
    weight.assign(N_decay_channel, 0.);
    double total = 0.;
    for (int i = 0; i < N_decay_channel; i++) {
        const decay_channel_info *ch = mother_decay_info->decay_channels[i];
        if (std::abs(ch->decay_Npart) < 2) continue;
        if (channel_reachable(kin[i], M)) {
            weight[i] = ch->branching_ratio;
            total += weight[i];
        }
    }
    if (total <= 0.) {
        // counted by check_particle_stable()
        daughter_list.push_back(*mother);
        return;
    }
    double target = ran.rand_uniform();
    if (std::abs(total - 1.) > 1e-12) target *= total;
    int picked = -1;
    double cumulated_branching_ratio = 0.;
    for (int i = 0; i < N_decay_channel; i++) {
        if (weight[i] <= 0.) continue;
        picked = i;  // the last reachable channel takes rounding leftovers
        cumulated_branching_ratio += weight[i];
        if (cumulated_branching_ratio > target) break;
    }
    const decay_channel_info *picked_channel = (
                                mother_decay_info->decay_channels[picked]);
    const ChannelKinematics &channel = kin[picked];
    const int N_decay_part = std::abs(picked_channel->decay_Npart);
    masses.resize(N_decay_part);
    for (int k = 0; k < N_decay_part; k++) {
        masses[k] = static_cast<float>(
                        resonance_table[channel.daughter_index[k]]->mass);
    }
    if (!(M > channel.threshold_pole)) {
        sample_daughter_masses(M, channel, masses, ran);
    }
    daughters.assign(N_decay_part, iSS_Hadron());
    for (int k = 0; k < N_decay_part; k++) {
        daughters[k].pid = picked_channel->decay_part[k];
        daughters[k].mass = masses[k];
    }
    if (N_decay_part == 2) {
        perform_two_body_decay(mother, &daughters[0], &daughters[1], ran);
    } else if (N_decay_part == 3) {
        perform_three_body_decay(mother, &daughters[0], &daughters[1],
                                 &daughters[2], ran);
    } else {
        perform_n_body_decay(mother, daughters, ran);
    }
    for (int k = 0; k < N_decay_part; k++) {
        daughter_list.push_back(daughters[k]);
    }
}


//! Uniform n-body phase space by the Raubold-Lynch method (GENBOD, as in
//! ROOT's TGenPhaseSpace): the invariant masses M_k of the first k
//! daughters are drawn from ordered uniform numbers and accepted with
//! probability prod_k p*(M_k -> M_{k-1} m_k) over its maximum; each step is
//! a two-body decay isotropic in the rest frame of M_k.  The daughters
//! share the mother's decay point, as in the two- and three-body routines.
void particle_decay::perform_n_body_decay(iSS_Hadron *mother,
                                          vector<iSS_Hadron> &daughters,
                                          RandomUtil::Random &ran) {
    const int n = daughters.size();
    vector<double> m(n);
    double m_sum = 0.;
    for (int k = 0; k < n; k++) {
        m[k] = daughters[k].mass;
        m_sum += m[k];
    }
    if (mother->mass <= m_sum) {
        cout << "Error:particle_decay::perform_n_body_decay: "
             << "the decay of " << mother->pid << " (M = " << mother->mass
             << ") into " << n << " daughters of total mass " << m_sum
             << " is closed" << endl;
        exit(1);
    }
    const double M = rest_frame_mass(*mother, m_sum);
    const double T = M - m_sum;
    // largest possible weight
    double w_max = 1.;
    double e_min = 0., e_max = T + m[0];
    for (int k = 1; k < n; k++) {
        e_min += m[k - 1];
        e_max += m[k];
        w_max *= breakup_momentum(e_max, e_min, m[k]);
    }
    vector<double> r(n), M_inv(n), p_star(n);
    double w;
    do {
        r[0] = 0.;
        r[n - 1] = 1.;
        for (int k = 1; k < n - 1; k++) r[k] = ran.rand_uniform();
        std::sort(r.begin() + 1, r.end() - 1);
        double sum = 0.;
        for (int k = 0; k < n; k++) {
            sum += m[k];
            M_inv[k] = r[k]*T + sum;
        }
        w = 1.;
        for (int k = 1; k < n; k++) {
            p_star[k] = breakup_momentum(M_inv[k], M_inv[k - 1], m[k]);
            w *= p_star[k];
        }
    } while (ran.rand_uniform()*w_max > w);

    // build up the daughters in the rest frame of the first k + 1 of them
    vector<double> E(n), px(n), py(n), pz(n);
    for (int k = 1; k < n; k++) {
        double cos_theta = 2.*ran.rand_uniform() - 1.;
        double sin_theta = sqrt(1. - cos_theta*cos_theta);
        double phi = 2.*M_PI*ran.rand_uniform();
        double ux = sin_theta*cos(phi);
        double uy = sin_theta*sin(phi);
        double uz = cos_theta;
        double p = p_star[k];
        if (k == 1) {
            E[0] = sqrt(p*p + m[0]*m[0]);
            px[0] = -p*ux;
            py[0] = -p*uy;
            pz[0] = -p*uz;
        } else {
            // daughters 0..k-1 move together against daughter k
            double E_sub = sqrt(p*p + M_inv[k - 1]*M_inv[k - 1]);
            double vx = -p*ux/E_sub, vy = -p*uy/E_sub, vz = -p*uz/E_sub;
            for (int i = 0; i < k; i++) boost(vx, vy, vz, E[i], px[i], py[i],
                                              pz[i]);
        }
        E[k] = sqrt(p*p + m[k]*m[k]);
        px[k] = p*ux;
        py[k] = p*uy;
        pz[k] = p*uz;
    }

    double life_time = 1e10;
    double M_width = get_particle_width(mother);
    if (M_width > 1e-10) {
        double tau = mother->E/M*1./M_width;
        life_time = -tau*log(ran.rand_uniform())*0.19733;  // fm
    }
    double vx = mother->px/mother->E;
    double vy = mother->py/mother->E;
    double vz = mother->pz/mother->E;
    for (int k = 0; k < n; k++) {
        boost_from_rest_frame(*mother, M, E[k], px[k], py[k], pz[k]);
        daughters[k].E = E[k];
        daughters[k].px = px[k];
        daughters[k].py = py[k];
        daughters[k].pz = pz[k];
        daughters[k].t = mother->t + life_time;
        daughters[k].x = mother->x + vx*life_time;
        daughters[k].y = mother->y + vy*life_time;
        daughters[k].z = mother->z + vz*life_time;
    }
}

//! This function returns the particle mass for a given particle id
double particle_decay::get_particle_mass(int POI_monval) {
    double mass = 0.0;
    unsigned int ipart = 0;
    for (ipart = 0; ipart < resonance_table.size(); ipart++) {
        if (POI_monval == resonance_table[ipart]->monval) {
            mass = resonance_table[ipart]->mass;
            break;
        }
    }
    if (ipart == resonance_table.size() && mass < 1e-10) {
        cout << "[Error]:particle_decay::get_particle_mass(): "
             << "Can not find particle with monval = " << POI_monval
             << endl;
        exit(1);
    }
    return(mass);
}

//! This function perform two body decay
void particle_decay::perform_two_body_decay(iSS_Hadron *mother,
                                            iSS_Hadron *daughter1,
                                            iSS_Hadron *daughter2) {
    perform_two_body_decay(mother, daughter1, daughter2, *ran_gen_ptr);
}


void particle_decay::perform_two_body_decay(iSS_Hadron *mother,
                                            iSS_Hadron *daughter1,
                                            iSS_Hadron *daughter2,
                                            RandomUtil::Random &ran) {
    double M_pole = mother->mass;
    double M_width = get_particle_width(mother);
    double m1 = daughter1->mass;
    double m2 = daughter2->mass;
    double M_min = m1 + m2;
    if (M_pole < M_min) {
        cout << "Error:particleSamples::perform_two_body_decay:"
             << "can not found decays!" << endl;
        cout << "M = " << M_pole << ", m1 = " << m1 << ", m2 = " << m2 << endl;
        exit(1);
    }
    //double M_sampled = sample_breit_wigner(M_pole, M_width, M_min);
    const double M_sampled = rest_frame_mass(*mother, M_min);
    double temp = M_sampled*M_sampled - m1*m1 - m2*m2;
    double p_lrf = sqrt(temp*temp - 4*m1*m1*m2*m2)/(2*M_sampled);

    // randomly pick emission angle
    double phi = ran.rand_uniform()*2*M_PI;
    double cos_theta = 2.*(ran.rand_uniform() - 0.5);
    double sin_theta = sqrt(1. - cos_theta*cos_theta);

    // compute daughter particles' energy and momentum in the rest frame
    double E1_lrf = sqrt(p_lrf*p_lrf + m1*m1);
    double p1_lrf_x = p_lrf*sin_theta*cos(phi);
    double p1_lrf_y = p_lrf*sin_theta*sin(phi);
    double p1_lrf_z = p_lrf*cos_theta;
    double E2_lrf = sqrt(p_lrf*p_lrf + m2*m2);
    double p2_lrf_x = -p1_lrf_x;
    double p2_lrf_y = -p1_lrf_y;
    double p2_lrf_z = -p1_lrf_z;
    
    // boost to the frame of the mother's four-momentum
    boost_from_rest_frame(*mother, M_sampled, E1_lrf, p1_lrf_x, p1_lrf_y,
                          p1_lrf_z);
    boost_from_rest_frame(*mother, M_sampled, E2_lrf, p2_lrf_x, p2_lrf_y,
                          p2_lrf_z);
    daughter1->E = E1_lrf;
    daughter1->px = p1_lrf_x;
    daughter1->py = p1_lrf_y;
    daughter1->pz = p1_lrf_z;
    daughter2->E = E2_lrf;
    daughter2->px = p2_lrf_x;
    daughter2->py = p2_lrf_y;
    daughter2->pz = p2_lrf_z;

    double life_time = 1e10;
    if (M_width > 1e-10) {
        // compute life-time = gamma*1/\Gamma
        double tau0 = mother->E/(M_sampled)*1./(M_width);
        life_time = -tau0*log(ran.rand_uniform())*0.19733;  // convert to fm
    }

    daughter1->t = mother->t + life_time;
    daughter1->x = mother->x + mother->px/mother->E*life_time;
    daughter1->y = mother->y + mother->py/mother->E*life_time;
    daughter1->z = mother->z + mother->pz/mother->E*life_time;
    daughter2->t = mother->t + life_time;
    daughter2->x = mother->x + mother->px/mother->E*life_time;
    daughter2->y = mother->y + mother->py/mother->E*life_time;
    daughter2->z = mother->z + mother->pz/mother->E*life_time;

    return;
}


//! This function perform 3 body decays
void particle_decay::perform_three_body_decay(iSS_Hadron *mother,
                                              iSS_Hadron *daughter1,
                                              iSS_Hadron *daughter2,
                                              iSS_Hadron *daughter3) {
    perform_three_body_decay(mother, daughter1, daughter2, daughter3,
                             *ran_gen_ptr);
}


void particle_decay::perform_three_body_decay(iSS_Hadron *mother,
                                              iSS_Hadron *daughter1,
                                              iSS_Hadron *daughter2,
                                              iSS_Hadron *daughter3,
                                              RandomUtil::Random &ran) {
    double M_pole = mother->mass;
    double M_width = get_particle_width(mother);
    double m1 = daughter1->mass;
    double m2 = daughter2->mass;
    double m3 = daughter3->mass;
    double M_min = m1 + m2 + m3;
    if (M_pole < M_min) {
        cout << "Error:particleSamples::perform_three_body_decay:"
             << "can not found decays!" << endl;
        cout << "M = " << M_pole << ", m1 = " << m1 << ", m2 = " << m2
             << ", m3 = " << m3 << endl;
        cout << "reso: " << mother->pid << ", m1: " << daughter1->pid
             << ", m2: " << daughter2->pid << ", m3: " << daughter3->pid
             << endl;
        exit(1);
    }
    //double M_sampled = sample_breit_wigner(M_pole, M_width, M_min);
    const double M_sampled = rest_frame_mass(*mother, M_min);
    // generate lrf E1, E2, and theta12 using accept and reject method
    double E1_lrf, E2_lrf, E3_lrf, p1_lrf, p2_lrf, cos12_lrf;
    do {
        do {
            E1_lrf = ran.rand_uniform()*(M_sampled - m1 - m2 - m3) + m1;
            E2_lrf = ran.rand_uniform()*(M_sampled - m1 - m2 - m3) + m2;
        } while (E1_lrf + E2_lrf > M_sampled);
        p1_lrf = sqrt(E1_lrf*E1_lrf - m1*m1);
        p2_lrf = sqrt(E2_lrf*E2_lrf - m2*m2);
        E3_lrf = M_sampled - E1_lrf - E2_lrf;
        cos12_lrf = (E3_lrf*E3_lrf - p1_lrf*p1_lrf - p2_lrf*p2_lrf - m3*m3)
                     /(2.*p1_lrf*p2_lrf);
    } while (cos12_lrf < - 1.0 || cos12_lrf > 1.0);
    // now we get the a good sample

    // sample the lifetime
    double life_time = 1e10;
    if (M_width > 1e-10) {
        double tau = mother->E/(M_sampled)*1./M_width;
        life_time = -tau*log(ran.rand_uniform())*0.19733;  // convert unit to fm
    }
    // compute the decay position
    double decay_time = mother->t + life_time;
    double decay_x = mother->x + mother->px/mother->E*life_time;
    double decay_y = mother->y + mother->py/mother->E*life_time;
    double decay_z = mother->z + mother->pz/mother->E*life_time;

    // compute the momentum of decay daughters
    double tp2_lrf_x = p2_lrf*sqrt(1. - cos12_lrf*cos12_lrf);
    double tp2_lrf_z = p2_lrf*cos12_lrf;
    double tp3_lrf_x = - tp2_lrf_x;
    double tp3_lrf_z = - (p1_lrf + tp2_lrf_z); 
    double phi = 2.*M_PI*ran.rand_uniform();
    double ksi = 2.*M_PI*ran.rand_uniform();
    double cos_theta = 2.*ran.rand_uniform() - 1.0;

    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double sin_ksi = sin(ksi);
    double cos_ksi = cos(ksi);
    double sin_theta = sqrt(1. - cos_theta*cos_theta);

    double p1_lrf_x = - p1_lrf*sin_theta*cos_ksi;
    double p1_lrf_y = p1_lrf*sin_theta*sin_ksi;
    double p1_lrf_z = p1_lrf*cos_theta;
    E1_lrf = sqrt(m1*m1 + p1_lrf_x*p1_lrf_x
                  + p1_lrf_y*p1_lrf_y + p1_lrf_z*p1_lrf_z);
    double p2_lrf_x = (tp2_lrf_x*(cos_phi*cos_theta*cos_ksi - sin_phi*sin_ksi)
                       - tp2_lrf_z*sin_theta*cos_ksi);
    double p2_lrf_y = (tp2_lrf_x*(-cos_phi*cos_theta*sin_ksi - sin_phi*cos_ksi)
                       + tp2_lrf_z*sin_theta*sin_ksi);
    double p2_lrf_z = tp2_lrf_x*(cos_phi*sin_theta) + tp2_lrf_z*cos_theta;
    E2_lrf = sqrt(m2*m2 + p2_lrf_x*p2_lrf_x
                  + p2_lrf_y*p2_lrf_y + p2_lrf_z*p2_lrf_z);
    double p3_lrf_x = (tp3_lrf_x*(cos_phi*cos_theta*cos_ksi - sin_phi*sin_ksi)
                       - tp3_lrf_z*sin_theta*cos_ksi);
    double p3_lrf_y = (tp3_lrf_x*(-cos_phi*cos_theta*sin_ksi - sin_phi*cos_ksi)
                       + tp3_lrf_z*(sin_theta*sin_ksi));
    double p3_lrf_z = tp3_lrf_x*cos_phi*sin_theta + tp3_lrf_z*cos_theta;
    E3_lrf = sqrt(m3*m3 + p3_lrf_x*p3_lrf_x
                  + p3_lrf_y*p3_lrf_y + p3_lrf_z*p3_lrf_z);

    boost_from_rest_frame(*mother, M_sampled, E1_lrf, p1_lrf_x, p1_lrf_y,
                          p1_lrf_z);
    boost_from_rest_frame(*mother, M_sampled, E2_lrf, p2_lrf_x, p2_lrf_y,
                          p2_lrf_z);
    boost_from_rest_frame(*mother, M_sampled, E3_lrf, p3_lrf_x, p3_lrf_y,
                          p3_lrf_z);

    daughter1->E = E1_lrf;
    daughter1->px = p1_lrf_x;
    daughter1->py = p1_lrf_y;
    daughter1->pz = p1_lrf_z;
    daughter1->t = decay_time;
    daughter1->x = decay_x;
    daughter1->y = decay_y;
    daughter1->z = decay_z;

    daughter2->E = E2_lrf;
    daughter2->px = p2_lrf_x;
    daughter2->py = p2_lrf_y;
    daughter2->pz = p2_lrf_z;
    daughter2->t = decay_time;
    daughter2->x = decay_x;
    daughter2->y = decay_y;
    daughter2->z = decay_z;

    daughter3->E = E3_lrf;
    daughter3->px = p3_lrf_x;
    daughter3->py = p3_lrf_y;
    daughter3->pz = p3_lrf_z;
    daughter3->t = decay_time;
    daughter3->x = decay_x;
    daughter3->y = decay_y;
    daughter3->z = decay_z;

    return;
}

//! This function sample mother particle mass according to breit-wigner
//! distribution
double particle_decay::sample_breit_wigner(double mass, double width,
                                           double M_min) {
    return sample_breit_wigner(mass, width, M_min, *ran_gen_ptr);
}


double particle_decay::sample_breit_wigner(double mass, double width,
                                           double M_min,
                                           RandomUtil::Random &ran) {
    // this function sample the Breit Wigner distribution for the mass
    // particle mass sampled is >= M_min
    // From Wiki: https://en.wikipedia.org/wiki/Cauchy_distribution
    // The CDF of BreitWigner (or Cauchy) distribution is
    // y = 1/pi arctan((x-x0)/gamma) + 1/2
    // The inverse CDF is 
    // x = x0 + gamma*tan(pi*(y - 1/2))

    // compute the minimum probability given by M_min
    double p_min = atan2((M_min - mass), width/2.)/M_PI + 0.5;
    // generate a random number from p_min to 1
    double y = (1. - p_min)*ran.rand_uniform() + p_min;
    // compute the corresponding sampled mass
    double mass_sampled = mass + width/2.*tan(M_PI*(y - 0.5));
    return(mass_sampled);
}
