// Copyright Chun Shen @ 2016

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
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
}

particle_decay::~particle_decay() {
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

//! This function checks whether the particle is stable
int particle_decay::check_particle_stable(iSS_Hadron *part) {
    int stable = 0;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (part->pid == resonance_table[i]->monval) {
            stable = resonance_table[i]->stable;
            break;
        }
    }
    return(stable);
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


//! This is a shell function to perform resonance decays
void particle_decay::perform_decays(
            iSS_Hadron *mother, vector<iSS_Hadron>* daughter_list) {
    particle_info* mother_decay_info = NULL;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (mother->pid == resonance_table[i]->monval) {
            mother_decay_info = resonance_table[i];
            break;
        }
    }
    if (mother_decay_info->stable == 1) {
        // the particle is a stable particle
        daughter_list->push_back(*mother);
        return;
    }
    int N_decay_channel = mother_decay_info->decays;
    double random_local = ran_gen_ptr->rand_uniform();
    double cumulated_branching_ratio = 0.0;
    decay_channel_info* picked_channel = NULL;
    for (int i_channel = 0; i_channel < N_decay_channel; i_channel++) {
        cumulated_branching_ratio += (
            mother_decay_info->decay_channels[i_channel]->branching_ratio);
        if (cumulated_branching_ratio > random_local) {
            picked_channel = mother_decay_info->decay_channels[i_channel];
            break;
        }
    }
    int N_decay_part = picked_channel->decay_Npart;
    if (N_decay_part == 2) {
        iSS_Hadron *daughter1 = new iSS_Hadron;
        iSS_Hadron *daughter2 = new iSS_Hadron;
        int decay_part1_monval = picked_channel->decay_part[0];
        int decay_part2_monval = picked_channel->decay_part[1];
        daughter1->pid = decay_part1_monval;
        daughter2->pid = decay_part2_monval;
        daughter1->mass = get_particle_mass(decay_part1_monval);
        daughter2->mass = get_particle_mass(decay_part2_monval);
        perform_two_body_decay(mother, daughter1, daughter2);
        daughter_list->push_back(*daughter1);
        daughter_list->push_back(*daughter2);
        delete daughter1;
        delete daughter2;
    } else if (N_decay_part == 3) {
        iSS_Hadron *daughter1 = new iSS_Hadron;
        iSS_Hadron *daughter2 = new iSS_Hadron;
        iSS_Hadron *daughter3 = new iSS_Hadron;
        int decay_part1_monval = picked_channel->decay_part[0];
        int decay_part2_monval = picked_channel->decay_part[1];
        int decay_part3_monval = picked_channel->decay_part[2];
        daughter1->pid = decay_part1_monval;
        daughter2->pid = decay_part2_monval;
        daughter3->pid = decay_part3_monval;
        daughter1->mass = get_particle_mass(decay_part1_monval);
        daughter2->mass = get_particle_mass(decay_part2_monval);
        daughter3->mass = get_particle_mass(decay_part3_monval);
        perform_three_body_decay(mother, daughter1, daughter2, daughter3);
        daughter_list->push_back(*daughter1);
        daughter_list->push_back(*daughter2);
        daughter_list->push_back(*daughter3);
        delete daughter1;
        delete daughter2;
        delete daughter3;
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
    double M_sampled = M_pole;
    double temp = M_sampled*M_sampled - m1*m1 - m2*m2;
    double p_lrf = sqrt(temp*temp - 4*m1*m1*m2*m2)/(2*M_sampled);

    // randomly pick emission angle
    double phi = ran_gen_ptr->rand_uniform()*2*M_PI;
    double cos_theta = 2.*(ran_gen_ptr->rand_uniform() - 0.5);
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
    
    // compute mother velocity
    double vx = mother->px/mother->E;
    double vy = mother->py/mother->E;
    double vz = mother->pz/mother->E;

    // perform the boost
    double v2 = vx*vx + vy*vy + vz*vz;
    double gamma = 1./sqrt(1. - v2);
    double gamma_m_1 = gamma - 1.;
    double vp1 = vx*p1_lrf_x + vy*p1_lrf_y + vz*p1_lrf_z;
    double vp2 = vx*p2_lrf_x + vy*p2_lrf_y + vz*p2_lrf_z;
    daughter1->E = gamma*(E1_lrf + vp1);
    daughter1->px = p1_lrf_x + (gamma_m_1*vp1/v2 + gamma*E1_lrf)*vx;
    daughter1->py = p1_lrf_y + (gamma_m_1*vp1/v2 + gamma*E1_lrf)*vy;
    daughter1->pz = p1_lrf_z + (gamma_m_1*vp1/v2 + gamma*E1_lrf)*vz;
    daughter2->E = gamma*(E2_lrf + vp2);
    daughter2->px = p2_lrf_x + (gamma_m_1*vp2/v2 + gamma*E2_lrf)*vx;
    daughter2->py = p2_lrf_y + (gamma_m_1*vp2/v2 + gamma*E2_lrf)*vy;
    daughter2->pz = p2_lrf_z + (gamma_m_1*vp2/v2 + gamma*E2_lrf)*vz;

    double life_time = 1e10;
    if (M_width > 1e-10) {
        // compute life-time = gamma*1/\Gamma
        double tau0 = mother->E/(M_sampled)*1./(M_width);
        life_time = -tau0*log(ran_gen_ptr->rand_uniform())*0.19733;  // convert to fm
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
    double M_sampled = M_pole;
    // generate lrf E1, E2, and theta12 using accept and reject method
    double E1_lrf, E2_lrf, E3_lrf, p1_lrf, p2_lrf, cos12_lrf;
    do {
        do {
            E1_lrf = ran_gen_ptr->rand_uniform()*(M_sampled - m1 - m2 - m3) + m1;
            E2_lrf = ran_gen_ptr->rand_uniform()*(M_sampled - m1 - m2 - m3) + m2;
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
        life_time = -tau*log(ran_gen_ptr->rand_uniform())*0.19733;  // convert unit to fm
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
    double phi = 2.*M_PI*ran_gen_ptr->rand_uniform();
    double ksi = 2.*M_PI*ran_gen_ptr->rand_uniform();
    double cos_theta = 2.*ran_gen_ptr->rand_uniform() - 1.0;

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

    double vx = mother->px/mother->E;
    double vy = mother->py/mother->E;
    double vz = mother->pz/mother->E;
    double v2 = vx*vx + vy*vy + vz*vz;
    double gamma = 1./sqrt(1. - v2);
    double gamma_m_1 = gamma - 1.;
    double vp1 = vx*p1_lrf_x + vy*p1_lrf_y + vz*p1_lrf_z;
    double vp2 = vx*p2_lrf_x + vy*p2_lrf_y + vz*p2_lrf_z;
    double vp3 = vx*p3_lrf_x + vy*p3_lrf_y + vz*p3_lrf_z;

    daughter1->E = gamma*(E1_lrf + vp1);
    daughter1->px = p1_lrf_x + (gamma_m_1*vp1/v2 + gamma*E1_lrf)*vx;
    daughter1->py = p1_lrf_y + (gamma_m_1*vp1/v2 + gamma*E1_lrf)*vy;
    daughter1->pz = p1_lrf_z + (gamma_m_1*vp1/v2 + gamma*E1_lrf)*vz;
    daughter1->t = decay_time;
    daughter1->x = decay_x;
    daughter1->y = decay_y;
    daughter1->z = decay_z;

    daughter2->E = gamma*(E2_lrf + vp2);
    daughter2->px = p2_lrf_x + (gamma_m_1*vp2/v2 + gamma*E2_lrf)*vx;
    daughter2->py = p2_lrf_y + (gamma_m_1*vp2/v2 + gamma*E2_lrf)*vy;
    daughter2->pz = p2_lrf_z + (gamma_m_1*vp2/v2 + gamma*E2_lrf)*vz;
    daughter2->t = decay_time;
    daughter2->x = decay_x;
    daughter2->y = decay_y;
    daughter2->z = decay_z;

    daughter3->E = gamma*(E3_lrf + vp3);
    daughter3->px = p3_lrf_x + (gamma_m_1*vp3/v2 + gamma*E3_lrf)*vx;
    daughter3->py = p3_lrf_y + (gamma_m_1*vp3/v2 + gamma*E3_lrf)*vy;
    daughter3->pz = p3_lrf_z + (gamma_m_1*vp3/v2 + gamma*E3_lrf)*vz;
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
    double y = (1. - p_min)*ran_gen_ptr->rand_uniform() + p_min;
    // compute the corresponding sampled mass
    double mass_sampled = mass + width/2.*tan(M_PI*(y - 0.5));
    return(mass_sampled);
}
