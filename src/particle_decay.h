// Copyright Chun Shen @ 2016
#ifndef SRC_particle_decay_h_
#define SRC_particle_decay_h_

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "data_struct.h"
#include "Random.h"

class particle_decay {
 private:
    const std::string table_path_;
    const AfterburnerType afterburner_type_;
    std::vector<particle_info*> resonance_table;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

    //! What a decay channel needs at decay time, fixed when the table is
    //! read: the table index of each daughter, the sum of the daughters'
    //! pole masses (the channel is open at or above it) and the sum of the
    //! lowest masses they can take (below the pole sum the channel can
    //! still be reached by giving its broad daughters masses from their
    //! Breit-Wigner distributions, see sample_daughter_masses()).
    struct ChannelKinematics {
        std::vector<int> daughter_index;
        double threshold_pole;
        double threshold_min;
        bool has_broad_daughter;
    };
    std::vector<std::vector<ChannelKinematics>> channel_kinematics_;
    std::vector<double> lowest_mass_;
    std::unordered_map<int, int> table_index_;
    long n_unreachable_ = 0;
    long n_mass_fallback_ = 0;

    //! fills table_index_, lowest_mass_ and channel_kinematics_
    void prepare_channel_kinematics();
    double compute_lowest_mass(int idx, std::vector<int> &state);
    bool is_broad(int idx) const;
    int find_index(int monval) const;
    void sample_daughter_masses(double M, const ChannelKinematics &channel,
                                std::vector<double> &masses);

 public:
    particle_decay(std::shared_ptr<RandomUtil::Random> ran_gen,
                   AfterburnerType afterburner_type,
                   std::string table_path="iSS_tables");
    ~particle_decay();

    //! This function reads in resonance decay table
    int read_resonances_list();

    //! This is a test function to check whether the resonance table is
    //! read in correctly
    void check_resonance_table();

    //! This function returns particle width in GeV
    double get_particle_width(iSS_Hadron *part);

    //! This function checks whether the particle is stable
    int check_particle_stable(iSS_Hadron *part);

    //! This function returns the electric charge of particle
    int get_particle_charge(int monval);

    //! This function returns the baryon number of particle
    int get_particle_baryon_number(int monval);

    //! This function returns the strange number of particle
    int get_particle_strange_number(int monval);

    //! This function returns the particle mass for a given particle id
    double get_particle_mass(int POI_monval);

    //! This is a shell function to perform resonance decays
    void perform_decays(iSS_Hadron *mother,
                        std::vector<iSS_Hadron>* daughter_list);


    //! This function perform two body decay
    void perform_two_body_decay(iSS_Hadron *mother,
                                iSS_Hadron *daughter1,
                                iSS_Hadron *daughter2);

    //! This function perform 3 body decays
    void perform_three_body_decay(iSS_Hadron *mother,
                                  iSS_Hadron *daughter1,
                                  iSS_Hadron *daughter2,
                                  iSS_Hadron *daughter3);

    //! Uniform n-body phase-space decay (Raubold-Lynch / GENBOD); used for
    //! channels with four or five daughters
    void perform_n_body_decay(iSS_Hadron *mother,
                              std::vector<iSS_Hadron> &daughters);

    //! This function sample mother particle mass according to breit-wigner
    //! distribution
    double sample_breit_wigner(double mass, double width, double M_min);
};

#endif  // SRC_particle_decay_h_
