// Copyright Chun Shen @ 2016
#ifndef SRC_particle_decay_h_
#define SRC_particle_decay_h_

#include <vector>
#include <string>
#include <memory>

#include "data_struct.h"
#include "Random.h"

class particle_decay {
 private:
    const std::string table_path_;
    const AfterburnerType afterburner_type_;
    std::vector<particle_info*> resonance_table;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

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

    //! This function sample mother particle mass according to breit-wigner
    //! distribution
    double sample_breit_wigner(double mass, double width, double M_min);
};

#endif  // SRC_particle_decay_h_
