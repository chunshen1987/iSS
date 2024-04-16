// Copyright @ 2020 Chun Shen

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_lambert.h>

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <array>

#include "zlib.h"
#include "data_struct.h"
#include "readindata.h"
#include "FSSW.h"
#include "RandomVariable1DArray.h"
#include "RandomVariable2DArray.h"
#include "ParameterReader.h"
#include "arsenal.h"
#include "Stopwatch.h"

using iSS_data::AMOUNT_OF_OUTPUT;
using iSS_data::VecD4;
using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::setw;
using iSS_data::hbarC;


//***************************************************************************
FSSW::FSSW(std::shared_ptr<RandomUtil::Random> ran_gen,
           Table* chosen_particles_in,
           std::vector<particle_info> particles_in,
           const std::vector<FO_surf_LRF> &FOsurf_ptr_in, int flag_PCE,
           ParameterReader* paraRdr_in, string path, string table_path,
           AfterburnerType afterburner_type) :
               path_(path), table_path_(table_path),
               afterburner_type_(afterburner_type), FOsurf_ptr(FOsurf_ptr_in) {

    ran_gen_ptr = ran_gen;
    momentum_sampler_ptr_ = std::shared_ptr<MomentumSamplerShell> (
            new MomentumSamplerShell(ran_gen_ptr));

    // get info
    flag_PCE_ = flag_PCE;

    particles = particles_in;
    int Nparticles = particles.size();

    //FOsurf_ptr = FOsurf_ptr_in;
    FO_length = FOsurf_ptr.size();

    paraRdr = paraRdr_in;
    echoLevel_ = paraRdr->getVal("JSechoLevel", 1);

    hydro_mode = paraRdr->getVal("hydro_mode");

    USE_OSCAR_FORMAT         = paraRdr->getVal("use_OSCAR_format");
    USE_GZIP_FORMAT          = paraRdr->getVal("use_gzip_format");
    USE_BINARY_FORMAT        = paraRdr->getVal("use_binary_format");
    INCLUDE_DELTAF           = paraRdr->getVal("include_deltaf_shear");
    INCLUDE_BULK_DELTAF      = paraRdr->getVal("include_deltaf_bulk");
    bulk_deltaf_kind_        = paraRdr->getVal("bulk_deltaf_kind");
    INCLUDE_DIFFUSION_DELTAF = paraRdr->getVal("include_deltaf_diffusion");
    if (paraRdr->getVal("RegVisYield", 0) == 1) {
        flagRegVisYield_ = true;
    } else {
        flagRegVisYield_ = false;
    }

    if (bulk_deltaf_kind_ == 21) {
        NEoS_deltaf_kind_ = 1;
    } else if (bulk_deltaf_kind_ == 20) {
        NEoS_deltaf_kind_ = 0;
    } else {
        NEoS_deltaf_kind_ = -1;
    }

    local_charge_conservation = paraRdr->getVal("local_charge_conservation");
    number_of_repeated_sampling = (
            paraRdr->getVal("number_of_repeated_sampling"));
    maximumSamplingEvents_ = paraRdr->getVal("maximum_sampling_events");

    flag_output_samples_into_files = (
                            paraRdr->getVal("output_samples_into_files"));
    flag_store_samples_in_memory = 1;
    Hadron_list = new vector< vector<iSS_Hadron>* >;

    if (paraRdr->getVal("perform_decays") == 1) {
        flag_perform_decays_ = true;
        decayer_ptr_ = std::unique_ptr<particle_decay>(new particle_decay(
                                ran_gen_ptr, afterburner_type_, table_path_));
    } else {
        flag_perform_decays_ = false;
    }

    int flag_include_spectators = paraRdr->getVal("include_spectators", 0);
    if (flag_include_spectators != 0) {
        flag_spectators_ = true;
        spectators_ptr_ = std::unique_ptr<Spectators>(new Spectators(
                                                flag_include_spectators));
        std::ostringstream spectatorFilename;
        spectatorFilename << path_ << "/spectators.dat";
        spectators_ptr_->readInSpectatorsFromFile(spectatorFilename.str());
    } else {
        flag_spectators_ = false;
    }

    // deal with chosen_particle_xxx tables
    number_of_chosen_particles = chosen_particles_in->getNumberOfRows();
    chosen_particles_sampling_table.resize(number_of_chosen_particles, 0);
    std::vector<int> unidentifiedPid_table;
    // first copy the chosen_particles table, but now using indecies 
    // instead of monval
    int current_idx = 0;
    for (int m = 0; m < number_of_chosen_particles; m++) {
        int monval = chosen_particles_in->get(1, m+1);
        for(int n = 0; n < Nparticles; n++) {
            if (particles[n].monval == monval) {
                chosen_particles_sampling_table[current_idx] = n;
                current_idx ++;
                break;
            } else if (n == Nparticles - 1) {
                unidentifiedPid_table.push_back(monval);
            }
        }
    }
    // check whether all chosen particles are in the particle list
    if (number_of_chosen_particles != current_idx) {
        messager_ << "not all chosen particles are "
                 << "in the pdg particle list!";
        messager_.flush("warning");
        messager_ << "There are " << number_of_chosen_particles - current_idx
                 << " particles can not be found in the pdg particle list!";
        messager_.flush("warning");
        messager_ << "Their monte carlo numbers are:";
        messager_.flush("warning");
        for (auto const &ipart : unidentifiedPid_table){
            messager_ << ipart;
            messager_.flush("warning");
        }
        number_of_chosen_particles = current_idx;
    }

    // sort particle list by its mass
    for (int m = 0; m < number_of_chosen_particles; m++) {
        for (int n = 0; n < number_of_chosen_particles-m-1; n++)
            if (particles[chosen_particles_sampling_table[n]].mass 
                > particles[chosen_particles_sampling_table[n+1]].mass) {
            // swap them
            int particle_idx = chosen_particles_sampling_table[n+1];
            chosen_particles_sampling_table[n+1] = (
                                        chosen_particles_sampling_table[n]);
            chosen_particles_sampling_table[n] = particle_idx;
        }
    }

    OSCAR_header_filename_ = table_path_ + "/OSCAR_header.txt";
    OSCAR_output_filename_ = "OSCAR.DAT";

    dN_dxtdy_for_one_particle_species.resize(FO_length, 0.);

    gsl_rng_env_setup();
    gsl_type_random_number = gsl_rng_default;
    gsl_random_r = gsl_rng_alloc(gsl_type_random_number);
    gsl_rng_set(gsl_random_r, ran_gen_ptr->get_seed());

    // arrays for bulk delta f coefficients
    if (INCLUDE_BULK_DELTAF == 1) {
        if (bulk_deltaf_kind_ == 0) {
            bulkdf_coeff_ = std::unique_ptr<Table> (new Table (
                table_path_
                + "/deltaf_tables/BulkDf_Coefficients_Hadrons_s95p-v0-PCE.dat"));
        } else if (bulk_deltaf_kind_ == 11) {
            load_bulk_deltaf_14mom_table(table_path_ + "/deltaf_tables");
        }
    }

    if (NEoS_deltaf_kind_ == 1) {
        load_CE_deltaf_NEOSBQS_table(table_path_ + "/deltaf_tables");
    } else if (NEoS_deltaf_kind_ == 0) {
        load_22mom_deltaf_NEOSBQS_table(table_path_ + "/deltaf_tables");
    }

    // load table for diffusion delta f coeffient
    if (INCLUDE_DIFFUSION_DELTAF == 1) {
        load_deltaf_qmu_coeff_table(
            table_path_ + "/deltaf_tables/Coefficients_RTA_diffusion.dat");
    }

    // create arrays for special functions who are needed to compute 
    // particle yields
    initialize_special_function_arrays();

}
//***************************************************************************


//***************************************************************************
FSSW::~FSSW() {
    if (INCLUDE_BULK_DELTAF == 1) {
        if (bulk_deltaf_kind_ == 11) {
            for (int i = 0; i < deltaf_bulk_coeff_14mom_table_length_T_; i++) {
                delete[] deltaf_bulk_coeff_14mom_c0_tb_[i];
                delete[] deltaf_bulk_coeff_14mom_c1_tb_[i];
                delete[] deltaf_bulk_coeff_14mom_c2_tb_[i];
            }
            delete deltaf_bulk_coeff_14mom_c0_tb_;
            delete deltaf_bulk_coeff_14mom_c1_tb_;
            delete deltaf_bulk_coeff_14mom_c2_tb_;
        }
    }

    if (INCLUDE_DIFFUSION_DELTAF == 1) {
        for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++) {
            delete [] deltaf_qmu_coeff_tb[i];
        }
        delete [] deltaf_qmu_coeff_tb;
    }

    gsl_rng_free(gsl_random_r);

    // clean arrays for special functions
    for (int i = 0; i < sf_tb_length; i++) {
        delete [] sf_bessel_Kn[i];
        if (INCLUDE_DIFFUSION_DELTAF == 1) {
            delete [] sf_expint_En[i];
        }
    }
    delete [] sf_bessel_Kn;
    if (INCLUDE_DIFFUSION_DELTAF == 1) {
        delete [] sf_expint_En;
    }

    for (unsigned int i = 0; i < Hadron_list->size(); i++) {
        (*Hadron_list)[i]->clear();
    }
    Hadron_list->clear();
}
//***************************************************************************


//***************************************************************************
inline long FSSW::determine_number_to_sample(
    double dN_dy_in, int model, double para1) {
// From a non-integer averaged particles number dN, return an actual interger
// particle number that can be used in sampling.
    if (dN_dy_in < 0) {
        cout << "FSSW::"
             << "determine_number_to_sample error: "
             << "dN_dy should be positive but receives " << dN_dy_in << endl;
        exit(-1);
    }
    double dN_dy = dN_dy_in;
    long number_to_sample;
    long dN_dy_int = long(dN_dy); // integer part
    double dN_dy_fraction = dN_dy - dN_dy_int; // fractional part
    double p,k; // parameters in NBD

    // determine actual number of particles
    switch (model) // explained in parameters.dat
    {
        case 1: // with possibly 1 more particle
            number_to_sample = dN_dy_int;
            // lucky! 1 more particle...
            if (ran_gen_ptr->rand_uniform() < dN_dy_fraction)
                number_to_sample++; 
            break;
        case 10: // use NBD
            k = para1*dN_dy_fraction;
            p = 1.0/(1.0+para1);
            if (k<1e-15)
                number_to_sample = dN_dy_int;
            else 
                number_to_sample = (
                    dN_dy_int + gsl_ran_negative_binomial(gsl_random_r, p, k));
            break;
        case 20:
            k = para1*dN_dy;
            p = 1.0/(1.0+para1);
            if (k<1e-15) 
                number_to_sample = dN_dy_int;
            else 
                number_to_sample = (
                                gsl_ran_negative_binomial(gsl_random_r, p, k));
            break;
        case 30:  // use Poisson distribution
            if (dN_dy < 1e-15)
                number_to_sample = 0;
            else
                number_to_sample = gsl_ran_poisson(gsl_random_r, dN_dy);
            break;
        default:
            cout << "FSSW::"
                 << "determine_number_to_sample error: "
                 << "the model specified (" << model << ") "
                 << "to determine the number of particles is not found" 
                 << endl;
            exit(-1);
    }

    return number_to_sample;
}
//***************************************************************************


//***************************************************************************
bool FSSW::particles_are_the_same(int idx1, int idx2) {
    if (particles[idx1].sign != particles[idx2].sign)
        return false;
    if (particles[idx1].gspin != particles[idx2].gspin)
        return false;
    if (particles[idx1].baryon != particles[idx2].baryon)
        return false;
    if (particles[idx1].strange != particles[idx2].strange)
        return false;
    if (particles[idx1].charge != particles[idx2].charge)
        return false;
    double tolerance = paraRdr->getVal("grouping_tolerance");
    if (std::abs((particles[idx1].mass-particles[idx2].mass)
            /(particles[idx2].mass+1e-30)) > tolerance)
        return false;
    if (flag_PCE_ == 1) {
        for (long l = 0; l < FO_length; l++) {
            double chem1 = FOsurf_ptr[l].particle_mu_PCE[idx1];
            double chem2 = FOsurf_ptr[l].particle_mu_PCE[idx2];
            if (std::abs((chem1-chem2)/(chem2+1e-30)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}



//***************************************************************************
void FSSW::shell() {
    sample_using_dN_dxtdy_4all_particles_conventional();
    if (flag_perform_decays_ && afterburner_type_ != AfterburnerType::SMASH) {
        perform_resonance_feed_down(Hadron_list);
    }
    if (flag_spectators_)
        addSpectatorsToHadronList();

    computeAvgTotalEnergyMomentum();

    if (USE_OSCAR_FORMAT) {
        combine_samples_to_OSCAR();
    } else if (USE_GZIP_FORMAT) {
        combine_samples_to_gzip_file();
    } else if (USE_BINARY_FORMAT) {
        combine_samples_to_binary_file();
    }
}


//***************************************************************************
void FSSW::combine_samples_to_OSCAR() {
    Stopwatch sw;
    sw.tic();
    messager_.info(" -- Now combine sample files to OSCAR file...");

    char line_buffer[500];

    // open file for output
    remove(OSCAR_output_filename_.c_str());
    ofstream oscar(OSCAR_output_filename_.c_str());

    // write header first
    ifstream header(OSCAR_header_filename_.c_str());
    if (!header.is_open()) {
        cout << endl 
            << "combine_samples_to_OSCAR error: OSCAR header file " 
            << OSCAR_header_filename_.c_str() << " not found." << endl;
        exit(-1);
    }
    while (true) {
        header.getline(line_buffer, 500);
        if (!header.eof()) {
            oscar << line_buffer << endl;
        } else {
            break;
        }
    }
    header.close();

    if (flag_store_samples_in_memory == 0) {
        // open control and sample files
        vector<ifstream*> controls(number_of_chosen_particles);
        vector<ifstream*> samples(number_of_chosen_particles);
        for (int m = 0; m < number_of_chosen_particles; m++) {
            int monval = particles[chosen_particles_sampling_table[m]].monval;
            // control files first
            std::stringstream samples_control_filename;
            samples_control_filename << path_ << "/samples_control_"
                                     << monval << ".dat";
            controls[m] = new ifstream ;
            controls[m]->open(samples_control_filename.str().c_str());
            if (!controls[m]->is_open()) {
                cout << endl
                     << "combine_samples_to_OSCAR error: control file "
                     << samples_control_filename.str() << " not found."
                     << endl;
                exit(-1);
            }
            std::stringstream samples_filename;
            samples_filename << path_ << "/samples_" << monval << ".dat";
            samples[m] = new ifstream ;
            samples[m]->open(samples_filename.str().c_str());
            if (!samples[m]->is_open()) {
                cout << endl
                     << "combine_samples_to_OSCAR error: sample file "
                     << samples_filename.str() << " not found." << endl;
                exit(-1);
            }
        }

        // big loop for generating OSCAR
        if (AMOUNT_OF_OUTPUT > 7) print_progressbar(-1);
        for (long sample_idx=1; sample_idx<=number_of_repeated_sampling; 
             sample_idx++) {
            // read-in number of particles for each species
            int number_of_particles[number_of_chosen_particles];
            long total_number_of_particles = 0;
            for (int m=0; m<number_of_chosen_particles; m++) {
                (*controls[m]) >> number_of_particles[m];
                total_number_of_particles += number_of_particles[m];
            }
            // sub-header for each event
            oscar << setw(10) << sample_idx << "  "
                  << setw(10) << total_number_of_particles << "  "
                  << setw(8) << 0.0 << "  " << setw(8) << 0.0 << endl;

            // now copy each line from samples file to OSCAR file
            long ipart = 1;
            for (int m = 0; m < number_of_chosen_particles; m++) {
                int monval = (
                        particles[chosen_particles_sampling_table[m]].monval);
                for (long ii = 1; ii <= number_of_particles[m]; ii++) {
                    oscar << setw(10) << ipart << "  "
                          << setw(10) << monval << "  ";
                    samples[m]->getline(line_buffer, 500);
                    oscar << line_buffer << endl;
                    ipart ++;
                }
            }
            if (AMOUNT_OF_OUTPUT > 7) {
                print_progressbar(static_cast<double>(
                            sample_idx/number_of_repeated_sampling));
            }
        }
    } else {
        for (unsigned int iev = 0; iev < Hadron_list->size(); iev++) {
            int total_number_of_particles = (*Hadron_list)[iev]->size();
            if (total_number_of_particles > 0) {
                // sub-header for each event
                oscar << setw(10) << iev << "  "
                      << setw(10) << total_number_of_particles << "  "
                      << setw(8) << 0.0 << "  " << setw(8) << 0.0 << endl;
                for (int ipart = 0; ipart < total_number_of_particles;
                     ipart++) {
                    oscar << setw(10) << ipart + 1 << "  "
                          << setw(10) << (*(*Hadron_list)[iev])[ipart].pid
                          << "  ";
                    sprintf(line_buffer,
                            "%24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e",
                            (*(*Hadron_list)[iev])[ipart].px,
                            (*(*Hadron_list)[iev])[ipart].py,
                            (*(*Hadron_list)[iev])[ipart].pz,
                            (*(*Hadron_list)[iev])[ipart].E,
                            (*(*Hadron_list)[iev])[ipart].mass,
                            (*(*Hadron_list)[iev])[ipart].x,
                            (*(*Hadron_list)[iev])[ipart].y,
                            (*(*Hadron_list)[iev])[ipart].z,
                            (*(*Hadron_list)[iev])[ipart].t);
                    oscar << line_buffer << endl;
                }
            }
        }
    }

    sw.toc();
    cout << endl
         << " -- combine_samples_to_OSCAR samples finishes " 
         << sw.takeTime() << " seconds."
         << endl;
}


void FSSW::combine_samples_to_gzip_file() {
    Stopwatch sw;
    sw.tic();
    messager_.info(" -- Now combine sample files to a gzip file...");

    // open file for output
    std::string gzip_output_filename = "particle_samples.gz";
    remove(gzip_output_filename.c_str());
    gzFile fp_gz = gzopen(gzip_output_filename.c_str(), "wb");

    for (auto const &ev_i: (*Hadron_list)) {
        int total_number_of_particles = ev_i->size();
        gzprintf(fp_gz, "%d \n", total_number_of_particles);
        for (auto &part_i: (*ev_i)) {
            gzprintf(fp_gz, "%d ", part_i.pid);
            gzprintf(fp_gz,
                     "%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
                     part_i.mass,
                     part_i.t, part_i.x, part_i.y, part_i.z,
                     part_i.E, part_i.px, part_i.py, part_i.pz);
        }
    }
    gzclose(fp_gz);

    sw.toc();
    cout << endl
         << " -- combine_samples_to_gzip_file finishes "
         << sw.takeTime() << " seconds."
         << endl;
}


void FSSW::combine_samples_to_binary_file() {
    Stopwatch sw;
    sw.tic();
    messager_.info(" -- Now combine sample files to a binary file...");

    // open file for output
    std::string binary_output_filename = "particle_samples.bin";
    remove(binary_output_filename.c_str());
    FILE *outbin = NULL;
    outbin = fopen(binary_output_filename.c_str(), "wb");
    for (auto const &ev_i: (*Hadron_list)) {
        int total_number_of_particles = ev_i->size();
        fwrite(&total_number_of_particles, sizeof(int), 1, outbin);
        for (auto &part_i: (*ev_i)) {
            fwrite(&part_i.pid, sizeof(int), 1, outbin);
            float array[] = {
                static_cast<float>(part_i.mass),
                static_cast<float>(part_i.t), static_cast<float>(part_i.x),
                static_cast<float>(part_i.y), static_cast<float>(part_i.z),
                static_cast<float>(part_i.E), static_cast<float>(part_i.px),
                static_cast<float>(part_i.py), static_cast<float>(part_i.pz),
            };
            fwrite(array, sizeof(float), 9, outbin);
        }
    }
    fclose(outbin);

    sw.toc();
    cout << endl
         << " -- combine_samples_to_binary_file finishes "
         << sw.takeTime() << " seconds."
         << endl;
}


//***************************************************************************
void FSSW::calculate_dN_dxtdy_for_one_particle_species(
                                                const int real_particle_idx) {
/*
   The p_integral_table is a table stores results for the integral
     int( p^2 / ( exp(sqrt(p^2+m^2)-mu) + 1), p=0..inf)
   and m_integral_table is a table stores results for the integral
     int( p^2 / ( exp(sqrt(p^2+m^2)-mu) - 1), p=0..inf)
   for different m and mu values.
   The m values changes with row indices and it starts with
   integral_table_m0 with step integral_table_dm.
   The m-mu values changes with column indices and it starts with
   integral_table_m_minus_mu0 with step integral_table_dm_minus_mu.
   Note that here m and mu represent in fact m/T and mu/T so they are
   unitless.
*/
    Stopwatch sw;
    sw.tic();

    std::vector<double> bulkvisCoefficients;

    // now loop over all freeze-out cells and particles
    const double unit_factor = 1.0/pow(hbarC, 3);  // unit: convert to unitless

    // loop over all the fluid cells
    for (long l = 0; l < FO_length; l++) {
        const FO_surf_LRF* surf = &FOsurf_ptr[l];
        const double temp = surf->Tdec;
        const double mu_B = surf->muB;

        double dsigma_dot_u = surf->da_mu_LRF[0];

        // bulk delta f contribution
        double bulkPi = 0.0;
        if (NEoS_deltaf_kind_ == 1) {
            getCENEOSBQSCoefficients(surf->Edec, surf->Bn,
                                     bulkvisCoefficients);
        } else if (NEoS_deltaf_kind_ == 0) {
            get22momNEOSBQSCoefficients(surf->Edec, surf->Bn,
                                        bulkvisCoefficients);
        }

        if (INCLUDE_BULK_DELTAF == 1) {
            if (bulk_deltaf_kind_ == 21 || bulk_deltaf_kind_ == 20) {
                bulkPi = surf->bulkPi;    // GeV/fm^3
            } else if (bulk_deltaf_kind_ == 11) {
                bulkPi = surf->bulkPi;    // GeV/fm^3
                getbulkvisCoefficients(temp, mu_B, bulkvisCoefficients);
            } else {
                if (bulk_deltaf_kind_ == 0) {
                    bulkPi = surf->bulkPi;        // unit in GeV/fm^3
                } else {
                    bulkPi = surf->bulkPi/hbarC;  // unit in fm^-4
                }
                getbulkvisCoefficients(temp, bulkvisCoefficients);
            }
        }

        // diffusion delta f
        double dsigma_dot_q = 0.0;
        double deltaf_qmu_coeff = 1.0;
        double prefactor_qmu = 0.0;
        if (INCLUDE_DIFFUSION_DELTAF == 1) {
            dsigma_dot_q = (  surf->qmuLRF_x*surf->da_mu_LRF[1]
                            + surf->qmuLRF_y*surf->da_mu_LRF[2]
                            + surf->qmuLRF_z*surf->da_mu_LRF[3]);

            deltaf_qmu_coeff = get_deltaf_qmu_coeff(temp, mu_B);

            double rho_B = surf->Bn;
            double Edec = surf->Edec;
            double Pdec = surf->Pdec;
            prefactor_qmu = rho_B/(Edec + Pdec);  // 1/GeV
        }

        // calculate dN / (dxt dy) for all particles
        double total_N = 0;
        const particle_info *particle = &particles[real_particle_idx];

        //int sign = particle->sign;
        const int degen = particle->gspin;
        //double mass = particle->mass;
        const int baryon  = particle->baryon;
        const int strange = particle->strange;
        const int charge  = particle->charge;
        const double mass = particle->mass;
        double mu = baryon*surf->muB + strange*surf->muS + charge*surf->muQ;
        if (flag_PCE_ == 1) {
            double mu_PCE = surf->particle_mu_PCE[real_particle_idx];
            mu += mu_PCE;
        }

        double prefactor = degen/(2.*M_PI*M_PI);

        // calculate dN / (dxt dy)
        std::array<double, 6> results_ptr = {0.0};
        calculate_dN_analytic(particle, mu, temp, results_ptr);

        double N_eq = unit_factor*prefactor*dsigma_dot_u*results_ptr[0];

        double deltaN_bulk = 0.0;
        if (INCLUDE_BULK_DELTAF == 1) {
            if (bulk_deltaf_kind_ == 1) {
                deltaN_bulk = (unit_factor*prefactor*dsigma_dot_u
                               *(- bulkPi*bulkvisCoefficients[0])
                               *(- bulkvisCoefficients[1]*results_ptr[1]
                                 + results_ptr[2]));
            } else if (bulk_deltaf_kind_ == 11) {
                deltaN_bulk = (
                    unit_factor*prefactor*dsigma_dot_u*bulkPi
                    *(  results_ptr[1]*mass*mass*bulkvisCoefficients[0]
                      + results_ptr[2]*baryon*bulkvisCoefficients[1]
                      + results_ptr[3]*bulkvisCoefficients[2]));
            } else if (bulk_deltaf_kind_ == 21) {
                // WSU Chapman-Enskog for NEoS BQS
                deltaN_bulk = (unit_factor*prefactor*dsigma_dot_u
                               *(- bulkPi*bulkvisCoefficients[0])
                               *(- bulkvisCoefficients[1]*results_ptr[1]
                                 + results_ptr[2]));
            } else if (bulk_deltaf_kind_ == 20) {
                // WSU 22-mom for NEoS BQS
                deltaN_bulk = (unit_factor*prefactor*dsigma_dot_u*bulkPi
                    *(  results_ptr[1]*mass*mass*bulkvisCoefficients[2]
                      + results_ptr[2]*(  baryon*bulkvisCoefficients[3]
                                        + strange*bulkvisCoefficients[4]
                                        + charge*bulkvisCoefficients[5])
                      + results_ptr[3]*(  bulkvisCoefficients[1]
                                        - bulkvisCoefficients[2]))
                );
            }
        }

        double deltaN_qmu = 0.0;
        if (INCLUDE_DIFFUSION_DELTAF == 1) {
            deltaN_qmu = (unit_factor*prefactor
                          *dsigma_dot_q/deltaf_qmu_coeff
                          *(- prefactor_qmu*results_ptr[4]
                            - baryon*results_ptr[5]));
        }

        total_N = N_eq + deltaN_bulk + deltaN_qmu;
        if (flagRegVisYield_) {
            total_N = std::min(2.*N_eq, total_N);
        }

        dN_dxtdy_for_one_particle_species[l] = std::max(0., total_N);
    }

    sw.toc();
    if (AMOUNT_OF_OUTPUT > 2) {
        cout << endl
             << " -- calculate_dN_dxtdy_for_one_particle_species finished in "
             << sw.takeTime() << " seconds." << endl;
    }
}


//***************************************************************************
void FSSW::calculate_dN_analytic(
        const particle_info* particle, double mu, double Temperature,
        std::array<double, 6> &results) {
/* calculate particle yield using analytic formula as a series sum of Bessel 
   functions. The particle yield emitted from a given fluid cell is 
   \dsigma_mu u^\mu times the return value from this function
*/
    double N_eq = 0.0;                  // equilibrium contribution
    double deltaN_bulk_term1 = 0.0;     // contribution from bulk delta f
    double deltaN_bulk_term2 = 0.0;     // contribution from bulk delta f
    double deltaN_bulk_term3 = 0.0;     // contribution from bulk delta f
    double deltaN_qmu_term1 = 0.0;      // contribution from baryon diffusion
    double deltaN_qmu_term2 = 0.0;      // contribution from baryon diffusion


    const int sign    = particle->sign;
    const double mass = particle->mass;
    const double beta = 1./Temperature;

    double lambda = exp(beta*mu);  // fugacity factor

    int truncate_order;
    if (mass < 0.7 && Temperature > 0.05) {
        truncate_order = 10.;
    } else {
        truncate_order = 1.;
    }
    // compute the sum in the series
    for (int n = 1; n <= truncate_order; n++) {
        double arg = n*mass*beta;  // argument inside bessel functions

        double theta = pow(-sign, n-1);
        double fugacity = pow(lambda, n);
        double K_2 = get_special_function_K2(arg);

        N_eq += theta/n*fugacity*K_2;
        if (std::isnan(N_eq)) {
            cout << "N_eq is nan"
                 << ", theta = " << theta << ", n = " << n
                 << ", T = " << Temperature << "GeV, mu = " << mu << " GeV"
                 << ", exp(mu/T) = " << lambda
                 << ", exp(n*mu/T) = " << fugacity
                 << ", n*m/T = " << arg
                 << ", K_2 = " << K_2 << endl;
            exit(1);
        }

        if (INCLUDE_BULK_DELTAF == 1) {
            double K_1 = get_special_function_K1(arg);
            if (bulk_deltaf_kind_ == 1 || bulk_deltaf_kind_ == 21) {
                // Chapman-Enskog
                deltaN_bulk_term1 += theta*fugacity*(mass*beta*K_1 + 3*K_2/n);
                deltaN_bulk_term2 += theta*fugacity*K_1;
            } else if (bulk_deltaf_kind_ == 11 || bulk_deltaf_kind_ == 20) {
                // 14-moment/22-moment
                double K_3 = get_special_function_K3(arg);
                deltaN_bulk_term1 += theta*fugacity*K_2;  // mass
                deltaN_bulk_term2 += theta*fugacity*(mass*beta*K_1 + 3*K_2/n);  // E
                deltaN_bulk_term3 += theta*fugacity*(mass*beta*K_2 + 3*K_3/n);  // E^2
            }
        }

        if (INCLUDE_DIFFUSION_DELTAF == 1) {
            deltaN_qmu_term1 += theta/n*fugacity*K_2;

            std::vector<double> sf_expint_En_ptr(
                                            sf_expint_truncate_order-1, 0.);
            get_special_function_En(arg, sf_expint_En_ptr);

            double mbeta = mass*beta;
            double I_1_n = 0.0;
            double I_1_1 = exp(-arg)/arg*(2./(arg*arg) + 2./arg - 1./2.);
            double I_1_2 = 3./8.*sf_expint_En_ptr[0];
            I_1_n = I_1_1 + I_1_2;

            double double_factorial = 1.;  // record (2k-5)!!
            double factorial = 2.;         // record k! start with 2!
            double factor_2_to_k_power = 4.;      // record 2^k start with 2^2
            for (int k = 3; k <= sf_expint_truncate_order; k++) {
                double_factorial *= (2*k - 5);
                factorial *= k;
                factor_2_to_k_power *= 2;

                double I_1_k = (
                    3.*double_factorial/factor_2_to_k_power/factorial
                    *sf_expint_En_ptr[k-2]);
                I_1_n += I_1_k;
            }
            I_1_n = -(mbeta*mbeta*mbeta)*I_1_n;
            deltaN_qmu_term2 += n*theta*fugacity*I_1_n;
        }
    }

    // equilibrium contribution
    double prefactor_Neq = mass*mass*Temperature;
    N_eq = prefactor_Neq*N_eq;

    // contribution from bulk viscosity
    if (INCLUDE_BULK_DELTAF == 1) {
        if (bulk_deltaf_kind_ == 1 || bulk_deltaf_kind_ == 21) {
            // Chapman-Enskog
            deltaN_bulk_term1 = mass*mass/beta*deltaN_bulk_term1;
            deltaN_bulk_term2 = mass*mass*mass/3.*deltaN_bulk_term2;
            deltaN_bulk_term3 = 0.0;
        } else if (bulk_deltaf_kind_ == 11 || bulk_deltaf_kind_ == 20) {
            // 14-moment/22-moment
            deltaN_bulk_term1 = mass*mass/beta*deltaN_bulk_term1;
            deltaN_bulk_term2 = mass*mass/(beta*beta)*deltaN_bulk_term2;
            deltaN_bulk_term3 = mass*mass*mass/(beta*beta)*deltaN_bulk_term3;
        }
    } else {
        deltaN_bulk_term1 = 0.0;
        deltaN_bulk_term2 = 0.0;
        deltaN_bulk_term3 = 0.0;
    }

    // contribution from baryon diffusion
    if (INCLUDE_DIFFUSION_DELTAF == 1) {
        deltaN_qmu_term1 = mass*mass/(beta*beta)*deltaN_qmu_term1;
        deltaN_qmu_term2 = 1./(3.*beta*beta*beta)*deltaN_qmu_term2;
    }

    // final results
    results[0] = N_eq;
    results[1] = deltaN_bulk_term1;
    results[2] = deltaN_bulk_term2;
    results[3] = deltaN_bulk_term3;
    results[4] = deltaN_qmu_term1;
    results[5] = deltaN_qmu_term2;
}


int FSSW::compute_number_of_sampling_needed(
                                            int number_of_particles_needed) {
    double dNdy_thermal_pion = 0.;
    calculate_dN_dxtdy_for_one_particle_species(1);  // 1 for thermal pion^+
    for (long l = 0; l < FO_length; l++) {
        dNdy_thermal_pion += dN_dxtdy_for_one_particle_species[l];
    }
    // particle_info *particle = &particles[1];
    // cout << "check  Name: " << particle->name
    //      << ", Monte-carlo index: " << particle->monval
    //      << ", dN/dy = " << dNdy_thermal_pion << endl;
    int nev_needed = static_cast<int>(number_of_particles_needed
                                      /(6.*dNdy_thermal_pion));
    if (hydro_mode == 2) {
        nev_needed *= 10;
    }
    nev_needed = std::max(1, std::min(maximumSamplingEvents_, nev_needed));
    return(nev_needed);
}


//***************************************************************************
void FSSW::sample_using_dN_dxtdy_4all_particles_conventional() {
/*
 Sample using the already calculated dN_dxtdy_4all array.
 Different model can be used in the sampling.
*/
    Stopwatch sw_total;
    sw_total.tic();
    if (echoLevel_ > 0) {
        messager_.info(
                " Function sample_using_dN_dxtdy_4all_particles started...");
    }

    std::vector<double> bulkvisCoefficients;

    // control variables
    int sampling_model    = paraRdr->getVal("dN_dy_sampling_model");
    double sampling_para1 = paraRdr->getVal("dN_dy_sampling_para1");

    int flag_sample_upto_desired_particle_number = paraRdr->getVal(
                                    "sample_upto_desired_particle_number");
    if (flag_sample_upto_desired_particle_number == 1) {
        long number_of_particles_needed = paraRdr->getVal(
                                                "number_of_particles_needed");
        number_of_repeated_sampling = compute_number_of_sampling_needed(
                                                number_of_particles_needed);
    }

    if (sampling_model > 100) {
        cout << "FSSW::"
             << "sample_using_dN_dxtdy_4all_particles error: sampling model " 
             << sampling_model << " is not supported." << endl;
        exit(-1);
    }

    if (flag_store_samples_in_memory == 1) {
        for (int isample = 0; isample < number_of_repeated_sampling;
                isample++) {
            Hadron_list->push_back(new vector<iSS_Hadron> );
        }
    }

    // use "conventional" sampling
    // reusable variables
    // dN_dxtdy for 1 particle
    //vector<double> dN_dxtdy_single_particle(FO_length, 0);
    if (echoLevel_ > 0) {
        messager_ << "Sampling using dN/dy with "
                  << "sample_using_dN_dxtdy_4all_particles function.";
        messager_.flush("info");
        messager_ << "number of repeated sampling = "
                  << number_of_repeated_sampling;
        messager_.flush("info");
    }
    for (int n = 0; n < number_of_chosen_particles; n++) {
        int real_particle_idx = chosen_particles_sampling_table[n];
        const particle_info *particle = &particles[real_particle_idx];
        const double mass = particle->mass;
        const int sign    = particle->sign;
        const int baryon  = particle->baryon;
        const int strange = particle->strange;
        const int charge  = particle->charge;
        if (echoLevel_ > 0) {
            messager_ << "Index: " << n << ", Name: " << particle->name
                 << ", Monte-carlo index: " << particle->monval;
            messager_.flush("info");
        }
        if (local_charge_conservation == 1) {
            if (particle->charge < 0) {
                if (echoLevel_ > 0) {
                    messager_ << "local charge conservation is turn on~ "
                              << "Skip the negative charge particles.";
                    messager_.flush("info");
                }
                continue;
            }
        }

        calculate_dN_dxtdy_for_one_particle_species(real_particle_idx);

        // prepare the inverse CDF
        RandomVariable1DArray rand1D(&dN_dxtdy_for_one_particle_species,
                                     ran_gen_ptr, 0);

        // first get total number of particles
        double dN_dy = rand1D.return_sum();

        // get y range for sampling
        const double y_LB = paraRdr->getVal("y_LB");
        const double y_RB = paraRdr->getVal("y_RB");
        double dN;
        if (hydro_mode != 2) {
            dN = (y_RB - y_LB)*dN_dy;
        } else {
            // for (3+1)-d case, dN_dy is total N (summing over all etas)
            dN = dN_dy;
        }
        if (echoLevel_ > 0) {
            messager_ << " -- Sampling using dN_dy=" << dN_dy << ", "
                      << "dN=" << dN << "...";
            messager_.flush("info");
        }

        Stopwatch sw;
        sw.tic();

        if (AMOUNT_OF_OUTPUT > 7) {
            print_progressbar(-1);
        }

        for (long repeated_sampling_idx = 1; 
             repeated_sampling_idx <= number_of_repeated_sampling;
             repeated_sampling_idx++) {
            long number_to_sample = determine_number_to_sample(
                dN, sampling_model, sampling_para1);

            long sample_idx = 1;
            while (sample_idx <= number_to_sample) {
                // first, sample eta and freeze-out cell index
                long FO_idx = rand1D.rand();
                const FO_surf_LRF *surf = &FOsurf_ptr[FO_idx];

                if (NEoS_deltaf_kind_ == 1) {
                    getCENEOSBQSCoefficients(surf->Edec, surf->Bn,
                                             bulkvisCoefficients);
                } else if (NEoS_deltaf_kind_ == 0) {
                    get22momNEOSBQSCoefficients(surf->Edec, surf->Bn,
                                                bulkvisCoefficients);
                }

                if (INCLUDE_BULK_DELTAF == 1) {
                    if (NEoS_deltaf_kind_ == -1) {
                        if (bulk_deltaf_kind_ == 11) {
                            // OSU 14-moment
                            getbulkvisCoefficients(surf->Tdec, surf->muB,
                                                   bulkvisCoefficients);
                        } else {
                            // bulk delta f at mu_B = 0
                            getbulkvisCoefficients(surf->Tdec,
                                                   bulkvisCoefficients);
                        }
                    }
                }

                // diffusion delta f
                double deltaf_qmu_coeff = 1.0;
                if (INCLUDE_DIFFUSION_DELTAF == 1)
                    deltaf_qmu_coeff = get_deltaf_qmu_coeff(surf->Tdec,
                                                            surf->muB);

                // next sample pt and phi
                double pT, phi, y_minus_eta_s;
                int status = sample_momemtum_from_a_fluid_cell(
                                mass, sign, baryon, strange, charge,
                                surf, bulkvisCoefficients, deltaf_qmu_coeff,
                                pT, phi, y_minus_eta_s);

                if (status == 0) {
                    continue;
                } else {
                    sample_idx++;
                }
                //number_of_success++; // to track success rate

                double eta_s = surf->eta;
                if (hydro_mode != 2) {
                    double rap = (y_LB + (y_RB - y_LB)
                                         *ran_gen_ptr->rand_uniform());
                    eta_s = rap - y_minus_eta_s;
                }

                add_one_sampled_particle(
                    repeated_sampling_idx, surf, particle->monval,
                    mass, pT, phi, y_minus_eta_s, eta_s);

                if (local_charge_conservation == 1 && particle->charge > 0) {
                    // sample a negative particle from the same fluid cell
                    double pT2, phi2, y_minus_eta_s2;
                    int status2 = 0;
                    do {
                        status2 = sample_momemtum_from_a_fluid_cell(
                                mass, sign, -baryon, -strange, -charge,
                                surf, bulkvisCoefficients, deltaf_qmu_coeff,
                                pT2, phi2, y_minus_eta_s2);
                    } while (status2 == 0);
                    add_one_sampled_particle(
                        repeated_sampling_idx, surf, -particle->monval,
                        mass, pT2, phi2, y_minus_eta_s2, eta_s);
                }
            }

            if (AMOUNT_OF_OUTPUT > 7) {
                print_progressbar(
                        static_cast<double>(repeated_sampling_idx)
                                            /number_of_repeated_sampling);
            }
        }
        if (AMOUNT_OF_OUTPUT > 7) print_progressbar(1);

        sw.toc();
        if (AMOUNT_OF_OUTPUT > 2) {
            cout << endl << "Sampling finished in " 
                 << sw.takeTime() << " seconds." << endl;
        }

    }   // n; particle loop

    sw_total.toc();
    if (echoLevel_ > 0) {
        cout << endl 
             << "sample_using_dN_dxtdy_4all_particles finished in " 
             << sw_total.takeTime() << " seconds." << endl;
    }
}


void FSSW::getbulkvisCoefficients(
            const double Tdec, std::vector<double> &bulkvisCoefficients) {
    bulkvisCoefficients.resize(3, 0.);
   const double Tdec_fm = Tdec/hbarC;  // [1/fm]
   double Tdec_fm_power[11];    // cache the polynomial power of Tdec_fm
   Tdec_fm_power[1] = Tdec_fm;
   for (int ipower = 2; ipower < 11; ipower++) {
       Tdec_fm_power[ipower] = Tdec_fm_power[ipower-1]*Tdec_fm;
   }
   if (bulk_deltaf_kind_ == 0) {
       // 14 moment expansion
       // load from file

       //B0 [fm^3/GeV^3]
       bulkvisCoefficients[0] = (
                       bulkdf_coeff_->interp(1, 2, Tdec_fm, 5)/pow(hbarC, 3));

       // D0 [fm^3/GeV^2]
       bulkvisCoefficients[1] = (
                       bulkdf_coeff_->interp(1, 3, Tdec_fm, 5)/pow(hbarC, 2));

       // E0 [fm^3/GeV^3]
       bulkvisCoefficients[2] = (
                       bulkdf_coeff_->interp(1, 4, Tdec_fm, 5)/pow(hbarC, 3));

       // parameterization for mu = 0
       // B0[fm^3/GeV^3]
       //bulkvisCoefficients[0] = (
       //              exp(-15.04512474*Tdec_fm + 11.76194266)/pow(hbarC, 3)); 
       // D0 [fm^3/GeV^2]
       //bulkvisCoefficients[1] = (
       //              exp( -12.45699277*Tdec_fm + 11.4949293)/hbarC/hbarC);  
       // E0 [fm^3/GeV^3]
       //bulkvisCoefficients[2] = (
       //            -exp(-14.45087586*Tdec_fm + 11.62716548)/pow(hbarC, 3));
   } else if(bulk_deltaf_kind_ == 1) {  // relaxation type
       // parameterization from JF
       // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV, 
       // do not trust it beyond
       bulkvisCoefficients[0] = (  642096.624265727 
                                 - 8163329.49562861*Tdec_fm_power[1] 
                                 + 47162768.4292073*Tdec_fm_power[2] 
                                 - 162590040.002683*Tdec_fm_power[3] 
                                 + 369637951.096896*Tdec_fm_power[4] 
                                 - 578181331.809836*Tdec_fm_power[5] 
                                 + 629434830.225675*Tdec_fm_power[6] 
                                 - 470493661.096657*Tdec_fm_power[7] 
                                 + 230936465.421*Tdec_fm_power[8] 
                                 - 67175218.4629078*Tdec_fm_power[9] 
                                 + 8789472.32652964*Tdec_fm_power[10]);

       bulkvisCoefficients[1] = (  1.18171174036192 
                                 - 17.6740645873717*Tdec_fm_power[1]
                                 + 136.298469057177*Tdec_fm_power[2] 
                                 - 635.999435106846*Tdec_fm_power[3] 
                                 + 1918.77100633321*Tdec_fm_power[4] 
                                 - 3836.32258307711*Tdec_fm_power[5] 
                                 + 5136.35746882372*Tdec_fm_power[6] 
                                 - 4566.22991441914*Tdec_fm_power[7] 
                                 + 2593.45375240886*Tdec_fm_power[8] 
                                 - 853.908199724349*Tdec_fm_power[9]
                                 + 124.260460450113*Tdec_fm_power[10]);
   } else if (bulk_deltaf_kind_ == 2) {
       // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV
       // do not trust it beyond
       bulkvisCoefficients[0] = (  21091365.1182649 
                                 - 290482229.281782*Tdec_fm_power[1]
                                 + 1800423055.01882*Tdec_fm_power[2]
                                 - 6608608560.99887*Tdec_fm_power[3]
                                 + 15900800422.7138*Tdec_fm_power[4]
                                 - 26194517161.8205*Tdec_fm_power[5]
                                 + 29912485360.2916*Tdec_fm_power[6]
                                 - 23375101221.2855*Tdec_fm_power[7]
                                 + 11960898238.0134*Tdec_fm_power[8]
                                 - 3618358144.18576*Tdec_fm_power[9]
                                 + 491369134.205902*Tdec_fm_power[10]);

       bulkvisCoefficients[1] = (  4007863.29316896 
                                 - 55199395.3534188*Tdec_fm_power[1]
                                 + 342115196.396492*Tdec_fm_power[2]
                                 - 1255681487.77798*Tdec_fm_power[3]
                                 + 3021026280.08401*Tdec_fm_power[4]
                                 - 4976331606.85766*Tdec_fm_power[5]
                                 + 5682163732.74188*Tdec_fm_power[6]
                                 - 4439937810.57449*Tdec_fm_power[7]
                                 + 2271692965.05568*Tdec_fm_power[8]
                                 - 687164038.128814*Tdec_fm_power[9]
                                 + 93308348.3137008*Tdec_fm_power[10]);
   } else if (bulk_deltaf_kind_ == 3) {
       bulkvisCoefficients[0] = (  160421664.93603
                                 - 2212807124.97991*Tdec_fm_power[1]
                                 + 13707913981.1425*Tdec_fm_power[2]
                                 - 50204536518.1767*Tdec_fm_power[3]
                                 + 120354649094.362*Tdec_fm_power[4]
                                 - 197298426823.223*Tdec_fm_power[5]
                                 + 223953760788.288*Tdec_fm_power[6]
                                 - 173790947240.829*Tdec_fm_power[7]
                                 + 88231322888.0423*Tdec_fm_power[8]
                                 - 26461154892.6963*Tdec_fm_power[9]
                                 + 3559805050.19592*Tdec_fm_power[10]);
       bulkvisCoefficients[1] = (  33369186.2536556 
                                 - 460293490.420478*Tdec_fm_power[1]
                                 + 2851449676.09981*Tdec_fm_power[2]
                                 - 10443297927.601*Tdec_fm_power[3]
                                 + 25035517099.7809*Tdec_fm_power[4]
                                 - 41040777943.4963*Tdec_fm_power[5]
                                 + 46585225878.8723*Tdec_fm_power[6]
                                 - 36150531001.3718*Tdec_fm_power[7]
                                 + 18353035766.9323*Tdec_fm_power[8]
                                 - 5504165325.05431*Tdec_fm_power[9]
                                 + 740468257.784873*Tdec_fm_power[10]);
   } else if (bulk_deltaf_kind_ == 4) {
       bulkvisCoefficients[0] = (  1167272041.90731 
                                 - 16378866444.6842*Tdec_fm_power[1]
                                 + 103037615761.617*Tdec_fm_power[2]
                                 - 382670727905.111*Tdec_fm_power[3]
                                 + 929111866739.436*Tdec_fm_power[4]
                                 - 1540948583116.54*Tdec_fm_power[5]
                                 + 1767975890298.1*Tdec_fm_power[6]
                                 - 1385606389545*Tdec_fm_power[7]
                                 + 709922576963.213*Tdec_fm_power[8]
                                 - 214726945096.326*Tdec_fm_power[9]
                                 + 29116298091.9219*Tdec_fm_power[10]);
       bulkvisCoefficients[1] = (  5103633637.7213 
                                 - 71612903872.8163*Tdec_fm_power[1]
                                 + 450509014334.964*Tdec_fm_power[2]
                                 - 1673143669281.46*Tdec_fm_power[3]
                                 + 4062340452589.89*Tdec_fm_power[4]
                                 - 6737468792456.4*Tdec_fm_power[5]
                                 + 7730102407679.65*Tdec_fm_power[6]
                                 - 6058276038129.83*Tdec_fm_power[7]
                                 + 3103990764357.81*Tdec_fm_power[8]
                                 - 938850005883.612*Tdec_fm_power[9]
                                 + 127305171097.249*Tdec_fm_power[10]);
   }
   return;
}


void FSSW::load_bulk_deltaf_14mom_table(string filepath) {
    std::string folder_name = "/urqmd";
    if (afterburner_type_ == AfterburnerType::SMASH) {
        folder_name = "/smash_box";
    }
    std::stringstream file_c0;
    file_c0 << filepath << folder_name << "/c0.dat";
    ifstream c0_table(file_c0.str().c_str());
    if (!c0_table) {
        messager_ << "Can not found file: " << file_c0.str();
        messager_.flush("error");
        exit(1);
    }
    std::stringstream file_c1;
    file_c1 << filepath << folder_name << "/c1.dat";
    ifstream c1_table(file_c1.str().c_str());
    if (!c1_table) {
        messager_ << "Can not found file: " << file_c1.str();
        messager_.flush("error");
        exit(1);
    }
    std::stringstream file_c2;
    file_c2 << filepath << folder_name << "/c2.dat";
    ifstream c2_table(file_c2.str().c_str());
    if (!c2_table) {
        messager_ << "Can not found file: " << file_c2.str();
        messager_.flush("error");
        exit(1);
    }

    // get table size
    c0_table >> deltaf_bulk_coeff_14mom_table_length_T_
             >> deltaf_bulk_coeff_14mom_table_length_mu_;

    string dummy;
    // read header
    std::getline(c0_table, dummy);

    for (int i = 0; i < 3; i++) {
        std::getline(c1_table, dummy);
        std::getline(c2_table, dummy);
    }

    // allocate 2D tables
    deltaf_bulk_coeff_14mom_c0_tb_ = (
            new double* [deltaf_bulk_coeff_14mom_table_length_T_]);
    deltaf_bulk_coeff_14mom_c1_tb_ = (
            new double* [deltaf_bulk_coeff_14mom_table_length_T_]);
    deltaf_bulk_coeff_14mom_c2_tb_ = (
            new double* [deltaf_bulk_coeff_14mom_table_length_T_]);
    for (int i = 0; i < deltaf_bulk_coeff_14mom_table_length_T_; i++) {
        deltaf_bulk_coeff_14mom_c0_tb_[i] = (
                new double [deltaf_bulk_coeff_14mom_table_length_mu_]);
        deltaf_bulk_coeff_14mom_c1_tb_[i] = (
                new double [deltaf_bulk_coeff_14mom_table_length_mu_]);
        deltaf_bulk_coeff_14mom_c2_tb_[i] = (
                new double [deltaf_bulk_coeff_14mom_table_length_mu_]);
    }

    // load 2D table
    double temp1, temp2;
    for (int j = 0; j < deltaf_bulk_coeff_14mom_table_length_mu_; j++) {
        for (int i = 0; i < deltaf_bulk_coeff_14mom_table_length_T_; i++) {
            c0_table >> temp1 >> temp2 >> deltaf_bulk_coeff_14mom_c0_tb_[i][j];
            c1_table >> temp1 >> temp2 >> deltaf_bulk_coeff_14mom_c1_tb_[i][j];
            c2_table >> temp1 >> temp2 >> deltaf_bulk_coeff_14mom_c2_tb_[i][j];
            if (j == 0) {
                if (i == 0) {
                    deltaf_bulk_coeff_14mom_table_T0_ = temp1;
                    deltaf_bulk_coeff_14mom_table_mu0_ = temp2;
                } else if (i == 1) {
                    deltaf_bulk_coeff_14mom_table_dT_ = (
                            temp1 - deltaf_bulk_coeff_14mom_table_T0_);
                }
            } else if (j == 1) {
                if (i == 0) {
                    deltaf_bulk_coeff_14mom_table_dmu_ = (
                            temp2 - deltaf_bulk_coeff_14mom_table_mu0_);
                }
            }
        }
    }
    c0_table.close();
    c1_table.close();
    c2_table.close();
}


void FSSW::load_CE_deltaf_NEOSBQS_table(string filepath) {
    std::string folder_name = "/urqmd";
    if (afterburner_type_ == AfterburnerType::SMASH) {
        folder_name = "/smash";
    }
    std::stringstream filename;
    filename << filepath << folder_name << "/NEoSBQS_CE_deltafCoeff.dat";
    ifstream NEoS_table(filename.str().c_str());
    if (!NEoS_table) {
        messager_ << "Can not found file: " << filename.str();
        messager_.flush("error");
        exit(1);
    }

    // define table size
    deltaf_coeff_NEOSBQS_table_length_e_ = 200;
    deltaf_coeff_NEOSBQS_table_length_nB_ = 200;

    string dummy;
    // read header
    std::getline(NEoS_table, dummy);

    // load 2D table
    for (int i = 0; i < deltaf_coeff_NEOSBQS_table_length_e_; i++) {
        for (int j = 0; j < deltaf_coeff_NEOSBQS_table_length_nB_; j++) {
            std::vector<double> temp(5, 0);
            for (int k = 0; k < 5; k++) {
                NEoS_table >> temp[k];
            }
            deltaf_coeff_CE_NEOSBQS_tb_.push_back(temp);
        }
    }
    deltaf_coeff_NEOSBQS_table_de_ = (  deltaf_coeff_CE_NEOSBQS_tb_[200][0]
                                      - deltaf_coeff_CE_NEOSBQS_tb_[0][0]);
    NEoS_table.close();
}


void FSSW::load_22mom_deltaf_NEOSBQS_table(string filepath) {
    std::string folder_name = "/urqmd";
    if (afterburner_type_ == AfterburnerType::SMASH) {
        folder_name = "/smash";
    }
    std::stringstream filename;
    filename << filepath << folder_name << "/NEoSBQS_22mom_deltafCoeff.dat";
    ifstream NEoS_table(filename.str().c_str());
    if (!NEoS_table) {
        messager_ << "Can not found file: " << filename.str();
        messager_.flush("error");
        exit(1);
    }

    // define table size
    deltaf_coeff_NEOSBQS_table_length_e_ = 200;
    deltaf_coeff_NEOSBQS_table_length_nB_ = 200;

    string dummy;
    // read header
    std::getline(NEoS_table, dummy);

    // load 2D table
    for (int i = 0; i < deltaf_coeff_NEOSBQS_table_length_e_; i++) {
        for (int j = 0; j < deltaf_coeff_NEOSBQS_table_length_nB_; j++) {
            std::vector<double> temp(8, 0);
            for (int k = 0; k < 8; k++) {
                NEoS_table >> temp[k];
            }
            deltaf_coeff_22mom_NEOSBQS_tb_.push_back(temp);
        }
    }
    deltaf_coeff_NEOSBQS_table_de_ = (  deltaf_coeff_22mom_NEOSBQS_tb_[200][0]
                                      - deltaf_coeff_22mom_NEOSBQS_tb_[0][0]);
    NEoS_table.close();
}


void FSSW::getbulkvisCoefficients(const double Tdec, const double mu_B,
                                  std::vector<double> &bulkvisCoefficients) {
    bulkvisCoefficients.resize(3, 0.);
    int idx_T = static_cast<int>((Tdec - deltaf_bulk_coeff_14mom_table_T0_)
                                 /deltaf_bulk_coeff_14mom_table_dT_);
    int idx_mu = static_cast<int>((mu_B - deltaf_bulk_coeff_14mom_table_mu0_)
                                  /deltaf_bulk_coeff_14mom_table_dmu_);
    double x_fraction = ((Tdec - deltaf_bulk_coeff_14mom_table_T0_)
                         /deltaf_bulk_coeff_14mom_table_dT_ - idx_T);
    double y_fraction = ((mu_B - deltaf_bulk_coeff_14mom_table_mu0_)
                         /deltaf_bulk_coeff_14mom_table_dmu_ - idx_mu);

    // avoid overflow and underflow
    int idx_T1 = std::max(
        0, std::min(deltaf_bulk_coeff_14mom_table_length_T_ - 1, idx_T));
    int idx_mu1 = std::max(
        0, std::min(deltaf_bulk_coeff_14mom_table_length_mu_ - 1, idx_mu));
    int idx_T2 = std::max(
        0, std::min(deltaf_bulk_coeff_14mom_table_length_T_ - 1, idx_T + 1));
    int idx_mu2 = std::max(
        0, std::min(deltaf_bulk_coeff_14mom_table_length_mu_ - 1, idx_mu + 1));

    double f1 = deltaf_bulk_coeff_14mom_c0_tb_[idx_T1][idx_mu1];
    double f2 = deltaf_bulk_coeff_14mom_c0_tb_[idx_T1][idx_mu2];
    double f3 = deltaf_bulk_coeff_14mom_c0_tb_[idx_T2][idx_mu2];
    double f4 = deltaf_bulk_coeff_14mom_c0_tb_[idx_T2][idx_mu1];
    double c0_T4 = (  f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));
    f1 = deltaf_bulk_coeff_14mom_c1_tb_[idx_T1][idx_mu1];
    f2 = deltaf_bulk_coeff_14mom_c1_tb_[idx_T1][idx_mu2];
    f3 = deltaf_bulk_coeff_14mom_c1_tb_[idx_T2][idx_mu2];
    f4 = deltaf_bulk_coeff_14mom_c1_tb_[idx_T2][idx_mu1];
    double c1_T3 = (  f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));
    f1 = deltaf_bulk_coeff_14mom_c2_tb_[idx_T1][idx_mu1];
    f2 = deltaf_bulk_coeff_14mom_c2_tb_[idx_T1][idx_mu2];
    f3 = deltaf_bulk_coeff_14mom_c2_tb_[idx_T2][idx_mu2];
    f4 = deltaf_bulk_coeff_14mom_c2_tb_[idx_T2][idx_mu1];
    double c2_T4 = (  f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));
    double T3 = Tdec*Tdec*Tdec;
    double T4 = T3*Tdec;
    bulkvisCoefficients[0] = (c0_T4 - c2_T4)/T4;
    bulkvisCoefficients[1] = c1_T3/T3;
    bulkvisCoefficients[2] = (4.*c2_T4 - c0_T4)/T4;
}


void FSSW::getCENEOSBQSCoefficients(const double Edec, const double nB,
                                    std::vector<double> &visCoefficients) {
    visCoefficients.resize(3, 0.);
    int idx_e = static_cast<int>((Edec - deltaf_coeff_CE_NEOSBQS_tb_[0][0])
                                 /deltaf_coeff_NEOSBQS_table_de_);
    // avoid overflow and underflow
    idx_e = std::max(
        0, std::min(deltaf_coeff_NEOSBQS_table_length_e_ - 2, idx_e));
    int Ne1 = idx_e*deltaf_coeff_NEOSBQS_table_length_e_;
    int Ne2 = (idx_e+1)*deltaf_coeff_NEOSBQS_table_length_e_;
    double e_frac = ((Edec - deltaf_coeff_CE_NEOSBQS_tb_[Ne1][0])
                     /deltaf_coeff_NEOSBQS_table_de_);

    double dnB1 = deltaf_coeff_CE_NEOSBQS_tb_[Ne1+1][1];
    double dnB2 = deltaf_coeff_CE_NEOSBQS_tb_[Ne2+1][1];
    int idx_nB1 = static_cast<int>((nB)/dnB1);
    int idx_nB2 = static_cast<int>((nB)/dnB2);
    // avoid overflow and underflow
    idx_nB1 = std::min(deltaf_coeff_NEOSBQS_table_length_nB_ - 2, idx_nB1);
    idx_nB2 = std::min(deltaf_coeff_NEOSBQS_table_length_nB_ - 2, idx_nB2);
    double nB_frac1 = std::min(
            1., (nB - deltaf_coeff_CE_NEOSBQS_tb_[Ne1+idx_nB1][1])/dnB1);
    double nB_frac2 = std::min(
            1., (nB - deltaf_coeff_CE_NEOSBQS_tb_[Ne2+idx_nB2][1])/dnB2);

    // check
    //cout << "ed = " << Edec
    //     << ", e0 = " << deltaf_coeff_CE_NEOSBQS_tb_[Ne1][0]
    //     << ", e1 = " << deltaf_coeff_CE_NEOSBQS_tb_[Ne2][0]
    //     << ", efrac = " << e_frac
    //     << ", nB = " << nB
    //     << ", nB1_0 = " << deltaf_coeff_CE_NEOSBQS_tb_[Ne1+idx_nB1][1]
    //     << ", nB1_1 = " << deltaf_coeff_CE_NEOSBQS_tb_[Ne1+idx_nB1+1][1]
    //     << ", nB_frac1 = " << nB_frac1
    //     << ", nB2_0 = " << deltaf_coeff_CE_NEOSBQS_tb_[Ne2+idx_nB2][1]
    //     << ", nB2_1 = " << deltaf_coeff_CE_NEOSBQS_tb_[Ne2+idx_nB2+1][1]
    //     << ", nB_frac2 = " << nB_frac2
    //     << endl;

    std::vector<double> interpCoeff(3, 0);
    for (int i = 2; i < 5; i++) {
        // c2_hat, zeta_hat, eta_hat
        double temp1 = (
              deltaf_coeff_CE_NEOSBQS_tb_[Ne1+idx_nB1  ][i]*(1 - nB_frac1)
            + deltaf_coeff_CE_NEOSBQS_tb_[Ne1+idx_nB1+1][i]*nB_frac1);
        double temp2 = (
              deltaf_coeff_CE_NEOSBQS_tb_[Ne2+idx_nB2  ][i]*(1 - nB_frac2)
            + deltaf_coeff_CE_NEOSBQS_tb_[Ne2+idx_nB2+1][i]*nB_frac2);

        interpCoeff[i-2] = temp1*(1. - e_frac) + temp2*e_frac;
    }
    visCoefficients[0] = 1./interpCoeff[1];
    visCoefficients[1] = 1./3. - interpCoeff[0];
    visCoefficients[2] = interpCoeff[2];
}


void FSSW::get22momNEOSBQSCoefficients(const double Edec, const double nB,
                                       std::vector<double> &visCoefficients) {
    visCoefficients.resize(6, 0.);
    int idx_e = static_cast<int>((Edec - deltaf_coeff_22mom_NEOSBQS_tb_[0][0])
                                 /deltaf_coeff_NEOSBQS_table_de_);
    // avoid overflow and underflow
    idx_e = std::max(
        0, std::min(deltaf_coeff_NEOSBQS_table_length_e_ - 2, idx_e));
    int Ne1 = idx_e*deltaf_coeff_NEOSBQS_table_length_e_;
    int Ne2 = (idx_e+1)*deltaf_coeff_NEOSBQS_table_length_e_;
    double e_frac = ((Edec - deltaf_coeff_22mom_NEOSBQS_tb_[Ne1][0])
                     /deltaf_coeff_NEOSBQS_table_de_);

    double dnB1 = deltaf_coeff_22mom_NEOSBQS_tb_[Ne1+1][1];
    double dnB2 = deltaf_coeff_22mom_NEOSBQS_tb_[Ne2+1][1];
    int idx_nB1 = static_cast<int>((nB)/dnB1);
    int idx_nB2 = static_cast<int>((nB)/dnB2);
    // avoid overflow and underflow
    idx_nB1 = std::min(deltaf_coeff_NEOSBQS_table_length_nB_ - 2, idx_nB1);
    idx_nB2 = std::min(deltaf_coeff_NEOSBQS_table_length_nB_ - 2, idx_nB2);
    double nB_frac1 = std::min(
            1., (nB - deltaf_coeff_22mom_NEOSBQS_tb_[Ne1+idx_nB1][1])/dnB1);
    double nB_frac2 = std::min(
            1., (nB - deltaf_coeff_22mom_NEOSBQS_tb_[Ne2+idx_nB2][1])/dnB2);

    // check
    //cout << "ed = " << Edec
    //     << ", e0 = " << deltaf_coeff_22mom_NEOSBQS_tb_[Ne1][0]
    //     << ", e1 = " << deltaf_coeff_22mom_NEOSBQS_tb_[Ne2][0]
    //     << ", efrac = " << e_frac
    //     << ", nB = " << nB
    //     << ", nB1_0 = " << deltaf_coeff_22mom_NEOSBQS_tb_[Ne1+idx_nB1][1]
    //     << ", nB1_1 = " << deltaf_coeff_22mom_NEOSBQS_tb_[Ne1+idx_nB1+1][1]
    //     << ", nB_frac1 = " << nB_frac1
    //     << ", nB2_0 = " << deltaf_coeff_22mom_NEOSBQS_tb_[Ne2+idx_nB2][1]
    //     << ", nB2_1 = " << deltaf_coeff_22mom_NEOSBQS_tb_[Ne2+idx_nB2+1][1]
    //     << ", nB_frac2 = " << nB_frac2
    //     << endl;

    std::vector<double> interpCoeff(6, 0);
    for (int i = 2; i < 8; i++) {
        double temp1 = (
              deltaf_coeff_22mom_NEOSBQS_tb_[Ne1+idx_nB1  ][i]*(1 - nB_frac1)
            + deltaf_coeff_22mom_NEOSBQS_tb_[Ne1+idx_nB1+1][i]*nB_frac1);
        double temp2 = (
              deltaf_coeff_22mom_NEOSBQS_tb_[Ne2+idx_nB2  ][i]*(1 - nB_frac2)
            + deltaf_coeff_22mom_NEOSBQS_tb_[Ne2+idx_nB2+1][i]*nB_frac2);

        interpCoeff[i-2] = temp1*(1. - e_frac) + temp2*e_frac;
    }

    for (int i = 0; i < 6; i++)
        visCoefficients[i] = interpCoeff[i];
}


void FSSW::load_deltaf_qmu_coeff_table(string filename) {
    ifstream table(filename.c_str());
    deltaf_qmu_coeff_table_length_T = 150;
    deltaf_qmu_coeff_table_length_mu = 100;
    delta_qmu_coeff_table_T0 = 0.05;
    delta_qmu_coeff_table_mu0 = 0.0;
    delta_qmu_coeff_table_dT = 0.001;
    delta_qmu_coeff_table_dmu = 0.007892;

    deltaf_qmu_coeff_tb = new double* [deltaf_qmu_coeff_table_length_T];
    for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++) {
        deltaf_qmu_coeff_tb[i] = new double [deltaf_qmu_coeff_table_length_mu];
    }

    // load 2D table
    double dummy;
    for (int j = 0; j < deltaf_qmu_coeff_table_length_mu; j++) {
        for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++) {
            table >> dummy >> dummy >> deltaf_qmu_coeff_tb[i][j];
        }
    }
    table.close();
}


double FSSW::get_deltaf_qmu_coeff(double T, double muB) {
    int idx_T = static_cast<int>((T - delta_qmu_coeff_table_T0)
                                 /delta_qmu_coeff_table_dT);
    int idx_mu = static_cast<int>((muB - delta_qmu_coeff_table_mu0)
                                  /delta_qmu_coeff_table_dmu);
    double x_fraction = ((T - delta_qmu_coeff_table_T0)
                         /delta_qmu_coeff_table_dT - idx_T);
    double y_fraction = ((muB - delta_qmu_coeff_table_mu0)
                         /delta_qmu_coeff_table_dmu - idx_mu); 

    //avoid overflow
    if (idx_mu > deltaf_qmu_coeff_table_length_mu - 2) {
        return(1e30);
    }
    if (idx_T > deltaf_qmu_coeff_table_length_T - 2) {
        return(1e30);
    }

    // avoid underflow
    if (idx_mu < 0) {
        return(1e30);
    }
    if (idx_T < 0) {
        return(1e30);
    }

    double f1 = deltaf_qmu_coeff_tb[idx_T][idx_mu];
    double f2 = deltaf_qmu_coeff_tb[idx_T][idx_mu+1];
    double f3 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu+1];
    double f4 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu];
    double coeff = (  f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction 
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));  // 1/fm^3
    return(coeff);
}


void FSSW::initialize_special_function_arrays() {
    if (echoLevel_ > 0) {
        messager_.info("Initializing special function arrays ... ");
    }
    sf_expint_truncate_order = 10;
    sf_x_min = 0.5;
    sf_x_max = 400;
    sf_dx = 0.05;
    sf_tb_length = static_cast<int>((sf_x_max - sf_x_min)/sf_dx) + 1;
    sf_bessel_Kn = new double* [sf_tb_length];
    if (INCLUDE_DIFFUSION_DELTAF == 1) {
        sf_expint_En = new double* [sf_tb_length];
    }
    for (int i = 0; i < sf_tb_length; i++) {
        sf_bessel_Kn[i] = new double [3];

        double sf_x = sf_x_min + i*sf_dx;

        if (INCLUDE_BULK_DELTAF == 1) {
            sf_bessel_Kn[i][0] = gsl_sf_bessel_K1(sf_x);     // store K_1
            sf_bessel_Kn[i][2] = gsl_sf_bessel_Kn(3, sf_x);     // store K_3
        } else {
            sf_bessel_Kn[i][0] = 0.0;
            sf_bessel_Kn[i][2] = 0.0;
        }

        sf_bessel_Kn[i][1] = gsl_sf_bessel_Kn(2, sf_x);  // store K_2

        if (INCLUDE_DIFFUSION_DELTAF == 1) {
            sf_expint_En[i] = new double [sf_expint_truncate_order-1];
            sf_expint_En[i][0] = gsl_sf_expint_E2(sf_x);     // store E_2
            for (int k = 1; k < sf_expint_truncate_order-1; k++) {
                sf_expint_En[i][k] = gsl_sf_expint_En(2*k+2, sf_x);
            }
        }
    }
}


double FSSW::get_special_function_K3(double arg) {
    double results;
    if (arg < sf_x_min || arg > sf_x_max-sf_dx) {
        if (AMOUNT_OF_OUTPUT > 5) {
            cout << "FSSW::get_special_function_K3: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results = gsl_sf_bessel_Kn(3, arg);
    } else {
        int idx = static_cast<int>((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        results = ((1. - fraction)*sf_bessel_Kn[idx][2]
                    + fraction*sf_bessel_Kn[idx+1][2]);
    }
    return(results);
}


double FSSW::get_special_function_K2(double arg) {
    double results;
    if (arg < sf_x_min || arg > sf_x_max-sf_dx) {
        if (AMOUNT_OF_OUTPUT > 5) {
            cout << "FSSW::get_special_function_K2: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results = gsl_sf_bessel_Kn(2, arg);
    } else {
        int idx = static_cast<int>((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        results = ((1. - fraction)*sf_bessel_Kn[idx][1]
                    + fraction*sf_bessel_Kn[idx+1][1]);
    }
    return(results);
}


double FSSW::get_special_function_K1(double arg) {
    double results;
    if (arg < sf_x_min || arg > sf_x_max-sf_dx) {
        if (AMOUNT_OF_OUTPUT > 5) {
            cout << "FSSW::get_special_function_K1: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results = gsl_sf_bessel_K1(arg);
    } else {
        int idx = static_cast<int>((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        results = ((1. - fraction)*sf_bessel_Kn[idx][0]
                    + fraction*sf_bessel_Kn[idx+1][0]);
    }
    return(results);
}


void FSSW::get_special_function_En(
                                double arg, std::vector<double> &results) {
    if (arg < sf_x_min || arg > sf_x_max-sf_dx) {
        if (AMOUNT_OF_OUTPUT > 5) {
            cout << "FSSW::get_special_function_En: "
                 << "out of the table bound!" << endl;
            cout << "sf_x_min = " << sf_x_min << ", sf_x_max = " << sf_x_max
                 << ", arg = " << arg << endl;
        }
        results[0] = gsl_sf_expint_E2(arg);
        for (int i = 1; i < sf_expint_truncate_order-1; i++) {
            results[i] = gsl_sf_expint_En(2*i+2, arg);
        }
    } else {
        int idx = static_cast<int>((arg - sf_x_min)/sf_dx);
        double fraction = (arg - sf_x_min - idx*sf_dx)/sf_dx;
        for (int i = 0; i < sf_expint_truncate_order-1; i++) {
            results[i] = ((1. - fraction)*sf_expint_En[idx][i] 
                          + fraction*sf_expint_En[idx+1][i]);
        }
    }
}


void FSSW::check_samples_in_memory() {
    cout << "check samples in memory ... " << endl;
    cout << "number of events: " << Hadron_list->size() << endl;
    for (unsigned int i = 0; i < Hadron_list->size(); i++) {
        cout << "event " << i << ": number of particles = "
             << (*Hadron_list)[i]->size() << endl;
        for (unsigned int j = 0; j < (*Hadron_list)[i]->size(); j++) {
            cout << "particle " << j
                 << " pid = " << (*(*Hadron_list)[i])[j].pid
                 << " mass = " << (*(*Hadron_list)[i])[j].mass
                 << endl;
        }
    }
}


void FSSW::perform_resonance_feed_down(
                    vector< vector<iSS_Hadron>* >* input_particle_list) {
    if (echoLevel_ > 0) {
        cout << "perform resonance decays... " << endl;
    }
    // loop over events
    unsigned int nev = input_particle_list->size();
    for (unsigned int ievent = 0; ievent < nev; ievent++) {
        // create a temporary particle list
        vector<iSS_Hadron> temp_list;
        // copy all particles into the temp list
        unsigned int Npart = (*input_particle_list)[ievent]->size();
        for (unsigned int ipart = 0; ipart < Npart; ipart++) {
            temp_list.push_back((*(*input_particle_list)[ievent])[ipart]);
        }
        (*input_particle_list)[ievent]->clear();
        // perform resonance decays
        for (unsigned int ipart = 0; ipart < temp_list.size(); ipart++) {
            vector<iSS_Hadron> *daughter_list = new vector<iSS_Hadron>;
            decayer_ptr_->perform_decays(&temp_list[ipart], daughter_list);
            for (unsigned int idaughter = 0; idaughter < daughter_list->size();
                    idaughter++) {
                if (decayer_ptr_->check_particle_stable(
                                        &(*daughter_list)[idaughter]) == 1) {
                    (*input_particle_list)[ievent]->push_back(
                                            (*daughter_list)[idaughter]);
                } else {
                    temp_list.push_back((*daughter_list)[idaughter]);
                }
            }
            daughter_list->clear();
            delete daughter_list;
        }
        temp_list.clear();
    }
}


void FSSW::addSpectatorsToHadronList() {
    cout << "Add spectators to the hadron list... " << endl;
    // loop over events
    unsigned int nev = Hadron_list->size();
    for (unsigned int ievent = 0; ievent < nev; ievent++) {
        auto spectatorList = spectators_ptr_->getSpectatorList();
        for (auto ispec: spectatorList) {
            (*Hadron_list)[ievent]->push_back(ispec);
        }
    }
}


double FSSW::get_deltaf_bulk(
        const double mass, const double pdotu, const double bulkPi,
        const double Tdec, const int sign, const int baryon,
        const int strange, const int charge,
        const double f0, const std::vector<double> bulkvisCoefficients) {
    if (INCLUDE_BULK_DELTAF== 0) return(0.0);
    double delta_f_bulk = 0.0;
    if (bulk_deltaf_kind_ == 0) {
        delta_f_bulk = (-(1. - sign*f0)*bulkPi
                        *(  bulkvisCoefficients[0]*mass*mass
                          + bulkvisCoefficients[1]*pdotu
                          + bulkvisCoefficients[2]*pdotu*pdotu));
    } else if (bulk_deltaf_kind_ == 1) {
        double E_over_T = pdotu/Tdec;
        double mass_over_T = mass/Tdec;
        delta_f_bulk = (- 1.0*(1. - sign*f0)*bulkvisCoefficients[0]
                        *(mass_over_T*mass_over_T/(3.*E_over_T)
                          - bulkvisCoefficients[1]*E_over_T)*bulkPi);
    } else if (bulk_deltaf_kind_ == 2) {
        double E_over_T = pdotu/Tdec;
        delta_f_bulk = (- 1.*(1. - sign*f0)*bulkPi
                *(- bulkvisCoefficients[0] + bulkvisCoefficients[1]*E_over_T));
    } else if (bulk_deltaf_kind_ == 3) {
        double E_over_T = pdotu/Tdec;
        delta_f_bulk = (- 1.*(1. - sign*f0)*bulkPi/sqrt(E_over_T)
                *(- bulkvisCoefficients[0] + bulkvisCoefficients[1]*E_over_T));
    } else if (bulk_deltaf_kind_ == 4) {
        double E_over_T = pdotu/Tdec;
        delta_f_bulk = (- 1.*(1. - sign*f0)*bulkPi
                *(bulkvisCoefficients[0] - bulkvisCoefficients[1]/E_over_T));
    } else if (bulk_deltaf_kind_ == 11) {
        // OSU 14-moment
        delta_f_bulk = ((1. - sign*f0)*bulkPi*(
                          bulkvisCoefficients[0]*mass*mass
                        + bulkvisCoefficients[1]*baryon*pdotu
                        + bulkvisCoefficients[2]*pdotu*pdotu));
    } else if (bulk_deltaf_kind_ == 20) {
        // WSU 22-moment NEoS BQS
        delta_f_bulk = ((1. - sign*f0)*bulkPi*(
                          mass*mass*bulkvisCoefficients[2]
                        + pdotu*(  baryon*bulkvisCoefficients[3]
                                 + strange*bulkvisCoefficients[4]
                                 + charge*bulkvisCoefficients[5])
                        + pdotu*pdotu*(  bulkvisCoefficients[1]
                                       - bulkvisCoefficients[2])));
    } else if (bulk_deltaf_kind_ == 21) {
        // WSU Chapman-Enskog NEoS BQS
        double E_over_T = pdotu/Tdec;
        double mass_over_T = mass/Tdec;
        delta_f_bulk = (- 1.0*(1. - sign*f0)*bulkvisCoefficients[0]
                        *(mass_over_T*mass_over_T/(3.*E_over_T)
                          - bulkvisCoefficients[1]*E_over_T)*bulkPi);
    }
    return(delta_f_bulk);
}


int FSSW::sample_momemtum_from_a_fluid_cell(
        const double mass, const int sign,
        const int baryon, const int strange, const int charge,
        const FO_surf_LRF *surf,
        const std::vector<double> bulkvisCoefficients,
        const double deltaf_qmu_coeff,
        double &pT, double &phi, double &y_minus_eta_s
        ) {
    const double Tdec = surf->Tdec;
    const double mu = std::min(mass, static_cast<double>(baryon*surf->muB
                                                         + strange*surf->muS
                                                         + charge*surf->muQ));
    const double deltaf_prefactor = (
                            1.0/(2.0*Tdec*Tdec*(surf->Edec + surf->Pdec)));
    const double prefactor_qmu = surf->Bn/(surf->Edec + surf->Pdec);
    const double dsigam_fac = (std::abs(surf->da_mu_LRF[0])
                            + sqrt(  surf->da_mu_LRF[1]*surf->da_mu_LRF[1]
                                   + surf->da_mu_LRF[2]*surf->da_mu_LRF[2]
                                   + surf->da_mu_LRF[3]*surf->da_mu_LRF[3]));
    int tries = 1;
    const int maximum_impatience = 5000;
    while (tries < maximum_impatience) {
        // refer to calculate_dNArrays function to see how
        // the rate is calculated
        // Basically it is "just Cooper-Frye"

        double p_mag = momentum_sampler_ptr_->Sample_a_momentum(mass, Tdec, mu,
                                                                sign);
        double phi_candidate = 2*M_PI*ran_gen_ptr->rand_uniform();
        double cos_theta = 2.*ran_gen_ptr->rand_uniform() - 1.;
        double sin_theta = sqrt(1. - cos_theta*cos_theta);
        double pT_candidate = p_mag*sin_theta;

        double px = pT_candidate*cos(phi_candidate);
        double py = pT_candidate*sin(phi_candidate);

        double p0 = sqrt(mass*mass + p_mag*p_mag);
        double pz = p_mag*cos_theta;

        double pdsigma = (  p0*surf->da_mu_LRF[0] + px*surf->da_mu_LRF[1]
                          + py*surf->da_mu_LRF[2] + pz*surf->da_mu_LRF[3]);

        double f0 = 1./(exp((p0 - mu)/Tdec) + sign);

        // delta f for shear viscosity
        double delta_f_shear = 0.;
        if (INCLUDE_DELTAF == 1) {
            double Wfactor = (
                px*px*surf->piLRF_xx + 2.*px*py*surf->piLRF_xy
                + 2.*px*pz*surf->piLRF_xz
                + py*py*surf->piLRF_yy + 2.*py*pz*surf->piLRF_yz
                + pz*pz*(- surf->piLRF_xx - surf->piLRF_yy));
            if (NEoS_deltaf_kind_ == 1) {
                delta_f_shear = (
                    (1. - sign*f0)*Wfactor/(2.*bulkvisCoefficients[2])
                    /(p0*Tdec));
            } else if (NEoS_deltaf_kind_ == 0) {
                delta_f_shear = (1. - sign*f0)*Wfactor*bulkvisCoefficients[0];
            } else {
                delta_f_shear = (1. - sign*f0)*Wfactor*deltaf_prefactor;
            }
        }

        // delta f for bulk viscosity
        double delta_f_bulk = 0.;
        if (INCLUDE_BULK_DELTAF == 1) {
            double bulkPi = 0.;
            if (bulk_deltaf_kind_ == 11 || bulk_deltaf_kind_ == 21
                || bulk_deltaf_kind_ == 20) {
                bulkPi = surf->bulkPi;         // GeV/fm^3
            } else if (bulk_deltaf_kind_ == 1) {
                bulkPi = surf->bulkPi/hbarC;   // 1/fm^4
            }
            delta_f_bulk = get_deltaf_bulk(mass, p0, bulkPi, Tdec, sign,
                                           baryon, strange, charge,
                                           f0, bulkvisCoefficients);
        }

        // delta f for diffusion
        double delta_f_qmu = 0.0;
        if (INCLUDE_DIFFUSION_DELTAF == 1) {
            double qmufactor = (- px*surf->qmuLRF_x - py*surf->qmuLRF_y
                                - pz*surf->qmuLRF_z);
            delta_f_qmu = ((1. - sign*f0)*(prefactor_qmu - baryon/p0)
                           *qmufactor/deltaf_qmu_coeff);
        }

        double fact1 = pdsigma/p0/dsigam_fac;
        double fact2 = (1. + delta_f_shear + delta_f_bulk + delta_f_qmu)/2.;
        fact1 = std::max(0., std::min(1., fact1));
        fact2 = std::max(0., std::min(1., fact2));
        double accept_prob = fact1*fact2;

        if (ran_gen_ptr->rand_uniform() < accept_prob) {
            // accept the sample
            // now we need to boost the momentum to the lab frame
            VecD4 pLRF = {p0, px, py, pz};
            VecD4 umu = {surf->u_tz[0], surf->u_tz[1],
                         surf->u_tz[2], surf->u_tz[3]};
            VecD4 pLab = {0., 0., 0., 0.};
            boost_vector_back_to_lab_frame(pLRF, pLab, umu);

            // assigned to the return variables
            // if E and p^z are too close, resample one
            double y = 0.5*log((pLab[0] + pLab[3])/(pLab[0] - pLab[3]));
            if (!std::isnan(y) && !std::isinf(y)) {
                pT = sqrt(pLab[1]*pLab[1] + pLab[2]*pLab[2]);
                phi = atan2(pLab[2], pLab[1]);
                y_minus_eta_s = y - surf->eta;
                return(1);
            }
        }
        tries++;
    }
    return(0);
}


void FSSW::add_one_sampled_particle(
                const int repeated_sampling_idx, const FO_surf_LRF *surf,
                const int particle_monval, const double mass,
                const double pT, const double phi,
                const double y_minus_eta_s, const double eta_s) {
    double rapidity_y = y_minus_eta_s + eta_s;

    const double px = pT*cos(phi);
    const double py = pT*sin(phi);
    const double mT = sqrt(mass*mass + pT*pT);
    const double p_z = mT*sinh(rapidity_y);
    const double E = mT*cosh(rapidity_y);
    const double z = surf->tau*sinh(eta_s);
    const double t = surf->tau*cosh(eta_s);

    iSS_Hadron *temp_hadron = new iSS_Hadron;
    temp_hadron->pid        = particle_monval;
    temp_hadron->mass       = mass;
    temp_hadron->E          = E;
    temp_hadron->px         = px;
    temp_hadron->py         = py;
    temp_hadron->pz         = p_z;
    temp_hadron->t          = t;
    temp_hadron->x          = surf->xpt;
    temp_hadron->y          = surf->ypt;
    temp_hadron->z          = z;
    (*Hadron_list)[repeated_sampling_idx-1]->push_back(*temp_hadron);
}


// this function boost a 4-vector from the fluid local rest frame to
// the lab frame. The fluid velocity is u^\mu.
void FSSW::boost_vector_back_to_lab_frame(
                            VecD4 &p_LRF, VecD4 &p_lab, VecD4 &umu) const {
    double p_dot_u = 0.;
    for (int i = 1; i < 4; i++)
        p_dot_u += p_LRF[i]*umu[i];

    p_lab[0] = p_LRF[0]*umu[0] + p_dot_u;
    for (int i = 1; i < 4; i++)
        p_lab[i] = p_LRF[i] + (p_dot_u/(umu[0] + 1) + p_LRF[0])*umu[i];
}


double FSSW::bilinearInterp(std::vector<std::vector<double>>&mat,
                            int idx_e1, int idx_e2, int idx_nB1, int idx_nB2,
                            double x_fraction, double y_fraction) {
    double f1 = mat[idx_e1][idx_nB1];
    double f2 = mat[idx_e1][idx_nB2];
    double f3 = mat[idx_e2][idx_nB2];
    double f4 = mat[idx_e2][idx_nB1];
    double f = (  f1*(1. - x_fraction)*(1. - y_fraction)
                + f2*(1. - x_fraction)*y_fraction
                + f3*x_fraction*y_fraction
                + f4*x_fraction*(1. - y_fraction));
    return(f);
}


void FSSW::computeAvgTotalEnergyMomentum() {
    std::vector<std::vector<double>> PmuList;
    for (auto const &ev_i: (*Hadron_list)) {
        std::vector<double> totalPmu(4, 0.);
        for (auto &part_i: (*ev_i)) {
            totalPmu[0] += part_i.E;
            totalPmu[1] += part_i.px;
            totalPmu[2] += part_i.py;
            totalPmu[3] += part_i.pz;
        }
        PmuList.push_back(totalPmu);
    }
    int nev = 0;
    std::vector<double> Pmu_avg(4, 0);
    std::vector<double> Pmu_err(4, 0);
    for (auto const &Pmu_i : PmuList) {
        for (int j = 0; j < 4; j++) {
            Pmu_avg[j] += Pmu_i[j];
            Pmu_err[j] += Pmu_i[j]*Pmu_i[j];
        }
        nev++;
    }
    if (echoLevel_ > 0) {
        messager_.info("Averaged total energy and momentum:");
    }
    for (int i = 0; i < 4; i++) {
        Pmu_avg[i] = Pmu_avg[i]/nev;
        Pmu_err[i] = sqrt((Pmu_err[i]/nev - Pmu_avg[i]*Pmu_avg[i])/nev);
        if (echoLevel_ > 0) {
            messager_ << "<P[" << i << "]> = " << Pmu_avg[i] << " +/- "
                      << Pmu_err[i] << " GeV.";
            messager_.flush("info");
        }
    }

}
