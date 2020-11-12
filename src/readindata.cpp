// Copyright @ 2012 Chun Shen and Zhi Qiu
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "data_struct.h"
#include "readindata.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "Table.h"

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using iSS_data::hbarC;
using iSS_data::Vec4;

read_FOdata::read_FOdata(ParameterReader* paraRdr_in, string path,
                         string table_path, string particle_table_path) :
        path_(path), table_path_(table_path),
        particle_table_path_(particle_table_path) {
    paraRdr = paraRdr_in;

    mode = paraRdr->getVal("hydro_mode");
    turn_on_bulk = paraRdr->getVal("turn_on_bulk");
    turn_on_rhob = paraRdr->getVal("turn_on_rhob");
    turn_on_diff = paraRdr->getVal("turn_on_diff");
    surface_in_binary = true;

    if (paraRdr->getVal("quantum_statistics") == 1) {
        quantum_statistics_ = true;
    } else {
        quantum_statistics_ = false;
    }

    if (mode == 1 || mode == 2) {
        // determine read in format in surface.dat from MUSIC simulation
        messager.info("read in hyper-surface from MUSIC simulations ...");
        ostringstream config_file;
        config_file << path_ << "/music_input";
        std::ifstream configuration(config_file.str().c_str());
        if (!configuration.is_open()) {
            messager << "read_FOdata::read_FOdata: "
                     << "can not find configuration file "
                     << config_file.str();
            messager.flush("error");
            exit(1);
        }
        string temp1;
        string temp_name;
        while (!configuration.eof()) {
            getline(configuration, temp1);
            std::stringstream ss(temp1);
            ss >> temp_name;
            if (temp_name == "Include_Bulk_Visc_Yes_1_No_0") {
                ss >> turn_on_bulk;
            }
            if (temp_name == "Include_Rhob_Yes_1_No_0") {
                ss >> turn_on_rhob;
            }
            if (temp_name == "turn_on_baryon_diffusion") {
                ss >> turn_on_diff;
            }
            if (temp_name == "freeze_surface_in_binary") {
                int flag_b;
                ss >> flag_b;
                if (flag_b == 1) {
                    surface_in_binary = true;
                } else {
                    surface_in_binary = false;
                }
            }
        }
        configuration.close();
        if (surface_in_binary)
            messager.info("the hyper-surface surface is in the binary format.");
        if (turn_on_bulk == 1)
            messager.info("the hyper-surface includes bulk viscosity.");
        if (turn_on_rhob == 1)
            messager.info("the hyper-surface includes net baryon density.");
        if (turn_on_diff == 1)
            messager.info("the hyper-surface includes baryon diffusion.");
    }
    n_eta_skip = 0;
    iEOS_MUSIC_ = 0;
    int afterburner_id = paraRdr->getVal("afterburner_type");
    if (afterburner_id == 1) {
        afterburner_type_ = AfterburnerType::UrQMD;
    } else if (afterburner_id == 2) {
        afterburner_type_ = AfterburnerType::SMASH;
    } else {
        afterburner_type_ = AfterburnerType::PDG_Decay;
    }
}

read_FOdata::~read_FOdata() {}

int read_FOdata::get_number_of_freezeout_cells(std::string surfaceFilename) {
    int number_of_cells = 0;
    if (mode == 0) {  // outputs from VISH2+1
        ostringstream decdatfile;
        decdatfile << path_ << "/decdat2.dat";
        Table block_file(decdatfile.str().c_str());
        number_of_cells = block_file.getNumberOfRows();
    } else if (mode == 1) {  // outputs from MUSIC boost-invariant
        ostringstream surface_file;
        surface_file << path_ << "/" << surfaceFilename;
        if (surface_in_binary) {
            number_of_cells = get_number_of_lines_of_binary_surface_file(
                                                        surface_file.str());
            n_eta_skip = 1;
        } else {
            Table block_file(surface_file.str().c_str());

            // determine number of the eta slides that are output
            double eta_target = block_file.get(4, 1);
            int num_temp = 0;
            for (int i = 0; i < block_file.getNumberOfRows(); i++) {
                if (block_file.get(4, i+1) == eta_target)
                    num_temp++;
            }
            number_of_cells = num_temp;
            n_eta_skip = block_file.getNumberOfRows()/number_of_cells;
        }
    } else if (mode == 2) {  // outputs from MUSIC full (3+1)-d
        int number_of_lines = 0;
        ostringstream surface_filename;
        surface_filename << path_ << "/" << surfaceFilename;
        if (surface_in_binary) {
            number_of_cells = get_number_of_lines_of_binary_surface_file(
                                                    surface_filename.str());
        } else {
            std::string temp_line;
            std::ifstream surface_file(surface_filename.str().c_str());
            while (std::getline(surface_file, temp_line)) {
                ++number_of_lines;
            }
            surface_file.close();
            number_of_cells = number_of_lines;
        }
    } else if (mode == 10) {  // outputs from hydro analysis
        ostringstream surface_file;
        surface_file << path_ << "/hyper_surface_2+1d.dat";
        Table block_file(surface_file.str().c_str());
        number_of_cells = block_file.getNumberOfRows();
    }
    return(number_of_cells);
}

int read_FOdata::get_number_of_lines_of_binary_surface_file(string filename) {
    std::ifstream surface_file(filename.c_str(), std::ios::binary);
    int count = 0;
    float temp = 0.;
    while(surface_file) {
        surface_file.read(reinterpret_cast<char*>(&temp), sizeof(float));
        count++;
    }
    int counted_line = count/34;
    surface_file.close();
    return(counted_line);
}


void read_FOdata::read_in_freeze_out_data(std::vector<FO_surf> &surf_ptr,
                                          std::string surface_filename) {
    if (mode == 0)     // VISH2+1 outputs
        read_FOsurfdat_VISH2p1(surf_ptr);
    else if (mode == 1)   // MUSIC boost invariant outputs
        read_FOsurfdat_MUSIC_boost_invariant(surf_ptr, surface_filename);
    else if (mode == 2)   // MUSIC full (3+1)-d outputs
        read_FOsurfdat_MUSIC(surf_ptr, surface_filename);
    else if (mode == 10)   // MUSIC boost invariant outputs
        read_FOsurfdat_hydro_analysis_boost_invariant(surf_ptr);
    regulate_surface_cells(surf_ptr);
}


void read_FOdata::read_in_chemical_potentials(
    std::vector<FO_surf> &surf_ptr, std::vector<particle_info> &particle_ptr) {
    int Nparticle = 0;
    int N_stableparticle = 0;
    Table mu_table;
    if (mode == 0) {      // VISH2+1 output
        std::ifstream particletable(table_path_ + "/EOS_particletable.dat");
        particletable >> N_stableparticle;
        particletable.close();
    } else if (mode == 1 || mode == 2) {   // music output
        // determine the type of the EOS
        ostringstream config_file;
        config_file << path_ << "/music_input";
        std::ifstream configuration(config_file.str().c_str());
        string temp1;
        string temp_name;
        while (!configuration.eof()) {
            getline(configuration, temp1);
            std::stringstream ss(temp1);
            ss >> temp_name;
            if (temp_name == "EOS_to_use") {
                ss >> iEOS_MUSIC_;
                break;
            }
        }
        configuration.close();
        std::ifstream particletable;
        if (iEOS_MUSIC_ == 2) {       // s95p-v1
            N_stableparticle = 0;
        } else if (iEOS_MUSIC_ == 3) {  // s95p-v1-PCE150
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE150/EOS_particletable.dat");
            particletable >> N_stableparticle;
            particletable.close();
        } else if (iEOS_MUSIC_ == 4) {  // s95p-v1-PCE155
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE155/EOS_particletable.dat");
            particletable >> N_stableparticle;
            particletable.close();
        } else if (iEOS_MUSIC_ == 5) {  // s95p-v1-PCE160
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE160/EOS_particletable.dat");
            particletable >> N_stableparticle;
            particletable.close();
        } else if (iEOS_MUSIC_ == 6) {  // s95p-v0-PCE165
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v0-PCE165/EOS_particletable.dat");
            particletable >> N_stableparticle;
            particletable.close();
        } else if (iEOS_MUSIC_ == 7) {        // s95p-v1.2 for UrQMD
            N_stableparticle = 0;
        } else if (iEOS_MUSIC_ == 8) {        // WB
            N_stableparticle = 0;
        } else if (iEOS_MUSIC_ == 9) {        // hotQCD + HRG(UrQMD)
            N_stableparticle = 0;
        } else if (iEOS_MUSIC_ == 91) {       // hotQCD + HRG(SMASH)
            N_stableparticle = 0;
        } else if (iEOS_MUSIC_ >= 10 && iEOS_MUSIC_ <=14) {   // NEoS
            N_stableparticle = 0;
        } else if (iEOS_MUSIC_ == 17) {       // BEST
            N_stableparticle = 0;
        } else {
            messager << "invalid iEOS_MUSIC_: " << iEOS_MUSIC_;
            messager.flush("error");
            exit(-1);
        }
    }
    if (mode == 10) {     // hydro_analysis output
        std::ifstream particletable(table_path_ + "/EOS_particletable.dat");
        particletable >> N_stableparticle;
        particletable.close();
    }

    // read particle resonance decay table
    Nparticle = read_resonances_list(particle_ptr);
    messager << "total number of particle species: " << Nparticle;
    messager.flush("info");

    if (N_stableparticle > 0) {
        messager.info(" -- EOS is partially chemical equilibrium ");
        flag_PCE_ = 1;
        int FO_length = surf_ptr.size();
        double** particle_mu = new double* [N_stableparticle];
        for (int i = 0; i < N_stableparticle; i++)
            particle_mu[i] = new double[FO_length];

        if (mode == 0)
            read_decdat_mu(FO_length, N_stableparticle, particle_mu);
        else if (mode == 1 || mode == 2)
            read_chemical_potentials_music(FO_length, surf_ptr,
                                           N_stableparticle, particle_mu);

        calculate_particle_mu_PCE(Nparticle, surf_ptr, FO_length,
                                  particle_ptr, particle_mu);
        for (int i = 0; i < N_stableparticle; i++)
            delete [] particle_mu[i];
        delete [] particle_mu;
    } else {
        messager.info(" -- EOS is chemical equilibrium. ");
        flag_PCE_ = 0;
    }
}

void read_FOdata::read_decdat(std::vector<FO_surf> &surf_ptr) {
    double temp, temp_vx, temp_vy;
    cout <<" -- Read in information on freeze out surface...";
    ostringstream decdat_stream;
    decdat_stream << path_ << "/decdat2.dat";
    std::ifstream decdat(decdat_stream.str().c_str());
    string input;
    getline(decdat, input, '\n');
    while (!decdat.eof()) {
        std::stringstream ss(input);

        FO_surf surf_elem;

        ss >> surf_elem.tau;
        surf_elem.eta = 0.0;

        ss >> surf_elem.da0;
        ss >> surf_elem.da1;
        ss >> surf_elem.da2;
        surf_elem.da3 = 0.0;

        ss >> temp_vx;
        ss >> temp_vy;
        surf_elem.u0 = 1./sqrt(1. - temp_vx*temp_vx - temp_vy*temp_vy);
        surf_elem.u1 = surf_elem.u0*temp_vx;
        surf_elem.u2 = surf_elem.u0*temp_vy;
        surf_elem.u3 = 0.0;

        ss >> surf_elem.Edec;
        ss >> surf_elem.Bn;
        ss >> surf_elem.Tdec;
        ss >> surf_elem.muB;
        ss >> surf_elem.muS;
        ss >> surf_elem.Pdec;

        ss >> surf_elem.pi33;
        ss >> surf_elem.pi00;
        ss >> surf_elem.pi01;
        ss >> surf_elem.pi02;
        surf_elem.pi03 = 0.0;
        ss >> surf_elem.pi11;
        ss >> surf_elem.pi12;
        surf_elem.pi13 = 0.0;
        ss >> surf_elem.pi22;
        surf_elem.pi23 = 0.0;

        ss >> temp;
        if (turn_on_bulk == 1)
            surf_elem.bulkPi = temp;
        else
            surf_elem.bulkPi = 0.0;

        surf_elem.qmu0 = 0.0e0;
        surf_elem.qmu1 = 0.0e0;
        surf_elem.qmu2 = 0.0e0;
        surf_elem.qmu3 = 0.0e0;

        surf_ptr.push_back(surf_elem);
        getline(decdat, input, '\n');
    }
    decdat.close();
    cout << "done" << endl;
    return;
}

void read_FOdata::read_surfdat(std::vector<FO_surf> &surf_ptr) {
    cout<<" -- Read spatial positions of freeze out surface...";
    ostringstream surfdat_stream;
    double dummy;
    char rest_dummy[512];
    surfdat_stream << path_ << "/surface.dat";
    std::ifstream surfdat(surfdat_stream.str().c_str());
    for (auto &surf_i: surf_ptr) {
        surfdat >> dummy >> dummy;
        surfdat >> surf_i.xpt;
        surfdat >> surf_i.ypt;
        surfdat.getline(rest_dummy, 512);
    }
    surfdat.close();
    cout << "done" << endl;
    return;
}

void read_FOdata::read_FOsurfdat_VISH2p1(std::vector<FO_surf> &surf_ptr) {
    cout << " -- Loading the decoupling data from VISH2+1 ...." << endl;
    // read the data arrays for the decoupling information
    read_decdat(surf_ptr);
    // read the positions of the freeze out surface
    read_surfdat(surf_ptr);
    return;
}

void read_FOdata::read_FOsurfdat_MUSIC_boost_invariant(
                std::vector<FO_surf> &surf_ptr, std::string surface_filename) {
    cout << " -- Read spatial positions of freeze out surface from MUSIC "
         << "(boost-invariant) ...";
    ostringstream surfdat_stream;
    double dummy;
    string input;
    double temp_tau, temp_xpt, temp_ypt, temp_eta;
    int idx = 0;
    surfdat_stream << path_ << "/" << surface_filename;
    std::ifstream surfdat;
    if (surface_in_binary) {
        surfdat.open(surfdat_stream.str().c_str(), std::ios::binary);
    } else {
        surfdat.open(surfdat_stream.str().c_str());
    }
    if (!surfdat.good()) {
        cout << "[Error] Surface file is not found! " << surfdat_stream.str()
             << endl;
        exit(1);
    }
    while (!surfdat.eof()) {
        FO_surf surf_elem;
        if (surface_in_binary) {
            float array[34];
            for (int ii = 0; ii < 34; ii++) {
                float temp = 0.;
                surfdat.read(reinterpret_cast<char*>(&temp), sizeof(float));
                array[ii] = temp;
            }
            surf_elem.tau = array[0];
            surf_elem.xpt = array[1];
            surf_elem.ypt = array[2];
            surf_elem.eta = 0.0;
            surf_elem.da0 = array[4];
            surf_elem.da1 = array[5];
            surf_elem.da2 = array[6];
            surf_elem.da3 = 0.0;
            surf_elem.u0  = array[8];
            surf_elem.u1  = array[9];
            surf_elem.u2  = array[10];
            surf_elem.u3  = array[11];

            surf_elem.Edec = array[12]*hbarC;
            surf_elem.Tdec = array[13]*hbarC;
            surf_elem.muB  = array[14]*hbarC;
            surf_elem.muS  = array[15]*hbarC;
            surf_elem.muC  = array[16]*hbarC;
            surf_elem.Pdec = array[17]*surf_elem.Tdec - surf_elem.Edec;

            surf_elem.pi00 = array[18]*hbarC;  // GeV/fm^3
            surf_elem.pi01 = array[19]*hbarC;  // GeV/fm^3
            surf_elem.pi02 = array[20]*hbarC;  // GeV/fm^3
            surf_elem.pi03 = array[21]*hbarC;  // GeV/fm^3
            surf_elem.pi11 = array[22]*hbarC;  // GeV/fm^3
            surf_elem.pi12 = array[23]*hbarC;  // GeV/fm^3
            surf_elem.pi13 = array[24]*hbarC;  // GeV/fm^3
            surf_elem.pi22 = array[25]*hbarC;  // GeV/fm^3
            surf_elem.pi23 = array[26]*hbarC;  // GeV/fm^3
            surf_elem.pi33 = array[27]*hbarC;  // GeV/fm^3

            surf_elem.bulkPi = array[28]*hbarC;   // GeV/fm^3

            surf_elem.Bn   = array[29];             // 1/fm^3
            surf_elem.qmu0 = array[30];
            surf_elem.qmu1 = array[31];
            surf_elem.qmu2 = array[32];
            surf_elem.qmu3 = array[33];
        } else {
            getline(surfdat, input, '\n');
            std::stringstream ss(input);

            ss >> temp_tau >> temp_xpt >> temp_ypt >> temp_eta;
            // freeze out position
            surf_elem.tau = temp_tau;
            surf_elem.xpt = temp_xpt;
            surf_elem.ypt = temp_ypt;
            surf_elem.eta = 0.0;

            // freeze out normal vectors
            ss >> surf_elem.da0;
            ss >> surf_elem.da1;
            ss >> surf_elem.da2;
            ss >> surf_elem.da3;
            surf_elem.da3 = 0.0;

            // flow velocity
            ss >> surf_elem.u0;
            ss >> surf_elem.u1;
            ss >> surf_elem.u2;
            ss >> surf_elem.u3;

            // thermodynamic quantities at freeze out
            ss >> dummy; surf_elem.Edec = dummy*hbarC;   
            ss >> dummy; surf_elem.Tdec = dummy*hbarC;
            ss >> dummy; surf_elem.muB = dummy*hbarC;
            ss >> dummy; surf_elem.muS = dummy*hbarC;
            ss >> dummy; surf_elem.muC = dummy*hbarC;
            ss >> dummy;              // (e+P)/T
            surf_elem.Pdec = dummy*surf_elem.Tdec - surf_elem.Edec;

            // dissipative quantities at freeze out
            ss >> dummy;                       // 1/fm^4
            surf_elem.pi00 = dummy*hbarC;  // GeV/fm^3
            ss >> dummy; surf_elem.pi01 = dummy*hbarC;
            ss >> dummy; surf_elem.pi02 = dummy*hbarC;
            ss >> dummy; surf_elem.pi03 = dummy*hbarC;
            ss >> dummy; surf_elem.pi11 = dummy*hbarC;
            ss >> dummy; surf_elem.pi12 = dummy*hbarC;
            ss >> dummy; surf_elem.pi13 = dummy*hbarC;
            ss >> dummy; surf_elem.pi22 = dummy*hbarC;
            ss >> dummy; surf_elem.pi23 = dummy*hbarC;
            ss >> dummy; surf_elem.pi33 = dummy*hbarC;

            if (turn_on_bulk == 1) {
                ss >> dummy; surf_elem.bulkPi = dummy*hbarC;  // GeV/fm^3
            } else {
                surf_elem.bulkPi = 0.0;
            }
            if (turn_on_rhob == 1) {
                ss >> dummy; surf_elem.Bn = dummy;   // 1/fm^3
            } else {
                surf_elem.Bn = 0.0;
            }
            if (turn_on_diff == 1) {
                ss >> surf_elem.qmu0;
                ss >> surf_elem.qmu1;
                ss >> surf_elem.qmu2;
                ss >> surf_elem.qmu3;
            } else {
                surf_elem.qmu0 = 0.0e0;
                surf_elem.qmu1 = 0.0e0;
                surf_elem.qmu2 = 0.0e0;
                surf_elem.qmu3 = 0.0e0;
            }
        }
        idx++;
        if (!surfdat.eof()) {
            surf_ptr.push_back(surf_elem);
        }
    }
    surfdat.close();
    cout << "done" << endl;
    return;
}

void read_FOdata::read_FOsurfdat_hydro_analysis_boost_invariant(
                                        std::vector<FO_surf> &surf_ptr) {
    cout << " -- Read spatial positions of freeze out surface from "
         << "hydro_analysis (boost-invariant) ...";
    ostringstream surfdat_stream;
    string input;
    double temp_tau, temp_xpt, temp_ypt;
    double temp_vx, temp_vy;
    int idx = 0;
    surfdat_stream << path_ << "/hyper_surface_2+1d.dat";
    std::ifstream surfdat(surfdat_stream.str().c_str());
    getline(surfdat, input, '\n' );
    while (!surfdat.eof()) {
        std::stringstream ss(input);
        ss >> temp_tau >> temp_xpt >> temp_ypt;

        FO_surf surf_elem;

        // freeze out position
        surf_elem.tau = temp_tau;
        surf_elem.xpt = temp_xpt;
        surf_elem.ypt = temp_ypt;
        surf_elem.eta = 0.0;

        // freeze out normal vectors
        ss >> surf_elem.da0;
        ss >> surf_elem.da1;
        ss >> surf_elem.da2;
        surf_elem.da3 = 0.0;

        // thermodynamic quantities at freeze out
        ss >> surf_elem.Tdec;

        // flow velocity
        ss >> temp_vx >> temp_vy;

        surf_elem.u0 = 1./sqrt(1. - temp_vx*temp_vx - temp_vy*temp_vy);
        surf_elem.u1 = surf_elem.u0*temp_vx;
        surf_elem.u2 = surf_elem.u0*temp_vy;
        surf_elem.u3 = 0.0;

        surf_elem.Edec = 0.0;   
        surf_elem.muB = 0.0;
        surf_elem.Pdec = 0.0;
        surf_elem.muS = 0.0;

        // dissipative quantities at freeze out
        surf_elem.pi00 = 0.0;  // GeV/fm^3
        surf_elem.pi01 = 0.0;
        surf_elem.pi02 = 0.0;
        surf_elem.pi03 = 0.0;
        surf_elem.pi11 = 0.0;
        surf_elem.pi12 = 0.0;
        surf_elem.pi13 = 0.0;
        surf_elem.pi22 = 0.0;
        surf_elem.pi23 = 0.0;
        surf_elem.pi33 = 0.0;

        surf_elem.bulkPi = 0.0;
        surf_elem.Bn = 0.0;

        surf_elem.qmu0 = 0.0e0;
        surf_elem.qmu1 = 0.0e0;
        surf_elem.qmu2 = 0.0e0;
        surf_elem.qmu3 = 0.0e0;

        surf_ptr.push_back(surf_elem);

        idx++;
        getline(surfdat, input, '\n' );
    }
    surfdat.close();

    cout << "done" << endl;
    return;
}

void read_FOdata::read_FOsurfdat_MUSIC(std::vector<FO_surf> &surf_ptr,
                                       std::string surface_filename) {
    cout << " -- Read spatial positions of freeze out surface from MUSIC...";
    ostringstream surfdat_stream;
    double dummy;
    surfdat_stream << path_ << "/" << surface_filename;
    std::ifstream surfdat;
    if (surface_in_binary) {
        surfdat.open(surfdat_stream.str().c_str(), std::ios::binary);
    } else {
        surfdat.open(surfdat_stream.str().c_str());
    }
    if (!surfdat.good()) {
        cout << "[Error] Surface file is not found! " << surfdat_stream.str()
             << endl;
        exit(1);
    }
    while (!surfdat.eof()) {
        FO_surf surf_elem;
        if (surface_in_binary) {
            float array[34];
            for (int i = 0; i < 34; i++) {
                float temp = 0.;
                surfdat.read(reinterpret_cast<char*>(&temp), sizeof(float));
                array[i] = temp;
            }
            surf_elem.tau = array[0];
            surf_elem.xpt = array[1];
            surf_elem.ypt = array[2];
            surf_elem.eta = array[3];
            surf_elem.da0 = array[4];
            surf_elem.da1 = array[5];
            surf_elem.da2 = array[6];
            surf_elem.da3 = array[7];
            surf_elem.u0  = array[8];
            surf_elem.u1  = array[9];
            surf_elem.u2  = array[10];
            surf_elem.u3  = array[11];

            surf_elem.Edec = array[12]*hbarC;
            surf_elem.Tdec = array[13]*hbarC;
            surf_elem.muB  = array[14]*hbarC;
            surf_elem.muS  = array[15]*hbarC;
            surf_elem.muC  = array[16]*hbarC;
            surf_elem.Pdec = array[17]*surf_elem.Tdec - surf_elem.Edec;

            surf_elem.pi00 = array[18]*hbarC;  // GeV/fm^3
            surf_elem.pi01 = array[19]*hbarC;  // GeV/fm^3
            surf_elem.pi02 = array[20]*hbarC;  // GeV/fm^3
            surf_elem.pi03 = array[21]*hbarC;  // GeV/fm^3
            surf_elem.pi11 = array[22]*hbarC;  // GeV/fm^3
            surf_elem.pi12 = array[23]*hbarC;  // GeV/fm^3
            surf_elem.pi13 = array[24]*hbarC;  // GeV/fm^3
            surf_elem.pi22 = array[25]*hbarC;  // GeV/fm^3
            surf_elem.pi23 = array[26]*hbarC;  // GeV/fm^3
            surf_elem.pi33 = array[27]*hbarC;  // GeV/fm^3

            surf_elem.bulkPi = array[28]*hbarC;   // GeV/fm^3

            surf_elem.Bn   = array[29];             // 1/fm^3
            surf_elem.qmu0 = array[30];
            surf_elem.qmu1 = array[31];
            surf_elem.qmu2 = array[32];
            surf_elem.qmu3 = array[33];
        } else {
            // freeze out position
            surfdat >> surf_elem.tau;
            surfdat >> surf_elem.xpt;
            surfdat >> surf_elem.ypt;
            surfdat >> surf_elem.eta;

            // freeze out normal vectors
            surfdat >> surf_elem.da0;
            surfdat >> surf_elem.da1;
            surfdat >> surf_elem.da2;
            surfdat >> surf_elem.da3;

            // flow velocity
            surfdat >> surf_elem.u0;
            surfdat >> surf_elem.u1;
            surfdat >> surf_elem.u2;
            surfdat >> surf_elem.u3;

            // thermodynamic quantities at freeze out
            surfdat >> dummy; surf_elem.Edec = dummy*hbarC;
            surfdat >> dummy; surf_elem.Tdec = dummy*hbarC;
            surfdat >> dummy; surf_elem.muB = dummy*hbarC;
            surfdat >> dummy; surf_elem.muS = dummy*hbarC;
            surfdat >> dummy; surf_elem.muC = dummy*hbarC;
            surfdat >> dummy;                    //(e+p)/T
            surf_elem.Pdec = dummy*surf_elem.Tdec - surf_elem.Edec;

            // dissipative quantities at freeze out
            surfdat >> dummy; surf_elem.pi00 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi01 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi02 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi03 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi11 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi12 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi13 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi22 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi23 = dummy*hbarC;
            surfdat >> dummy; surf_elem.pi33 = dummy*hbarC;
            if (turn_on_bulk == 1) {
                surfdat >> dummy; surf_elem.bulkPi = dummy*hbarC;
            } else {
                surf_elem.bulkPi = 0.0;
            }
            if (turn_on_rhob == 1) {
                surfdat >> dummy; surf_elem.Bn = dummy;   // 1/fm^3
            } else {
                surf_elem.Bn = 0.0;
            }
            if (turn_on_diff == 1) {
                surfdat >> surf_elem.qmu0;
                surfdat >> surf_elem.qmu1;
                surfdat >> surf_elem.qmu2;
                surfdat >> surf_elem.qmu3;
            } else {
                surf_elem.qmu0 = 0.0e0;
                surf_elem.qmu1 = 0.0e0;
                surf_elem.qmu2 = 0.0e0;
                surf_elem.qmu3 = 0.0e0;
            }
        }
        if (!surfdat.eof()) {
            if (surf_elem.Tdec > 0.01)
                surf_ptr.push_back(surf_elem);
        }
    }
    surfdat.close();
    cout << "done" << endl;
}


void read_FOdata::regulate_surface_cells(std::vector<FO_surf> &surf_ptr) {
    double pi_init[4][4];
    double pi_reg[4][4];
    double u_flow[4];
    for (auto &surf_i: surf_ptr) {
        surf_i.u0 = sqrt(1. + surf_i.u1*surf_i.u1
                            + surf_i.u2*surf_i.u2
                            + surf_i.u3*surf_i.u3);
        surf_i.qmu0 = ((  surf_i.u1*surf_i.qmu1
                        + surf_i.u2*surf_i.qmu2
                        + surf_i.u3*surf_i.qmu3)/surf_i.u0);
        u_flow[0] = surf_i.u0;
        u_flow[1] = surf_i.u1;
        u_flow[2] = surf_i.u2;
        u_flow[3] = surf_i.u3;
        pi_init[0][0] = surf_i.pi00;
        pi_init[0][1] = surf_i.pi01;
        pi_init[0][2] = surf_i.pi02;
        pi_init[0][3] = surf_i.pi03;
        pi_init[1][0] = surf_i.pi01;
        pi_init[1][1] = surf_i.pi11;
        pi_init[1][2] = surf_i.pi12;
        pi_init[1][3] = surf_i.pi13;
        pi_init[2][0] = surf_i.pi02;
        pi_init[2][1] = surf_i.pi12;
        pi_init[2][2] = surf_i.pi22;
        pi_init[2][3] = surf_i.pi23;
        pi_init[3][0] = surf_i.pi03;
        pi_init[3][1] = surf_i.pi13;
        pi_init[3][2] = surf_i.pi23;
        pi_init[3][3] = surf_i.pi33;

        regulate_Wmunu(u_flow, pi_init, pi_reg);

        surf_i.pi00 = pi_reg[0][0];
        surf_i.pi01 = pi_reg[0][1];
        surf_i.pi02 = pi_reg[0][2];
        surf_i.pi03 = pi_reg[0][3];
        surf_i.pi11 = pi_reg[1][1];
        surf_i.pi12 = pi_reg[1][2];
        surf_i.pi13 = pi_reg[1][3];
        surf_i.pi22 = pi_reg[2][2];
        surf_i.pi23 = pi_reg[2][3];
        surf_i.pi33 = pi_reg[3][3];
    }
}


void read_FOdata::read_decdat_mu(int FO_length, int N_stable, 
                                 double** particle_mu) {
    cout << " -- Read chemical potential for stable particles...";
    ostringstream decdat_mu_stream;
    double dummy;
    decdat_mu_stream << path_ << "/decdat_mu.dat";
    std::ifstream decdat_mu(decdat_mu_stream.str().c_str());

    //For backward compatibility: decdat_mu.dat can be one line or FO_length lines
    for (int j = 0; j < FO_length; j++) {
        decdat_mu >> dummy;  //not used in the code plz ignore it

        if (decdat_mu.eof()) {
            for(int k=j; k<FO_length; k++)
                for(int i=0; i<N_stable; i++)
                    particle_mu[i][k]=particle_mu[i][j-1];
            break;
        }

        for(int i=0; i<N_stable; i++) {
            decdat_mu >> particle_mu[i][j];
        }
    }

    cout<<"done" << endl;
    return;
}

void read_FOdata::read_chemical_potentials_music(
    int FO_length, std::vector<FO_surf> &FOsurf_ptr, int N_stable, double** particle_mu) {
    cout << " -- Interpolating chemical potentials for stable particles "
         << "(MUSIC IEOS = " << iEOS_MUSIC_ << ") ...";

    Table mu_table;
    if (iEOS_MUSIC_ == 3) {
        mu_table.loadTableFromFile(
                    table_path_ + "/EOS_tables/s95p-v1-PCE150/EOS_Mu.dat");
    } else if (iEOS_MUSIC_ == 4) {
        mu_table.loadTableFromFile(
                    table_path_ + "/EOS_tables/s95p-v1-PCE155/EOS_Mu.dat");
    } else if (iEOS_MUSIC_ == 5) {
        mu_table.loadTableFromFile(
                    table_path_ + "/EOS_tables/s95p-v1-PCE160/EOS_Mu.dat");
    } else if (iEOS_MUSIC_ == 6) {
        mu_table.loadTableFromFile(
                    table_path_ + "/EOS_tables/s95p-v1-PCE165/EOS_Mu.dat");
    }

    double edec_pre = 0.0e0;
    for (int j = 0; j < FO_length; j++) {
        double edec = FOsurf_ptr[j].Edec;
        if (std::abs(edec - edec_pre) > 1e-15) {
            edec_pre = edec;
            for (int i = 0; i < N_stable; i++) {
                particle_mu[i][j] = mu_table.interp(1, i+2, edec);
            }
        } else {
            for (int i = 0; i < N_stable; i++)
                particle_mu[i][j] = particle_mu[i][j-1];
        }
    }

    cout << "done" << endl;
    return;
}

int read_FOdata::read_resonances_list(std::vector<particle_info> &particle) {
    double eps = 1e-15;
    cout << " -- Read in particle resonance decay table...";
    std::string reso_filename;
    if (afterburner_type_ == AfterburnerType::SMASH) {
        reso_filename = particle_table_path_ + "/pdg-SMASH.dat";
    } else if (afterburner_type_ == AfterburnerType::UrQMD) {
        reso_filename = particle_table_path_ + "/pdg-urqmd_v3.3+.dat";
    } else {
        reso_filename = particle_table_path_ + "/pdg-s95pv1.dat";
    }
    std::ifstream resofile(reso_filename.c_str());
    if (!resofile.good()) {
        cout << "[Error] Can not found pdg file: " << reso_filename << endl;
        exit(1);
    }

    int local_i = 0;
    int dummy_int;
    while (!resofile.eof()) {
        particle_info particle_i;

        resofile >> particle_i.monval;
        resofile >> particle_i.name;
        resofile >> particle_i.mass;
        resofile >> particle_i.width;
        resofile >> particle_i.gspin;        //spin degeneracy
        resofile >> particle_i.baryon;
        resofile >> particle_i.strange;
        resofile >> particle_i.charm;
        resofile >> particle_i.bottom;
        resofile >> particle_i.gisospin;     //isospin degeneracy
        resofile >> particle_i.charge;
        resofile >> particle_i.decays;
        for (int j = 0; j < particle_i.decays; j++) {
            decay_channel_info *temp_decay_channel = new decay_channel_info;
            resofile >> dummy_int;
            resofile >> temp_decay_channel->decay_Npart;
            resofile >> temp_decay_channel->branching_ratio;
            resofile >> temp_decay_channel->decay_part[0];
            resofile >> temp_decay_channel->decay_part[1];
            resofile >> temp_decay_channel->decay_part[2];
            resofile >> temp_decay_channel->decay_part[3];
            resofile >> temp_decay_channel->decay_part[4];
            std::string rest_of_line;
            std::getline(resofile, rest_of_line);
            particle_i.decay_channels.push_back(temp_decay_channel);
        }

        //decide whether particle is stable under strong interactions
        if (particle_i.decay_channels[0]->decay_Npart == 1) {
            particle_i.stable = 1;
        } else {
            particle_i.stable = 0;
        }

        if (!resofile.eof()) {
            particle.push_back(particle_i);
        } else {
            particle_i.baryon = 0;
        }

        //add anti-particle entry
        if (particle_i.baryon > 0) {
            local_i++;
            particle_info particle_j;
            particle_j.monval = -particle_i.monval;
            ostringstream antiname;
            antiname << "Anti-" << particle_i.name;
            particle_j.name = antiname.str();
            particle_j.mass = particle_i.mass;
            particle_j.width = particle_i.width;
            particle_j.gspin = particle_i.gspin;
            particle_j.baryon = -particle_i.baryon;
            particle_j.strange = -particle_i.strange;
            particle_j.charm = -particle_i.charm;
            particle_j.bottom = -particle_i.bottom;
            particle_j.gisospin = particle_i.gisospin;
            particle_j.charge = -particle_i.charge;
            particle_j.decays = particle_i.decays;
            particle_j.stable = particle_i.stable;
            for (int j = 0; j < particle_j.decays; j++) {
                decay_channel_info *temp_anti_decay_channel = (
                                                    new decay_channel_info);
                temp_anti_decay_channel->decay_Npart = (
                        particle_i.decay_channels[j]->decay_Npart);
                temp_anti_decay_channel->branching_ratio =
                        particle_i.decay_channels[j]->branching_ratio;
                for (int k = 0; k < 5; k++) {
                    int decay_part_monval = (
                            particle_i.decay_channels[j]->decay_part[k]);
                    if (decay_part_monval == 0) {
                        // a null entry
                        temp_anti_decay_channel->decay_part[k] = 0;
                    } else {
                        // find the index for decay particle in the
                        // current resonance table
                        int idx;
                        for (idx = 0; idx < local_i; idx++) {
                            if (particle[idx].monval == decay_part_monval) {
                                break;
                            }
                        }
                        double temp_br = (
                            particle_i.decay_channels[j]->branching_ratio);
                        if (idx == local_i && particle_i.stable == 0
                            && temp_br > eps) {
                            messager << "Can not find decay particle index for "
                                     << "anti-baryon!";
                            messager.flush("error");
                            messager << "particle monval : "
                                     << decay_part_monval;
                            messager.flush("error");
                            exit(1);
                        }
                        if (particle[idx].baryon == 0
                            && particle[idx].charge == 0
                            && particle[idx].strange == 0) {
                            temp_anti_decay_channel->decay_part[k] = (
                                particle_i.decay_channels[j]->decay_part[k]);
                        } else {
                            temp_anti_decay_channel->decay_part[k] = (
                                - particle_i.decay_channels[j]->decay_part[k]);
                        }
                    }
                }
                particle_j.decay_channels.push_back(temp_anti_decay_channel);
            }
            particle.push_back(particle_j);
        }
        local_i++;   // Add one to the counting variable "i" for the meson/baryon
    }
    for (auto &particle_i: particle) {
        if (quantum_statistics_) {
            if (particle_i.baryon == 0) {
                particle_i.sign = -1;
            } else {
                particle_i.sign = 1;
            }
        } else {
            particle_i.sign = 0;
        }
    }
    cout << "done." << endl;
    return(particle.size());
}

void read_FOdata::calculate_particle_mu_PCE(int Nparticle,
                                            std::vector<FO_surf> &FOsurf_ptr,
                                            int FO_length,
                                            std::vector<particle_info> &particle,
                                            double** particle_mu) {
    int Nstable_particle;
    int Idummy;
    char cdummy[256];
    cout << " -- Read particle table and calculating chemical potential "
         << "for particles..." << endl;
    std::ifstream particletable;
    if (mode == 0) {
        particletable.open(table_path_ + "/EOS_particletable.dat");
    } else if (mode == 1 || mode == 2) {
        if (iEOS_MUSIC_ == 3) {
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE150/EOS_particletable.dat");
        } else if (iEOS_MUSIC_ == 4) {
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE155/EOS_particletable.dat");
        } else if (iEOS_MUSIC_ == 5) {
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE160/EOS_particletable.dat");
        } else if (iEOS_MUSIC_ == 6) {
            particletable.open(
                table_path_
                + "/EOS_tables/s95p-v1-PCE165/EOS_particletable.dat");
        } else {
            messager << "invalid EOS option for MUSIC: " << iEOS_MUSIC_;
            messager.flush("error");
            exit(-1);
        }
    } else {
        messager << "invalid hydro mode: " << mode;
        messager.flush("error");
        exit(-1);
    }

    particletable >> Nstable_particle;
    std::vector<int> stable_particle_monval(Nstable_particle, 0);
    for (int i = 0; i < Nstable_particle; i++) {
        particletable >> Idummy >> stable_particle_monval[i];
        particletable.getline(cdummy, 256);
    }
    particletable.close();

    for (int k = 0; k < FO_length; k++) {
        FOsurf_ptr[k].particle_mu_PCE.resize(Nparticle, 0.0);
    }

    // assign chemical potentials for stable particles
    for (int i = 0; i < Nstable_particle; i++) {
        for (int j = 0; j < Nparticle; j++) {
            if (particle[j].monval == stable_particle_monval[i]) {
                particle[j].stable = 1;
                for (int k=0; k < FO_length; k++)
                    FOsurf_ptr[k].particle_mu_PCE[j] = particle_mu[i][k];
                break;
            }
        }
    }

    // calculating chemical potentials for unstable resonances
    print_progressbar(-1);
    for (int i = 0; i < Nparticle; i++) {
        if (particle[i].stable == 0) {
            for (int j = 0; j < particle[i].decays; j++) {
                for (int k = 0;
                     k < std::abs(particle[i].decay_channels[j]->decay_Npart);
                     k++) {
                    for (int l = 0; l < Nparticle; l++) {
                        if (particle[i].decay_channels[j]->decay_part[k]
                            == particle[l].monval) {
                            for (int m = 0; m < FO_length; m++)
                                FOsurf_ptr[m].particle_mu_PCE[i] += (
                                    particle[i].decay_channels[j]->branching_ratio
                                    *FOsurf_ptr[m].particle_mu_PCE[l]);
                            break;
                        }
                        if (l == Nparticle-1) {
                            messager << "warning: can not find particle"
                                     <<  particle[i].name;
                            messager.flush("error");
                        }
                    }
                }
            }
        }
        print_progressbar(static_cast<double>(i)/Nparticle);
    }
    print_progressbar(1);
}


void read_FOdata::regulate_Wmunu(double u[4], double Wmunu[4][4],
                                 double Wmunu_regulated[4][4]) {
    double gmunu[4][4] = {
        {-1, 0, 0, 0},
        { 0, 1, 0, 0},
        { 0, 0, 1, 0},
        { 0, 0, 0, 1}
    };
    double u_dot_pi[4];
    double u_mu[4];
    for (int i = 0; i < 4; i++) {
        u_dot_pi[i] = (- u[0]*Wmunu[0][i] + u[1]*Wmunu[1][i]
                       + u[2]*Wmunu[2][i] + u[3]*Wmunu[3][i]);
        u_mu[i] = gmunu[i][i]*u[i];
    }
    double tr_pi = - Wmunu[0][0] + Wmunu[1][1] + Wmunu[2][2] + Wmunu[3][3];
    double u_dot_pi_dot_u = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            u_dot_pi_dot_u += u_mu[i]*Wmunu[i][j]*u_mu[j];
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Wmunu_regulated[i][j] = (
                Wmunu[i][j] + u[i]*u_dot_pi[j] + u[j]*u_dot_pi[i]
                + u[i]*u[j]*u_dot_pi_dot_u 
                - 1./3.*(gmunu[i][j] + u[i]*u[j])*(tr_pi + u_dot_pi_dot_u));
        }
    }
}


