// Copyright @ Chun Shen 2021

#include "Spectators.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

Spectators::Spectators(int mode) {
    mode_ = mode;
}

Spectators::~Spectators() {}

void Spectators::readInSpectatorsFromFile(std::string filename) {
    std::ifstream spectatorFile(filename.c_str());
    if (!spectatorFile.good()) {
        std::cout << "[Error]: Spectator file : " << filename << " not found!"
                  << std::endl;
        exit(1);
    }
    std::string input;
    getline(spectatorFile, input, '\n');  // read in comment line
    double t, x, y, z;
    double mass, px, py, rap;
    int e_charge;
    getline(spectatorFile, input, '\n');
    while (!spectatorFile.eof()) {
        std::stringstream ss(input);
        ss >> t >> x >> y >> z >> mass >> px >> py >> rap >> e_charge;
        getline(spectatorFile, input, '\n');
        iSS_Hadron hadron;
        if (e_charge == 0) {
            hadron.pid = 2112;
        } else {
            hadron.pid = 2212;
        }
        hadron.mass = mass;
        hadron.t = t;
        hadron.x = x;
        hadron.y = y;
        hadron.z = z;
        double mT = sqrt(mass*mass + px*px + py*py);
        hadron.E = mT*cosh(rap);
        hadron.px = px;
        hadron.py = py;
        hadron.pz = mT*sinh(rap);
        spectator_list_.push_back(hadron);
    }
}


std::vector<iSS_Hadron> Spectators::getSpectatorList() const {
    return(spectator_list_);
}
