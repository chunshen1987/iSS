// Copyright @ Chun Shen 2021

#include "Spectators.h"

#include <fstream>
#include <iostream>
#include <sstream>

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
        std::cout << t << " " << x << " " << y << " " << z << std::endl;
        getline(spectatorFile, input, '\n');
    }
}
