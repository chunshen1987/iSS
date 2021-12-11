// Copyright @ Chun Shen 2021

#ifndef SPECTATORS_H_
#define SPECTATORS_H_

#include <vector>
#include <string>

#include "data_struct.h"

class Spectators {
 private:
    int mode_;
    std::vector<iSS_Hadron> spectator_list_;

 public:
    Spectators(int mode);
    ~Spectators();

    int getMode() const {return(mode_);}

    void readInSpectatorsFromFile(std::string filename);

};

#endif  // SPECTATORS_H_
