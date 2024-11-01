// Copyright @ Chun Shen 2021

#ifndef SPECTATORS_H_
#define SPECTATORS_H_

#include <string>
#include <vector>

#include "data_struct.h"

class Spectators {
  private:
    int mode_;
    std::vector<iSS_Hadron> spectator_list_;

  public:
    Spectators(int mode);
    ~Spectators();

    int getMode() const { return (mode_); }

    int getNumberOfSpectators() const { return (spectator_list_.size()); }

    void readInSpectatorsFromFile(std::string filename);

    std::vector<iSS_Hadron> getSpectatorList() const;
};

#endif  // SPECTATORS_H_
