#include "Histogram.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

Histogram::Histogram(double x_min, double x_max, int nbins) {
    bin_min_ = x_min;
    bin_max_ = x_max;
    nbins_ = nbins;
    nev_ = 0;
    bin_width_ = (bin_max_ - bin_min_)/std::max(eps, (nbins_ - 1));
    bin_val_.resize(nbins_, 0.);
    bin_counts_.resize(nbins_, 0);
    bin_y_.resize(nbins_, 0.);
    bin_y_curr_event_.resize(nbins_, 0.);
    bin_y_variance_.resize(nbins_, 0.);
}


Histogram::~Histogram() {
    bin_val_.clear();
    bin_y_.clear();
    bin_counts_.clear();
}


void Histogram::fill(double x, double y) {
    int idx = static_cast<int>((x - bin_min_)/bin_width_);
    if (idx >= 0 && idx < nbins_) {
        bin_val_[idx] += x;
        //bin_y_[idx] += y;
        bin_y_curr_event_[idx] += y;
        bin_counts_[idx]++;
    }
}


void Histogram::add_an_event() {
    nev_++;
    for (int i = 0; i < nbins_; i++) {
        bin_y_[i] += bin_y_curr_event_[i];
        bin_y_variance_[i] += bin_y_curr_event_[i]*bin_y_curr_event_[i];
        bin_y_curr_event_[i] = 0.;
    }
}


void Histogram::output_histogram(std::string filename) {
    auto bin_x = get_x();
    std::ofstream of(filename.c_str());
    of << "# x  y  y_err  bin_counts" << std::endl;
    for (int i = 0; i < nbins_; i++) {
        double y_err = sqrt(bin_y_variance_[i]/nev_
                            - bin_y_[i]*bin_y_[i]/(nev_*nev_))/sqrt(nev_);
        of << std::scientific << std::setprecision(6) << std::setw(10)
           << bin_x[i] << "  " << bin_y_[i]/nev_ << "  "
           << y_err << "  " << bin_counts_[i] << std::endl;
    }
    of.close();
}
