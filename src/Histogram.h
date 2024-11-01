// Copyright @ Chun Shen 2020

#ifndef HISOTGRAM_H_
#define HISOTGRAM_H_

#include <string>
#include <vector>

class Histogram {
  private:
    std::vector<double> bin_val_;
    std::vector<double> bin_y_;
    std::vector<double> bin_y_curr_event_;
    std::vector<double> bin_y_variance_;
    std::vector<int> bin_counts_;
    const double eps = 1e-16;
    double bin_min_;
    double bin_max_;
    double nbins_;
    double bin_width_;
    int nev_;

  public:
    Histogram(double x_min, double x_max, int nbins);
    ~Histogram();

    void fill(double x, double y = 0.);

    void add_an_event();

    double get_bin_width() const { return (bin_width_); }

    std::vector<double> get_x() {
        std::vector<double> bin_x(nbins_, 0.);
        for (int i = 0; i < nbins_; i++) {
            if (bin_counts_[i] > 0) {
                bin_x[i] = bin_val_[i] / bin_counts_[i];
            } else {
                bin_x[i] = bin_min_ + (i + 0.5) * bin_width_;
            }
        }
        return (bin_x);
    };

    std::vector<int> get_bin_counts() const { return (bin_counts_); }

    std::vector<double> get_y() {
        std::vector<double> bin_y(nbins_, 0.);
        for (int i = 0; i < nbins_; i++)
            bin_y[i] =
                (bin_y_[i]
                 / std::max(eps, static_cast<double>(bin_counts_[i])));
        return (bin_y);
    };

    void output_histogram(std::string filename);
};

#endif  // HISOTGRAM_H_
