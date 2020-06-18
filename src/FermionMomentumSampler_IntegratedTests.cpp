// Copyright @ Chun Shen 2020

#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <memory>
#include <iomanip>
#include "FermionMomentumSampler.h"
#include "Histogram.h"
#include "Random.h"

const int NSAMPLES = 1000000;

double FermiDiracDistribution(double p, double m, double T, double mu);
double GetNormFDDis(double m, double T, double mu);
double PDFFDDis(double p1, double p2, double m, double T, double mu);
void TestFermionSampler(double m, double T, double mu, std::shared_ptr<RandomUtil::Random> ran_gen_ptr);


int main() {
    int randomSeed = -1;
    auto ran_gen_ptr = std::shared_ptr<RandomUtil::Random>(
            new RandomUtil::Random(randomSeed));
    std::vector<double> mass_list = {0.938, 1.5, 2.0, 2.5};
    std::vector<double> T_list = {0.05, 0.1, 0.15, 0.2};
    std::vector<double> mu_list = {-0.9, -0.5, 0.0, 0.5, 0.9};
    for (auto &mass_i: mass_list) {
        for (auto &T_i: T_list) {
            for (auto &mu_i: mu_list) {
                TestFermionSampler(mass_i, T_i, mu_i, ran_gen_ptr);
            }
        }
    }
    return(0);
}


double FermiDiracDistribution(double p, double m, double T, double mu) {
    double E = sqrt(p*p + m*m);
    double res = 1./(exp((E - mu)/T) + 1.);
    return(res);
}


double GetNormFDDis(double m, double T, double mu) {
    double p_max = 5.;
    double dp = 0.1;
    int Np = p_max/dp;
    double res = 0.;
    for (int i = 0; i < Np; i++) {
        double p_local = (i + 0.5)*dp;
        double pdf = FermiDiracDistribution(p_local, m, T, mu);
        res += p_local*p_local*pdf;
    }
    res *= dp;
    return(res);
}


double PDFFDDis(double p1, double p2, double m, double T, double mu) {
    double norm = GetNormFDDis(m, T, mu);
    int npT = 20;
    double dpT = (p2 - p1)/npT;
    double pdf = 0.;
    for (int ipT = 0; ipT < npT; ipT++) {
        double p = p1 + (ipT + 0.5)*dpT;
        pdf += p*p*FermiDiracDistribution(p, m, T, mu)*dpT;
    }
    return(pdf/norm);
}


void TestFermionSampler(double m, double T, double mu,
                      std::shared_ptr<RandomUtil::Random> ran_gen_ptr) {
    std::cout << "Testing the Sampler with m = " << m << " GeV, T = " << T
              << " GeV, mu = " << mu << " GeV." << std::endl;
    FermionMomentumSampler *fermion_sampler;
    double m0tilde = m/T - mu/T;
    if (m0tilde > 50.) {
        fermion_sampler = new FermionMomentumSampler(ran_gen_ptr, 50., 1);
    } else if (m0tilde > 30.) {
        fermion_sampler = new FermionMomentumSampler(ran_gen_ptr, 30., 2);
    } else {
        fermion_sampler = new FermionMomentumSampler(ran_gen_ptr, 0., 10);
    }
    Histogram hist(0., 5., 100);
    for (int i = 0; i < NSAMPLES; i++) {
        double p = fermion_sampler->Sample_a_momentum(m, T, mu);
        hist.fill(p);
    }
    auto hist_dx = hist.get_bin_width();
    auto hist_x = hist.get_bin_x();
    std::vector<double> pdf_bin(hist_x.size() + 1, 0.);
    for (unsigned int i = 0; i < pdf_bin.size(); i++) {
        pdf_bin[i] = i*hist_dx;
    }
    auto hist_y = hist.get_bin_count();
    std::vector<double> hist_y_err(hist_x.size(), 0.);
    std::vector<double> pdf(hist_x.size(), 0.);
    for (unsigned int i = 0; i < hist_x.size(); i++) {
        pdf[i] = PDFFDDis(pdf_bin[i], pdf_bin[i+1], m, T, mu)/hist_dx;
        hist_y_err[i] = sqrt(hist_y[i])/NSAMPLES/hist_dx;
        hist_y[i] = hist_y[i]/NSAMPLES/hist_dx;
    }
    std::stringstream filename;
    filename << "TestFermionSampler_m_" << m << "_T_" << T << "_mu_" << mu
             << ".dat";
    std::ofstream of(filename.str().c_str());
    of << "# pT [GeV]  pdf  sample  sample_err" << std::endl;
    for (unsigned int i = 0; i < hist_x.size(); i++) {
        of << std::scientific << std::setprecision(6)
           << hist_x[i] << "  " << pdf[i] << "  "
           << hist_y[i] << "  " << hist_y_err[i] << std::endl;
    }
    of.close();
    delete fermion_sampler;
}

