// Ver 1.1
// Note that all calculations are done at a given particle rapidity y; and all
// "y_minus_y_minus_eta_s" appearences in the code are y-y_minus_eta_s.

#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include "arsenal.h"
#include "Poisson.h"
#include "Pearson_distribution.h"
using namespace std;

int main() {
    srand48(time(NULL));

    double chi_1 = 0.;
    double chi_2 = 0.1;
    double chi_3 = 0.;
    double chi_4 = 0.1;
    PearsonDistribution test_Pearson(chi_1, chi_2, chi_3, chi_4);

    double lambda = 5.0;
    Poisson test_Poisson(lambda);

    double hist[11];
    for (int i = 0; i < 11; i++) {
        hist[i] = 0.0;
    }
    long n_sample = 50000000;
    double mu_1 = 0.0;
    double mu_2 = 0.0;
    double mu_3 = 0.0;
    double mu_4 = 0.0;
    for (long i = 0; i < n_sample; i++) {
        double net_baryons = test_Pearson.rand();
        double total_baryons = 20.;
        //do {
        //    total_baryons = test_Poisson.rand();
        //} while (total_baryons < fabs(net_baryons));
        //double d_baryons = (total_baryons + net_baryons)/2.;
        //double d_b_bars = (total_baryons - net_baryons)/2.;
        //double frac_baryon = d_baryons - floor(d_baryons);
        //double frac_b_bar = d_b_bars - floor(d_b_bars);
        //int n_baryons = static_cast<int>(d_baryons);
        //if (drand48() < frac_baryon) {
        //    n_baryons++;
        //}
        //int n_b_bar = static_cast<int>(d_b_bars);
        //if (drand48() < frac_b_bar) {
        //    n_b_bar++;
        //}
        //int n_Netbaryons = n_baryons - n_b_bar;
        //double n_Netbaryons = d_baryons - d_b_bars;
        double p_1 = test_Pearson.pdf(floor(net_baryons));
        double p_2 = test_Pearson.pdf(ceil(net_baryons));
        double p = test_Pearson.pdf(net_baryons);
        int n_Netbaryons;
        //if (drand48() < fabs(p_1 - p)/fabs(p_1 - p_2)) {
        if (drand48() < fabs(log(p_1/p)/log(p_1/p_2))) {
            n_Netbaryons = static_cast<int>(ceil(net_baryons));
        } else {
            n_Netbaryons = static_cast<int>(floor(net_baryons));
        }
        if (n_Netbaryons > -6 && n_Netbaryons < 6) {
            hist[n_Netbaryons + 5]++;
        }
        mu_1 += n_Netbaryons;
        mu_2 += n_Netbaryons*n_Netbaryons;
        mu_3 += n_Netbaryons*n_Netbaryons*n_Netbaryons;
        mu_4 += n_Netbaryons*n_Netbaryons*n_Netbaryons*n_Netbaryons;
    }

    for (int i = 0; i < 11; i++) {
        hist[i] = hist[i]/n_sample;
    }

    ofstream check("check_net_baryon_hist.dat");
    for (int i = 0; i < 11; i++) {
        check << scientific << setw(18) << setprecision(8)
              << i - 5 << "  " << hist[i] << endl;
    }
    check.close();

    mu_1 = mu_1/n_sample;
    mu_2 = mu_2/n_sample;
    mu_3 = mu_3/n_sample;
    mu_4 = mu_4/n_sample;

    double sample_chi_1 = mu_1;
    double sample_chi_2 = mu_2 - mu_1*mu_1;
    double sample_chi_3 = (mu_3 - 3.*mu_1*mu_2 + 2.*mu_1*mu_1*mu_1);
    double sample_chi_4 = (mu_4 - 4.*mu_1*mu_3 + 6.*mu_1*mu_1*mu_2
                           - 3.*mu_1*mu_1*mu_1*mu_1);

    cout << scientific << setw(18) << setprecision(8)
         << "check: mean = " << sample_chi_1
         << ", chi_2 = " << sample_chi_2
         << ", chi_3 = " << sample_chi_3
         << ", chi_4 = " << sample_chi_4
         << endl;
}


 
