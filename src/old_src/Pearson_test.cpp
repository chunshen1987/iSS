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

#include "arsenal.h"
#include "Pearson_distribution.h"
using namespace std;

int main() {
    double hist_size = 10.;
    int hist_length = 101;
    double *hist = new double[hist_length];
    double *hist_x = new double[hist_length];
    double dx = hist_size/(hist_length - 1.);
    for (int i = 0; i < hist_length; i++) {
        hist_x[i] = -hist_size/2. + i*dx;
        hist[i] = 0.0;
    }

    double chi_1 = 0.;
    double chi_2 = 0.1;
    double chi_3 = 0.0;
    double chi_4 = 0.08;
    PearsonDistribution test(chi_1, chi_2, chi_3, chi_4);

    // sw.tic();
    // sum = 0;
    // for (long i=1; i<=10000; i++) sum += binomial_coefficient(i, i/2);
    // cout << "sum=" << sum << endl;
    // sw.toc();
    // cout << "1 takes time: " << sw.takeTime() << endl;

    // sw.tic();
    // sum = 0;
    // for (long i=1; i<=10000; i++) sum += 1.0/(i+1)/beta_function(i-i/2+1, i/2+1);
    // cout << "sum=" << sum << endl;
    // sw.toc();
    // cout << "2 takes time: " << sw.takeTime() << endl;


    // for (long i=0; i<10; i++) cout << i << "  " << nbd.pdf(i) << endl;
    // for (double i=0; i<10; i+=0.25) cout << setw(10) << i << "  " << nbd.envelopPdfTab->map(i) << endl;
    // for (double i=0; i<0.1; i+=0.003641) cout << setw(10) << i << "  " << nbd.envelopInvCDFTab->map(i) << endl;
    // nbd.envelopPdfTab->printFunction();
    
    long n_sample = 5000000;
    for (long i = 0; i < n_sample; i++) {
        double random_sample = test.rand();
        int idx = static_cast<int>((random_sample + (hist_size + dx)/2.)/dx);
        if (idx >= 0 && idx < hist_length) {
            hist[idx] += 1.;
        }
    }

    ofstream checkof("check_Pearson.dat");
    for (int i = 0; i < hist_length; i++) {
        hist[i] /= (n_sample*dx);
        hist[i] += 1e-30;
    }

    double norm = hist[50]/test.pdf(hist_x[50]);
    for (int i = 0; i < hist_length; i++) {
        checkof << scientific << setw(18) << setprecision(6)
                << hist_x[i] << "  "
                << test.pdf(hist_x[i])*norm << "  "
                << hist[i] << endl;
    }

    checkof.close();
    delete[] hist_x;
    delete[] hist;
}


 
