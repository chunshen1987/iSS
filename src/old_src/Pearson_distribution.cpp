#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "Pearson_distribution.h"

#define ZERO 1e-15

using namespace std;

PearsonDistribution::PearsonDistribution(double chi_1_in, double chi_2_in,
                                         double chi_3_in, double chi_4_in) {  
    resetDistribution(chi_1_in, chi_2_in, chi_3_in, chi_4_in);
}

PearsonDistribution::~PearsonDistribution() {  
}


//! define the probability distribution function
double PearsonDistribution::pdf(double x) {
    double temp = (x - lambda)/a;
    double temp1 = 1. + temp*temp;
    double result = pow(temp1, -m)*exp(-nu*atan(temp));
    return(result);
}


//! reset parameter for pdf with lambda_in.
//! calculate mean, std, and mode for the given pdf
//! it then construct the envelop functions for better sampling efficiencies.
void PearsonDistribution::resetDistribution(double chi_1_in, double chi_2_in,
                                            double chi_3_in, double chi_4_in) {
    chi_1 = chi_1_in;
    chi_2 = chi_2_in;
    chi_3 = chi_3_in;
    chi_4 = chi_4_in;

    double mu_1 = chi_1;
    double mu_2 = chi_2;
    double mu_3 = chi_3;
    double mu_4 = chi_4 + 3.*chi_2*chi_2;

    double beta1 = mu_3*mu_3/(mu_2*mu_2*mu_2);
    double beta2 = mu_4/(mu_2*mu_2);

    r = 6.*(beta2 - beta1 - 1.)/(2.*beta2 - 3.*beta1 - 6.);
    m = r/2. + 1.;
    nu = (- r*(r - 2.)*(mu_3/pow(mu_2, 1.5))
          /sqrt(16.*(r - 1.) - beta1*(r - 2.)*(r - 2.)));
    a = sqrt(mu_2*(16.*(r - 1.) - beta1*(r - 2.)*(r - 2.)))/4.;
    lambda = chi_1 - (r - 2.)*(mu_3/pow(mu_2, 1.5))*sqrt(mu_2)/4.;

    mean = chi_1;
    std = sqrt(chi_2);
    
    mode = 0.;

    // fill in envelop pdf and inverse CDF
    // first test left boundary
    int nstep = 50;
    int step_left = nstep;
    int step_right = nstep;
    double step_width = std/2.;
    // then fill envelop functions
    constructEnvelopTab(mode, step_width, step_left, step_right);
}

double PearsonDistribution::rand() {
    return(sampleUsingPDFAndEnvelopFunc());
}

double PearsonDistribution::rand(double chi_1, double chi_2,
                               double chi_3, double chi_4) {
    resetDistribution(chi_1, chi_2, chi_3, chi_4);
    return(rand());
}
