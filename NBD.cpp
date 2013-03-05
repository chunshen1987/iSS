// Ver. 1.4.3
// Use rand(p,r) to sample.
// Can use rand() if p and r are the same as last one (~4 times faster).
// p: success probability; r: number of failure

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "NBD.h"

#define ZERO 1e-15

using namespace std;


//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------
//----------------------------------------------------------------------
NBD::NBD(double p, double r)
{
  recalculateMode(p,r);
}

//----------------------------------------------------------------------
double NBD::pdf(double k_in)
// NBD pdf. Parameters p and r are passed through the bridge_XX private
// variables.
{
  if (k_in<0) return 0;
  int k = (floor)(k_in);
  double prefactor = binomial_coefficient(k+bridge_r-1, k);
  double answer = prefactor*pow(1-bridge_p,bridge_r)*pow(bridge_p,k);
  return answer;
}


//----------------------------------------------------------------------
void NBD::recalculateMode(double p, double r)
// Recalculate mode and maximun corresponding to the given p and r.
// Also copy p and r to the bridge_xx variables.
// Furthermore, it construct the envelop functions for better sampling
// efficiencies.
{
  bridge_p = p; bridge_r = r;
  if (r<=1)
  {
    mode = 1e-30;
    maximum = pow(1-p,r);
  } else
  {
    mode = p*(r-1)/(1-p); // note it is not an integer here
    maximum = pdf(mode);
  }
  mean = p*r/(1-p);
  std = sqrt(p*r)/(1-p);

  //std /=10; // finer

  // fill in envelop pdf and inverse CDF
  // first test left boundary
  int step_left;
  for (step_left=6; step_left>0; step_left--) if (mode-std*step_left>=0) break;
  // then fill envelop functions
  constructEnvelopTab(mode, std, step_left, 6);
}


//----------------------------------------------------------------------
long NBD::rand()
// Sample NBD with last used p and r.
{
  if (bridge_p<ZERO) return 0;
  if (bridge_p+ZERO>1.0) return 0;
  return (long)sampleUsingPDFAndEnvelopFunc();
}

//----------------------------------------------------------------------
long NBD::rand(double p, double r)
// Sample NBD with parameters p and r.
{
  if (p<ZERO) return 0;
  if (p+ZERO>1.0) return 0;
  recalculateMode(p,r);
  return (long)rand();
}
