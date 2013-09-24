#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "Poisson.h"

#define ZERO 1e-15

using namespace std;

Poisson::Poisson(double lambda_in)
{  
   if(lambda_in > 0)
      resetDistribution(lambda_in);
   else
   {
      cout << "Error: parameter of Poisson distribution is negative!" << endl;
      exit(1);
   }
}

Poisson::~Poisson()
{  
}

double Poisson::pdf(double k)
{
  if (k<0) return 0;
  int kint = (floor)(k);
  double result = exp(kint*log(lambda) - lambda - lgamma(kint+1));
  return(result);
}

void Poisson::resetDistribution(double lambda_in)
// reset parameter for pdf with lambda_in.
// calculate mean, std, and mode for the given pdf
// it then construct the envelop functions for better sampling efficiencies.
{
   if(lambda_in <= 0)
   {
      cout << "Error: parameter of Poisson distribution is negative!" << endl;
      exit(1);
   }
   lambda = lambda_in;
   mean = lambda;
   std = sqrt(lambda);

   mode = floor(lambda);

   // fill in envelop pdf and inverse CDF
   // first test left boundary
   int nstep = 20;
   int step_left, step_right;
   double step_width = std/2.;
   step_right = nstep;  // no constrain boundary on the right side
   //find the boundary at 0 on the left side
   for (step_left=nstep; step_left>0; step_left--) if (mode-step_width*step_left>=0) break;
   // then fill envelop functions
   constructEnvelopTab(mode, step_width, step_left, step_right);
   return;
}

long Poisson::rand()
{
  return (long)sampleUsingPDFAndEnvelopFunc();
}

long Poisson::rand(double lambda)
{
   resetDistribution(lambda);
   return (long)rand();
}
