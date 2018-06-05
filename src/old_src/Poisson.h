#include <vector>
#include "RandomVariable.h"

#ifndef Poisson_h
#define Poisson_h

class Poisson: public RandomVariable
{
    private:
       double lambda;  // parameter for Poisson distribution

    public:
       Poisson(double lambda_in = 1);
       ~Poisson();

       double mean, std;      // mean and standard deviation
       int mode;              // maximum position of pdf
       double pdf(double k);
       void resetDistribution(double lambda_in);
       
       long rand();           // return random integer according to Poisson distribution
       long rand(double lambda);    // return random integer according to Poisson distribution with given parameter lambda


};

#endif
