// Ver. 1.0.2
// Less general than RandomVariable class: supports only 1d arrays.
// After initialization use rand() to sample.
// Note that if use sampleAccToInvCDF function to sample, then due to
// the rounding down problem the last point is always inaccessible.
// Note that the CDF array constructed has one more element than the
// corresponding PDF. The rule is that if the pdf f(x) on x_i is
// f_i=f(x_i), i=1...N, then the CDF F(x) takes value
// F_j=F(x_j)=sum(f(x_i), i<j). Thus F(x_0)=0, F(x_1)=f_0, ..., F(x_N)=1.
// The rule is that if a random number is between F(x_n) and F(x_(n+1)),
// the associated sample point is x_n.

#include <iostream>
#include <stdlib.h>
#include "arsenal.h"
#include "RandomVariable1DArray.h"

using namespace std;


//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------
//----------------------------------------------------------------------
RandomVariable1DArray::RandomVariable1DArray(vector<double> *data_in, double zero)
// Initilize from a given double* array data_in.
// Values smaller than "zero" will be set to "zero".
{
  // get data info
  data_size = data_in->size();

  // allocate memory
  invCDF = new vector<double>(data_size+1, zero);

  // generate inverse CDF
  data_sum = 0.0;
  (*invCDF)[0] = 0.0;
  for (long l=0; l<data_size; l++)
  {
    double val = (*data_in)[l];
    
    // enforce positiveness
    if (val<zero) val = zero;

    data_sum += val;

    // copy to inverse CDF
    (*invCDF)[l+1] = data_sum;
  }

}


//----------------------------------------------------------------------
RandomVariable1DArray::~RandomVariable1DArray()
{
  delete invCDF;
}


//----------------------------------------------------------------------
long RandomVariable1DArray::rand()
// Sample according to inverse CDF, which gives a (left, right) pair.
{
  return binarySearch(invCDF, (*invCDF)[0]+(data_sum-(*invCDF)[0]-1e-15)*drand48());
}

//----------------------------------------------------------------------
double RandomVariable1DArray::return_sum()
// Return the maximum value in the inverse CDF.
{
  return data_sum;
}

