NOT FULLY TESTED

// Ver. 1.1.1
// Less general than RandomVariable class: supports only arrays.
// After initialization using sample*** to get col and row indices
// samples.
// Note that if due to the rounding down problem the last point 
// is always inaccessible.
// Note that the CDF array constructed has one more element than the
// corresponding PDF. The rule is that if the pdf f(x) on x_i is
// f_i=f(x_i), i=1...N, then the CDF F(x) takes value
// F_j=F(x_j)=sum(f(x_i), i<j). Thus F(x_0)=0, F(x_1)=f_0, ..., F(x_N)=1.
// The rule is that if a random number is between F(x_n) and F(x_(n+1)),
// the associated sample point is x_n.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include "arsenal.h"
#include "RandomVariableNDArray.h"


using namespace std;


//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------
//----------------------------------------------------------------------
RandomVariableNDArray::RandomVariableNDArray(vector<double>* data_in, vector<long>& sizes_in, double zero, bool destructive_in)
// Initilize from a given table data_in.
// Values smaller than "zero" will be set to "zero".
// If destructive is set to true, the original data will be over-written
// to be the inverse CDF; otherwise a new memory space will be allocated.
{
  // get data dimension and sizes
  shape = new vector<long>(sizes_in);
  dimension = shape->size();
  product_of_sizes = new vector<long>(dimension);
  long product = 1;
  for (int i=0; i<dimension; i++)
  {
    product *= (*shape)[i];
    (*product_of_sizes)[i] = product;
  }

  // allocating memory
  destructive = destructive_in;
  if (destructive)
  {
    invCDF = data_in;
  }
  else
  {
    invCDF = new vector<double>(*data_in);
  }

  // add one more element
  invCDF->insert(invCDF->begin(),0.0);

  // summing data; find max and sum at the same time
  // also flatten it into 2 1d arrays: pdf and cdf
  data_sum = 0.0;
  for (long current_index=1; current_index<=(*product_of_sizes)[dimension-1]; current_index++)
  {
    double val = (*invCDF)[current_index];

    // enforce positiveness
    if (val<zero) val = zero;

    // copy to 1d invCDF
    data_sum += val;
    (*invCDF)[current_index] = data_sum;

  }

}

//----------------------------------------------------------------------
RandomVariableNDArray::~RandomVariableNDArray()
{
  delete shape;
  delete product_of_sizes;
  if (!destructive) delete invCDF;
}



//----------------------------------------------------------------------
long RandomVariableNDArray::sampleAccToInvCDF(long* indices)
// Sample according to inverse CDF, which gives an array of indices.
// All indices start with 0.
// Also returns the 1d index (starting from 0).
{
  long l = binarySearch(invCDF, (*invCDF)[0]+(data_sum-(*invCDF)[0])*drand48());
  idx_1_to_n(l, indices);
  return l;
}


//----------------------------------------------------------------------
inline void RandomVariableNDArray::idx_1_to_n(long l, long* indices)
{
  long reduced_digit = l;
  for (int i=0; i<dimension; i++)
  {
    indices[i] = reduced_digit % (*shape)[i];
    reduced_digit = (reduced_digit-indices[i])/(*shape)[i];
  }
}

//----------------------------------------------------------------------
inline long RandomVariableNDArray::idx_n_to_1(long* indices)
{
  long l = 0;
  for (int i=0; i<dimension; i++)
  {
    l += indices[i]*(*product_of_sizes)[i];
  }
  return l;
}