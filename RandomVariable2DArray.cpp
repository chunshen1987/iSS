// Ver. 1.1.3
// Less general than RandomVariable class: supports only arrays.
// After initialization using sample*** to get col and row indices
// samples.
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
#include "RandomVariable2DArray.h"

using namespace std;


//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------
//----------------------------------------------------------------------
RandomVariable2DArray::RandomVariable2DArray(double **data_in, long size_right_in, long size_left_in, double zero)
// Initilize from a given double** array data_in.
// Values smaller than "zero" will be set to "zero".
// If elements of data_in are data_in[i][j], then size_right_in
// is the allowed max range of j and size_left_in is for i.
{
  // get data info
  size_right = size_right_in;
  size_left = size_left_in;
  data_size = size_right*size_left;

  // allocate memory
  invCDF = new vector<double>(data_size+1, zero);

  // copy data; find max and sum at the same time
  // also flatten it into 2 1d arrays: pdf and cdf
  data_max = zero; data_sum = zero;
  long current_index = 0;
  (*invCDF)[0] = 0.0;
  current_index ++;
  
  for (long i=0; i<size_left; i++)
  for (long j=0; j<size_right; j++)
  {
    double val = data_in[i][j];
    
    // enforce positiveness
    if (val<zero) val = zero;

    data_sum += val;

    // copy to data
    (*invCDF)[current_index] = data_sum;

    // index increment
    current_index ++;
  }

}


//----------------------------------------------------------------------
RandomVariable2DArray::~RandomVariable2DArray()
{
  delete invCDF;
}


//----------------------------------------------------------------------
void RandomVariable2DArray::sampleAccToInvCDF(long* idx_left, long* idx_right)
// Sample according to inverse CDF, which gives a (left, right) pair.
// indices start with 0.
{
  long l = binarySearch(invCDF, (*invCDF)[0]+(data_sum-(*invCDF)[0]-1e-15)*drand48());
  *idx_right = l % size_right;
  *idx_left = (l-*idx_right) / size_right;
}


//----------------------------------------------------------------------
inline long RandomVariable2DArray::idx_1d_to_2d_left(long l)
// l -> l/data_row+1
{
  return l / size_right;
}

//----------------------------------------------------------------------
inline long RandomVariable2DArray::idx_1d_to_2d_right(long l)
// l -> l%data_row+1
{
  return l % size_right;
}

//----------------------------------------------------------------------
inline long RandomVariable2DArray::idx_2d_to_1d(long left, long right)
// (left, right) -> left*size_right + right
{
  return left*size_right + right;
}
