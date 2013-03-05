// Ver 1.1.3

#ifndef RandomVariable2DArray_h
#define RandomVariable2DArray_h


class RandomVariable2DArray
{
  private:
    vector<double>* invCDF;
    long size_left, size_right, data_size; // col, row, and col*row of data
    double data_zero, data_max, data_sum; // "tolerance zero value", max and sum
  public:
    RandomVariable2DArray(double**, long, long, double=1e-30);
    ~RandomVariable2DArray();
    void sampleAccToInvCDF(long*, long*);
    inline long idx_1d_to_2d_left(long);
    inline long idx_1d_to_2d_right(long);
    inline long idx_2d_to_1d(long, long);
};

#endif

/*----------------------------------------------------------------------
 Change logs:

 03-15-2012:
 -- Ver 1.0:
    First version.
 03-27-2012:
 -- Ver 1.1:
    Restricted to only allowing sampleAccToInvCDF function; use internal
    vector data to speed up.
 04-10-2012:
 -- Ver 1.1.1:
    Bug fix: The current_index variable in the constructor should start
    with 0 instead of 1.
 04-15-2012:
 -- Ver 1.1.2:
    Bug fix: The sampling should use a random number between the
    smallest (first) and largest (last) element of the inverse CDF, instead
    of between 0 and the largest element.
 07-31-2012:
 -- Ver 1.1.3:
    Bug fix: correct inverse CDF implementd.
    Note that the CDF array constructed has one more element than the
    corresponding PDF. The rule is that if the pdf f(x) on x_i is
    f_i=f(x_i), i=1...N, then the CDF F(x) takes value
    F_j=F(x_j)=sum(f(x_i), i<j). Thus F(x_0)=0, F(x_1)=f_0, ..., F(x_N)=1.
    The rule is that if a random number is between F(x_n) and F(x_(n+1)),
    the associated sample point is x_n.
-----------------------------------------------------------------------*/
