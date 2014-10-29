// Ver 1.0.2

#ifndef RandomVariable1DArray_h
#define RandomVariable1DArray_h


class RandomVariable1DArray
{
  private:
    vector<double>* invCDF;
    long data_size;
    double data_zero, data_sum; // "tolerance zero value", max and sum
  public:
    RandomVariable1DArray(vector<double>*, double=1e-30);
    ~RandomVariable1DArray();
    long rand();
    double return_sum();
};

#endif

/*----------------------------------------------------------------------
 Change logs:

 04-11-2012:
 -- Ver 1.0:
    First version.
 04-16-2012:   
 -- Ver 1.0.1:
    Bug fix: The sampling should use a random number between the
    smallest (first) and largest (last) element of the inverse CDF, instead
    of between 0 and the largest element.
 07-31-2012:
 -- Ver 1.0.2:
    Bug fix: correct inverse CDF implementd.
    Note that the CDF array constructed has one more element than the
    corresponding PDF. The rule is that if the pdf f(x) on x_i is
    f_i=f(x_i), i=1...N, then the CDF F(x) takes value
    F_j=F(x_j)=sum(f(x_i), i<j). Thus F(x_0)=0, F(x_1)=f_0, ..., F(x_N)=1.
    The rule is that if a random number is between F(x_n) and F(x_(n+1)),
    the associated sample point is x_n.
-----------------------------------------------------------------------*/
