// Ver 1.1.2
#include <vector>

#ifndef RandomVariableNDArray_h
#define RandomVariableNDArray_h

using namespace std;

class RandomVariableNDArray
{
  private:
    vector<double>* invCDF; // holds the inverse CDF of the data
    vector<long>* shape; // shape of the data; lower dimension goes first: A[3][2]->shape[0]=2, shape[1]=3
    int dimension; // dimension of data
    vector<long> *product_of_sizes; // list of length for each dimension and their products
    double data_sum; // "tolerance zero value", max and sum
    bool destructive; // true: the array passed in will be re-used to store inverse CDF; false: a new memory area will be allocated
  public:
    RandomVariableNDArray(vector<double>*, vector<long>&, double zero=1e-30, bool destructive=true);
    ~RandomVariableNDArray();
    inline double get_sum() {return data_sum;};
    long sampleAccToInvCDF(long*);
    inline void idx_1_to_n(long, long*);
    inline long idx_n_to_1(long*);
};

#endif

/*----------------------------------------------------------------------
 Change logs:

 03-15-2012:
 -- Ver 1.0:
    First version.
 04-16-2012:
 -- Ver 1.1.1:
    Bug fix: The sampling should use a random number between the
    smallest (first) and largest (last) element of the inverse CDF, instead
    of between 0 and the largest element.
 07-31-2012:
 -- Ver 1.1.2:
    Bug fix: correct inverse CDF implementd.
    Note that the CDF array constructed has one more element than the
    corresponding PDF. The rule is that if the pdf f(x) on x_i is
    f_i=f(x_i), i=1...N, then the CDF F(x) takes value
    F_j=F(x_j)=sum(f(x_i), i<j). Thus F(x_0)=0, F(x_1)=f_0, ..., F(x_N)=1.
    The rule is that if a random number is between F(x_n) and F(x_(n+1)),
    the associated sample point is x_n.
-----------------------------------------------------------------------*/
