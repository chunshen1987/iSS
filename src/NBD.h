// Ver 1.4.3

#include <string>
#include <vector>
#include "RandomVariable.h"

#ifndef NBD_h
#define NBD_h


class NBD: public RandomVariable
{
  private:
    double bridge_p, bridge_r; // used to transfer parameters to pdf function
  public:
    NBD(double p=0.5, double r=2);
    double mode, maximum; // maximum corresponding to the "mode". Refreshed by recalculateMode call
    double mean, std; // standard deviation. Refreshed by recalculateMode call
    double pdf(double);
    void recalculateMode(double p, double r); // here "mode" is in statistical language
    long rand();
    long rand(double p, double r);
};

#endif

/*-----------------------------------------------------------------------
 Change logs:

 02-04-2012:
 -- Ver 1.0:
    First version.
 02-07-2012:
 -- Ver 1.4
    Rewritten upon the updated RandomVariable (Ver 1.4) class.
 04-04-2012:
 -- Ver 1.4.1:
    Bug fix: NBD should only generate discrete samples; its continuous
    extension is not properly normalized to 1 and the expectation value
    it generates changes too.
 04-15-2012:
 -- Ver 1.4.2:
    Bug fix: In pdf function, the (long) cast is changed to floor to
    be more explicit, and k is declared to be "int" for better
    compatibility.
 07-06-2012:
 -- Ver 1.4.3:
    When p is close to 1, the returned value will be set to 0.
-----------------------------------------------------------------------*/
