// Ver 1.5

#include "TableFunction.h"

#ifndef RandomVariable_h
#define RandomVariable_h


class RandomVariable
{
  public:
    RandomVariable();
    ~RandomVariable();
    TableFunction *pdfTab, *invCDFTab, *envelopPdfTab, *envelopInvCDFTab;
    double drand(double, double);
    virtual double pdf(double);
    virtual double invCDF(double);
    virtual double envelopPdf(double);
    virtual double envelopInvCDF(double);
    double sampleUsingInvCDF(double, double);
    double sampleUsingInvCDF();
    double sampleUsingPDFDirect(double, double, double);
    double sampleUsingPDFDirect(double);
    double sampleUsingPDFAndEnvelopFunc(double, double, double);
    double sampleUsingPDFAndEnvelopFunc(double factor=1);
    double calculateMoments(long, double, double);
    void constructEnvelopTab(double, double, int, int);
    void calculateInvCDFFromPdf(double, double, long number_of_points=30);
    void calculateInvCDFFromPdf(long number_of_points=30);
};

#endif

/*----------------------------------------------------------------------
 Change logs:

 02-04-2012:
 -- Ver 1.0:
    First version.
 02-07-2012:
 -- Ver 1.4:
    Totally rewritten using TableFunction class.
 03-15-2012:
 -- Ver 1.5:
    Functions added: calculateInvCDFFromPdf.
    Bug fix: each table should be deleted separately in the destruction
    function.
-----------------------------------------------------------------------*/
