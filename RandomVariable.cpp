// Ver. 1.5
// There are 4 functions: pdf, invCDF, envelopPdf, envelopInvCDF.
// The pdf and invCDF functions are for the actual variable.
// The envelopPdf and envelopInvCDF functions are functions used in
// sampling process. The envelopPdf is assumed that when multiplied
// with a factor, it is pointwisely bigger than acutal pdf function.
//
// If any of these functions are not overloaded, then the corresponding
// table function (e.g. pdfTab) is used. These table function should
// be pre-loaded with mapping tables before they can be used, see
// TableFunction class documentation. The interpolation method can be set
// though the "interpolation_model" variable in the TableFunction class.
//
// When either overloaded or loaded with mapping table, the function is
// considered to be "ready".
//
// There are 3 ways to sample a random variables:
//
// (1) Use sampleUsingInvCDF to sample directly according to inverser CDF.
//     Require the invCDF function to be ready.
// (2) Use sampleUsingPDFDirect to sample according to pdf directly,
//     when its maximum is known. Require the pdf function to be ready.
// (3) Use sampleUsingPDFAndEnvelop to sample pdf, using envelop
//     distributions. See function comments for details.
//     Require functions envelopPdf and envelopInvCDF to be ready.
//
// Note that the CDF array should have one more element than the
// corresponding PDF. The rule is that if the pdf f(x) on x_i is
// f_i=f(x_i), i=1...N, then the CDF F(x) takes value
// F_j=F(x_j)=sum(f(x_i), i<j). Thus F(x_0)=0, F(x_1)=f_0, ..., F(x_N)=1.
// The rule is that if a random number is between F(x_n) and F(x_(n+1)),
// the associated sample point is x_n.

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "arsenal.h"
#include "TableFunction.h"
#include "RandomVariable.h"

// used in sampleXXXXXX functions
#define MAXITER 1000
#define DEFAULTSAMPLERANGE 300 // sample in [-#,#]

// used in calculateMoments function
#define EPSILON 1e-8
#define MAXRECURSION 400

RandomVariable* zq_global_randomVariable = NULL;
int zq_global_moments_order = 0;

using namespace std;


//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------
//----------------------------------------------------------------------
RandomVariable::RandomVariable()
{
  pdfTab = new TableFunction;
  invCDFTab = new TableFunction;
  envelopPdfTab = new TableFunction;
  envelopInvCDFTab = new TableFunction;
}

//----------------------------------------------------------------------
RandomVariable::~RandomVariable()
{
  delete pdfTab;
  delete invCDFTab;
  delete envelopPdfTab;
  delete envelopInvCDFTab;
}

//----------------------------------------------------------------------
double RandomVariable::drand(double LB, double RB)
// Uniform distribution between LB and RB with boundary-protection.
{
  double width = RB-LB;
  double dw = width*1e-15;
  return LB+dw+(width-2*dw)*drand48();
}

//----------------------------------------------------------------------
double RandomVariable::pdf(double x)
// The pdf function. The parameter "model" is used in interpolation.
{
  if (pdfTab->isMappingTableLoaded())
  {
    return pdfTab->map(x);
  }
  else
  {
    cout << "RandomVariable::pdf error: no mapping table loaded." << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------
double RandomVariable::invCDF(double x)
// The inverse CDF function. The parameter "model" is used in interpolation.
{
  if (invCDFTab->isMappingTableLoaded())
  {
    return invCDFTab->map(x);
  }
  else
  {
    cout << "RandomVariable::invCDF error: no mapping table loaded." << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------
double RandomVariable::envelopPdf(double x)
// The envelop pdf function. The parameter "model" is used in interpolation.
{
  if (envelopPdfTab->isMappingTableLoaded())
  {
    return envelopPdfTab->map(x);
  }
  else
  {
    cout << "RandomVariable::envelopPdfTab error: no mapping table loaded." << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------
double RandomVariable::envelopInvCDF(double x)
// The envelop inverse CDF function.
{
  if (envelopInvCDFTab->isMappingTableLoaded())
  {
    return envelopInvCDFTab->map(x);
  }
  else
  {
    cout << "RandomVariable::envelopInvCDF error: no mapping table loaded." << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------
double RandomVariable::sampleUsingInvCDF(double LB, double RB)
// Return a random number according to the inverse CDF. The CDF value is
// sampled between LB and RB. (Thus 0<=LB<=RB<=1.)
{
  return invCDF(drand(LB,RB));
}
//----------------------------------------------------------------------
double RandomVariable::sampleUsingInvCDF()
// Short version, only usable if invCDFTab is loaded.
{
  return sampleUsingInvCDF(invCDFTab->getXMin(), invCDFTab->getXMax());
}


//----------------------------------------------------------------------
double RandomVariable::sampleUsingPDFDirect(double normalization_factor, double LB, double RB)
// Sample according to pdf. The function first sample a number between
// LB and RB, then use
// pdf / normalization_factor
// to determine if keeping (rand()<it) a sample or not (rand()>it).
// The normalization_factor factor is chosen so that the pdf function
// is less than normalization_factor.
// -- normalization_factor: see above
// -- LB, RB: only numbers in [LB,RB] are sampled
{
  double x;
  long numberOfIteration = 0;
  while (1)
  {
    x = drand(LB,RB);
    if (drand48()<pdf(x)/(normalization_factor+1e-60)) return x;
    numberOfIteration++;
    if (numberOfIteration>MAXITER)
    {
      cout << "RandomVariable::sampleUsingPDFDirect warning: max number of iteration MAXITER reached." << endl;
      return x;
    }
  }
  return 0;
}
//----------------------------------------------------------------------
double RandomVariable::sampleUsingPDFDirect(double normalization_factor)
// Short version, only usable if pdfTab is loaded.
{
  return sampleUsingPDFDirect(normalization_factor, pdfTab->getXMin(), pdfTab->getXMax());
}



//----------------------------------------------------------------------
double RandomVariable::sampleUsingPDFAndEnvelopFunc(double normalization_factor, double LB, double RB)
// Sample according to pdf. The function first sample according to
// the envelop inverse CDF, then use the ratio
// pdf / (normalization_factor*envelopPdf)
// to determine if keeping (rand()<it) a sample or not (rand()>it).
// The normalization_factor factor is chosen so that the pdf function
// is always less than normalization_factor*envelopPdf.
// If the envelop functions are constructed by constructEnvelopTab function
// then the factor is 1.
// -- normalization_factor: see above
// -- LB, RB: only numbers in [LB,RB] are sampled
{
  double x;
  long numberOfIteration = 0;
  while (1)
  {
    x = envelopInvCDF(drand(LB,RB));
    if (drand48()<pdf(x)/(normalization_factor*envelopPdf(x)+1e-60)) return x;
    numberOfIteration++;
    if (numberOfIteration>MAXITER)
    {
      cout << "RandomVariable::sampleUsingPDFAndEnvelop warning: max number of iteration MAXITER reached." << endl;
      return x;
    }
  }
  return 0;
}
//----------------------------------------------------------------------
double RandomVariable::sampleUsingPDFAndEnvelopFunc(double normalization_factor)
// Short version, only usable if envelopInvCDFTab and envelopPdfTab are loaded.
{
  return sampleUsingPDFAndEnvelopFunc(normalization_factor, envelopInvCDFTab->getXMin(), envelopInvCDFTab->getXMax());
}



//----------------------------------------------------------------------
double zq_moments_integrad_hook(double x) // "pdf", hook version
{
  return pow(x,zq_global_moments_order)*zq_global_randomVariable->pdf(x);
}
//----------------------------------------------------------------------
double RandomVariable::calculateMoments(long order, double LB, double RB)
// Calcualte the moment
// \int x^order pdf(x) dx
{
  zq_global_randomVariable = this;
  zq_global_moments_order = order;
  return qiu_simpsons(&zq_moments_integrad_hook, LB, RB, EPSILON, MAXRECURSION);
}


//----------------------------------------------------------------------
void RandomVariable::constructEnvelopTab(double center, double width, int step_left, int step_right)
// This function construct envelop functions tables using step functions.
// It scan the range from "center - step_left*width" to "center + step_right*width".
// For each interval of size "width", it assumes that the pdf is monotonic
// and the larger of the values at the two ends will be used as the value
// of the envelop pdf in this region, from which envelop CDF will be
// calculated.
// Note that if you overloaded the envelopX functions then this function
// will not work since it only manipulates the envelopXTab tables.
// Note also that the envelop pdf table has width/2 wider range on both
// sides for easy implementation, compared with the envelop inverse
// CDF function, which should not be a problem if the constructed
// table is used by functions like sampleUsingPDFAndEnvelopTab.
{
  envelopPdfTab->deleteMappingTable();
  envelopInvCDFTab->deleteMappingTable();
  double sum = 0;
  double LB = center-step_left*width;
  double RB = LB + width;
  double pdfLB = pdf(LB), pdfRB = pdf(RB), pdfLarger; // x and y values at the two ends
  int idx = 1;
  // initial values
  envelopPdfTab->setMappingTable(idx, LB-width/2, 0);
  envelopInvCDFTab->setMappingTable(idx, 0, LB);
  idx ++;
  // fill in the loop
  for (int ii=0; ii<step_left+step_right; ii++)
  {
    pdfLarger = pdfLB>pdfRB? pdfLB : pdfRB;
    envelopPdfTab->setMappingTable(idx, LB+width/2, pdfLarger);
    sum += width*pdfLarger;
    envelopInvCDFTab->setMappingTable(idx, sum, RB);
    idx ++;
    LB = RB;
    RB += width;
    pdfLB = pdf(LB);
    pdfRB = pdf(RB);
  }
  // final values
  envelopPdfTab->setMappingTable(idx, LB+width/2, 0);
  // set interpolation scheme
  envelopPdfTab->interpolation_model = 10; // nearest neighbor direct interpolation
  envelopInvCDFTab->interpolation_model = 2; // linear mono
}

//----------------------------------------------------------------------
void RandomVariable::calculateInvCDFFromPdf(double LB, double RB, long number_of_points)
// Sum up and pdf from LB to RB using step (RB-LB)/(number_of_points-1), 
// and store the resulted table to the invCDFTab.
// The previous invCDFTab will be deleted.
{
  invCDFTab->resetMappingTable(number_of_points);
  double dx = (double)(RB-LB)/(number_of_points-1);
  /*
  // Trapezoid method
  double pre_x = LB, x = LB, sum = 0.0, pre_val = 0.0, val = 0.0;
  invCDFTab->setMappingTable(1, sum, x);
  for (long i=1; i<number_of_points; i++)
  {
    x += dx;
    val = pdf(x);
    sum += (pre_val + val)*dx/2.0;
    invCDFTab->setMappingTable(i+1, sum, x);
    pre_x = x; pre_val = val;
  }
  */
  // Riemann sum
  double x = LB, sum = 0.0, val = 0.0;
  invCDFTab->setMappingTable(1, sum, x);
  for (long i=1; i<number_of_points; i++)
  {
    x += dx;
    val = pdf(x);
    sum += val*dx;
    invCDFTab->setMappingTable(i+1, sum, x);
  }
}


//----------------------------------------------------------------------
void RandomVariable::calculateInvCDFFromPdf(long number_of_points)
// Short version; usable only if pdfTab is loaded.
{
  double LB = pdfTab->getXMin();
  double RB = pdfTab->getXMax();
  RB = RB-(RB-LB)*1e-10;
  calculateInvCDFFromPdf(LB, RB, number_of_points);
}