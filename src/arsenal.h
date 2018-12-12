// Version 1.8.0
// Zhi Qiu

#ifndef arsenal_h
#define arsenal_h

#include <cmath>
#include "stdlib.h"
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::istream;
using std::ostream;

double sixPoint2dInterp(double x, double y,
    double v00, double v01, double v02, double v10, double v11, double v20);

double interpCubicDirect(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpCubicMono(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpLinearDirect(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpLinearMono(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpNearestDirect(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpNearestMono(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);

inline double interpCubic4Points(double A0, double A1, double A2, double A3, double dx, double xx)
// Assuming that A0, A1, A2, A3 are the values located at 0, dx, 2*dx, and 3*dx, interpolate the value
// at xx using cubic polynomials.
{
    double deltaX = xx - dx;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
}

double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy=1e-10);

double invertTableDirect_hook(double xx);
double invertTableDirect(vector<double>* x, vector<double>* y, double y0, double x0, double relative_accuracy=1e-10);

double stringToDouble(string);
vector<double> stringToDoubles(string);
string toLower(string str);
string trim(string str);

vector< vector<double>* >* readBlockData(istream &stream_in);
void releaseBlockData(vector< vector<double>* >* data);

double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom);
double adaptiveSimpsons(double (*f)(double), double a, double b,  double epsilon=1e-15, int maxRecursionDepth=50);

double qiu_simpsons(double (*f)(double), double a, double b, double epsilon=1e-15, int maxRecursionDepth=50);

long binarySearch(vector<double>* A, double value, bool skip_out_of_range=false);

void formatedPrint(ostream&, int count, ...);
long double gamma_function(long double x);
long double log_gamma_function(long double x);
double beta_function(double, double);
double binomial_coefficient(double n, double k);

void print_progressbar(double percentage, int length=50, string symbol="#");
void display_logo(int which=1);

inline bool is_integer(double x, double tolerance=1e-30)
// Check if a double number is close to an integer
{
  if (std::abs(x-round(x))<tolerance) return true;
  else return false;
}

void GaussLegendre_getWeight(int npts,double* x,double* w, double A, double B, int opt);

void get_bin_average_and_count(istream& is, ostream& os, vector<double>* bins, long col_to_bin=0, void (*func)(vector<double>*)=NULL, long wanted_data_columns=-1, bool silence=false); // Note that col_to_bin starts with 1, and bins is assumed to be monotonically increasing

#endif
