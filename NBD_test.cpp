// Ver 1.1
// Note that all calculations are done at a given particle rapidity y; and all
// "y_minus_y_minus_eta_s" appearences in the code are y-y_minus_eta_s.

#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "arsenal.h"
#include "NBD.h"
#include "Stopwatch.h"

using namespace std;

int main()
{
    double r = 0.16;
    double p = 1.0/(1.0+0.16);
    NBD nbd(p,r);
    Stopwatch sw;
    long sum;


    // sw.tic();
    // sum = 0;
    // for (long i=1; i<=10000; i++) sum += binomial_coefficient(i, i/2);
    // cout << "sum=" << sum << endl;
    // sw.toc();
    // cout << "1 takes time: " << sw.takeTime() << endl;

    // sw.tic();
    // sum = 0;
    // for (long i=1; i<=10000; i++) sum += 1.0/(i+1)/beta_function(i-i/2+1, i/2+1);
    // cout << "sum=" << sum << endl;
    // sw.toc();
    // cout << "2 takes time: " << sw.takeTime() << endl;


    // for (long i=0; i<10; i++) cout << i << "  " << nbd.pdf(i) << endl;
    // for (double i=0; i<10; i+=0.25) cout << setw(10) << i << "  " << nbd.envelopPdfTab->map(i) << endl;
    for (double i=0; i<0.1; i+=0.003641) cout << setw(10) << i << "  " << nbd.envelopInvCDFTab->map(i) << endl;
    // nbd.envelopPdfTab->printFunction();

    // for (long i=1; i<=10000000; i++)
    // {
    //     cout << nbd.rand() << endl;
    // }
}


 