#include "arsenal.h"
#include <vector>
#include <iostream>

using namespace std;

int main()
{
  vector<double> A(5);
  vector<double> B(5);
  for (int i=0; i<5; i++)
  {
    A[i] = i;
    B[i] = i*i;
  }
  cout << "--- Nearest direct method:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpNearestDirect(&A,&B,xx,true) << endl;
  }
  
  cout << "--- Linear direct method:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpLinearDirect(&A,&B,xx,true) << endl;
  }

  cout << "--- Cubic direct method:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpCubicDirect(&A,&B,xx,true) << endl;
  }

  cout << "--- Nearest mono method:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpNearestMono(&A,&B,xx,true) << endl;
  }  

  cout << "--- Linear mono method:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpLinearMono(&A,&B,xx,true) << endl;
  }

  cout << "--- Cubic mono method:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpCubicMono(&A,&B,xx,true) << endl;
  }

  cout << "--- Inline cubic direct:" << endl;
  for (double xx=-1; xx<7; xx+=0.5)
  {
    cout << "xx=" << xx <<" yy="
      << interpCubic4points(B[0],B[1],B[2],B[3],1,xx) << endl;
  }
    
}

