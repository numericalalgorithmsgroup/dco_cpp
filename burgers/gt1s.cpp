#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "burgers.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_T;

int main() {
  vector<DCO_T> y(n); 
  for (int i=1;i<n-1;i++) y[i]=sin(2*PI*i/n);
  dco::derivative(y[24])=1;
  burgers(y);
  vector<double> yp(dco::value(y)), dydy_v(dco::derivative(y));
  for(int i=1;i<n-1;i++) 
    cout << i << ": " << yp[i] << " " << dydy_v[i] << endl;
  return 0;
}
