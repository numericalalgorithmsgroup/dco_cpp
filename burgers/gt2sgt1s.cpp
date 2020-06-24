#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "burgers.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_BT;
typedef dco::gt1s<DCO_BT>::type DCO_T;

int main() {
  vector<DCO_T> y(n);
  for (int i=1;i<n-1;i++) y[i]=sin(2*PI*i/n);
  dco::value(dco::derivative(y[24]))=1;
  dco::derivative(dco::value(y[25]))=1;
  burgers(y);
  vector<double> ddydyy_v_v(dco::derivative(dco::derivative(y)));
  for(int i=1;i<n-1;i++) 
    cout << i << ": " << dco::passive_value(y[i]) << " "
         << ddydyy_v_v[i] << endl;
  return 0;
}
