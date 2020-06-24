#include <iostream>
#include <vector>
using namespace std;

#include "burgers.h"

#include "dco.hpp"
typedef double DCO_BT;
typedef dco::gt1s<DCO_BT>::type DCO_T;

int main() {
  vector<double> y_in(n,0);
  for (int i=1;i<n-1;i++) y_in[i]=sin(2*PI*i/n);
  vector<double> dymdy(n,0);
  vector<DCO_T> y(n);
  for (int i=1;i<n-1;i++) {
    for(int j=0;j<n;j++) y[j]=y_in[j];
    dco::derivative(y[i])=1;
    burgers(y);
    dymdy[i]=dco::derivative(y[n/2]);
  }
  for (int i=1;i<n-1;i++) cout << i << ": " << dco::value(y[i]) << " " << dymdy[i] << endl;
  return 0;
}
