#include <iostream>
#include <vector>
using namespace std;

#include "burgers.h"

#include "dco.hpp"
typedef double DCO_BT;
typedef dco::gt1s<DCO_BT>::type DCO_T;

#include <omp.h>

int main() {
  vector<double> y_in(n,0),y_out(n,0);
  for (int i=1;i<n-1;i++) y_in[i]=sin(2*PI*i/n);
  omp_set_num_threads(NT);
  vector<double> dymdy(n,0);
#pragma omp parallel 
{
  int tid=omp_get_thread_num();
  int nt=omp_get_max_threads();
  vector<DCO_T> y(n);
  for (int i=1;i<n-1;i++) {
    if((i-1)%nt!=tid) continue;
    for(int j=0;j<n;j++) y[j]=y_in[j];
    dco::derivative(y[i])=1;
    burgers(y);
    if (i==1) y_out=dco::value(y);
    dymdy[i]=dco::derivative(y[n/2]);
  }
}
  for (int i=1;i<n-1;i++) 
    cout << i << ": " << y_out[i] << " " << dymdy[i] << endl;

  return 0;
}
