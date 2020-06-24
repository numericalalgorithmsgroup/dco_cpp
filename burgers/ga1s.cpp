#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "burgers.h"

#include "dco.hpp"
typedef dco::ga1s<double> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;

int main() {
  vector<DCO_T> y(n,0);
  for (int i=1;i<n-1;i++) y[i]=sin(2*PI*i/n);
  DCO_M::global_tape=DCO_TAPE_T::create();
  DCO_M::global_tape->register_variable(y);
  vector<DCO_T> yc(y);
  burgers(yc);  
  DCO_M::global_tape->register_output_variable(yc[25]);
  dco::derivative(yc[25])=1.;
  DCO_M::global_tape->interpret_adjoint();
  vector<double> v_dydy(dco::derivative(y));
  cerr << dco::size_of(DCO_M::global_tape) << "B" << endl;
  DCO_TAPE_T::remove(DCO_M::global_tape);
  for(int i=1;i<n-1;i++)
    cout << i << ": " << dco::value(yc[i]) 
              << " " << v_dydy[i] << endl;
  return 0;
}

