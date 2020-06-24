#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "dco.hpp"
typedef double DCO_BT;
typedef dco::ga1s<DCO_BT> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;
typedef DCO_M::external_adjoint_object_t DCO_EAO_T;

void newton(const vector<DCO_T>& yp, vector<DCO_T>& y);
#include "burgers.h"

#include "burgers_symbolic_nls.h"

int main() {
  vector<DCO_T> y(n,0);
  for (int i=1;i<n-1;i++) y[i]=sin(2*PI*i/n);
  vector<double> v_dydy(n,0);
  DCO_M::global_tape=DCO_TAPE_T::create();
  DCO_M::global_tape->register_variable(y);
  vector<DCO_T> yc(y);
  burgers(yc);  
  DCO_M::global_tape->register_output_variable(yc[25]);
  dco::derivative(yc[25])=1.;
  DCO_M::global_tape->interpret_adjoint();
  for(int i=1;i<n-1;i++) v_dydy[i]=dco::derivative(y[i]);
  cerr << dco::size_of(DCO_M::global_tape,DCO_TAPE_T::size_of_stack|DCO_TAPE_T::size_of_internal_adjoint_vector|DCO_TAPE_T::size_of_checkpoints) << "b" << endl;
  DCO_TAPE_T::remove(DCO_M::global_tape);
  for(int i=1;i<n-1;i++)
    cout << i << ": " << dco::value(yc[i]) 
              << " " << v_dydy[i] << endl;
  return 0;
}

