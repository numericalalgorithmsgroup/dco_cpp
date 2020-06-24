#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "burgers.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_BT;
typedef dco::ga1s<DCO_BT> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;
typedef DCO_M::external_adjoint_object_t DCO_EAO_T;

#include "burgers_symbolic_nls.h"

int main() {
  vector<DCO_T> y(n,0);
  for (int i=1;i<n-1;i++) y[i]=sin(2*PI*i/n);
  DCO_M::global_tape=DCO_TAPE_T::create();
  DCO_M::global_tape->register_variable(y);
  dco::derivative(dco::value(y[24]))=1;
  vector<DCO_T> yc(y);
  burgers(yc);  
  DCO_M::global_tape->register_output_variable(yc[25]);
  dco::value(dco::derivative(yc[25]))=1.;
  DCO_M::global_tape->interpret_adjoint();
  vector<double> v_ddydyy_v(n,0);
  for(int i=1;i<n-1;i++) 
    v_ddydyy_v[i]=dco::derivative(dco::derivative(y[i]));
  cerr << dco::size_of(DCO_M::global_tape) << "b" << endl;
  DCO_TAPE_T::remove(DCO_M::global_tape);
  for(int i=1;i<n-1;i++)
    cout << i << ": " << dco::passive_value(yc[i]) << " "
         << v_ddydyy_v[i] << endl;
  return 0;
}

