#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_BT;
typedef dco::ga1s<DCO_BT> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;

int main() {
  vector<DCO_T> L(n,0.05); DCO_T P=0;
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  DCO_M::global_tape=DCO_TAPE_T::create();
  DCO_M::global_tape->register_variable(L);
  DCO_TAPE_T::position_t tpos=DCO_M::global_tape->get_position();
  vector<vector<double> > ddPdLL(n,vector<double>(n,0));
  for(int j=0;j<n;j++) {
    dco::derivative(dco::value(L[j]))=1;
    libor(L,P,Z);
    DCO_M::global_tape->register_output_variable(P);
    dco::value(dco::derivative(P))=1;
    DCO_M::global_tape->interpret_adjoint_to(tpos);
    for(int i=0;i<n;i++) {
      ddPdLL[i][j]=dco::derivative(dco::derivative(L[i]));
      dco::derivative(L[i])=0; 
    }
    dco::derivative(dco::value(L[j]))=0;
    DCO_M::global_tape->reset_to(tpos);
  }
  DCO_TAPE_T::remove(DCO_M::global_tape);
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      cout << "ddP/dL[" << i << "]dL[" << j << "]=" << ddPdLL[i][j] << endl;
  return 0;
}
