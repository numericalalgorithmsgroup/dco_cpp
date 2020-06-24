#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::ga1s<double> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;

void libor(vector<double>& L, double& P, const vector<vector<double>>& Z, vector<double>& dPdL) {
  P=0;
  DCO_M::global_tape=DCO_TAPE_T::create();
  vector<DCO_T> L_t(n,0); for (int i=0;i<n;i++) L_t[i]=L[i];
  DCO_M::global_tape->register_variable(L_t);
  DCO_TAPE_T::position_t tpos=DCO_M::global_tape->get_position();
  DCO_T P_t=0; 
  for (int j=0;j<p;j++) {
    vector<DCO_T> Lc_t(L_t);
    path_calc(j,Lc_t,Z);
    portfolio(Lc_t,P_t);
    DCO_M::global_tape->register_output_variable(P_t);
    dco::derivative(P_t)=1./p;
    DCO_M::global_tape->interpret_adjoint();
    P+=dco::value(P_t);
    DCO_M::global_tape->reset_to(tpos);
  }
  P/=p;
  for (int i=0;i<n;i++) dPdL[i]+=dco::derivative(L_t[i]);
  cerr << dco::size_of(DCO_M::global_tape) << "b" << endl;
  DCO_TAPE_T::remove(DCO_M::global_tape);
}

int main() {
  vector<double> L(n,0.05); double P=0;
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  vector<double> dPdL(n,0); 
  libor(L,P,Z,dPdL);
  cout << "P=" << dco::value(P) << endl;
  for(int i=0;i<n;i++) 
    cout << "dPdL[" << i << "]=" << dPdL[i] << endl;
  return 0;
}
