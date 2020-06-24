#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::ga1s<double> DCO_M;
typedef typename DCO_M::type DCO_T;
typedef typename DCO_M::tape_t DCO_TAPE_T;

void libor(const vector<DCO_T>& L, DCO_T& P, const vector<vector<double>> &Z) {
  DCO_T Ps=0;
  DCO_M::jacobian_preaccumulator_t jp(dco::tape(L));
  for (int j=0;j<p;j++) {
    jp.start();
    vector<DCO_T> Lc(L);
    path_calc(j,Lc,Z);
    portfolio(Lc,P);
    jp.register_output(P);
    jp.finish();
    Ps+=P;
  }
  P=Ps/p;
}

int main() {
  vector<DCO_T> L(n,0.05); DCO_T P=0;
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  vector<double> dPdL(n,0); 
  DCO_M::global_tape=DCO_TAPE_T::create();
  DCO_M::global_tape->register_variable(L);
  libor(L,P,Z);
  DCO_M::global_tape->register_output_variable(P);
  dco::derivative(P)=1;
  DCO_M::global_tape->interpret_adjoint();
  for(int i=0;i<n;i++) dPdL[i]=dco::derivative(L[i]);
  cerr << dco::size_of(DCO_M::global_tape) << "b" << endl;
  DCO_TAPE_T::remove(DCO_M::global_tape);
  cout << "P=" << dco::value(P) << endl;
  for(int i=0;i<n;i++)
    cout << "dPdL[" << i << "]=" << dPdL[i] << endl;
  return 0;
}
