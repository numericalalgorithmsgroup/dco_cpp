#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::ga1sm<double> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;

#include <omp.h>

void libor(vector<double>& L, double& P, const vector<vector<double>> &Z, vector<double>& dPdL) {
  int nt=omp_get_max_threads();
  P=0;
  #pragma omp parallel reduction(+:P)
  {
    int tid=omp_get_thread_num();
    DCO_TAPE_T *tape=DCO_TAPE_T::create();
    vector<DCO_T> L_t(n,0); for (int i=0;i<n;i++) L_t[i]=L[i];
    tape->register_variable(L_t);
    DCO_TAPE_T::position_t tpos=tape->get_position();
    DCO_T P_t=0; double Ps_t=0;
    for (int j=0;j<p;j++) {
      if(j%nt!=tid) continue;
      vector<DCO_T> Lc_t(L_t);
      path_calc(j,Lc_t,Z);
      portfolio(Lc_t,P_t);
      tape->register_output_variable(P_t);
      dco::derivative(P_t)=1./p;
      tape->interpret_adjoint();
      Ps_t+=dco::value(P_t);
      tape->reset_to(tpos);
    }
    Ps_t/=p;
    P+=Ps_t;
    for (int i=0;i<n;i++) {
      #pragma omp atomic
      dPdL[i]+=dco::derivative(L_t[i]);
    }
    cerr << dco::size_of(tape) << "b" << endl;
    DCO_TAPE_T::remove(tape);
  }
}

int main() {
  omp_set_num_threads(4);
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
