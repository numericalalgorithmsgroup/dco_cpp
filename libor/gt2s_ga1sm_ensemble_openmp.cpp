#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_BT;
typedef dco::ga1sm<DCO_BT> DCO_M;
typedef DCO_M::type DCO_T;
typedef DCO_M::tape_t DCO_TAPE_T;

#include <omp.h>

int main() {
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  vector<vector<double> > ddPdLL(n,vector<double>(n,0));
  omp_set_num_threads(4); int nt=omp_get_max_threads();
  #pragma omp parallel
  {
    vector<DCO_T> L(n,0.05); DCO_T P=0;
    DCO_TAPE_T *tape=DCO_TAPE_T::create();
    tape->register_variable(L);
    DCO_TAPE_T::position_t tpos=tape->get_position();
    int tid=omp_get_thread_num();
    for(int j=0;j<n;j++) {
      if (j%nt!=tid) continue;
      dco::derivative(dco::value(L[j]))=1;
      libor(L,P,Z);
      tape->register_output_variable(P);
      dco::derivative(P)=1;
      tape->interpret_adjoint_to(tpos);
      for(int i=0;i<n;i++) {
        ddPdLL[i][j]=dco::derivative(dco::derivative(L[i]));
        dco::derivative(L[i])=0; 
      }
      dco::derivative(dco::value(L[j]))=0;
      tape->reset_to(tpos);
    }
    DCO_TAPE_T::remove(tape);
  }
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      cout << "ddP/dL[" << i << "]dL[" << j << "]=" << ddPdLL[i][j] << endl;
  return 0;
}
