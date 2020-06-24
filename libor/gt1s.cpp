#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_T;

int main() {
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  vector<double> dPdL(n,0); 
  DCO_T P=0;
  for (int i=0;i<n;i++) {
    vector<DCO_T> L(n,0.05); 
    dco::derivative(L[i])=1;
    libor(L,P,Z);
    dPdL[i]=dco::derivative(P);
  }
  cout << "P=" << dco::value(P) << endl;
  for (int i=0;i<n;i++)
    cout << "dPdL[" << i << "]=" << dPdL[i] << endl;
  return 0;
}
