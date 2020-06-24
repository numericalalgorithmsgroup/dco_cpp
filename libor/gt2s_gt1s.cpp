#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

#include "dco.hpp"
typedef dco::gt1s<double>::type DCO_BT;
typedef dco::gt1s<DCO_BT>::type DCO_T;

int main() {
  vector<DCO_T> L(n,0.05); DCO_T P=0;
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  vector<vector<double> > ddPdLL(n,vector<double>(n,0));
  for (int i=0;i<n;i++) {
    dco::value(dco::derivative(L[i]))=1;
    for (int j=0;j<=i;j++) {
      dco::derivative(dco::value(L[j]))=1;
      libor(L,P,Z);
      ddPdLL[i][j]=ddPdLL[j][i]=dco::derivative(dco::derivative(P));
      dco::derivative(dco::value(L[j]))=0;
    }
    dco::value(dco::derivative(L[i]))=0;
  }
  for (int i=0;i<n;i++) 
    for (int j=0;j<n;j++) 
      cout << "ddP/dL[" << i << "]dL[" << j << "]=" << ddPdLL[i][j] << endl;
  return 0;
}
