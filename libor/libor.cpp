#include <vector>
#include <iostream>
using namespace std;

#include "libor.h"

int main() {
  vector<double> L(n,0.05); double P=0;
  srand(0); default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<vector<double>> Z(p,vector<double>(m));
  for (int j=0; j<p;j++)
    for (int i=0;i<m;i++)
      Z[j][i]=0.3+distribution(generator);
  libor(L,P,Z);
  cout << "P=" << P << endl;
  return 0;
}
