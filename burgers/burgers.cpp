#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "burgers.h"

int main() {
  vector<double> y(n); 
  for (int i=1;i<n-1;i++) y[i]=sin(2*PI*i/n);
  burgers(y);
  for(int i=1;i<n-1;i++) cout << y[i] << endl;
  return 0;
}
