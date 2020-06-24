/*
  Gaussian elimination on row-major serialization
  of three-diagonal band matrix
  Uwe Naumann (Aug 2017)
*/

#ifndef _GAUSS_INCLUDED_
#define _GAUSS_INCLUDED_

#include <vector>
using namespace std;

// A:=L+U-I s.t. A=L*U
template <typename T>
inline void LU(vector<T>& A) {
  int n=(A.size()-4)/3+2;
  for (int k=1;k<n-1;k++) {
    A[2+k*3]/=A[k*3];
    A[(k+1)*3]-=A[2+k*3]*A[(k+1)*3-2];
  }
}

// y:=L^-1*b
template <typename T>
inline void FS(const vector<T>& LU, vector<T>& b) {
  int n=b.size();
  for (int i=1;i<n-1;i++)
    b[i]-=LU[i*3-1]*b[i-1];
}

// x:=U^-1*y
template <typename T>
inline void BS(const vector<T>& LU, vector<T>& y) {
  int n=y.size();
  for (int k=n-1;k>1;k--) {
    y[k-1]-=LU[(k-1)*3+1]*y[k];
    y[k-1]/=LU[(k-1)*3];
  }
}

#endif

