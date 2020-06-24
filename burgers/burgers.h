#pragma once

#include <vector>
#include <cmath>
using namespace std;

// making passive read-only variables global for 
// more compact code listings
const int n=50, m=100, bw=32;
const double v=1e-2;

#include "utils.h"
#include "gauss.h"

//*** right-hand side of ODE dy/dt = v d^2 y / d x^2 - y dy / dx
//***                                `---,---------´   `---,---´
//***                                    diffusion         advection
template <typename T>
inline void g(const vector<T>& y, vector<T>& r) {
  int n=y.size();
  for (int i=1;i<n-1;i++) {
    //*** central finite difference scheme for diffusion term on
    //*** equidistant 1D grid
    T diffusion = v*(n-1)*(n-1)*(y[i-1]-2*y[i]+y[i+1]);
    //*** first-order upwind scheme for advection term on equidistant
    //*** 1D grid
    T advection = -y[i]*(n-1);
    if (advection < 0) { advection *= y[i]-y[i-1]; }
    else               { advection *= y[i+1]-y[i]; }
    r[i] = diffusion + advection;
  }
}

// tangent of g (hand-written)
template <typename T>
inline void g_t(const vector<T>& y, const vector<T>& y_t, vector<T>& r_t) {
  int n=y.size();
  for (int i=1;i<n-1;i++)  {
    T diffusion_t = v*(n-1)*(n-1)*(y_t[i-1]-2*y_t[i]+y_t[i+1]);
    T advection = -y[i]*(n-1);
    T advection_t = -y_t[i]*(n-1);
    if (advection < 0) { advection_t = advection * (y_t[i]-y_t[i-1]) + advection_t * (y[i]-y[i-1]); }
    else               { advection_t = advection * (y_t[i+1]-y_t[i]) + advection_t * (y[i+1]-y[i]); }
    r_t[i] = diffusion_t + advection_t;
  }
}

//*** Jacobian of r with respect to y (i.e. Jacobian of implementation g)
//***  - for the computation, the tangent of g (i.e. g_t) is used
//***  - A is known to be a tridiagonal matrix and uses custom data format, it stores row-wise
//***       a b
//***       c d e
//***         f g h
//***           i j k
//***             l m
//***    where in our case (a, b, c, k, l, m) = 0
//***  - Curtis/Powell/Reid(CPR)-algorithm is used to compute compressed Jacobian
//***  - if transpose=true, the transposed Jacobian is returned in A using the same data format.
template <typename T>
inline void dgdy(const vector<T>& y, vector<T>& A,
                 bool transpose=false) {
  int n=y.size();
  for (int i=0;i<3;i++) {
    vector<T> y_t(n,0),r_t(n);
    //*** CPR seeding
    for (int j=i+1;j<n;j+=3) y_t[j]=1;
    
    //*** computes compressed column in r_t
    g_t(y,y_t,r_t);

    for (int j=i+1;j<n-1;j+=3) {
      A[j*3]=r_t[j]; // diagonal elements (e.g. g)
      if (transpose) {
        //*** store column r_t in row
        A[j*3-1]=r_t[j-1]; // left of diagonal (i.e. f)
        A[j*3+1]=r_t[j+1]; // right of diagonal (i.e. h)
      } else {
        A[(j-1)*3+1]=r_t[j-1]; // above diagonal (i.e. e)
        A[(j+1)*3-1]=r_t[j+1]; // below diagonal (i.e. i)
      }
    }
  }
}

//*** residual of nonlinear system
//*** - yp is solution of previous timestep
template <typename T>
inline void f(const vector<T>& y, const vector<T>& yp, vector<T>& r) {
  int n=y.size();
  g(y,r);
  r[0]=y[0];
  for (int i=1;i<n-1;i++) r[i]=y[i]-yp[i]-r[i]/m;
  r[n-1]=y[n-1];
}

//*** Jacobian of residual of nls wrt. state (i.e. y)
//*** - see description of dgdy
template <typename T>
inline void dfdy(const vector<T>& y,
                 vector<T>& A, bool transpose=false) {
  int n=y.size();
  dgdy(y,A,transpose);
  for (auto& e : A) e/=-m;
  for (int i=1;i<n-1;i++) A[i*3]+=1; 
}

//*** Newton solver for nls
template <typename T>
inline void newton(const vector<T>& yp, vector<T>& y) {
  int n=y.size();
  const double eps=1e-14;
  //*** see comment of dgdy for description of data structure A
  vector<T> A((n-2)*3+4,0), r(n,0);
  f(y,yp,r); 
  while (norm(r)>eps) {
    dfdy(y,A);
    LU(A); FS(A,r); BS(A,r);
    for (int i=1;i<n-1;i++) y[i]-=r[i]; 
    f(y,yp,r); 
  }
}

//*** implicit Euler integration
template <typename T>
void burgers(vector<T>& y) {
  for (int j=0;j<m;j++) {
    vector<T> yp=y;
    newton(yp,y); 
  }
}
