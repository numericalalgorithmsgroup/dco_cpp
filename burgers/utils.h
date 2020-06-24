#pragma once

template <typename T>
inline T norm(vector<T>& v) {
  int n=v.size();
  T r=0;
  for (int i=1;i<n-1;i++) r+=v[i]*v[i];
  return sqrt(r);
}

const double PI=3.141592653589793;
