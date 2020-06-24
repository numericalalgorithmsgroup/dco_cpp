#include <cmath>
#include <vector>
#include <random>
using namespace std;

const int p=10000;
const int m=40;
const int n=m+40;
const int no=15;

const double delta=0.25;
const vector<int> maturities({4,4,4,8,8,8,20,20,20,28,28,28,40,40,40});
const vector<double> swaprates({.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055 });
const vector<double> sigma(n,0.2);

template <typename T>
inline void path_calc(
    const int path,
    vector<T>& L,
    const vector<vector<double>>& Z
) {
  for(int j=0;j<m;j++) {
    double aux1=sqrt(delta)*Z[path][j];
    T S=0.0;
    for (int i=j+1;i<n;i++) {
      double aux2=delta*sigma[i-j-1];
      S+=(aux2*L[i])/(1.0+delta*L[i]);
      L[i]=L[i]*exp(aux2*S+sigma[i-j-1]*(aux1-0.5*aux2));
    }
  }
}

template <typename T>
inline void portfolio( const vector<T>& L, T& P ) {
  vector<T> B(n),S(n);
  T b=1.0;
  T s=0.0;
  for (int j=m;j<n;j++) {
    b=b/(1.0+delta*L[j]); B[j]=b;
    s=s+delta*b; S[j]=s;
  }
  P=0;
  for (int i=0;i<no;i++){
    int j=maturities[i]+m-1;
    T swapval=B[j]+swaprates[i]*S[j]-1.0;
    if (swapval<0) P+=-100.0*swapval;
  }
  for (int i=0;i<m;i++) P=P/(1.0+delta*L[i]);
}

template<typename T>
inline void libor(const vector<T>& L, T& P, const vector<vector<double>>& Z) {
  T Ps=0;
  for (int j=0;j<p;j++) {
    vector<T> Lc(L);
    path_calc(j,Lc,Z);
    portfolio(Lc,P);
    Ps+=P;
  }
  P=Ps/p;
}
