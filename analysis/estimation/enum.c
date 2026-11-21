#include <stdio.h>
#include <stdlib.h>
#include "inst.h"
static int L[N]; static long long hist[4*NF+2]; static long long Z=0;
static int nb_cnt[N], nb_j[N][N], nb_s[N][N];
static int alt(const int *l){ /* positively alternating: distinct magnitudes, sorted by magnitude signs + - + ... */
  int a[M]; for(int i=0;i<M;i++)a[i]=l[i];
  for(int i=0;i<M;i++)for(int j=i+1;j<M;j++) if(abs(a[j])<abs(a[i])){int t=a[i];a[i]=a[j];a[j]=t;}
  for(int i=0;i+1<M;i++) if(abs(a[i])==abs(a[i+1])) return 0;
  for(int i=0;i<M;i++) if((a[i]>0)!=(i%2==0)) return 0;
  return 1; }
static void leaf(void){ int A=0; int l[M];
  for(int f=0;f<NF;f++){ for(int i=0;i<M;i++) l[i]=FS[f][i]*L[FI[f][i]]; A+=alt(l);
                         for(int i=0;i<M;i++) l[i]=-l[i]; A+=alt(l); }
  Z++; hist[A]++; }
static void rec(int v){ if(v==N){leaf();return;}
  for(int q=0;q<Q;q++){ int a=LAB[q], ok=1;
    for(int t=0;t<nb_cnt[v];t++){ if(a==-nb_s[v][t]*L[nb_j[v][t]]){ok=0;break;} }
    if(ok){ L[v]=a; rec(v+1); } } }
int main(void){
  for(int e=0;e<NE;e++){ int a=EA[e],b=EB[e]; /* b>a: check when placing b against a */
    nb_j[b][nb_cnt[b]]=a; nb_s[b][nb_cnt[b]]=ES[e]; nb_cnt[b]++; }
  rec(0);
  printf("Z %lld\n",Z); for(int A=0;A<4*NF+2;A++) if(hist[A]) printf("A %d %lld\n",A,hist[A]);
  return 0; }
