#include <bits/stdc++.h>
#define FOR(i, n) for (int i = 0; i < (n); ++i)
#define ROF(i, n) for (int i = (n) - 1; i >= 0; --i)
#define REP(i, n) for (int i = 1; i <= (n); ++i)
using namespace std;

const int MAXN = 100001;
char S[MAXN];

int N, SA[MAXN], C[MAXN], *RA, *RA2, RA_tmp[2][MAXN], k;

bool cmp(const int& i, const int& j) {
  if (RA[i] != RA[j]) return RA[i] < RA[j];
  return (i + k < N && j + k < N) ? RA[i + k] < RA[j + k] : i > j;
}

void suffix_array() {
  RA = RA_tmp[0]; RA2 = RA_tmp[1];
  iota(SA, SA + N, 0);
  FOR (i, N) RA[i] = S[i];
  for (k = 1; k < N; k <<= 1) {
    sort(SA, SA + N, cmp);
    RA2[SA[0]] = 1;
    for (int i = 1; i < N; ++i) RA2[SA[i]] = RA2[SA[i - 1]] + cmp(SA[i - 1], SA[i]);
    swap(RA, RA2);
    if (RA[SA[N - 1]] == N) break;
  }
}

int main() {
  scanf("%s", S);
  N = strlen(S);
  suffix_array();

  FOR (i, N) printf("%d\n", SA[i]);
  return 0;
}
