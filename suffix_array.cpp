#include <bits/stdc++.h>
#define FOR(i, n) for (int i = 0; i < (n); ++i)
using namespace std;

void counting_sort(int *a, int *b, int n, int *r, int dict) {
  int c[dict], sum = 0;
  memset(c, 0, dict * sizeof(int));
  FOR (i, n) ++c[r[a[i]]];
  FOR (i, dict) { int t = sum; sum += c[i]; c[i] = t; }
  FOR (i, n) b[c[r[a[i]]]++] = a[i];
}

void suffix_array(char *s, int n, int *sa) {
  int sa2[n], rank[2 * n], dict = 128;
  vector<bool> bits(n - 1);

  FOR (i, n) sa[i] = i, rank[i] = s[i], rank[i + n] = 0;
  for (int k = 1; k < n; k <<= 1) {
    counting_sort(sa, sa2, n, rank + k, dict);
    counting_sort(sa2, sa, n, rank    , dict);
    FOR (i, n - 1) bits[i] = rank[sa[i]] != rank[sa[i + 1]] ||
                             rank[sa[i] + k] != rank[sa[i + 1] + k];
    FOR (i, n) rank[sa[i]] = (i == 0) ? 1 : (rank[sa[i - 1]] + bits[i - 1]);
    dict = rank[sa[n - 1]] + 1;
    if (dict == n + 1) break;
  }
}

void permuted_lcp(char *s, int n, int *sa, int *plcp) {
  int inv_sa[n], l = 0;

  FOR (i, n) inv_sa[sa[i]] = i;
  FOR (i, n) {
    if (l > 0) --l;
    if (inv_sa[i] != 0)
      while (i + l < n && sa[inv_sa[i] - 1] + l < n &&
             s[i + l] == s[sa[inv_sa[i] - 1] + l]) ++l;
    plcp[inv_sa[i]] = l;
  }
}

int main() {
  char str[1024];
  scanf("%s", str);
  int n = strlen(str);

  int sa[n], plcp[n];
  suffix_array(str, n, sa);
  permuted_lcp(str, n, sa, plcp);

  FOR (i, n) printf("%2d | %2d | %s\n", sa[i], plcp[i], str + sa[i]);
  return 0;
}
