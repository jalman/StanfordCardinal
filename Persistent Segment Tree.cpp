struct node {
  node *l, *r;
  int sum;
} *tree[MAXN], pool[MAXN * MAXLGN * MAXLGN];

node *allocate(node * const l, node * const r, const int &sum) {
  static int pool_n = 0;
  pool[pool_n].l = l;
  pool[pool_n].r = r;
  pool[pool_n].sum = sum;
  return pool + (pool_n++);
}

node *update(const node * const old, const int &l, const int &r, const int &i, const int &k) {
  if (old != NULL && old->sum + k == 0) return NULL; // assumes there are no negative elements
  node *x;
  if (old == NULL) x = allocate(NULL, NULL, k);
  else x = allocate(old->l, old->r, old->sum + k);
  if (l < r) {
    int m = (l + r) / 2;
    if (i <= m) x->l = update(x->l, l, m, i, k);
    if (i >  m) x->r = update(x->r, m + 1, r, i, k);
  }
  return x;
}

int query(const node * const x, const int &l, const int &r, const int& i) {
  if (x == NULL) return 0;
  if (l < r) {
    int m = (l + r) / 2;
    if (i <= m) return query(x->l, l, m, i) + (x->r ? x->r->sum : 0);
    else return query(x->r, m + 1, r, i);
  } else return x->sum;
}
