+++
title = "uwu"
date = 2024-12-27
+++

Hello and welcome to this personal blog of mine! don't expect much except hopefully weekly brain farts, since that's all i can produce for now :3

im grinding for my national OI so i guess that's something to lock in for and write about, fingers crossed i don't waste this year too :pray:

my country does this silly thing where they make us take random contests (mostly from AtCoder and Codeforces) each Saturday and Sunday, tomorrow we have one for AtCoder and im looking for a good performance since i messed up the last Codeforces one, hopefully i'll make a writeup for this one

ehhhh yeah that's it, hopefully i don't give up on blogging :sob:

testing some features of zola (the SSG im using):

#### Code:

```c++
void apply(int p, int value) {
  t[p] += value;
  if (p < n) d[p] += value;
}

void build(int p) {
  while (p > 1) p >>= 1, t[p] = max(t[p<<1], t[p<<1|1]) + d[p];
}

void push(int p) {
  for (int s = h; s > 0; --s) {
    int i = p >> s;
    if (d[i] != 0) {
      apply(i<<1, d[i]);
      apply(i<<1|1, d[i]);
      d[i] = 0;
    }
  }
}

void inc(int l, int r, int value) {
  l += n, r += n;
  int l0 = l, r0 = r;
  for (; l < r; l >>= 1, r >>= 1) {
    if (l&1) apply(l++, value);
    if (r&1) apply(--r, value);
  }
  build(l0);
  build(r0 - 1);
}

int query(int l, int r) {
  l += n, r += n;
  push(l);
  push(r - 1);
  int res = -2e9;
  for (; l < r; l >>= 1, r >>= 1) {
    if (l&1) res = max(res, t[l++]);
    if (r&1) res = max(t[--r], res);
  }
  return res;
}
```

#### Math:
$$\int_{i=0}^n x^i = x^{i+1}$$
