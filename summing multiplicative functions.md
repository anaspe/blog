# Summing Multiplicative Functions
### Dirichlet's Convolution
we define $$(g*f)(n) = \sum_{d|n} g(d)f(\frac{n}{d})$$
we'll show in these notes that we can compute the sum of any dirichlet's convolution in $O(\sqrt{n})$ and some other cool stuff! 
### Dirichlet's Hyperbola method
notice that we can write a dirichlet convolution as: $$ \sum_{d|n} f(d)g(\frac{n}{d}) = \sum_{ab=n} f(i)g(j)$$
visually this means we'll iterate over all the lattice points that intersect the hyperbola $xy=n$, think about it for a moment and it'll make sense!
we also notice: (very cute abuse of notation!)$$\sum_{i=1}^n \sum_{ab=i} f(a)g(b) = \sum_{ab \le n} f(a)b(b)$$
or in other words, we iterate over all the lattice points under the hyperbola $xy=n$, i'll prolly elaborate more on this in the lecture itself so yes, the idea of *Dirichlet's Hyperbola method* is, choosing two points $\alpha$ and $\beta$ such that $\alpha \beta = n$, we can just compute the number of lattice point horizontally $\le \alpha$ and vertically while $\le \beta$ and remove the intersection we can compute the sum! so: $$ \sum_{ab \le n} f(a)b(b)$$
$$ \sum_{a \le \alpha} \Bigg(f(a) \cdot \sum_{b \le x/a} g(b) \Bigg) + \sum_{b \le \beta} \Bigg(g(a) \cdot \sum_{a \le x/b} f(a) \Bigg) - \Bigg(\sum_{a \le \alpha} f(a) \Bigg) \cdot \Bigg( \sum_{b \le \beta} g(b) \Bigg)$$
define $F(x) = \sum_{n \le x} f(n)$ and $G(x) = \sum_{n \le x} g(n)$, we write: $$ \sum_{a \le \alpha} f(a)G(x/a) + \sum_{b \le \beta} F(x/b)g(b) - F(\alpha)G(\beta)$$it's almost always fine to choose $\alpha = \sqrt{x}$ and $\beta = \sqrt{x}$, we can change them depends on how hard it is to compute sums and the functions themselves, so in conclusion, we can compute the sum of a dirichlet's convolution of two multiplicative functions f and g in $O(\sqrt{n})$ if we can compute $F(n), G(n), f(n), g(n)$ in $O(1)$, to elaborate on why *Dirichlet's Hyperbola method* works, i probably elaborated in the lecture so go there!
- exercices:
https://projecteuler.net/problem=401
https://projecteuler.net/problem=625 (requires a little bit of cleverness! read all the notes before attempting)
### ~~Flat~~ Almost-Linear Sieving
ever heard of sieve of eroitosneneess? ~~idk how to spell it~~ we'll we can use the idea of it and compute almost all multiplicative functions, here's one for phi:
```cpp
// very cute way to compute phi(n) for all x <= n in O(xloglogx)
vector<int> phi_sieve(int x) {
  vector<int> result(x+1);
  for(int i = 1 ; i <= x ; i++) result[i]=i;
  for(int p = 2; p <= x; p++) {
    if(result[p] == p) {
      for(int k = 1; k <= x/p ; k++) {
        result[p*k] = (result[p*k]/p)*(p-1);
      }
    }
  }

  return result;
}
```
and here is a possibly bad way to sieve mobius (it's bad cause i wrote it, although it seems to work surprisingly): 
```cpp
// non-linear mobius sieve since idc bout' linear sieves and this is O(nloglogn)
vector<int> mobius_sieve(int x) {
  vector<int> result(x+1);
  vector<bool> primes(x+1, 1);
  result[1] = 1;
  for(int p = 2 ; p <= x ; p++) {
    if(primes[p]) {
      result[p] = -1;
      for(int i = 2; i <= x/p ; i++) {
        primes[p*i] = 0;
        result[p*i] = -result[i];
        if(p*p <= p*i and p*i % p*p == 0) result[p*i] = 0;
      }
    }
  }

  return result;
}
```
and there is linear sieving, which sucks and i don't care about it since $O(nloglogn)$ is already linear to me and way easier to work with
### Summing $\mu$ and $\varphi$
let's start with summing $\mu$ ~~since we'll need it later since phi = (mu * 1)~~
since $(\mu * u)(n)$ = 1 for $n \neq 1$ and else 0, so for $v \ge 1$ we write: $$\sum_{n \leq \sqrt{v}} \mu(n)\left \lfloor \frac{v}{n}\right\rfloor + \sum_{n \leq \sqrt{v}} M\left(\frac{v}{n}\right) - \lfloor \sqrt{v} \rfloor M\left(\sqrt{v}\right) = 1$$
We can suppose we’ve sieved at least the first $\sqrt{x}$ values of $\mu$, and that we compute all $M(x/n)$ for $n > 1$, we write: $$M(x) = 1 - \sum_{n \leq \sqrt{x}} \mu(n)\left \lfloor \frac{x}{n}\right\rfloor - \sum_{2 \leq n \leq \sqrt{x}} M\left(\frac{x}{n}\right) + \lfloor \sqrt{x} \rfloor M\left(\sqrt{x}\right)$$
since $M(x/1) = M(x)$ in the second sum (easy to see for y'all math mature ppl!)
so we can compute $M(x)$ in $O(\sqrt{x})$ (assuming we have $M(x/2), M(x/3)...M(1)$ and already sieved $\mu$ up to $\sqrt{x}$), so we:
	- sieve $\mu$ up to $\sqrt{n}$
	- for each v in increasing order in an indexing of $1, 2, \cdots, \sqrt{x} \ div \ n \cdots$, set: $$M(v) = 1 - \sum_{n \leq \sqrt{v}} \mu(n)\left\lfloor \frac{v}{n}\right\rfloor - \sum_{2 \leq n \leq \sqrt{v}} M\left(\frac{v}{n}\right) + \lfloor \sqrt{v} \rfloor M\left(\sqrt{v}\right)$$
analyzing this we get $$O\left(\sum_{v \leq \sqrt{x}} \sqrt{v}\right) = O\left(\sqrt{x} \sqrt{\sqrt{x}}\right) = O\left(x^{3/4}\right)$$
so sieving with $O(\sqrt{n})$ memory and $O(n^{3/4})$ pretty neat if you ask me! we can optimize this to $O(x^{2/3})$ but it's not worth the headache for me so we'll stop here, if you are struggling to get it, i was too! i'll code it later, ask me on discord or something
#### summing $\varphi$
yes, $\varphi$ = $\mu * N$ and yes. one benefit of the method above is that it's dp-ish so the complexity of this one would be $O(x^{2/3})$ *can you see why?*, if you coded the last one this should be trivial

## Working Examples:
#### [Project Euler 401](https://projecteuler.net/problem=401)
basically, we have to compute $$ \sum_{i=1}^N \sum_{d|i} d^2$$
we can brainless-ly apply DHM here so: $$ \sum_{n\le x} (d*u)(n) $$ where $d(n)=n^2$ and $u(n)=1$, we know $D(n) = \frac{n(n+1)(2n + 1)}{6}$ and $U(n) = n$, so:
$$ \bigg\{\sum_{a\le \sqrt{x}} G\Bigl(\Big\lfloor \frac{x}{a} \Big\rfloor\Bigl) \bigg\} + \bigg\{\sum_{b\le \sqrt{x}} b^2\cdot \Big\lfloor \frac{x}{b} \Big\rfloor \bigg\} - \lfloor \sqrt{x} \rfloor \cdot G(\lfloor \sqrt{x} \rfloor) $$
and it does halt and return the right result in ~5 seconds for $n = 10^{15}$ using pypy, code:
```python
import math

def G(n):
    return n*(n+1)*(2*n + 1)/6

def summ(n):
    z = int(math.sqrt(n))
    p = 0
    for i in range(1, z+1):
        p += G(int(n/i))
        p += i**2 * int(n/i)
    return p - z*G(z)

print(summ(10**15) % 10**9)
```
#### [Project Euler 211](https://projecteuler.net/problem=211)
we can just sieve! really trivial, runs in 20s
```cpp
bool issquare(int x) {
  return (int)sqrt(x)*(int)sqrt(x) == x;
}

signed main() {
  int t = 64'000'000, r = 0;
  vector<int> vec(t+1);
  for(int i = 1; i <= t ; i++) {
    for(int j = 1; j <= t/i ; j++) {
      vec[i*j] += j*j;
    }

    if(issquare(vec[i])) r += i;
  }

  cout << r;
}
```
~~couldn't find an example that sieves mobius/phi~~
#### [Project Euler 625](https://projecteuler.net/problem=625)
notice that the number between $1\cdots n$ such that $\gcd(i, n) = a$ is $\varphi(\frac{n}{a})$ iff $a|n$ (i need to read more number theory proofs!) so we can write: $$ \sum_{i=1}^j \gcd(i, j) = \sum_{d|j} \varphi(\frac{j}{d}) \cdot d =(\varphi * N)(i)$$
where $N(n) = n$, writing this with the hyperbola method where $P(n) = \sum_{i=1}^n \varphi(i)$ and $G(n) = \sum_{i=1}^n i = \frac{n(n+1)}{2}$:
$$ \sum_{i=1}^n (\varphi*N)(i) = \sum_{a \le \alpha} \varphi(a)G(n/a) + \sum_{b \le \beta} bP(n/b) - G(\alpha)P(\beta)$$
we can compute $P(n)$ in $O(n^{2/3})$, and sieve $\varphi$ t'ill $\alpha$, one can derive that it'll take $O(n^{2/3} \beta^{5/2})$ to compute the $\beta$ sum, and since im (surprisingly) able to sieve $\varphi$ up to $2 \cdot 10^8$ using a non-linear sieve, i'll set $\alpha = 2\cdot 10^8$ and $\beta = 500$, the code is left as an exercise ~~since im too lazy to code it, although this should work, or at least barely pass the minute rule~~