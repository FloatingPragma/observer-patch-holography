# Local interacting spin dynamics with volume-independent error

On the same signed cubic source, a quasilocal unitary dressing of the
one-fermion-per-site inputs has a controlled interacting Heisenberg limit
at fixed lattice spacing. The error for a fixed local observable and time
is independent of the total volume and tends to zero when `U=c^2/J` with
fixed `J>0` and c tending to infinity. This construction uses the source
Hamiltonian alone. It does not identify the dressed state with an exact
source ground state or interleave the undressed repair instrument.

## Norm and integer penalty

In a finite open cubic box, retaining every occupation state, write

```
G = K+UV/2 = nu N + T,
nu=U/2,       N=sum_x (n_x-1)^2,       J=c^2/U.
```

Each onsite summand of N is a projection. Its zero space P is exactly
the tensor product of the two single-occupation spin states. The physical
particle number remains conserved. Use connected supports and even CAR
interactions with norm

```
||Phi||_k = sup_x sum_(X containing x) exp(k |X|) ||Phi_X||.
```

Terms on disjoint supports commute. The degree is at most six and each
two-spin hopping bond has norm at most c, so `||T||_1<=6 exp(2)c<48c`.
Set `B=144c`. All estimates below require

```
nu >= 2^20 B,             rho=4096 B/nu <= 1/256.       (1)
```

There is no bound on volume in (1). These conservative constants
establish a mathematical limit; they are not an executed resource budget.

For `k'>0`, `delta=k-k'>0`, counting overlapping supports gives

```
||[Phi,Psi]||_(k') <= (4/delta) ||Phi||_k ||Psi||_k.
```

Terms whose first support contains x contribute at most
`2||Psi||_k sum_(X containing x) |X| exp(k'|X|)||Phi_X||`;
the reversed assignment gives the other half. Use
`n exp(-delta n)<=1/delta`. Splitting the loss delta into m equal
parts and using `m! >= (m/e)^m` bounds the mth BCH term by
`q^m ||Z||_k`, where `q=12||S||_k/delta`. Its tail from order m is
at most `q^m/(1-q)||Z||_k` when `q<1`.

Average each local term separately under `exp(i theta N)`, and denote
this norm-contractive projection by A. Both A and `1-A` preserve support;
the latter has norm at most two. The inverse of `ad_N` on nonzero
integer frequencies is

```
I(V) = (1/(2pi)) integral_0^(2pi)
                  i(theta-pi) exp(i theta N) V exp(-i theta N) dtheta.
```

Its matrix entry of frequency m is `V_m/m`, its norm is at most
`pi/2<2`, and it is anti-Hermitian for Hermitian V. In particular
`[I(V),N]=-V` when A(V)=0. The operations preserve physical particle
number and even CAR parity.

## Five exact local transformations

Start with `D0=A(T)`, `V0=T-D0`. At step j set

```
Sj=I(Vj)/nu,
exp(Sj)(nu N+Dj+Vj)exp(-Sj)=nu N+Dj+Wj,
Wj=sum_(m>=1) ad_Sj^m(Dj)/m!
      +sum_(m>=1) m ad_Sj^m(Vj)/(m+1)!,
D_(j+1)=Dj+A(Wj),       V_(j+1)=(1-A)(Wj).             (2)
```

These are exact convergent local-interaction series, not a formal
truncation. The cancellation uses `[Sj,nu N]=-Vj`.
Take five steps, each with norm loss `delta=1/12`. Final decay is
`7/12>1/2`. Write `dj=||Dj||`, `vj=||Vj||` at its current decay.
Initially `d0+v0<=B`. As long as `dj+vj<=2B`,

```
qj<=288 vj/nu<=576 B/nu<1/2,
wj<=qj(dj+vj)/(1-qj),
v_(j+1)<=2wj<=rho vj.
```

The induction closes because
`dj+vj <= B+3 sum_(i<j) wi <= B+3B rho/[2(1-rho)]<2B`.
Consequently, for `Y=exp(S4)...exp(S0)`,

```
Y G Y* = nu N + D5 + V5,       [D5,N]=0,
||D5||_(1/2)<=2B,             ||V5||_(1/2)<=r:=B rho^5. (3)
```

Every diagonal interaction in D5 commutes with the penalty on its own
support. Its compression to P is a local spin interaction; no global
interacting spectral gap is used.

The first diagonal correction is `A([S0,V0])/2`:
`A([S0,D0])=0`, because D0 has zero penalty frequency and S0 has none.
A hop from P creates one hole and one doublon, with penalty frequency
two. Its P compression is therefore

```
P [S0,V0] P/2 = -P T^2 P/U
             = Hspin = J sum_edges(S_x.S_y-1/4).       (4)
```

First-step terms beyond these have norm at most `2(288 B/nu)^2 B`.
Subsequent diagonal increments sum to at most
`B rho^2/[2(1-rho)]<=B rho^2`. Compression is contractive support by
support. Thus

```
||P D5 P-Hspin||_(1/2) <= eta:=2^25 B^3/nu^2.          (5)
```

This bound does not use the additional bipartite fourth-order
cancellation in the finite spectral theorem.

## A local propagation estimate

Suppose an even connected interaction Phi has `||Phi||_(1/2)<=d`, and
a perturbation R has `||R||_(1/2)<=r`. For an observable A on s sites,
the commutator-chain expansion with
`M_xy=sum_(Z containing x,y)||Phi_Z||` has weighted row sum at most
`4d` at spatial decay `1/4`: connected Z has diameter at most
`|Z|-1`, and `|Z| exp(-|Z|/4)<=4`. Its matrix exponential bounds
commutators by `2||A||||R_Z||` times the double sum, over x in A and
y in Z, of `exp(8d|t|-dist(x,y)/4)`. Also retain the elementary cap
`2||A||||R_Z||`.

Choose `L=ceil(64d|t|)`. Terms intersecting the union of radius-L balls
about A contribute at most `16s||A||r(L+1)^3`. For the others, sum the
double-sum bound first over terms containing each site. In three
dimensions a distance-n sphere has at most `6(n+1)^2` sites, and

```
sum_(n>L) 6(n+1)^2 exp(-n/4) <= 2^14 exp(-L/8).
```

Use `(n+1)^2<=256 exp(n/8)` and `1/(1-exp(-1/8))<9`.
The far contribution is at most `2^15 s||A||r`, since
`exp(8d|t|-L/8)<=1`. Duhamel integration yields

```
||tau_(Phi+R)^t(A)-tau_Phi^t(A)||
 <= 2^16 s||A|| r |t| (1+64d|t|)^3.                  (6)
```

This also allows time-dependent R with the same norm. An onsite term
can be removed by an onsite interaction picture, preserving all supports
and norms. For (3), D5 and the spin observable commute with N. The large
onsite coefficient nu therefore does not enter d or the propagation speed.

## Uniform local dynamics on dressed inputs

Embed any spin state into P, and prepare the source state by Y*.
Compare the dressed observable `Y* A Y` with A. Equations (3), (5) and
(6), using `||Hspin||_(1/2)<=48J`, give the expectation error uniformly
over every spin input and untouched reference:

```
E_A(t) <= min(2||A||,
  2^16 s||A|| |t| [r(1+128B|t|)^3 + eta(1+3072J|t|)^3]). (7)
```

All virtual charge outputs remain present. This controls local
observables, not a volume-independent distance between entire states.
For a bare local source observable on the same dressed input, add
`8s||A||B/nu`. Each factor obeys
`||exp(Sj)A exp(-Sj)-A||<=2s||A||||Sj||_0`; telescope the five
isometric conjugations and use `sum ||Sj||_0<=2B/[nu(1-rho)]`.

At fixed J and `U=c^2/J`, `r=O(c^-4)`, `eta=O(c^-1)` and
`B/nu=O(c^-1)`. Thus (7) is `O(c^-1)` for each fixed local observable
and bounded time interval, without a total-volume factor. `local.py`
evaluates the constants rationally. At `c=2^256`, `U=c^2`, J=1 it
certifies small error at unit time for an eight-site observable;
this is a scale certificate, not a source execution.

The connected-support bounds also make finite-box conjugations and
dynamics Cauchy on fixed local observables as the boundary recedes.
They define quasilocal automorphisms and the interacting spin dynamics
on the infinite fixed-spacing lattice. Every spin state has the stated
dressed source realization; neighboring spin-bond terms do not commute.
This removes the extensive finite-box error on this state/readout domain.

Integer-spectrum local averaging is also treated in
[Abanin, De Roeck, Ho and Huveneers, Section 3](https://arxiv.org/html/1509.05386v3#S3),
and local effective Hamiltonians in
[Bravyi, DiVincenzo and Loss, Section 4](https://arxiv.org/html/1105.0675v1#S4).
The five-step constants, source compression and local estimate used here
are derived above.

## Boundary of the result

The lattice spacing is fixed. No sharp spin shell, Lorentz covariance,
surviving continuum connected interaction, physical particle dictionary
or Einstein law follows from (7). Its prepared dressed state is not
identified with the actual ground of G at finite c. The bound does not
cover arbitrary undressed inputs or interleaved native repair: bare repair
protects P, whereas Y*P contains virtual charge states. Those joins require
their own estimates on the continuing process.
