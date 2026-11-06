# Complete memory correlations and reduced dynamics

The specified family feedback service applies `U_i(h)=exp(-ih H_i)` after a
native event sends memory to address `i`, with `H_i=g X tensor P_i` and
`g>=0`. The `P_i` are the source's real rank-one family projectors. A disabled
pulse is identity and its record remains present. The pulse step and time are
dimensionless model coordinates. Neither the auxiliary two-component register
nor a transfer frequency is identified with physical chirality or a mass.

The result below sums every memory lag. It gives a completely positive
corrected reduced semigroup with a uniform finite-time error, including an
arbitrary untouched reference. Its preparation is a stationary memory state
tensored with the field/reference state, and its output discards memory and
its record. It does not identify complete output tapes with a semigroup.

## Source transitions and exact reduction

For independent flags `(C,E)` the probabilities are

```
(1-p)(1-q), (1-p)q, p(1-q), pq,      0<=p<1, 0<=q<=1.
```

The native branch has `p=1/12,q=1/6`. The separate current-gated source has
`p=sin(theta sqrt(cL))^2/12` and `q=sin(theta sqrt(cR))^2/6`; its supplied
Casimirs and interaction angles are additional premises. The symbolic result
applies throughout the stated probability domain. Native raw labels are
summed exactly without conditioning on enabled pulses.

Let `T_im` be the full transition and `E_im` its enabled part. Every column of
`T` sums to one; every column of `E` sums to the duty `d`, equal to `1,p,q,pq`
for always, C, E and CE policies. The native reset law gives

```
T=(1-p)R_reset+pS,
```

with stochastic `S`. Hence every zero-sum vector of trace-class memory blocks
contracts by at most `p` in the sum of trace norms, including a reference.
For `p<1` the stationary state `pi` is unique.

It is supported on `{2,3,9,10}`, and `P9=P10`. The stationary unread experiment
therefore admits exact lumping onto the three projector rays `{2,3,10}`.
This reduction does not retain the distinction between old-memory records 9
and 10; the separate four-address record package retains that distinction.

Put `a=p(1-q)`, `b=pq`, `Delta=1-a^2+ab>0`. The three-ray matrices are

```
T = [[1-a,1-b,1-a],[a,0,0],[0,b,a]],
pi = (1-a,a(1-a),ab)^T/Delta,
E_always = T,
E_C  = [[b,a,b],[a,0,0],[0,b,a]],
E_E  = [[q,q-b,q],[0,0,0],[0,b,0]],
E_CE = [[b,0,b],[0,0,0],[0,b,0]].
```

The producer checks the symbolic lumping identities against the actual
12-address destination rule for every duty. It also computes the native
12-address group inverse independently and verifies the resulting complete
lag sum. No transient input is deleted from that source check.

## Summing every lag

Let `Pi=pi 1^T`. The group inverse on zero-sum vectors is

```
R=(I-T+Pi)^(-1)-Pi = sum_(k>=0)(T^k-Pi).
```

The reset contraction proves convergence. On the three-ray source it reduces
to the exact formula `R=(I+T-2Pi)/Delta`. Both left/right inverse identities
and `R pi=0`, `1^T R=0` are checked symbolically.

Let `xi_i(n)` indicate that step `n` applies an enabled pulse of ray `i`,
and write `w=E pi`. In the stationary process,

```
C0=Cov(xi(n),xi(n))=diag(w)-w w^T,
Ck=E T^(k-1) diag(w)-w w^T,       k>=1.
```

The row index denotes the later pulse. Between the enabled endpoints there
are `k-1` unrestricted transitions. This temporal orientation matters for the
coherent correction. The complete one-sided covariance sum is

```
K=C0/2+sum_(k>=1) Ck
 =diag(w)/2-w w^T/2+E R diag(w).
```

Its symmetric part `D=K+K^T` is positive semidefinite: it is the limit of
`N^(-1) Cov(sum_(n=1)^N xi(n))`. Absolute convergence justifies the limit.
The source contraction also gives an explicit geometric bound for the
omitted tail; a finite scan of lag lengths is unnecessary.

## Coherent ordering and diffusion

Define `L_i=-i[H_i,.]`. The first two reduced generator coefficients are

```
G0=sum_i w_i L_i,
G1=sum_(i,j) K_ij L_i L_j.
```

For a direct finite-dimensional derivation, inject the stationary product
state with `J`, sum the memory blocks with `L`, and put `Q=I-JL`. Expand
`Psi_h=M0+hM1+h^2 M2+O(h^3)`, where

```
M0=T tensor id,
(M1)_im=E_im L_i,
(M2)_im=E_im L_i^2/2.
```

With the lifted group inverse `Rscript`, the first and second slow step
coefficients are `A=L M1 J=G0` and
`B=L M2 J+L M1 Rscript M1 J`. The logarithmic coefficient is
`B-G0^2/2=G1`. This is an expansion at `h=0` on the reduced slow block; it does
not take a logarithm of the full noninvertible memory channel.

Set

```
H_ord=-(i/2)sum_(i,j) K_ij[H_i,H_j].
```

Then the exact decomposition is

```
G1=-i[H_ord,.]-(1/2)sum_(i,j)D_ij[H_i,[H_j,.]].
```

The second term is a Lindblad dissipator because `D>=0`; its noise operators
are the Hermitian linear combinations of `H_i` from a square root of `D`.
Thus `exp(t(G0+hG1))` is a completely positive trace-preserving semigroup for
`h,t>=0`. Discarding the diffusion while retaining the ordering potential
changes this reduced second-order dynamics.

For always duty at gain one the entire correlation sum gives

```
H_ord=I2 tensor A,
A=i*a*b*(1-a)/(10 Delta^2)
  * [[0,3+sqrt(5),0],[-3-sqrt(5),0,-2],[0,2,0]].
```

At the native rates the prefactor is `603/6666724` and the norm is
`603 sqrt(18+6sqrt(5))/6666724`, approximately `0.0005069709792586407`.
The generator contains `h A`; at general gain the potential scales as `g^2`.
The zero `(0,2)` entry concerns the full stationary lag sum; a nonstationary
finite ordered word can have a different entry.

The potential multiplies the auxiliary identity, whereas the leading scalar
coefficient multiplies `X`. It is an imaginary antisymmetric family potential,
not an imaginary scalar Yukawa coefficient. The scalar/potential projection
in the finite-response package supplies the corresponding conditional
relativistic distinction. No physical charge-parity or mass identification
follows from the auxiliary grading.

## Uniform finite-time control

Use the sum of block trace norms. Set `u=1-p` and

```
alpha=2gd, beta=2g^2 d, c=2alpha/u,
b=beta+alpha c, gamma=b+alpha^2/2,
d2=(2beta+3alpha c)/u, mu=4g^3 d/3.
```

For `h=t/N`, define

```
e3=alpha gamma+h gamma^2/2+(alpha+h gamma)^3/6,
r3=2alpha d2+(beta+b)c+h(beta+b)d2
   +(mu+e3)(1+hc+h^2 d2).
```

The unread comparison, with an arbitrary untouched reference, satisfies

```
||L Psi_h^N J-exp(t(G0+hG1))||_diamond
 <= min(2, alpha h^2(c+h d2)/u + N h^3 r3).
```

This diamond norm is unhalved. At fixed `t` the bound is `O(h^2)`. The
executable evaluates it with exact rational parameters and retains the exact
fraction, rather than interpreting a numerical Choi test as a diamond bound.

To prove the estimate, complete probability and `||H_i||<=g` give
`||M1||<=alpha`, `||M2||<=beta`, `||Psi_h-M0||<=alpha h`, and cubic one-step
remainder at most `mu h^3`. The group inverse has norm at most `1/u` on
zero-sum blocks, while `||Q||<=2`.

Put `C=Rscript M1 J`,
`Dscript=Rscript(M2 J+M1 C-C G0)`, and `J_h=J+hC+h^2 Dscript`.
Then `LC=LDscript=0`, `||C||<=c`, `||Dscript||<=d2`, `||B||<=b`, and
`||G1||<=gamma`. The coefficients through `h^2` in
`Psi_h J_h` and `J_h exp(h(G0+hG1))` agree. Taylor expansion of the contraction
semigroup bounds the latter exponential's cubic remainder by `h^3 e3`.
Multiplication gives intertwining defect at most `h^3 r3`.

The actual preparation differs from `J_h` by a zero-sum memory vector of norm
at most `h c+h^2 d2`. For such a vector `y`, telescoping against `M0` gives
`||L Psi_h^N y||<=alpha h ||y||/u`. This yields the initial-slip term. Summing
intertwining defects gives `N h^3 r3`, since the exact channel and effective
semigroup are contractions. Finally `LJ_h=I`. Every estimate survives tensoring
an untouched reference; the cap two follows from comparison of two channels.

For the finite-horizon coefficient

```
K_N=C0/2+sum_(k=1)^(N-1)(1-k/N)Ck,
```

the same contraction gives

```
||G1-G1_N|| <= 8g^2 d^2[p^(N-1)/(1-p)+1/(N(1-p)^2)].
```

The estimate uses the geometric tail and the convergent weighted-lag sum.
Neither estimate is uniform for arbitrarily large time at fixed nonzero
pulse step. Nonstationary or memory/field-correlated preparations require
additional initial-state terms.

## Independent numerical controls

The producer constructs the full 432-dimensional native channel directly and
compares it to the actual source `joint_channel` implementation. For all four
native duties it compares the exact stationary unread channel at fixed time
one and `N=8,16,32,64,128` to both the first-order mean unitary and the corrected
semigroup. The complete quantum channel is retained in each calculation.
Normalized Choi trace distances are diagnostics on a particular reference
input; the rational estimate supplies the separate worst-case bound.

The declared current-gated alternatives satisfy the same symbolic theorem
when their actual rates and duties are substituted. Physical feedback
selection, clock calibration, interacting particles and a joint early-Universe
observation law are not supplied by this result.
