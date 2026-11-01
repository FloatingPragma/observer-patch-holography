# Contact resonances narrow; fixed wavepackets become transparent

The same K=0 Weyl/contact action now has a controlled continuum conclusion
on **both signs of energy**, including its negative-energy resonances.
For c=v/a, fixed v>0 and any real finite coupling U(a), its energy-normalized
bright scattering channel tends strongly to the identity on every fixed
L2 wavepacket. The convergence is uniform over the coupling choices for
bounded packet densities on a fixed bounded energy window. This strengthens
[the earlier positive-energy bound](../contact_continuum/README.md).

It is deliberately **not operator-norm convergence** on windows containing
negative energies: repulsive U can be tuned to keep S(E0)=-1 at one chosen
negative E0. The resonance narrows and does not retain a finite response
for a fixed packet. A packet whose bandwidth shrinks with a is a different
experiment and is not covered by the fixed-packet conclusion.

The original action has U>=0. The proof also covers the algebraic signed
rank-one extension U in R; this does not replace the source's repulsive law.
The model energy zero and total momentum K=0 stay fixed. A filled physical
vacuum, diverging energy recentering, additional interactions, nonzero-K
scaling and many-body continuum reconstruction remain outside this theorem.
No particle identity or complete interacting QFT is inferred.

## Source measure and a quantitative monotonicity theorem

Use the full eight-node contact measure mu at c=1 from
[the on-shell attachment](../contact_scattering/ONSHELL.md). It has mass one,
support in [-sqrt(12),sqrt(12)], no contact weight on the dark zero-energy
point subspace, and finite inverse second moment M>=1/12. On |t|<2,

```
rho_1(t)=t^2 sum_(n>=0) b_n t^(2n),
b_n=A_n/(4 pi^2 4^n),  0<A_n<=1.
```

Write a_1(x)=Re g_1(x+i0). The prior principal-value decomposition is

```
a_1(x) = -M x + 2 sum_n b_n x^3 I_n(x)
         + x^3 integral_(|t|>=1) dmu(t)/(t^2(x^2-t^2)),
I_n(x)=x^(2n-1) atanh(x)
       - sum_(k=0)^(n-1) x^(2k)/(2(n-1-k)+1).
```

For 0<x<=1/4 the derivative of every exterior integrand is negative:
x^2(x^2-3t^2)/[t^2(x^2-t^2)^2]<0. The negative polynomial terms in
(x^3 I_n)' may also be dropped for an upper bound. Since
atanh(x)<=x/(1-x^2), the remaining term is at most
(2n+3)x^(2n+2)/(1-x^2). Summing the geometric majorant, with pi>3, gives

```
a_1'(x) <= -1/12
 + x^2 (3-x^2/4) / [18(1-x^2)(1-x^2/4)^2]
 <= -1/12 + 6112/535815
 = -154157/2143260 < -1/16.
```

The positive majorant increases on this interval. The explicit monomial
formulas and the geometric b_n bound give locally uniform convergence of
the differentiated series; the exterior integral is separated from its
singularity. These facts justify differentiating, rather than treating a
singular integral's kernel derivative as an ordinary integral. At zero the
formula extends with a_1'(0)=-M. Evenness of a_1' proves the same bound on
[-1/4,1/4]. No broadened numerical grid supplies this monotonicity.

## Uniform integrated scattering bound

Let 0<Emax<=c/4. Define

```
A_c(E)=c^2 a_c(E)=c a_1(E/c),
B_c(E)=pi c^2 rho_c(E).
```

Then A_c'(E)<=-1/16 throughout the window, while the exact eight-node
angular series gives

```
0<=B_c(E)<=w=16 Emax^2/(63 pi c).
```

For U!=0 put tau=c^2/U. The inherited stationary rank-one formula yields

```
|S(E)-1|^2 = 4 B_c(E)^2 / [(tau-A_c(E))^2+B_c(E)^2]
          <= 4 w^2 / [(tau-A_c(E))^2+w^2].
```

At E=0 take the continuous value S=1; one point does not affect an L2
operator. Changing variable y=tau-A_c(E), with dy/dE>=1/16, and extending
the resulting integral to the real line gives

```
integral_(-Emax)^Emax |S(E)-1|^2 dE
  <= 64 pi w = 1024 Emax^2/(63 c).
```

The trivial bound is 8Emax. For U=0 the left side is zero. Neither bound
depends on the size or sign of U; arbitrarily large, oscillating and tuned
coupling sequences are included. The exact rational calculator evaluates
these inequalities and an explicit spacing for a requested squared-norm
error. An exceptional set where |S-1|>=epsilon has Lebesgue measure at most
the integral bound divided by epsilon^2.

## Which common Hilbert space and which limit

The contact-generated free ac channel has cyclic spectral measure
rho_c(E)dE. Multiplication by sqrt(rho_c) maps it to Lebesgue-normalized
energy amplitudes on the interior band. For each c extend the resulting
scattering multiplier by one outside that band to L2(R,dE). The other
reducing ac channels already scatter identically. Dark point states and
any isolated pair remain in the full model Hilbert space; they are not
incoming ac packets and are not included in this assertion.

For a fixed packet f supported in [-Emax,Emax] with |f(E)|^2<=D,

```
||(S_c,U-I)f||^2 <= D min(8Emax,1024 Emax^2/(63c)).
```

Bounded, compactly supported functions are dense in L2(R), and all these
multipliers obey ||S-I||<=2. Approximation therefore proves strong
convergence for every fixed L2 packet as c->infinity, for every U(c). This
is an explicit spectral comparison map, not an established physical
momentum-space identification or a many-body asymptotic theorem.

## Retained resonance and numerical cross-checks

For any fixed E0<0 and c>=4|E0|, the previous sign theorem gives a_c(E0)>0
and the density is positive. Thus the exact choice U_c=1/a_c(E0)>0 gives
S(E0)=-1. Continuity near E0 makes the essential operator-norm distance
from identity equal to two on a window containing E0. Consequently a
strong limit must not be replaced by uniform pointwise or norm convergence.
Also a_c(E0)~ -M E0/c^2, so U_c/c^2 tends to -1/(M E0).

The diagnostic computes the actual principal value without extrapolating
a finite broadening. Integrating one sine axis analytically leaves the
arcsine-convolution density of s=sin(k1)^2+sin(k2)^2. For positive x=E/c,
only s>x^2/4 contributes to the real part. The substitution
s=x^2/4+u^2 removes its square-root endpoint. Tests independently integrate
the original arcsine convolution and compare the resulting finite-eta
resolvent against the original full-torus source quadrature. The diagnostic
is floating arithmetic, separate from the rational proof certificate.
The retained E0=-1/16 rows keep the resonance distance two while the
response at E0/2 decreases. They do not identify this resonance as a mass.

The only imported scattering step is the previous trace-class/rank-one
attachment, whose premises were checked on this source. Its primary
references are [Behrndt, Malamud and Neidhardt](https://arxiv.org/abs/0902.3568)
and [Simon's account of Kato's work](https://arxiv.org/abs/1711.00528).
The all-energy derivative and packet bounds above are derived here.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.contact_wavepacket --check
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/contact_wavepacket
```
