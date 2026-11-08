# A source test for a scalar dispersion comparison

## 1. The actual joint energy-momentum measure

Let C be one of the bounded, positive-energy-transfer source creators in
LOCAL_CREATORS.md. The unchanged commuting source generators L,P give
the positive joint measure mu(E,k) of C Omega. It has bounded energy
support. Write sigma(dk)=integral mu(dE,dk), and disintegrate into
conditional probability measures mu_k(dE) for sigma-almost every k.
The mean e(k) and variance v(k) are finite measurable functions.

For every real measurable candidate omega(k) and compact packet fhat,
conditional expectation gives the exact identity

    integral |E-omega(k)|² |fhat(k)|² dmu
       =integral v(k)|fhat(k)|² dsigma
        +integral |e(k)-omega(k)|² |fhat(k)|² dsigma.       (1)

Thus the least scalar one-packet derivative defect is the source
conditional energy variance. It vanishes precisely when this prepared
measure is concentrated on the graph E=e(k), on the packet support.
The optimal mean need not be smooth or relativistic. A positive variance
is a failure of the full chosen preparation to be a single graph; it
neither excludes a spectral subbranch nor identifies a physical species.
A bounded local operator in a field theory can create several particles
and need not have a single-particle spectral measure.

There is also an exact strong-limit criterion for a fixed scalar
comparison. Put D=L-omega(P) on the joint spectral domain. If
exp(itD)phi has a strong limit psi as t tends to either infinity, then
exp(isD)psi=psi for every fixed s. Unitarity consequently gives
||phi-psi||=||exp(itD)phi-psi||, so phi=psi and Dphi=0 in the spectral
sense. Thus this one-packet phase comparison has a strong limit only
on its actual dispersion graph. Failure of (1) to vanish excludes that
particular comparison, without excluding a larger comparison space or
a source spectral subbranch.

One can introduce energy itself as an extra comparison coordinate and
use the joint spectral theorem to make the one-vector evolution exact.
That construction does not supply a tensor-product multipacket metric,
separated velocity supports or a scattering identification. Those are
independent requirements in the product-defect estimate.

## 2. Narrow filters retain their complete probability accounting

For a fixed candidate omega, define the spectral window
D_delta={|E-omega(k)|<=delta}. The complete ideal detector has this
success effect and its complement. Its original probability is
p_delta=mu(D_delta), and its unnormalized derivative defect is at most

    delta sqrt(p_delta).                                 (2)

As delta decreases to zero, continuity from above gives

    p_delta -> mu({E=omega(k)}).                          (3)

A graph carrying zero source weight therefore has vanishing original
success probability under arbitrary narrowing. Dividing by that success
probability can produce a small conditional defect; it does not produce
a fixed nonzero original preparation weight. Equations(2)-(3) use the
actual joint spectral measure. They neither assert an available native
implementation of every such joint window nor replace its complete
instrument by a normalized success branch.

## 3. Exact finite source control

The executable instantiates (1) on the same selected eight-site source
H_B at J=1 and B=1/4. Its coherent local inputs are
A_q Omega=Z_0(I-q epsilon_0 X_0)Omega/2, q=±1.
This is the local two-record instrument, without fine-site dephasing.

The even cube translations a in {000,011,101,110} act by permuting
sites x to x xor a. Each commutes exactly with H_B and Q. The four
character projections are

    P_k=(1/4)sum_a (-1)^(k dot a)T_a,
    k in {000,001,010,011}.

The code checks every group product, all source commutators, all
projector products and sum_k P_k=I as integer matrices. The complete
preparation and character instrument has eight records P_k A_q and
sum_(q,k) A_q* P_k A_q=I. A subsequent energy measurement is a
destructive spectral audit with all outcomes retained. It is separate
from the coherent initial preparation.

For phi=P_k A_q Omega set n=||phi||²,
s_1=<phi,(H_B+6I)phi>, s_2=||(H_B+6I)phi||². The optimal unnormalized
scalar defect is s_2-s_1²/n. Subtracting the ground energy or the shift6
does not change this variance. Twelve negative-projector bonds and
eight selection terms give -13I<=H_B<=I, hence ||H_B+6I||<=7.

The certified normalized ground proposal error is epsilon. Since
||P_k A_q||<=1, the errors in n,s_1,s_2 are at most
2epsilon,14epsilon,98epsilon. Exact rational moments of the retained
integer ground proposal, these error intervals, and outward division
certify every character weight and variance. All256 source rows enter
the contractions. No eigensolver output is treated as an error bound.

For the entire original preparation the least squared scalar
energy-character defect has the strict bounds

    q=+1: sum_k(s_2-s_1²/n)>49/1000,
    q=-1: sum_k(s_2-s_1²/n)>77/1000.                       (4)

The receipt gives narrower intervals and conditional means for every
character. They are excitation energies in the selected finite source
clock, not particle rest masses. This excludes fitting a single scalar
energy to the whole chosen finite preparation using these four common
translation labels. It does not exclude a different interpolator,
a spectral subbranch, or a sharp excitation in the zero-field continuum.
The two local finite charge inputs are not asserted to have identical
weights or energy measures.

The variance test specifies the spectral quantity in the asymptotic
estimate. It does not supply a Lorentz dispersion, boost covariance,
a waveoperator or an experimental cross section.
