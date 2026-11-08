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

The even cube reflections a in {000,011,101,110} act by permuting
sites x to x xor a. These are graph symmetries of the open cube, not
small translation momenta. Each commutes exactly with H_B and Q. The four
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
reflection-character labels. It does not exclude a different interpolator,
a spectral subbranch, or a sharp excitation in the zero-field continuum.
The two local finite charge inputs are not asserted to have identical
weights or energy measures.

The variance test specifies the spectral quantity in the asymptotic
estimate. It does not supply a Lorentz dispersion, boost covariance,
a waveoperator or an experimental cross section.

## 4. Selecting atomic subbranches without guessing a dispersion

Restrict the same measure mu to a compact momentum window and its bounded
positive transfer band, retaining all other records as complementary
outcomes. Its mass need not be one. Let K_l be nested finite dyadic
partitions generating the momentum Borel sets, and I_j nested dyadic
energy partitions of the transfer band. Use half-open bins, with the last
endpoint included, so atoms on boundaries are counted once. Define the
statistic of original complete joint probabilities

    A_(j,l)=sum_(K in K_l) max_(I in I_j) mu(I times K).    (5)

For fixed j this increases under momentum refinement and is at most
mu(total). Disintegrate mu=mu_k sigma(dk), with mu_k normalized only as
a mathematical conditional kernel. Conditional expectations with respect
to K_l converge in L1 to each of the finitely many functions mu_k(I).
The finite maximum is Lipschitz in their l1 distance. Consequently

    lim_l A_(j,l)=integral max_(I in I_j) mu_k(I) dsigma.  (6)

The maximum mass of nested energy bins decreases to the largest atom
of a finite measure. To see this, an atom gives the corresponding lower
bound at every depth. Conversely, bins carrying a fixed positive weight
have a convergent subsequence of centres on the compact energy band;
every neighbourhood of the limit eventually contains such a bin and
has at least that weight. Continuity from above gives an atom. The same
argument applies to the restriction to any half-open parent bin; an
excluded boundary cannot carry a positive mass of that restricted measure.
Dominated convergence gives

    A_* := lim_j lim_l A_(j,l)
         =integral m(k) dsigma,
    m(k)=max_E mu_k({E}), with max(empty)=0.               (7)

Thus A_*>0 is equivalent to a nonzero atomic dispersion subbranch in
this actual prepared measure. A continuous component and multiple atomic
branches are allowed. The order of limits is essential. For the measure
uniform on E=k, k in[0,1], aligned dyadic partitions give
A_(j,l)=min(1,2^(l-j)). Equation(7) is1, whereas reversing the two limits
is0. A finite-width momentum bin can conceal an exact moving branch.

The branch can be selected measurably from the source measure itself.
For every dyadic parent energy bin I, let m_I(k) be the limit of the
largest masses among its descendant bins. It is a measurable function
and equals the largest atom within I by the preceding argument. On
m(k)>0 start at the whole band. Choose the left child when its m_I equals
m(k), and otherwise choose the right child. At least one child has that
value. Each choice is measurable. At fixed k only finitely many atoms
can have weight m(k)>0, and the nested bins retain at least one of them;
their vanishing diameter selects exactly one such atom. Their endpoints
converge to a measurable value e_*(k), with
mu_k({e_*(k)})=m(k). No measurability of an unspecified maximizer is assumed.

The joint spectral projection onto
{(E,k):m(k)>0 and E=e_*(k)} is therefore a genuine reducing projection
for the same L and P. Its original weight in C Omega is A_*. Its
complement retains the rest of the response. Subtract this graph from
the conditional kernels, without renormalizing them, and repeat the same
construction. The resulting disjoint graph projections resolve every
conditional atom: an atom of weight w>0 can be preceded by only finitely
many atoms of weight at least w. The sum of all branch weights is at
most the original preparation probability. A zero first value A_*
means that this prepared measure has no atomic branch.

For any fixed finite partitions, (5) is also the maximum probability of
a record event choosing one energy bin in each momentum cell. Therefore
complete output half-trace error epsilon changes A_(j,l) by at most
epsilon. Sampling requires simultaneous probability guarantees before
maximization. The source joint-clock constructions supply particular
finite joint-cell readouts and their full error bounds; this argument
does not declare arbitrary continuum graph projections local or natively
implemented. Sharp projections can be nonlocal, and the double limit has
no derived effective source-size or observation-time modulus.

The source estimates used here prove neither A_*>0 nor a positive
shrinking-bin lower bound uniform in both required limits. They also do
not supply branch smoothness, separated velocities, a uniform residue
bound on a usable momentum region or the source-error estimates needed
for scattering. The moment and current identities in
`spin_shell_diagnostic` and `spin_current_spectrum` admit atom-free
positive-measure controls; those controls are not alternative source
ground states. The exact source's projected current/force estimate is
the additional physical premise. Finite observability and finite
spectral projectors alone do not prove it.
