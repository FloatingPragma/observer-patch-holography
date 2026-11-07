# A controlled local limit for weighted repair

Consider the declared connected weighted family in
[`weighted_source_family`](../weighted_source_family/PROOF.md). Each of B cells
has the weight word (25,25,8,8,8,8,8,8,8,8,8,8). The graph can be its retained
federation graph. The stationary population reference depends on these weights
and the conserved occupation, independently of the positive fixed seam clock.
The reference chooses a mean occupation r∈{2,4} per cell and total m=rB,
hence density one sixth or one third on 12B binary sites. The density is a
declared input. The native program admission
rejects the weighted connected extension; no admission rule is changed.

The proof supplies a bound for every integer B≥1 and every set of k cells,
0≤k≤B. It retains every local occupation sector, including zero-probability
sectors. The use of diagonal states identifies half trace distance with total
variation of their entire configuration distributions. No global state-space
enumeration is needed to evaluate the bound.

## Product reference and exact conditioning

For r=4 let

\[
 z=\frac{1+\sqrt{6401}}{1600},\qquad
 p_h=\frac{25z}{1+25z},\qquad p_l=\frac{8z}{1+8z}.
\]

The positive root satisfies 800z²−z−2=0. Multiplying by
(1+25z)(1+8z) shows that 2p_h+10p_l=4 is exactly that equation.
The polynomial is negative at z=1/20 and positive at z=3/50, and is increasing
between those values. Thus

\[
 \tfrac12<p_h<\tfrac35,\qquad \tfrac14<p_l<\tfrac13.              \tag{1}
\]

For r=2 use z=(sqrt(1256)−16)/1000. Its equation is
1000z²+32z−1=0, obtained by the same denominator multiplication from
2p_h+10p_l=2. This polynomial is negative at z=1/55 and positive at z=1/50;
its derivative is positive for z>0. The same occupation formulas therefore
give

\[
 \tfrac3{10}<p_h<\tfrac13,\qquad \tfrac18<p_l<\tfrac16.         \tag{1b}
\]

Each polynomial has exactly one positive root: its leading coefficient is
positive, its constant term is negative, and its discriminant is positive.
The two choices define different product references. No selection between
them follows from the limit theorem.

Let Q_B be the independent Bernoulli product with 2B high-weight sites and
10B low-weight sites. Equivalently Q_B=Q_1^{\otimes B}, grouping the same
twelve labeled ports into each cell. If S is the total occupation, then
E_Q S=rB. Fix the same declared r throughout each comparison. For every
occupied-site set T,

\[
 Q_B(T)=\frac{z^{|T|}\prod_{i\in T}w_i}
                {(1+25z)^{2B}(1+8z)^{10B}}.
\]

Conditioning on S=rB therefore gives exactly the canonical stationary state
P_B(T)=prod_{i∈T}w_i/e_{rB}(w). The event has positive probability q_B.
With natural logarithms and 0 log 0=0,

\[
 D_{\rm KL}(P_B\Vert Q_B)=\log(1/q_B).                         \tag{2}
\]

The conditional reference is diagonal on the global fixed-number sector.
Its restriction to a cell is a full number mixture. Conditioning that mixture
on local number two gives the native weighted H66 reference, as in the owning
family proof. Such conditioning is neither the full restriction nor a
deterministic preparation.

## A uniform conditioning-probability bound

First consider X∼Bin(n,p), n≥1 and 0<p<1. Its variance is np(1−p)≤n/4.
Chebyshev's inequality, obtained by applying Markov's inequality to
(X−np)², gives probability at least 3/4 in |X−np|<sqrt(n). This interval
contains at most 2sqrt(n)+1≤3sqrt(n) integers. Its largest mass, and hence a
binomial modal mass, is at least 1/(4sqrt(n)).

Write b_j=Pr(X=j). For 0≤j<n,

\[
 \frac{b_{j+1}}{b_j}=\frac{(n-j)p}{(j+1)(1-p)}.                 \tag{3}
\]

This ratio decreases with j, and crosses one at j+1=(n+1)p. A mode is
floor((n+1)p), with its neighboring lower index also a mode when (n+1)p is
an integer. Consequently floor(np) is a mode or one below a mode. In the
latter case, writing a=floor(np),

\[
 \frac{b_a}{b_{a+1}}
 =\frac{(a+1)(1-p)}{(n-a)p}
 \ge\frac{n(1-p)}{n(1-p)+1}.
\]

Similarly ceil(np) is a mode or one above a mode, with its mass at least
np/(np+1) times a modal mass. If np is an integer it is itself a mode, so the
same lower bounds apply. All these indices lie in {0,…,n}.

Under Q_B let X count the 2B high sites and Y count the 10B low sites. They
are independent. For r=4 choose

\[
 a=\lfloor2Bp_h\rfloor,\qquad b=4B-a=\lceil10Bp_l\rceil.
\]

The equality uses the exact mean identity, rather than rounding the two means
independently. By (1), the lower neighbor factor for X is at least 4/9 and the
upper neighbor factor for Y is at least 5/7. Thus

\[
 \begin{split}
 q_B&\ge\Pr(X=a)\Pr(Y=b)\\
 &\ge\frac{4/9}{4\sqrt{2B}}\,
       \frac{5/7}{4\sqrt{10B}}
 =\frac{5}{252\sqrt{20}\,B}
 >\frac1{240B}.                                               \tag{4}
 \end{split}
\]

For the final strict inequality, sqrt(20)<9/2 and 252·9/2=1134<1200.
This elementary bound holds for B=1 as well as arbitrarily large B. It does
not claim a central-limit approximation, an optimal constant, or the same
bound for an arbitrary total occupation under a fixed product reference.

For r=2 take a=floor(2Bp_h), b=2B−a=ceil(10Bp_l). Equation (1b) gives
lower-neighbor factor at least (4B/3)/(4B/3+1)≥4/7 for X and upper-neighbor
factor at least (5B/4)/(5B/4+1)≥5/9 for Y. Their product is 20/63, identical
to (4/9)(5/7). Thus (4) holds unchanged with its own fugacity and total 2B.
All subsequent bounds apply to either declared mean. At B=1 the r=2 canonical
state is exactly the native H66 reference; at larger B the one-cell restriction
has number fluctuations and is not confined to H66.

## Entropy distributes the conditioning cost

Let X_1,…,X_B denote the twelve-bit cells under P_B. Their joint law is
exchangeable: the total-number condition and the repeated weights are
unchanged by every permutation of whole cells. This statement does not
require a cell permutation to preserve the repair graph or its dynamics.

For a finite probability law, define H(X)=−sum_x p(x)log p(x), and conditional
entropy by averaging the entropies of the positive-probability conditional
laws. The finite log-sum inequality proves nonnegativity of relative entropy.
Applying it to each conditional joint law and the product of its marginals
gives H(U|V,W)≤H(U|V). The entropy chain rule follows by expanding
log p(u,v)=log p(v)+log p(u|v) on its support.

Set c_j=H(X_j|X_1,…,X_{j−1}). Cell exchangeability and conditioning imply

\[
 c_{j+1}\le H(X_{j+1}|X_2,\ldots,X_j)=c_j.
\]

Hence the average of the first k terms is at least the average of all B:

\[
 H(P_k)\ge\frac{k}{B}H(P_B),\qquad 1\le k\le B.               \tag{5}
\]

Here P_k is the full law on any k distinct cells with their local port labels
retained. Equation (5) applies when k does not divide B. It uses the entire
configuration entropy, rather than the entropy of the two occupation counts.

The cross-entropy against Q_1^{\otimes k} is the sum of k identical one-cell
expectations. Therefore it is k/B times the corresponding global
cross-entropy. The canonical one-cell mean need not equal its product-reference
mean for this identity. Subtracting (5) and using (2) and (4) gives

\[
 D_{\rm KL}(P_k\Vert Q_1^{\otimes k})
 \le\frac{k}{B}\log(1/q_B)
 \le\frac{k}{B}\log(240B).                                    \tag{6}
\]

The finite ingredient `logSum` is proved in
[`FiniteConditionalRepair.lean`](../../Lean/Thermodynamics/FiniteConditionalRepair.lean).
The cell attachment, entropy-chain argument and quantitative limit here are
prose proofs with executable controls; no additional Lean theorem is claimed.

For completeness, finite Pinsker follows from the same log-sum inequality.
For distributions P,Q with Q faithful, take A={x:P(x)>Q(x)} and
delta=P(A)−Q(A)=sum_x|P(x)−Q(x)|/2. Log-sum on A and its complement bounds
D_KL(P||Q) below by binary relative entropy d(P(A)||Q(A)). At fixed
v=Q(A)∈(0,1), the second derivative of d(u||v) with respect to u is
1/[u(1−u)]≥4, while d(v||v)=0 and its first derivative vanishes there.
Integration twice gives d(u||v)≥2(u−v)². Boundary values follow by continuity;
if delta=0 the conclusion is immediate. Thus the half trace distance obeys

\[
 \boxed{\quad
 \tfrac12\|P_k-Q_1^{\otimes k}\|_1
 \le\min\!\left(1,\sqrt{\frac{k\log(240B)}{2B}}\right).
 \quad}                                                       \tag{7}
\]

For fixed k this tends to zero. For growing k it suffices that
k log B/B tends to zero. Any subset of ports contained in k cells satisfies
the same bound by restriction. Two capacities B and B′, with a common local
port labeling, differ on that local state by at most the sum of their bounds.
This is a static reference limit; it does not construct a dynamical embedding.

## Complete local experiments and computational certificates

A common exactly trace-preserving completely positive quantum instrument
on the retained ports contracts trace distance, including its classical
outcomes and terminal quantum output.
One proof decomposes a Hermitian difference into its positive and negative
parts: positivity and trace preservation bound the output trace norm by the
sum of their traces. An independently appended reference has no effect on
the input trace norm. Thus (7) controls arbitrary specified finite local
instruments, including outcome-adaptive compositions whose success, failure,
hold and stop outcomes are all retained. No normalization conditioned on a
successful outcome is included. No bound is asserted for access to a
correlated exterior or for the continuing connected dynamics crossing the
window boundary. The latter has the separate complete-record bounds in the
owning weighted-family proof; its local-clock boundary hazard does not tend
to zero merely by increasing B.

Stored finite-precision primitive maps can instead have a certified completely
positive trace-gain bound Φ*(I)≤G I. For such a map the same positive/negative
decomposition gives ||Φ(X)||_1≤G||X||_1 for Hermitian X, including an untouched
reference. If an independent implementation bound compares Φ(σ) with a
target instrument Ψ(σ) by half norm at most eta, the complete output bound is
G epsilon+eta for an input half trace distance at most epsilon. This expression
is not capped at one when outputs can be supernormalized. Exact contractivity
is used only for an exactly trace-preserving common instrument; rounded
primitive maps require their own gain and implementation certificates.

The source cell-count grammar B=20·4^ell gives V=240·4^ell and m=20r·4^ell.
Evaluating (7) uses integers with O(log B) bits at fixed precision,
rather than enumerating choose(12B,rB) states. The certificate uses
ln n < (7/10)bit_length(n) for n>1. Indeed n<2^{bit_length(n)} and ln2<7/10;
the latter follows since the first five positive terms of the exponential
series at 7/10 sum to more than two. Square roots are enclosed by integer
square root after exact scaling. The number of precision bits is retained
as the bound becomes small.

The fugacity intervals use integer-square-root enclosures of sqrt(6401)
or sqrt(1256).
Finite examples compute the coefficient of t^m in
(1+25t)^{2B}(1+8t)^{10B}, grouped by the two occupation counts. Within each
count class both references are uniform, so summing the classwise absolute
differences gives the half trace distance of the full local state, with no
measurement information discarded. Every local count class is retained.
These examples check the formulas; the all-domain argument is (1)–(7).

The capacity index, density and fugacity are explicit mathematical choices.
V counts binary sites and is not a supplied cosmological capacity in nats.
The theorem supplies neither preparation time, an admitted weighted growing
observer, a physical spatial volume law, nor a continuum or particle spectrum.
