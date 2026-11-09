# A finite gauge preparation with controlled Wilson-gap weight

The construction uses the periodic cubic coordinate field in
[FIELD_ROUTE](../particle_observer_outputs/FIELD_ROUTE.md). Its compact
group, native M15 representation, full Dirac Fock space, Wilson kernel,
electric metric and Hamiltonian clock are retained. The preparation is
a specified program service, without a claim that it is a ground state
or that architecture choices are unique. No particle data enter it.
The executable preparation code `V_(EL)` and the evolution cutoff `V_K`
are finite. The coordinate space `L²(G^E,Haar)` is their analytic
comparison parent; no program allocating an infinite carrier is asserted.

## 1. Reused results and the precise gap requirement

The [free spatial-overlap calculation](../interacting_extension/chiral_lattice/PROOF.md#complete-free-spatial-overlap-control)
gives, for the actual spatial beta `gamma_tuple[3]`,

```
H_W = beta (aD_W-I),
R(k)^2 = 1+2 sum_(i<j) (1-cos k_i)(1-cos k_j) >= 1.
```

This is the free **periodic** Fourier calculation. It is not assumed for
an open Wilson box with boundary modes. The
[Euclidean chiral patch](../interacting_extension/dependencies/chiral_patch.py)
and [electric connection](../interacting_extension/chiral_electric/PROOF.md)
check supplied admissible backgrounds and fixed-rank bundles. Neither
prepares a fluctuating gauge state with a bound on small Wilson gaps.
The source [coordinate completion](../interacting_extension/chiral_lattice/PROOF.md#source-representation-and-conditional-coordinate-completion)
supplies Haar integration and the dense matrix-coefficient algebra.

The [polynomial field proof](../particle_observer_outputs/FIELD_ROUTE.md#a-bounded-overlap-polynomial-on-the-full-coordinate-hilbert-space)
gives a global norm bound without a gap. Its sign accuracy requires
`sigma_min(H_W)>=1/2`. The relevant additional premise is the prepared
state's weight outside that set, on the same Hamiltonian clock. The
construction below supplies such a finite preparation and an explicit
short-time bound. The [energy-confinement member](REFINEMENT.md) has an
all-time bound at an explicitly different plaquette coefficient.
Neither supplies a uniform continuum estimate.

## 2. A finite source-word amplitude

Let V be the number of sites and E=3V the number of oriented positive
links of a periodic cubic graph. A one-cell torus has three loop links.
For the actual unitary native representation R of dimension 15 put

```
s(U) = (1+Re Tr R(U)/15)/2
     = 1-||R(U)-I||_F²/60,                 0 <= s <= 1.
```

For an integer n>=1 choose L=n². Let P_G be normalized averaging over
the independent gauge transformations at every site and define

```
p_L(U) = product_edges s(U_e)^L,
psi_L = P_G p_L,
Z_L = integral |psi_L|² dHaar,
Psi_L = psi_L/sqrt(Z_L).                               (1)
```

The amplitude is nonnegative and gauge invariant. The gauge average acts
on amplitudes; it is not a classical mixture of background fields. Each
factor is a polynomial of native matrix coefficients and their conjugates of degree
at most L on its link. Gauge averaging preserves that degree. Thus
`Psi_L` lies in the actual source word space `V_(E L)`. It is a
normalizable finite polynomial, with no sharp spectral cutoff or
delta-function gauge configuration.
The bad-gap probability on a finite code is the expectation of the
compressed parent effect `P_K 1_Bad P_K`. It is not the spectral gap of
the truncated Wilson matrix `P_K H_W P_K`.

## 3. Explicit Haar normalization

Use the unnormalized Frobenius metric on the image R(G) in real ambient
dimension `2*15²=450`. Every matrix in that image has Frobenius norm
sqrt15, measured from the **zero matrix**. A maximal 1/n-separated set
has disjoint ambient balls of radius 1/(2n), all contained in the ambient
ball of radius `sqrt15+1/(2n)` centered at zero. Their number is at most

```
(1+2 sqrt15 n)^450 <= (1+8n)^450.
```

Maximality makes the balls of radius 1/n centered at those group points
cover R(G). The Frobenius distance and Haar measure are invariant under
left translation. Consequently the Haar probability of the identity
ball of radius 1/n is at least `(1+8n)^(-450)`.

On that ball, `s(U)>=1-1/(60n²)`. Bernoulli's inequality gives
`s(U)^(n²)>=59/60`. Since gauge averaging preserves the Haar integral,
and the constant normalized function is an allowed test vector,

```
integral s^L >= (59/60)(1+8n)^(-450),
Z_L >= |integral psi_L|²
    >= (59/60)^(2E) (1+8n)^(-900E) > 0.                 (2)
```

All powers are exported in factored form. No enormous integer dimension
or normalization denominator is allocated to evaluate this certificate.

## 4. Small-gap suppression in every gauge frame

If `max_e ||R(U_e)-I||<=delta`, each of the three covariant shift pairs
differs from its free pair by at most delta. Their orthogonal spin
projections preserve that bound. Therefore

```
||H_W(U)-H_W(I)|| <= 3 delta,
sigma_min(H_W(U)) >= 1-3 delta.                         (3)
```

The smallest singular value is globally 1-Lipschitz in operator norm:
for every unit vector v, `||Av||>=||Bv||-||A-B||`, and take the infimum
in both directions. No eigenvalue differentiability is assumed.

Let `Bad={U:sigma_min(H_W(U))<1/2}`. Its definition is gauge invariant.
For every U in Bad, and **every** site gauge transformation g, (3) shows
that some transformed link has operator-norm distance greater than 1/6
from identity. Frobenius distance is at least operator-norm distance,
so one factor satisfies `s(U_e^g)<=2159/2160`. All other factors are at
most one. Positivity of the integrand in (1) therefore gives

```
psi_L(U) <= (2159/2160)^L                    on Bad,
Prob_PsiL(Bad)
 <= (60/59)^(2E) (1+8n)^(900E) (2159/2160)^(2n²).       (4)
```

The positive exponentials can be enclosed with integers only. Since
`(1-x)^k<=exp(-kx)` and e>=2, a valid upper bound is

```
Prob_PsiL(Bad) <= min(1,2^b),
b = 2E+900E ceil(log2(1+8n))-floor(n²/1080).              (5)
```

For each fixed E and requested binary accuracy, doubling n eventually
makes b sufficiently negative. This is a finite algorithm and its
volume dependence is explicit. The proof does not infer a growing-volume
vacuum or an inexpensive preparation from existence of this algorithm.

## 5. Electric energy in the actual source metric

The coordinate completion uses the invariant metric
`g_ab=-Tr_M15(T_a T_b)` for anti-Hermitian generators, equivalently the
positive trace metric for the corresponding Hermitian generators. The
native quadratic Casimir packet verifies that conventional fundamental
SU(3) and SU(2) generators have metric 2, while `q=6Y` has metric 120.
In an orthonormal source basis their M15 operator norms satisfy

```
sum_A ||t_A||² <= 11/2+36/120 = 29/5.                   (6)
```

The eleven nonabelian conventional generators each have norm at most
one on the retained fundamental/antifundamental blocks, and the largest
integer hypercharge magnitude is six. The bound deliberately allows
more than their exact norm sum.

Specify electric energy as one half of this positive Casimir on each
link, a permitted fixed coefficient of the same field program. On tensor
words of length at most L, the differentiated generator has norm at most
`L ||t_A||`. Every irreducible component has the same Casimir under the
Peter-Weyl coordinate realization. Gauge averaging commutes with these
Casimirs and preserves the finite word space. Thus the normalized state
in (1) satisfies the electric-domain bound

```
<H_E> <= (29/10) E L².                                 (7)
```

This bound uses the finite representation content, not division of an
uncontrolled derivative norm by a small normalization. The state is
smooth and belongs to every finite-volume electric-operator domain.
The chosen concentration can have a very large electric energy. It is
not a variational proof of a gauge vacuum.

## 6. Finite preparation and a complete executable seed

The per-link word space has dimension at most 451^L: each of L slots can
hold identity, one of 225 native matrix entries or its conjugate. Hence
the complete coordinate preparation code has dimension at most 451^(EL),
requiring at most 9EL logical bits. Haar inner products and gauge averages
are finite representation computations: form the tensor generator
matrices, project to their common invariant nullspace and resolve the
finite Gram matrix. These exact finite operations give the coordinates
of (1). The coefficient list at the large certified L is not materialized
by this package.

Let `a=<1,Psi_L>`, real and positive. If a=1 use identity; otherwise the
Householder unitary

```
u=(Psi_L-1)/sqrt(2(1-a)),       W=I-2|u><u|
```

sends the specified Haar blank to Psi_L. Both vectors are gauge invariant,
so W is gauge invariant. Its support lies in the finite source word code;
it acts as identity on the entire orthogonal complement. It is a complete
unitary on every input and reference, not a reset of an unknown state.
It may be a declared preparation primitive or compiled under the existing
[finite apparatus contract](../particle_observer_outputs/READOUT_CONTRACT.md).
No retry or success-conditioned histogram is used. Compiler error and
preparation work require separate budgets if such a physical service is
claimed; no large compiled tape is asserted here.

For a one-cell torus at L=1, each native character is conjugation invariant,
so the gauge average is explicit. The five native irreps and their five
conjugates are inequivalent because their nonzero q labels are distinct.
Schur orthogonality gives the normalized one-link character vector

```
Psi_1 = (30 |trivial> + sum_(ten irreps) |character_r>)/sqrt(910),
integral s² = 91/360.
```

The code constructs the actual 11-dimensional Householder and verifies
its full matrix unitarity and prepared column. Three copies give the
1331-dimensional character code. The apparatus measures the eleven
listed electric irreps plus the outside-menu complement on each link.
Each listed projector retains its full Peter-Weyl multiplicity; it is
not only the character vector. Thus all 12³=1728 original records are
retained, including 397 zero-probability complement records.

The one-link probabilities are 90/91 for the trivial representation and
1/910 for each of the ten irreps. Their source-metric Casimirs are

| Native block | Source-metric Casimir |
| --- | ---: |
| `(3,2), q=1` | 21/20 |
| `(bar3,1), q=-4` | 4/5 |
| `(bar3,1), q=2` | 7/10 |
| `(1,2), q=-3` | 9/20 |
| `(1,1), q=6` | 3/10 |

Conjugates have the same Casimir. The exact three-link mean electric
energy is 99/9100 in source units. These are complete quantum-coordinate
preparation records with a finite source code. The L=1 seed does **not**
satisfy a useful bound from (5), and its code is not claimed invariant
under the full interacting Hamiltonian.

## 7. A same-clock accuracy interval

Keep the field preparation with all 60 fermion modes at one site filled
and the other sites empty. Its determinant character is trivial; tensoring
it with (1) is gauge invariant. Fermion number N=60 is conserved by the
kinetic/electric/plaquette program used here. Set additional onsite
interactions to zero in this declared control. On this invariant number
sector, the polynomial-field bounded part obeys

```
B=N(2+epsilon_m)/a+|kappa|E,
H_m+B >= 0,
<H_m+B> <= K=(29/10) E L²+2B.                          (8)
```

For any positive operator A and finite first energy moment,
`||(exp(-itA)-I)Psi||²<=2|t|<A>`. This follows from
`2(1-cos x)<=2|x|` in the spectral integral and needs no second energy
moment. Applying it to `H_m+B`, with the irrelevant scalar phase removed,
gives

```
sqrt(Prob_Psi(t)(Bad)) <= sqrt(eta0)+sqrt(2K|t|).         (9)
```

If (5) gives `eta0<=2^(-2b)`, the rational duration
`T=2^(-2b)/(2K)` guarantees bad weight at most `2^(2-2b)` throughout
the interval. It is the same Hamiltonian time as the field readout, not
a preparation-dependent rescaling.

The exact-sign coordinate Hamiltonian has the same electric and plaquette
terms and the same state. Its bounded difference from the polynomial
field, on this number sector, is at most
`(N/a)[epsilon_m+(2+epsilon_m)sqrt(eta(t))]` on the evolving vector.
Duhamel's formula therefore bounds their state-vector distance, and the
TV of every complete final detector, by

```
(N T/a)[epsilon_m+(2+epsilon_m)2^(1-b)].                 (10)
```

This comparison has one evolution interval followed by a complete final
instrument. It does not cover arbitrary intervening measurements that
could change the prepared gap weight. The owning source cutoff theorem
applies with initial coordinate degree EL; the certificate also computes
its separate finite-cutoff bound for degree-preserving final records.
The finite evolution uses `P_K H_m P_K`, including the compressed
polynomial of the parent kernel. Replacing that polynomial by one
evaluated after compressing `H_W` is a different operator and is not
covered by this comparison.
The same errors must be added if both preparations or both evolutions are
implemented approximately.

For E=3 and a 20-bit requested initial gap bound, the conservative planner
chooses n=8192, L=67108864 and total word degree 201326592. Its binary
bad-weight bound is far smaller than the target, but its code and
preparation scale are not practical. The energy bound makes (9)'s certified
duration extremely short. No nonzero time interval uniform in volume or
refinement is proved.
The [conserved-energy construction](REFINEMENT.md) supplies an all-time
bound for a different, explicitly declared high-plaquette-cost member
of this Hamiltonian family. It does not change the 1/100 control above.

## 8. An exact native gap-closure control

Gap control cannot be omitted solely because the native group is fixed.
On one periodic cell choose actual weak SU(2) links

```
U1=(I+i sqrt3 sigma1)/2,  U2=(I+i sqrt3 sigma2)/2,  U3=I.
```

All native weak doublets are unitarily equivalent to the fundamental
two-dimensional representation; this constant internal basis change
preserves the Wilson spectrum. On each weak-doublet block, the actual
eight-dimensional Dirac/weak
Wilson matrix has squared eigenvalues 0 and 3, each of rank 4. Native M15
contains four such weak doublets and seven weak singlets. The complete
60-dimensional source symbol therefore has squared eigenvalues
`0(rank 16),1(rank 28),3(rank 16)`, verified by exact matrix arithmetic.
Continuity gives a positive-Haar-measure neighborhood with gap below 1/2.
The Haar blank does not guarantee the admitted gap.

This calculation evaluates a coordinate symbol. A point-coordinate value
is not a normalizable L²(G) preparation, and its auxiliary trace weights
are not substituted for the quantum seed records in section 6.

## 9. Physical boundary

The construction supplies a finite gauge-invariant word preparation,
an explicit gap-probability and electric-energy bound, complete seed
records and a same-clock finite-time overlap comparison. The declared
energy-confinement member has an all-time gap estimate. Its very large
coupling and code costs enter the local-dynamics error bounds as well.
There is no interacting-vacuum, controlled refinement, chiral-regulator,
Lorentz-symmetry or physical-particle theorem here. A relativistic
interpretation requires common local quantum correlations at fixed
physical scales, with controlled normalization and surviving interactions
along the actual spacing, degree, volume and coupling family.
