# Endpoint variational energy and a controlled complement

## The actual Hamiltonian and the trial domain

Use the unchanged [primary member](../particle_common_contract/PROOF.md).
Write Hhat=aH=H0+V with the exact dimensionless split of the
[spectral-scale proof](../particle_primary_spectral_scale/PROOF.md).
All expressions below evaluate the interaction coefficient at one.
The physical product vector Omega fills both negative-beta Dirac rows
in every native block and has constant links. It has H0 energy
E00=-90 times the number of cells. It is an auxiliary variational vector,
not a replacement for the selected polynomial density or its idle
reference. An upper bound on its Rayleigh quotient bounds the ground
energy of the entire physical source, even if that ground lies in a
different occupation sector.

Split V into the five number-preserving native-block hopping operators
J_r and the five real character-loop operators M_r. Let B be the actual
number of interior bonds and P the number of plaquettes. Define

```
h_r=J_r Omega,     l_r=M_r Omega,
M_r=-(1/210) sum_p (chi_r(U_p)+conjugate(chi_r(U_p))).        (1)
```

Whole parent products are compressed only after multiplication. Each
displayed vector contains at most one matrix entry on any link, so it
is retained at every selected cutoff. The inherited packet calculation
gives mutually orthogonal vectors, orthogonal also to Omega:

```
||h_r||²=B d_r,       ||l_r||²=g=P/22050,
(H0-E00)h_r=Delta_r h_r,    Delta_r=6+(105/2)C_r,
(H0-E00)l_r=D_r l_r,        D_r=210 C_r.                   (2)
```

The native order is q=(-4,-3,1,2,6), dimensions (3,2,6,3,1), and
Casimirs (4/5,9/20,21/20,7/10,3/10). The trial space is the span of
Omega and these ten packets. Its dimension is eleven for every selected
box. No invariant-subspace claim is made.

Every trial vector is physical and has the same five occupations
N_q=2d_q times the number of cells. The physical cube action in the
[channel proof](../particle_primary_channel_labels/PROOF.md) fixes Omega
and each global packet sum, so these vectors have its A1 label. That
finite cubic label is not continuum spin, a protected topological defect,
or a single charged-particle identification. No assertion that the full
physical ground has these occupations or cubic label is required.

## Exact compressed matrix at interaction coefficient one

The source graph is bipartite. Multiplying matter by minus one on one
sublattice reverses every hopping term and preserves the product vector
up to a phase. Odd hopping vacuum expectations therefore vanish.
Pointwise in the link coordinates, the squared hopping norm of each
native block on Omega is B d_r. Different bonds or blocks have distinct
particle-hole occupations. Thus its contraction with any single
mean-zero magnetic loop vanishes. A magnetic operator commutes with
every hopping operator in the parent coordinate algebra. These facts give

```
<h_r,V h_s>=0,                <h_r,V l_s>=0.               (3)
```

Equation(3) is also true after the finite cutoff: both endpoint packet
vectors are retained, so their matrix element of the compressed V equals
the parent matrix element. This does not infer that arbitrary compressed
hopping and multiplication operators commute.

For loop packets a triple of distinct plaquettes has an edge in only one
loop: two distinct elementary plaquettes share at most one edge, so the
other two cannot cover all four edges of a third. The same reasoning
handles two repeated plaquettes and a distinct third. Independent Haar
integration then leaves only three loops on the same plaquette. Define

```
N_rs=integral (chi_r+conjugate chi_r)
              (sum_t chi_t+conjugate chi_t)
              (chi_s+conjugate chi_s) dg.
```

Each integer N_rs is computed by two independent exact character methods.
The first uses the actual fifteen exterior-basis weights and the SU3,
SU2 Weyl constant term, divided by Weyl group order12; the Q6 exponent
enforces the central U1 integral. The second counts singlets in the
ordered triples of native and conjugate irreducible blocks. The color
factors are 3, conjugate3 or1, and the weak factors are2 or1, so their
triple singlet multiplicities are explicit. The methods agree entrywise.
They also give

```
integral (Re chi_M15/15)²=1/90,
integral (Re chi_M15/15)³=1/500,
sum_(r,s) N_rs=54,
<l_r,V l_s>=-P N_rs/210³.                                (4)
```

In the orthonormal packet basis the exact compression of Hhat-E00 is

```
              [ 0       u_h*          u_l*       ]
K_trial =     [ u_h     diag(Delta)      0        ],
              [ u_l       0            A        ]
u_h,r=sqrt(B d_r),    u_l,r=sqrt(g),
A_rs=D_r delta_rs-N_rs/420.                              (5)
```

The rational matrix A is positive by strict diagonal dominance. Matrix(5)
is an endpoint compression, with no perturbative extrapolation. As a
separate check, the two-vector space spanned by Omega and v=V Omega has

```
z=||v||²=15B+P/4410,
<v,(Hhat-E00)v>=720B+(5389/171500)P,
<Omega,V³ Omega>=-P/171500.                              (6)
```

Its negative Ritz shift is `(b-sqrt(b²+4z))/2`, where b is the second
line of (6) divided by z. The eleven-vector compression contains that
space and resolves the different physical packet denominators.

## Rational endpoint upper bound

Matrix(5) has exactly one negative eigenvalue -d. Its Schur equation is

```
d = sum_r B d_r/(Delta_r+d) + g 1* (A+dI)^(-1) 1.           (7)
```

The right side is positive and strictly decreasing for d>=0. The root
lies between zero and sqrt(z). Exact rational bisection encloses it in
[d_lo,d_hi]. The code uses positive-pivot rational solves and retains a
root bracket of width at most2^-60. It uses no rounded eigenvalue as a
certificate.

Set `t_h,r=1/(Delta_r+d_lo)` and `t_l=(A+d_lo I)^(-1)1`.
These coefficients are positive. The normalized physical trial is

```
psi=(Omega-sum_r t_h,r h_r-sum_r t_l,r l_r)/sqrt(1+s),
s=sum_r B d_r t_h,r²+g sum_r t_l,r²,
q=sum_r B d_r t_h,r+g sum_r t_l,r,
e=(-q-d_lo s)/(1+s).                                    (8)
```

Its exact rational Rayleigh shift satisfies `-d_hi<=e<0` and
`0<=e+d<=e+d_hi`. Consequently `(E00+e)/a` is a strict improvement
over the physical product-vector ground upper bound, for the full
interacting source at the specified member. The source's inherited
spectral lower bound supplies a full ground enclosure; a tighter
independently proved endpoint lower bound may replace that lower end.
Ritz upper levels are not asserted to be exact levels of H.

## Explicit error from the omitted physical space

A residual bound is required even though the ground upper bound itself
uses only the variational principle. The following finite product-vector
lemma supplies it without a dense spatial allocation. Let
`X=sum_i X_i`, where the X_i are centered, bounded, even local operators
on a product vector, with norms bounded by b_i. Put

```
J=sum_i b_i,    D=max_i sum_(j whose support meets i) b_j,
U4(J,D)=3J²D²+16JD³.
```

Then `<X⁴><=U4(J,D)`. Expand into ordered quadruples. An isolated support
has a vanishing expectation and commutes past the other even factors.
A surviving support graph is either two connected pairs or one connected
four-vertex graph. The three pairings contribute at most3(JD)². Each
connected graph contains one of the sixteen labelled spanning trees;
successive leaf summation bounds each tree by JD³. Repeated term indices
are included as overlapping vertices. This is an upper bound, without
an assumption that the terms commute.

For the actual parent perturbation there is a sharper disconnected term.
Group each bond and each plaquette into one local summand. Their vectors
X_i Omega are pairwise orthogonal by the exact CAR/Haar packet argument.
For a support graph consisting of two disjoint pairs, factorization is
therefore zero unless each pair repeats its local summand. Summing the
three possible pairings gives at most3z², where z=||V Omega||² is the
exact variance in (6). The connected bound is unchanged. Thus the code
uses `3z²+16JD³`. The same argument for V_t uses its exact variance s
and gives `3s²+16J_t D_t³`. No Wick rule is imposed on overlapping terms.

For parent V, `J=30B+P/7` and `D<=3(1272/7)`, using the actual at-most
three-cell supports and the inherited six-bond/twelve-plaquette incidence
bound. Let V_t be the parent perturbation with the packet weights in
(8). Its bond and plaquette norm bounds are respectively

```
b_t=2 sum_r d_r t_h,r,       p_t=(sum_r d_r t_l,r)/105,
J_t=B b_t+P p_t,             D_t<=3(6b_t+12p_t).             (9)
```

For the bond bound each native block has rank-two Dirac hopping, with
second-quantized norm at most2d_r. These parent operators are bounded
on the fixed finite spatial graph with unrestricted compact link
coordinates. Omega is the same product vector used above.

In that parent algebra V and V_t commute: different native-block
quadratic hopping operators act on disjoint CAR modes, coordinate
coefficients commute, and all magnetic terms are scalar multiplication.
The finite cutoff need not preserve this commutation. It does retain
V_t Omega, and hence

```
||V_finite V_t Omega||²
 <= ||V_parent V_t Omega||²
 <= sqrt((3z²+16JD³)(3s²+16J_t D_t³)) =: M.               (10)
```

The second inequality is Cauchy-Schwarz applied to V²Omega and V_t²Omega
after using the parent commutation. The code rounds the final square
root outward and takes its minimum with the independent bound J²s from
`||V||<=J` and `||V_t Omega||²=s`. Let C=-N/420. The exactly known projection of V V_t Omega
onto the trial space has squared norm `q²+g||C t_l||²`. Moreover H0
maps every packet into that trial space. Thus an omitted-component
bound is `R_out²=M-q²-g||C t_l||²`, and the full dimensionless residual
of (8) is bounded by

```
||(Hhat-E00-e)psi||²
 <= [(q+e)²+(d_lo+e)²s+R_out²]/(1+s) =: rhat².             (11)
```

The spectral theorem places some full physical eigenvalue in
`[(E00+e-rhat)/a,(E00+e+rhat)/a]`. The same theorem bounds the trial
spectral weight outside a radius R by `rhat²/R²`. These estimates retain
the complement; they do not identify the enclosed level as the ground,
a first excitation or a pole. At n0 the connected-quadruple terms in (10)
dominate its exact disconnected covariances, and the global-norm
alternative is smaller. The reported n0
source-energy residual upper bound is approximately102.09, around
Rayleigh energy -723.77. No narrow or isolated branch is certified.
A sharper connected moment or resolvent estimate is required to turn
this trial into such a branch.

The trial spectral measure differs from the actual polynomial state's
instrument records. The ground upper bound is a property of the same H
and may be used in an independently controlled actual-state energy
readout. It supplies neither an original probe residue nor a substitute
preparation. The selected reference, its clock and every original record
remain those of that readout. Neither a neutral pair threshold nor a
single-particle rest mass is inferred from the variational enclosure.
