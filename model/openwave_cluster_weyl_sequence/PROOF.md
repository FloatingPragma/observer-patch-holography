# Actual breakup spectrum of the complete native strong source

This result concerns the declared infinite one-dimensional strong Hamiltonian,
with every colour/relative channel and the common correlated triquark/meson
rule retained. It proves an essential-spectrum inclusion at the native Lambda
rest energy. It does not prove an outgoing weak amplitude or a lifetime.

## Exact separated-cluster intertwiner

Let J_sep be the CAR product of a colour-singlet three-quark vector and a
colour-singlet quark-antiquark vector supported in disjoint intervals, with
more than three empty separating cells. The meson has even fermion parity.
The ordered CAR product is an isometry from the tensor product of these two
fixed-region sectors; the spatial regions determine the cluster assignment
uniquely, so no indistinguishable-quark norm or exchange term is discarded.

Every count, local square, gradient, pair term and local/bond completion
factorizes between the clusters. Single-particle hopping acts in its own
cluster and retains all successors. A compact triple annihilator can act
only in the baryon; a compact quark-antiquark annihilator can act only in the
meson. The full inverse-rest creation matrix and its shifted difference
therefore factorize too, including every spin channel and Pauli sign. The
nuclear three-triplet term vanishes at total occupation5. Between the two
singlets the total prefix colour generator is zero, so every electric gap
Casimir is zero; within either cluster the other cluster's total generator
annihilates the state. Hence, on these supports including their successors,

    H5 J_sep = J_sep (H3 tensor I + I tensor H2).

This identity is local to separated supports; H5 is not replaced by a free
operator on overlapping clusters. The executable independently applies every
full source component to two onsite/split native products and verifies exact
rational equality to the tensor action, including all outside outputs. The
term-by-term argument above proves the identity for arbitrary finite supports.

## Uniform confining-tail shortcut

The existing finite Schur theorem bounds the full proton and pion ground
energies E_p(k), E_pi(-k) on the certified shell bracket. At any fixed momentum,
the relative diagonal has finite-dimensional occupation blocks and grows
linearly with span. Its resolvent is compact; bounded hopping and the compact
correlated term preserve compact resolvent. Thus the bottom energy has a
normalized eigenvector. Choose any vector in its ground spectral subspace
at the selected k*. Ground simplicity or a unique eigenvector is unnecessary. For the proton,
the native Q_R subspace of occupied span at least R>=2 obeys

    Q_R H3(k) Q_R >= C_p + (4 kappa/3)(R-2).

This follows from the same split/three-site colour energies used in the
existing complete complement proof. Every extra span costs4 kappa/3; the
bounded hopping and local-spin lower bounds are unchanged. The pion's two
relative tails have the identical slope, with their existing floor C_pi.
Compact correlated propagation vanishes in either Q_R. Exact onsite
variational energies give uniform E_upper for each cluster over the entire
small shell bracket, with positive denominators C-E_upper. No root-finder
value is substituted for the unknown exact on-shell momentum.

Let psi(k) be a normalized isolated ground eigenvector and q_R=||Q_R psi||.
Only nearest-neighbour hopping crosses the span boundary. Its norm is bounded
by b=6t for the baryon and b=4t for the meson. Applying the exact Q_R resolvent
to the eigenvalue equation gives

    q_R <= b q_(R-1) / [C+(4 kappa/3)(R-2)-E_upper],   q_1<=1.

The products decay factorially. The executable rounds their upper bounds
outward at32R bits; at most R*2^(-32R) is added, which also tends to zero.
No high relative state is discarded from the Hamiltonian.

Let phi_R=P_(span<=R) psi / ||P_(span<=R) psi||. It is a finite-dimensional
relative vector, defined by the full isolated spectral projection rather than
an unverified finite-cutoff eigenvector. The complete commutator estimate is

    ||(H(k)-E(k))phi_R|| <= 2b q_R / sqrt(1-q_(R+1)^2).

Both sides of the hopping boundary are retained; all other terms commute
with this span projector. The exact projection defines the coefficients of
the mathematical family. The executable certifies its errors; it does not
numerically enumerate every coefficient or compile this spectral service
into a finite observer preparation. The existing isolated-spectrum proof and
finite-support core also give convergent finite-matrix constructions of the
same vectors; no physical preparation claim is needed for Weyl's criterion.

## Normalizable Weyl vectors and actual source spectrum

For each cluster, multiply its finite relative vector by
L^(-1/2) exp(ikx) on x=1,...,L, zero outside, using its leftmost-cell translation
coordinate. This is a normalized vector on the infinite chain. One lattice
shift changes the box envelope by squared norm2/L. The sum of the two
translation-block norms is bounded by6t+4||J_tri|| for the baryon and
4t+4||J_meson|| for the pion: each directional bare hopping row/column has at
most N entries t, and the complete compact inverse-rest blocks are retained.
The boundary residual is therefore at most T sqrt(2/L).

Use momenta k and -k, and translate the pion box by D=4(L+R+3). Both clusters
and every one-step successor remain separated. Their CAR product has exact
norm one and full H5 residual at energy E_p(k)+E_pi(-k) bounded by

    sum_clusters [2b q_R/sqrt(1-q_(R+1)^2) + T sqrt(2/L)].

Set L=2^(4R) and let R tend to infinity. The bound tends to zero. Meanwhile
the antiquark lies beyond D, so the overlap with every fixed finite occupation
vector is eventually zero. Density of finite-support vectors proves weak
convergence to zero. These are normalizable singular Weyl vectors of the
actual interacting H5, not finite-ring states or delta-momentum vectors.
Weyl's criterion therefore puts the complete pair-band energy in the essential
spectrum. Total momentum concentrates around zero; no fixed total-momentum
fibre is substituted for this normalizable full-space construction.

The previously certified strict endpoint inequalities place the entire native
Lambda rest enclosure between the pair-band energies at half-angle tangents
7651/40960 and11477/61440. Continuity of the full isolated bands supplies at
least one k* in this bracket with E_p(k*)+E_pi(-k*)=E_Lambda. Applying the
construction at that exact, implicitly specified k* proves

    E_Lambda belongs to the essential spectrum of the complete strong H5.

The executable evaluates exact uniform error certificates for R=4,8,16; it
does not claim that a numerical approximate root is the exact energy shell.
All GeV/momentum conversions and strong/strange coefficients are inherited
calibrations. No new independent observation is used.

For every normalized member, the spectral measure outside an energy window
of half-width w is bounded by residual²/w², with its entire complement
retained. This energy instrument commutes with H5 and has zero nonselective
work, including arbitrary idle references. The existence of this breakup
spectrum and the separate positive local-current response do not imply a
nonzero on-shell weak spectral density. An outgoing wave operator and its
actual weak overlap, complete weak dynamics and physical3+1 identification
remain open.
