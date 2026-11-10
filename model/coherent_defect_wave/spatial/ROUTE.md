# Clifford routing with a controlled one-particle limit

The exact square-root SWAP has a quadratic spectator band. The construction
below is a separate declared routing law. It reuses the repository's actual
four-component Clifford matrices and executes their spatial shifts in
position space. Its free one-particle continuum approximation has an explicit
operator-error bound. A finite local coupling to the native defect-wave
unitary is also specified. The bound does not cover the coupled theory.

## Finite routing program

`dirac_route.py` consumes `gamma_matrices` from
`model/common_source/overlap_regulator.py`. All consumed entries are checked
to be Gaussian integers before conversion to exact arithmetic. Set
`β=γ_4` and `α_i=iβγ_i`. The sixteen exact relations are

\[
\alpha_i^\dagger=\alpha_i,\quad\beta^\dagger=\beta,\quad
\{\alpha_i,\alpha_j\}=2\delta_{ij}I,\quad
\{\alpha_i,\beta\}=0,\quad\beta^2=I.
\]

The source of these matrices supplies their convention, not an identification
between the routing coin and native physical spin. No calibrated internal
mass kernel, Wilson sign function or continuum field action is consumed.

On a finite periodic cubic lattice let `T_i` move a position ket by one site
in direction `i`. Put

\[
P_i^\pm=\frac{I\pm\alpha_i}{2},\qquad
S_i=P_i^+T_i+P_i^-T_i^\dagger,\qquad
D_a=e^{-ia m\beta}S_3S_2S_1.
\]

The two projectors are orthogonal and sum to the identity, so each `S_i`
is unitary. The code executes all three shifts on sparse position/coin
amplitudes. An inverse reverses the shifts and the local mass coin. The
receipt checks the exact inverse and all four Bloch columns, including
non-real translation eigenvalues. After one tick the support moves at most
one site in each coordinate and three in Manhattan distance.

The four-dimensional coin, cubic shifts and mass angle `a m` are declared
architecture parameters. The original native stationary triplet has three
colors. Reserving one as a background leaves two excitation colors, so it
does not by itself implement this four-component routing coin. Even a
two-component one-particle encoding needs a collision rule when several
excitations occupy the same cell.

## Quantitative free relativistic limit

Assign the lattice spacing and tick duration the common value `a` in units
where `c=1`. For a plane wave of physical momentum `p`, the shift symbol is
`S_i(ap_i)=exp(−ia p_i α_i)`. Hence

\[
D_a(p)=e^{-ia m\beta}e^{-ia p_3\alpha_3}
e^{-ia p_2\alpha_2}e^{-ia p_1\alpha_1},\qquad
H_D(p)=\sum_i p_i\alpha_i+m\beta.
\]

The exact Clifford relations give

\[
H_D(p)^2=(|p|^2+m^2)I.
\]

Write the four Hermitian summands as `H_r=h_r A_r`, with
`h=(p_1,p_2,p_3,m)`. The unitary Duhamel identity bounds the error for two
factors by `a²||[H_r,H_s]||/2`. Iteration over four factors and
`||[H_r,H_s]||=2|h_rh_s|` give

\[
\|D_a(p)-e^{-ia H_D(p)}\|
\le a^2\sum_{r<s}|h_rh_s|.
\]

Telescoping `n` unitary factors at `t=na` gives

\[
\boxed{\quad
\|D_a(p)^n-e^{-itH_D(p)}\|
\le t a\sum_{r<s}|h_rh_s|
\le t a(3P^2+3|m|P),\quad |p_i|\le P.
\quad}
\]

The bound holds uniformly over a fixed finite momentum window. Finite
periodic Fourier transformation makes it the same operator bound on all
wave packets supported in that window. Increasing the number of sites at
fixed physical box size enlarges its represented momentum range; increasing
the physical box size refines the momentum grid. An infinite-volume comparison
requires the corresponding box and momentum-window embeddings. The
receipt checks the bound against independent 4 by 4 matrix exponentials
for 50, 100 and 200 ticks at fixed duration.

This is a free Dirac one-particle limit with a supplied mass coefficient.
The compact-window estimate does not eliminate modes elsewhere in the
Brillouin zone, supply a chiral anomaly cancellation theorem, or select the
observed number of species. A free relativistic symbol gives no nontrivial
interacting four-dimensional continuum theorem.

## Local many-particle implementation and native coupling

A specific finite extension introduces four fermionic routing modes per
cell. Their full local Fock space has dimension sixteen and includes every
occupation sector. The canonical anticommutation relations are a declared
register algebra. Second quantization of each `S_i` gives its exact finite
many-particle unitary; no one-particle truncation is used in this definition.

In the local basis diagonalizing `α_i`, pair the two positive and two
negative lanes. For each pair, an intercell fermionic SWAP exchanges the
positive lane at `x` with the negative lane at `x+e_i`; an onsite SWAP of
the lanes then shifts positive occupations right and negative occupations
left. Distinct lane pairs are disjoint modes. The fermionic signs are those
of the second-quantized permutation, including double occupations. This
supplies a finite-range implementation and a collision convention instead
of treating multiple colors as distinguishable classical walkers.

Attach the full native link/wave Hilbert space to each such cell. Let
`P_occ` project onto its nonempty routing occupations. The local gate

\[
C_x=(I-P_{\rm occ,x})\otimes I_{
\rm native,x}+P_{\rm occ,x}\otimes U_x
\]

is unitary and commutes with the native interior constraint. The tick
`Γ(D_a) product_x C_x` is a fully specified finite joint operator. It
couples the routing occupations to the literal native defect-wave program
on the same local state rather than averaging its outcomes into a classical
potential. In the receipt, one routing fermion coherently occupies the
tested cell or an untouched reference cell, keeping total number and parity
fixed. Occupation of the tested cell controls `U` on its physical native
preparation. The reduced two-position routing-state purity is strictly less
than one, proving that the gate is not a product operation. This witness
does not require coherence between different total fermion parities.

This extension introduces neutral routing modes with respect to the native
star constraint. Matching their gauge transformation, spin and statistics
to physical matter requires a further representation attachment. Declaring
CAR is not a derivation of physical spin-statistics from protected defects.
The finite coupling also does not implement the square-root SWAP family;
they are distinct candidate spatial programs.

A fixed controlled `U` is not close to the identity when `a→0`. Therefore
the free routing error bound cannot be applied to the coupled tick by
discarding `C_x`. One explicit weak-tick alternative is

\[
K_x=\frac{U_x-U_x^\dagger}{2i},\qquad
\widetilde C_x(a)=e^{-ia\lambda P_{\rm occ,x}\otimes K_x}.
\]

Here `K_x` is Hermitian, has norm at most one, and commutes with the native
constraint. It is a fixed algebraic function of the actual native gate and
requires no branch choice for a unitary logarithm. The coupling `λ` and
weak-tick replacement are explicit architecture choices. The stationary
native background is annihilated by `K_x`; that fact alone gives no
interacting low-energy vacuum or nontrivial continuum response. A local
many-particle refinement requires a vacuum/state family and uniform
correlation, renormalization and error estimates for the complete joint
generator. The free bound is only its kinetic input.

## Relation to the particle problem

The construction gives an executable local alternative with anticommuting
spatial generators and a controlled relativistic free window. It explains
which additional register and coupling decisions are required to replace
the quadratic spectator propagation. The parameter `m` is an input; the
native defect-wave law does not determine its value through this argument.
A GeV mass requires the interacting pole or threshold relative to a stated
state and a common physical clock. No Standard Model species or GeV value
is assigned by this routing test.
