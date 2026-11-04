# Protected compact-texture records and local band transport

## Declared source and complete record law

The source is the compact opposite-degree texture of
[compact_defect_pair](../compact_defect_pair/PROOF.md), on an open cubic lattice
of radius `L>=a+2`, with half-separation `a>=1`, hopping norm below two and mass
amplitude `M>=8`. Write `f_x` for its unnormalized direction and retain the two
zero cores. Noncore texture variables `v_x` are exact rational vectors with

```math
\|v_x-f_x\|\le\eta\|f_x\|,\qquad 0<\eta<1/4.
```

The texture is a classical control register. Outside its compact support it
equals north. Inside, the fermion Hamiltonian uses the actual normalized
direction `v_x/|v_x|`. A proposal replaces one noncore vertex. It is accepted
only if it satisfies the displayed budget, is nonzero, has nonnegative scalar
product with the old vertex, and has nonnegative scalar product with every
other vertex of every incident shell triangle. A proposal affecting the shared
face of the two unit shells is checked on both shells. Every rejected proposal
is retained as an identity transition. The cores and exterior are fixed by the
declared law.

For the displayed source `a=1`, `eta=1/512`, the four equiprobable proposals
replace the shared-face vertex `(0,0,0)` by `(eta,0,-1)`, `(-eta,0,-1)`, or
zero, or replace the upper core by north. The first two are accepted and the
last two are rejected. These proposals define a continuing three-state texture
process: the original texture, the positive perturbation and the negative
perturbation. Its row-stochastic transition matrix is

```math
T=\begin{pmatrix}1/2&1/4&1/4\\0&3/4&1/4\\0&1/4&3/4\end{pmatrix}.
```

All twelve source-state/proposal records enter the instrument. Its stationary
texture law is `(0,1/2,1/2)`; this statement concerns the classical texture
register and does not specify a stationary quantum field state. The proposal
probabilities, guards, gate controls and pulse services are simulator design.
No measured mass determines them.

## Exact winding preservation on every record

Every admitted shell triangle starts with three nonzero vectors whose pairwise
scalar products are nonnegative. During a one-vertex replacement, its normalized
interpolation is the normalization of a convex combination of the two unchanged
vertices, the old vector and the proposed vector. Every pair in that set has
nonnegative scalar product by the guard. If the nonnegative coefficients are
`t_i` and at least one is positive, then

```math
\left\|\sum_i t_i v_i\right\|^2
=\sum_i t_i^2\|v_i\|^2+2\sum_{i<j}t_it_j\langle v_i,v_j\rangle>0.
```

The normalized interpolation therefore gives a continuous homotopy of the
complete oriented shell. Its degree is invariant. Adjacent triangles use the
same updated vertex, so their restrictions agree throughout the homotopy.
Rejection changes no data. Induction gives degrees `(-1,+1)` on every complete
finite history, with a constant degree-zero exterior, regardless of the proposal
schedule. The certificate independently recomputes the signed rational
regular-value count on every triangle after every displayed record.

This is protection by an explicit guarded repair law. Removing its restrictions
permits setting every noncore direction to north, giving shell degrees `(0,0)`.
The zero core amplitudes and four-mode spectral count survive that erasure.
Spectral isolation therefore does not itself protect winding. Neither winding
degree is identified with a fermion-number operator, and fermion occupation
does not create a texture.

## Uniform band and uncontrolled leakage

For every admitted texture, the onsite potential has exactly four zero
coordinates and every other eigenvalue is `+M` or `-M`. Weyl bounds and the
compact-source Schur argument therefore apply without a change of constants:

```math
\dim\operatorname{ran}P_v=4,\quad
|E|\le r={b^2\over M-4},\quad
\operatorname{spec}(H_v|_{1-P_v})\cap(-g,g)=\varnothing,\quad g=M-2,
```

where `b^2=7/4` for `a=1` and `b^2=3/2` for `a>=2`. The core-probability bound
`[1+b^2/(M-2-r)^2]^-1` is unchanged. These statements hold in every containing
cube and on the infinite lattice by the same exterior-resolvent argument.
No angular topological assumption enters this band count.

Normalization obeys `|v/|v|-f/|f||<=2|v-f|/|f|<=2 eta`. Hence any two admitted
source Hamiltonians satisfy `|H_w-H_v|<=epsilon=4 M eta`, independent of volume.
Take the circular resolvent contour of radius `b=(g+r)/2` about zero and put
`d=(g-r)/2`. Every point on it is at distance at least `d` from both spectra.
The resolvent identity gives

```math
\|P_w-P_v\|
\le {1\over2\pi}\int_{|z|=b}
       \|(z-H_w)^{-1}(H_w-H_v)(z-H_v)^{-1}\|\,|dz|
\le \beta:={b\epsilon\over d^2}.
```

For `M=64`, `a=1`, `eta=1/512`, the exact bound is
`beta=3572880/221206129<1`. An instantaneous texture quench without transporting
the state has leakage probability at most `beta^2`. After `n` such updates and
arbitrary intervening evolution under each actual source Hamiltonian, the
triangle inequality bounds leakage amplitude by `n beta` and probability by
`min(1,n^2 beta^2)`. The bound holds path by path, so averaging or retaining all
classical outcomes preserves it. The numerical source controls exhibit nonzero
quench leakage; an unchanged field state is not exact band transport.

## Exact band transport and the same currents and pairing

For two consecutive band projectors `P` and `Q` with `|P-Q|<1`, define

```math
S=QP+(1-Q)(1-P),\qquad
U=S\,[S^*S]^{-1/2},\qquad S^*S=1-(P-Q)^2.
```

The positive factor is invertible and commutes with `P`. The identity `SP=QS`
then proves that `U` is unitary and `UPU*=Q`. Apply this supplied control after
an accepted texture update and the identity after a rejected update. If a one-body
state lies in the admitted band at one record endpoint, it lies in the target
band at every later endpoint, including arbitrary intervening source evolution.
This is an exact endpoint statement. Individual gates in a compiled control
need not preserve the instantaneous band.

The literal source has particle-hole identity `H C=-C conjugate(H)`, where
`C=D_stagger i sigma2` is fixed, antisymmetric and unitary. Its symmetric low
spectral window gives `P C=C conjugate(P)` and similarly for `Q`. Every operation
in the polar formula respects this antiunitary covariance; therefore

```math
U C=C\overline U,\qquad U C U^T=C.
```

Extend by the actual native fifteen-state matter module and its three declared
family copies. The spatial control is `U tensor I45`, so it commutes with all
twelve uniform native currents. If `W` is the original four-column band frame,
the transported frame `W'=UW` satisfies

```math
W'^*H_wW'=A_w,\qquad
W'^*C\overline{W'}=W^*C\overline W.
```

Thus the native-current inclusion and the local pairing
`Delta=C tensor conjugate(K)` use the same transported states. The coefficient
`A_w` is computed from the changed Hamiltonian; it is not dropped or assumed
equal to `A_v`. Second quantization preserves parity, carries
`N_P=dGamma(P tensor I45)` to `N_Q`, and consequently carries
`lambda(N_P-90)^2` to `lambda(N_Q-90)^2`. This is the operator attachment used by
the [compact interacting response](../compact_interacting_response/PROOF.md).
It does not turn the controlled deformation into free particle motion.

The spectral and leakage bounds concern the spatial one-body operator. On Fock
space the exact map carries the corresponding localized CAR algebra and its
complement, retaining all modes. It supplies no full interacting ground-state
selection or many-body spectral gap.

The field part of each full record has Kraus operator `sqrt(q_r) Gamma(U_r)`
and a classical transition from its old texture to its target. Since every
`U_r` is unitary and the complete proposal probabilities sum to one, every
source texture has `sum_r K_r* K_r=I`. All accepted, rejected and repeated-state
records are retained. There is no selection of a successful subensemble.

## Finite local pulse realization and its work

The lattice sites admit a snake ordering in which adjacent sites share a bond.
Interleave the two spin modes at each site. Adjacent coordinates of the resulting
`N=2(2L+1)^3` path are either onsite or nearest neighbors. Adjacent-row Givens
elimination reduces any `N`-dimensional unitary to a diagonal phase matrix using
at most `N(N-1)/2` two-mode gates. Reversing those gates and applying the phases
reconstructs the original unitary. Each two-mode unitary is generated by an
onsite/neighbor number-preserving quadratic Hamiltonian; uniform replication
over the 45 internal slots preserves every uniform native current throughout
the control. The producer replays every gate for all three nontrivial texture
pairs on the complete `N=686` source. The reverse transition uses the adjoint
sequence. No spatial or internal modes are discarded.

This compilation assumes individually addressed quadratic pulses, their timing
and cancellation of the source drift during those pulses. It is not a compiler
into the unchanged C/E repair alphabet. Pulse precision is an additional input:
if each completed one-body branch approximates its ideal `U_r` within operator
norm `zeta`, `n` controlled records have one-body band-leakage probability at
most `min(1,n^2 zeta^2)` and complete one-body record-state trace-distance error
at most `min(2,2n zeta)`. On a fixed `k`-particle sector, telescoping the `k`
unitary tensor factors before antisymmetric restriction gives the bound
`min(2,2nk zeta)`. A full-Fock bound with occupation-changing pairing uses
`k=45N` unless the supplied control error is bounded directly on Fock space.
The numerical replay residuals check an implementation; they are not interval
certificates for physical pulse errors.

The source Hamiltonian changes during a controlled texture deformation. Even
ideal transport need not conserve its expectation. On the one-body band both
endpoint Hamiltonians have norm at most `r`, so the absolute endpoint work is
at most `2r`. Its value is computed from `W^*(U^*H_wU-H_v)W`. The gate drive and
texture controller supply this work. No autonomous texture formation, moving-core
dispersion, physical species, local Gauss dressing or conserved gravitational
stress is inferred. The construction closes the finite guarded-protection and
controlled-state-transport requirements on this declared law.
