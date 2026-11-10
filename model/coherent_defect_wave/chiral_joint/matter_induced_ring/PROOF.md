# A matter-generated two-link ring on the complete periodic source

## Source and preparation

The Hamiltonian is the actual unit-parameter
[periodic lepton/Higgs source](../periodic_spectrum/PROOF.md), with two
six-state oriented quantum links, two five-state quantum scalars and twelve
fermion modes. Its full active dimension is 3686400. Each cell retains both
charged-lepton Yukawa tensors and scalar energy; the two edges carry the
actual coefficients `+i sigma3/4` and `-i sigma3/4`. Link energies have
coefficient one. Quarks are empty reducing spectators at `y_u=y_d=0`.
No coupling or clock is fitted to a measured quantity.

The declared preparation is one native spin-up e^c occupation at site 1,
with both scalars and both links blank. A static weak/color singlet reference
of hypercharge -1 at site 1 compensates the source charge. Its lepton number
and spin projection are opposite to the source values, and its fermion parity
matches the odd source parity. The complete source-plus-reference preparation
can therefore preserve these charges. The compensator has declared zero
energy. This is a specified pure preparation, with no ground-state claim.

The complete conserved source block has boundary charges `(0,-1)`, lepton
number -1 and `2Sz=1`. Its diagonal Gauss weight space has dimension 345;
the exact joint kernel has dimension 72. The independent character count,
exact kernel, full Hamiltonian closure and all 24 endpoint current equations
are checked. Other conserved blocks have zero probability in this preparation;
the complete readout assigns them to the other-states outcome.

## Complete retained and eliminated spaces

P fixes the matter occupation to the initial one-fermion mask 64 and retains
every compatible scalar and link configuration. Q is its orthogonal
complement in the complete 72-state physical block. The matter projector
commutes with every native current on the full CAR domain: at site 0 it is
the matter vacuum projector, and at site 1 it is the occupied spin-up e^c
projector with every other matter mode empty. Scalar and link factors are
identities. `projection_current_checks()` verifies all 24 local statements.

P has dimension five and Q has dimension 67. In the retained orthonormal
basis supplied in the receipt,

`A = PHP = diag(0,3,2,3,2)`.

The five states are the blank configuration, a scalar/link singlet on edge
10, the e^c closed loop, a scalar/link singlet on edge 01, and the L closed
loop. The scalar/link states retain the actual charged scalar alternatives.
Discarding these states changes P. Every hopping and every Yukawa component
has zero P-to-P block; A contains only electric and scalar number energies.

Write the full metric-self-adjoint block operator as

`H = [[A,C],[B,D]]`, `C=B^T G_Q`, `D^T G_Q=G_Q D`.

The physical basis has a nontrivial exact Gram matrix. The implementation
uses its complete P and Q blocks, normalizes P explicitly, and checks their
orthogonality. It never treats the rational coordinate basis as orthonormal.
The full physical embedding has independent columns, so `G_Q` is positive
definite. All 67 Q states enter the exact inverse.

## Gauge-invariant loop operators

For the actual e^c block of dimension one and L block of dimension two,
define

`R_B=sum_(a,b in B) U01_ab U10_ba`, `A_B=(R_B+R_B^dagger)/d_B`.

These are the original raising matrix units of the two six-state quantum
links, with all link alternatives retained. They are not unitary classical
link matrices. Acting on both links blank, `A_e` creates the normalized
e^c loop, and `A_L` creates the normalized sum of four L-loop configurations.
Each operator commutes with all twelve actual gauge currents at both
endpoints. The [complete loop operator checks](../periodic_boson_records/PROOF.md)
are replayed, including their norm-one identities and finite instruments.
The L loop carries the full native doublet representation, including its
hypercharge action. No identification from center labels is used.

The periodic graph has two vertices and two oriented edges. This closed
word is its noncontractible two-edge ring. A two-dimensional plaquette is
not present in this observation domain.

## Exact energy-dependent elimination

Eight unit-norm Yukawa monomials at each cell contribute a norm bound 16.
Four number projectors contribute 4. Forty actual quarter-norm hopping
monomials contribute 10. Thus `||H||<=30` on the full active carrier, for
all fermion occupations and scalar/link states. Compression gives the same
bound for D in the physical Gram metric.

For every real `z<-30`, `D-zI` is positive definite with spectral lower
bound `-z-30`. Exact elimination gives

`H_eff(z)=A-C(D-zI)^(-1)B`, `F(z)=H_eff(z)-zI`.

The code solves the symmetric rational system
`(G_Q D-z G_Q) X=G_Q B`, checks the residual exactly, and verifies on all
five P columns

`(H-zI) (I;-X) = (F(z);0)`.

Here the retained coordinates are embedded with their actual normalization.
This identity retains the complete eliminated state. The five-by-five
effective matrices at `z=-32` and `z=-64` are rational and included in the
receipt. Both blank-to-loop entries are strictly negative at both settings,
and each entry differs between the two settings. Their certified distances
from the eliminated spectrum are respectively 2 and 34. A has both entries
exactly zero. The loop couplings therefore exceed electric-only behavior
on this retained space and depend on the resolvent energy.

For `|z|>30`, the convergent expansion is

`H_eff(z)=A+sum_(n>=0) C D^n B / z^(n+1)`.

The relevant exact coefficients are:

| Q internal powers n | Blank to e^c loop | Blank to L loop |
| --- | --- | --- |
| 0 | 1/16 | 0 |
| 1 | 1/16 | 0 |
| 2 | 23/16 | 3/8 |
| 3 | 113/16 | 9/4 |

The leading entries are consequently `1/(16z)` and `3/(8z^3)`.
These asymptotic coefficients are distinct from the full rational values
evaluated at a specified z. Exact inversion, rather than a discarded
perturbative remainder, establishes their nonzero values at the two settings.
If the series is truncated after internal power n, the operator remainder
is bounded by

`900*30^(n+1)/(|z|^(n+2)*(1-30/|z|))`.

The bound follows from `||B||,||C||<=30`, in their physical metrics. It is
conservative and is not used to replace the exact rational calculation.

## Irreducible histories and open-edge terms

Each coefficient `C D^n B` inserts Q after every intermediate action.
Irreducibility means that the history has no intermediate return to P.
Same-edge return paths and all local self-energy terms are retained in the
complete matrices. They do not account for the off-diagonal blank-to-loop
entries by themselves: either loop changes both link registers.

`exact_parts()` decomposes every column of the unchanged Hamiltonian into
its diagonal term, its two actual oriented hops, and its two local Yukawa
terms. It proves exact closure of each part and equality of their sum to H.
The minimal e^c-loop history contains one hop on each edge. Removing either
edge eliminates its second-order matrix element. This is a product of the
two lower-order hopping operators; a primitive magnetic loop term is not
inserted into H. Their composition through Q supplies the closed-ring kernel.

For the L loop, every nonzero minimal four-factor history contains one hop
on each edge and two actual quantum Yukawa factors. Their exact signed sum
is 3/8. The lower Q-irreducible coefficients vanish. No intermediate return
to P accounts for this coefficient. The receipt lists every nonzero minimal
history, with chronology and coefficient. The quantum scalar and virtual
matter pairs are retained throughout this calculation.

## Complete finite-time records under the original Hamiltonian

The energy-dependent kernel is not used as a time-evolution generator.
Instead evolve the declared preparation by the unchanged `exp(-itH)` at
`t=1/16384` in its source clock. The complete three-outcome projective
instrument is the normalized e^c-loop state, the normalized L-loop state,
and their orthogonal complement. Their ranks on the prepared conserved
block are 1, 1 and 70. The pointer has three orthogonal records. The usual
controlled cyclic shift by these orthogonal projectors is a unitary
recording operation; all unnormalized postmeasurement states are retained.

With initial state psi and either final loop state ell, the exact moments
`<ell,H^n psi>` for `n=0,...,4` are

`e^c: (0,0,1/16,3/16,329/128)`,

`L: (0,0,0,0,3/8)`.

Every power uses the complete 72-state H. In particular,

`||H^4 psi||^2=238254041/8192 < 171^2`.

The degree-four Taylor state has full-domain error at most

`delta=30*171*t^5/(120*(1-30*t))`.

This follows by bounding all higher powers using `||H||<=30`; the bound
does not vanish after compression. Triangle inequalities give amplitude
intervals with strictly positive lower endpoints:

`e^c lower = t^2/32 - 329*t^4/3072 - delta`,

`e^c upper = t^2/32 + t^3/32 + 329*t^4/3072 + delta`,

`L lower = t^4/64-delta`, `L upper = t^4/64+delta`.

Squaring gives rational probability enclosures for both loop records.
The other-states probability is one minus their sum, with outward interval
endpoints. Its positive probability and its entire unnormalized state are
retained. No selection probability is divided out. The electric/scalar-only
diagonal A gives zero loop arrival from the initial state at every time;
the full H gives the certified positive records.

## Scientific boundary and relation to ring constructions

The result is a matter-generated gauge-invariant kernel and a finite loop
record on one declared interacting quantum source. P retains a charged
matter occupation and its reference. The Schur coefficient depends on z.
No autonomous energy-independent pure-gauge Hamiltonian is derived. The
two-edge ring does not supply a spatial plaquette action, a Coulomb phase,
a relativistic photon or a physical gauge coupling. No color spectrum is
included in this lepton-reduced source. Results do not transfer to a
119-state link with additional Clebsch-Gordan transitions.

[Hermele, Fisher and Balents](https://arxiv.org/abs/cond-mat/0305401) derive
ring dynamics by degenerate perturbation theory in an easy-axis constrained
spin model and analyze a three-dimensional phase with further hypotheses.
The shared method is elimination of intermediate states. Their constraint,
reference sector and phase assumptions are not supplied by the finite
periodic construction here, so their phase or photon conclusions are not
transferred. The [open-edge elimination](../field_limit/PROOF.md) on this
source family gives quartic virtual interactions on a different projection;
its one-edge domain has no closed ring. The source loop operators used here
contain no added magnetic term.

## Replay and custody

The receipt contains exact rationals only. Its check mode validates every
consumed receipt pin, computes the transitive Python input closure, binds
the package's own code and proof, and independently recomputes the complete
result without writing. Tests cover all retained/eliminated states, full
current equations, irreducible ancestry, exact reconstruction, rejected
resolvent settings and the complete positive probability records.
