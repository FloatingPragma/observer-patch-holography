# The actual periodic quantum reference and charged records

## Complete carrier and common Hamiltonian

This package calculates the two-oriented-link source defined by
`two_cell.symmetries.periodic_column`. The basis is
`(mask,link01,link10,H_left,H_right)`: twelve lepton CAR modes, two six-state
quantum links and two five-state Higgs memories. Its active dimension is
`2^12*6^2*5^2=3686400`. Empty quark occupations and unexcited colored link
alternatives remain reducing spectators at the declared zero quark Yukawa
coefficients. This is an active lepton/link/scalar control, not the complete
R15-family interacting ground.

The Hamiltonian contains both local quantum charged-lepton Yukawa terms,
both scalar number terms and both excited-link projectors, all with unit
coefficients. Its hopping is

`(i/4)c_left^dagger sigma3 U01 c_right
 -(i/4)c_right^dagger sigma3 U10 c_left + h.c.`.

The full quantum link matrices are raising operators, not classical
transporters. The source clock is `exp(-itH)` and `||H||<=30`. These are
precisely the existing periodic completion's terms and boundary signs.
Its spectrum is calculated independently of the open-edge spectrum.

## Full Gauss projection and translation reduction

The local constraints are the actual `periodic_current_column` generators.
At a cell the matter and scalar currents are supplemented by its endpoint
current from each incident quantum link, plus any declared static boundary
charge. Both opposite link orientations are included.

The basis construction first enumerates every simultaneous local
hypercharge and weak-Cartan weight-zero state. It then applies the actual
weak raising operators into their complete adjacent-weight image, without
compressing that image back into weight zero. Exact rational RREF gives
the simultaneous kernel. On weight zero, the SU(2) commutator and adjoint
identity imply that raising-annihilated vectors are also annihilated by
lowering. Hence these vectors satisfy every weak Gauss generator. Color
acts trivially on the specified active lepton module.

An independent SU(2) by SU(2) character difference checks every nullity:

`dim singlet = dim(0,0)-dim(2,0)-dim(0,2)+dim(2,2)`,

where the entries denote twice the local weak-Cartan weights. The exact
neutral physical dimension is 3424, from 18296 diagonal-weight states.
The right-negative and right-positive boundary sectors each have 2456
physical states, from 12846 diagonal-weight states. Every L and 2Sz sector
is included. The largest neutral physical block has dimension 502; the
largest charged block has dimension 330. The exact constraint components
used by RREF have at most 36 columns.
Here q is the dimensionless U(1) hypercharge of a weak-singlet apparatus
reference. This notation does not assume electroweak symmetry breaking or
identify the excitation with an observed electron.

The coordinate basis B has an explicit row-selection left inverse E and
positive Gram matrix `G=B^dagger B`. Use the gauge-invariant diagonal phase
`D=i^(N_left)` to make H real; the physical embedding is **D B**. Exact
identities verify `H_real B=B H_coordinate` and Gram-Hermiticity. Every
outgoing full-H column is retained.

For the neutral domain the actual periodic translation T preserves every
Gauss constraint and commutes with H. It obeys `T^2=(-1)^L`. In the real
frame, divide T by `i^(L mod2)` to obtain a rational involution S. The code
constructs both complete eigenspaces of S, checks their exact dimensions,
intertwines H into each, and retains their physical embeddings. The largest
502-dimensional neutral block splits into two 251-dimensional translation
blocks. This is an exact reduction, not a restriction to a selected state
or a guessed low-energy branch.

## Complete certified energy spectrum

Every generalized Hermitian energy block is further divided only into
exact connected components of its joint H/Gram support. Their direct sum
is the full physical space. The complete rational eigenbasis certificate
from `two_cell.spectrum.exact_dyadic_eigenbasis_error` is reused on each component:
its exact Gram estimate proves invertibility, and its full integer residual
bounds every ordered eigenvalue by Hermitian minmax. Numerical eigensolvers
only propose the basis during explicit `--generate-witnesses`. The compressed
integer artifact `eigenbasis_witnesses.npz` stores whole-block proposals for
detector consumers and complete proposals for every connected component.
Its manifest binds each proposal to exact rational K/G entries, dimension,
boundary and conserved labels, translation sign and component coordinates.
The receipt hashes the artifact. Both `--write` and `--check` read those
same dyadic numerators and repeat the exact checks without an eigensolve.
Missing, altered or mismatched witnesses fail replay. No ground conclusion relies on a Ritz value or on
a vacuum-generated Krylov subspace.

All exported energies have dyadic centers on the `2^-32` grid and outward
radius `2^-26`. Each interval is checked against its exact rational error
bound before export. The complete neutral ground is unique, has
`(L,2Sz,T)=(2,0,+1)`, and lies in

`[-18632299339/4294967296, -18632299211/4294967296]`,

approximately `-4.33817023295`. Its gap to every other neutral physical
state is at least

`1106815801/4294967296`, approximately `0.25770063535`.

The right-negative boundary ground window is in L=1, 2Sz=+/-1; the
right-positive window is in L=3, 2Sz=+/-1. Both lie in

`[-14465768193/4294967296, -14465768065/4294967296]`,

approximately `-3.36807410442`. Exact spin reversal extends to both
oriented links: complex conjugation reverses i and the spin half-turn
reverses sigma3, leaving both hopping terms and the real singlet Yukawa
terms invariant. It preserves each boundary sector and squares to
`(-1)^L`. The isolated odd-sector windows are therefore exact Kramers
pairs, rather than an inferred equality of nearby decimal eigenvalues.
These finite symmetries do not identify physical Lorentz spin.

## A same-boundary conjugation specific to this periodic source

Let A use the inherited complete fermionic map

`A c_site A^-1 = i c_other_site^dagger`,

conjugate coefficients, transpose each quantum link's two internal indices
separately, and exchange/conjugate the scalar memories. The raising link
stays raising. Complex conjugation and CAR reordering cancel the hopping
sign change on each oriented link. The local Yukawa term maps to its
adjoint at the other cell. Scalar and link energies are unchanged. Thus
A commutes with this H and reflects the boundary charge.

Compose with the actual translation: `C=T A`. Then

`C c_left C^-1=-i c_left^dagger`,

`C c_right C^-1=+i c_right^dagger`.

C commutes with H and sends a boundary pair `(q_left,q_right)` to its
negative at the **same** positions. It sends L to 4-L and 2Sz to -2Sz.
The code checks the actual link matrix-unit identities, occupied full-H
columns and all local Gauss generators. The inherited complete CAR
intertwiner and term identities prove the result on all occupations.

Because the neutral ground is unique, its density is C-invariant. Each
local e^c occupation therefore has mean one half. Same-site addition and
removal have identical complete spectral measures when their apparatus
charges and readout effects are transformed together. This conclusion
belongs to the periodic Hamiltonian. The open source lacked the required
translation and does not acquire this same-boundary equality.

## The actual ground preparation and charged instrument

Prepare `P/3424` on the entire neutral Gauss space with its retained graded
reference. A basis can be chosen with definite L, 2Sz, parity and T phase.
Give the reference opposite L and 2Sz, matching parity and inverse T phase.
Their purification is jointly Gauss neutral, has zero total conserved
charges, even total parity and translation eigenvalue one. A finite
controller inside that joint sector can prepare it. This is a declared
preparation service, not idle dynamics or a claim of efficiency.

The complete binary ground filter has success probability `1/3424` and
failure `3423/3424`. Its successful source state is the actual unique
periodic ground. The preparation reference and failure are retained.

Append the same fermionic reservoir and compensator as in the complete
open-source instrument. The swap acts on the twelve source CAR modes and
the apparatus, with both quantum links and scalars as spectators. Its
all-occupation current and CAR identities therefore tensor to this larger
source without modification. The even-ground Kraus branches are

`a^dagger/sqrt2, -n/sqrt2, -a/sqrt2, (I-n)/sqrt2`.

Each has original conditional-on-ground probability exactly one quarter,
by the proved half occupation. The two successful exchanges receive a
complete energy-window measurement in the appropriate recorded boundary
sector. Both unsuccessful exchanges are explicit **no-energy-readout**
records; their post-swap quantum states and apparatus records are retained.
This is a complete adaptive protocol, not omitted failure spectroscopy or
renormalization of the successful outputs.

Actual probe matrices obey `C_probe E_source=E_target T_probe` using the
returned phased physical embeddings on both ends. The add/remove maps
satisfy the exact CAR metric identity. Ground-vector errors are bounded
by the full residual and the proved complete gap. Target energy-window
errors use equal-rank spectral-projector bounds from the complete spectrum,
including Gram errors and all alternatives. The receipt reports outward
rational probability intervals; an interval containing zero is not evidence
that its channel is absent.

The first successful addition/removal transfer lies in
`[2083265509/2147483648,2083265637/2147483648]`, approximately
`0.97009612853`, in this source clock's units. Its original-swap weight lies
in `[29863037/134217728,14931519/67108864]`, approximately `0.22249697`.
The positive lower weight certifies an observable charged threshold on this
finite source. Each successful exchange has 328 complete energy windows;
the total original-swap success probability for each exchange is `1/4`.
These are finite response energies, not identified Standard Model rest
masses or GeV predictions.

## What the charged translation measurement means

For a charged success branch, retain a two-position register for its static
compensating charge and translate it together with the source. The joint
translation preserves the full Gauss kernel and commutes with H. Starting
with the reference at one cell, its translated position is orthogonal.
The declared reference translation is the ordinary position swap with
square I and no additional boundary phase. Its charge, fermion parity and
retained internal record are unchanged by that swap. This is an apparatus
convention, not a boundary condition derived for a physical moving particle.
The odd-source sector therefore obeys `T_joint^2=-I`, so the complete four phase
outcomes have conditional probabilities

`(0,1/2,0,1/2)`.

This remains true after **each energy projection**: H preserves the
reference position and commutes with the joint translation. Every
successful energy interval is therefore duplicated equally at phases
`+pi/2` and `-pi/2`, with the zero outcomes also retained.
For a zero-weight energy outcome this statement is the factorization of
its joint energy/phase weight, without conditioning on a zero-probability
event. The adaptive protocol has 2626 records conditional on ground
preparation, including the two no-energy-readout failures and all zero
phase records. Multiplying each by `1/3424` and retaining the ground-filter
failure `3423/3424` gives 2627 unconditional records.

This degeneracy is reference-position structure. It does not provide a
second squared momentum or an intrinsic dispersion measurement of the
charged source object. Both phases have squared momentum `pi^2/(4a^2)`.
The reference carries no kinetic Hamiltonian in this construction. A
source-relative kinetic interpretation needs a separate, controlled
reference and boundary argument. The same-H periodic energy and joint
phase records are meaningful finite outputs without satisfying that further
condition.

The common positive finite energy reference is `H-E0`. The complete
neutral and +/- boundary minima lie at or above E0, so it is nonnegative
on the corresponding source-plus-apparatus Gauss domain when the apparatus
Hamiltonian is zero. The driven swap's work supply is separate. A controlled
finite mass estimate requires additional same-source dispersion and
reference-error evidence; a full continuum is not a prerequisite for such
an estimate, but these two momenta alone cannot provide it.
