# A fixed-charge, two-doublet native core and its neutral bath

## Why a new interface is needed

The preceding primitive-spin impurity has a useful interacting radiating
branch, but that branch has exchange spin3/2. Further, primitive occupation
Q differs between the two orientations of its exchange spin1/2: [Q,S+]=S+.
It cannot be identified with the charged-lepton spin doublet of the current
row11 CAR/Gauss source, where Q is scalar on the entire spin doublet. An
isometry preserving Q and all spin matrices would preserve the commutator
and Casimir, so neither discrepancy disappears by changing basis. These are
restrictions on that direct identification, not on every reconstructed spin.

## Literal native fixed-charge spin code

Use M=2N+1 logical sites, each with two physical native rails. Define
E|0>_j=|10>_j and E|1>_j=|01>_j. Every logical site has exactly one primitive
occupation, independently of its spin orientation. Logical spin rotations
therefore commute with the primitive number on each pair. On the complete
four-state rail pair, extend SU(2) by spin zero on |00> and |11>; thus it
also commutes with charge on the entire local carrier. Physical spacetime
spin is not inferred solely from this internal representation.

A logical swap S_ij is the physical permutation R_ij swapping the two first
rails and the two second rails. Exactly R_ij E=E S_ij. A whole bath-rung
permutation swaps two logical sites and therefore consists of four disjoint
physical rail swaps. All these physical swaps are implemented by retained
native retries: on a requested seam, the recoded flag coin XOR n_first is
conditionally fair, identity at0 and the endpoint swap at1. No working-state
reset is used. The controller reconstructs loads from public preparation
and retained currents. The auxiliary zero pair supplies a retained fair
macro decision; decision0 holds, decision1 requests the physical swap word.

On the N2 census, the complete primitive total-five sector has252 binary
load configurations; all are included, with both passive counter values.
The32-dimensional one-unit-per-pair code is a selected invariant subspace
of the ideal full permutation source, not the whole carrier. The220 other
load configurations remain. Finite caps can enter them and do not discard
them. The arbitrary-N proof is on the whole binary carrier, before code
selection.

## Complete capped source and its positive transfer

For each requested physical swap cap retries at B and put delta=2^-B,
u=1−delta. A k-leg word has first-failure probabilities delta*u^j for
j=0,...,k−1 and successful probability u^k. Its partial physical permutation
R_j is the product of the first j disjoint swaps; R_0=I. The passive counter
is toggled only after COMPLETE success on an odd-type address. An even-type
address never toggles it. Every failed and partial outcome remains.

The complete one-macro transfer on all physical loads and memory is

  P=sum_(j=0..k−1) delta*u^j (I+R_j)/2
        +u^k (I+R_k tensor X^epsilon)/2.                    (1)

All coefficients sum to1. Each R_j and R_k tensor X^epsilon is a selfadjoint
involution. Thus (1) is positive, selfadjoint and stochastic on the full
source. It is not inferred from the positive ideal code restriction.
Its difference from the ideal transfer is at most1−u^k<=k delta.
Its mean primitive cost is1+u+...+u^k, including the initial decision and
all native retries. There are N odd two-leg internal-rung addresses, two
even two-leg impurity addresses, and N−1 even four-leg rung transports.
Hence W=2N+1, full per-address generator error at most6N delta, and ideal
expected primitive cost(8N+1)t for mean Wt macro marks. Finite cap errors,
random durations and all controller records remain separate from time units.

## Source core and full reference

In the counter-odd source sector H=W(I−P_odd). After taking the explicitly
controlled cap limit and restricting through E, it is exactly

  H=sum_r Pi_triplet,r + Pi_singlet,(s,a0)+Pi_singlet,(s,b0)
                       +(1/2)sum_r(I−R_(r,r+1)),           (2)

where the first terms come from odd addresses and the two impurity terms
from even addresses. No Hamiltonian coefficient is inserted in place of
the native count matrix. All scheduler weights are1.

The isolated core is Hc=Pi_triplet,ab+Pi_singlet,sa+Pi_singlet,sb.
The lower rung-singlet doublet has energy1/2. The rung-triplet J=1/2 doublet
has energy5/2. The J=3/2 quartet has energy1 and is retained. On the full
three-local-state permutation module (including all off-code rails), the
possible core energies are1/2,1,2,5/2: these are the trivial, sign and
standard S3 permutation representations. Equivalently Hc satisfies the
polynomial with these four roots and is Hermitian. Therefore Hc>=I/2 on
the complete carrier. All other terms in(2) are positive, and the impurity
spin times all rung singlets attains1/2. This is the actual ideal-source
reference energy in the executed total-charge sector, not a family-specific
subtraction. Finite caps need not fix it; their full reference error is
bounded and explicitly retained.

## Two spin-half families and all neutral products

With one rung triplet, the impurity couples to total spin1/2 or3/2. Relative
to the reference1/2, their exact restrictions are

  K_(1/2)=I+L_path/2+|0><0|,
  K_(3/2)=I+L_path/2−(1/2)|0><0|.                            (3)

Both follow from −S_s·T_0, with eigenvalues+1 and−1/2. The entire two sectors,
all magnetic substates and all remaining source states are retained.

Let l_s=|s>_imp tensor |singlet>_core. The higher doublet for spin up is
h_up=(|110>+|101>−2|011>)/sqrt6 in impurity/core-rung logical order. Choose
h_down by the actual total-spin lowering operation, fixing its phase.
Both l_s and h_s transform in the same fundamental SU(2) representation.
Their primitive core charge is exactly3−2=+1 in every spin component.
A bath rung contains two logical rail pairs, hence exactly two primitive
occupations in every singlet/triplet polarization. Its relative charge
is exactly0. Thus each outgoing polarization is individually neutral.

When the higher mode leaves the core, the source whole-rung permutation
produces the Clebsch isometry

  W|up> = (|up,t0>−sqrt2|down,t+>)/sqrt3,
  W|down> = (sqrt2|up,t−>−|down,t0>)/sqrt3,                  (4)

up to the fixed overall phase chosen for h_down. W†W=I and W intertwines
fundamental spin with the spin-half part of spin-half tensor spin-one.
Every branch has the same lower core charge+1 and a separately neutral bath
mode. The two polarizations in(4) are both retained, with weights1/3 and2/3.
They are spin channels, not the old charged-versus-neutral alternatives.
The vertex comes from the literal offdiagonal term−R/2 in(2); it is not an
independently installed weak contact or a discarded projection.

## Radiation, retained bound mode and finite errors

For the higher J=1/2 branch on the half-line, the endpoint diagonal is5/2,
bulk diagonal2 and offdiagonal−1/2. If U|r>=(−1)^r|r>, then
K_(1/2)=4I−U K_Neumann U, where K_Neumann is the original neutral bath
I+L_halfline/2. Consequently its endpoint survival is exactly
J0(t)^2+J1(t)^2 and tends to0. Complete probability remains1; every fixed
core region empties. All escaped products in(4) have the same lower charge
and neutral ejecta. The exact lower reference stays stationary.

The retained J=3/2 branch has endpoint diagonal1, bulk diagonal2. Its bound
mode has ratio r=1/2, energy3/4, normalized squared endpoint weight3/4,
and asymptotic endpoint survival9/16. It is not another spin-half family.
The continuum remainder is integrable on[1,3] and decays locally. The same
endpoint recurrence/resolvent argument as the prior impurity proof establishes
there are no other point eigenvalues. Finite chains retain recurrences.

After subtracting2I, both finite and infinite higher-branch operators have
norm at most1; finite-volume vector error is at most2|t|^N/N!. The bound
branch has norm at most3/2, and the conservative bound2 may be used. Full
cap vector error is at most6N|t|2^-B, effect error at most twice that. First
choose N for a fixed observation window and then B; only then take the
large-time local-relaxation limit. No primitive event time is assigned the
status of physical time by these bounds.

## Exact interface to the current row11 operators

The row11 tagged CAR/Gauss source transports flavor a=3s+f and its integer
charge independently. On two selected family coordinates f=0,1 define

  V|s,0>=E l_s,       V|s,1>=E h_s.                         (5)

This is an exact isometry intertwining charge+1 and all three spin matrices.
Its family matrix units are V(I_spin tensor |f><g|)V†; they commute with
primitive charge and transform as spin scalars. Core energy restricts to
I_spin tensor diag(1/2,5/2). The emitted vertex is the spin intertwiner(4),
so it preserves the full rotation/charge algebra instead of just a diagonal
spin component. This supplies a concrete operator dictionary for two source
families and their neutral-emission channel.

It does NOT yet intertwine the row11 mass analyzer A_N^r, its clock, or the
full three-family source. The source core contains a retained spin3/2 quartet,
not a third spin1/2 family. No basis map can turn these finite spectra into
three selected lepton masses. The row11 dressing/ghost operators can carry
the neutral family matrix units without changing their Gauss values because
the flavor charges coincide. Actual occupation-controlled composition on
that full CAR/Gauss chart is a separate source service, not supplied by a
formal tensor product or by(5). The declared coherent native source services
and physical spacetime-spin/field attachment keep their own exact scopes.
Absolute charged rest energy needs the common charged/vacuum reference and
propagation map, not just the core energy difference2. Physical μ/e/ν and
weak chirality remain unassigned.
