# Complete two-cell quantum energies

## 1. One source and one clock

The control retains twelve lepton CAR modes, the actual six-state active
quantum link, and both five-state quantum Higgs memories. Its basis is
`(mask, link, hleft, hright)` and its full dimension is
`2^12*6*5*5=614400`. The left modes are 0 through 5 and the right modes are
6 through 11. Each site's internal order is the native `(e^c,L_1,L_2)`
module with spin order `spin*3+internal`. The link order is the blank,
the e^c raising component and the four weak-doublet raising components.

The same declared clock is `exp(-it H)`, with

`H = H_hop + H_Y,left + H_Y,right + N_H,left + N_H,right + P_link`,

`H_hop = sum_(a,b,spin) (i/2)(-1)^spin
 c_left,spin,a^dagger U_ab c_right,spin,b + h.c.`.

The coefficient is the actual side-two spectral z-edge coefficient. Both
charged-lepton Yukawa coefficients are one; the quark Yukawa coefficients
are zero. Each local quantum term is the complete operator in the
[quantum Higgs proof](../../chiral_higgs/PROOF.md). No scalar expectation
or classical transporter is substituted. This is a two-cell control,
without a claimed three-dimensional spatial refinement.

The full native construction also contains the quark occupation factors
and colored link alternatives. Empty quarks and unexcited colored link
rows form a reducing spectator sector for these specified coefficients:
quark hopping conserves quark number and quark Yukawa creation is switched
off. The active 614400-state domain retains every lepton occupation,
every active link alternative and both entire scalar memories. It is not
a vacuum-generated Krylov subspace.

`hamiltonian.column` defines every full-domain column. Its complete
320-column local operator equals the primary Higgs operator exactly.
All 24576 fermion/link hopping columns retain the CAR signs and adjoints.
The local and path Ward identities prove commutation with every local
Gauss generator throughout the full active domain. The local scalar
number and excited-link projector also commute with their currents.
The scalar/link terms have combined norm at most three. Ten forward
hopping monomials and their adjoints, each of coefficient magnitude one
half, contribute at most ten. Four local pair monomials and their
adjoints contribute at most eight at each site. Thus `||H||<=29` on the
entire carrier. Hermiticity and this bound imply a complete finite unitary
clock independently of a selected state.

Two exactly conserved labels are

`L = sum_sites (N_L - N_ec)`,

`2Sz = N_spin_up - N_spin_down`.

The first is a specified source occupation charge, not an independently
identified measured lepton number. The scalar and link registers have
zero charge for these labels. Spin-diagonal hopping and spin-singlet
pair creation preserve both quantities.

## 2. Every physical state, including boundary charge

The [exact Gauss basis](BASIS_PROOF.md) retains the joint kernel of the
actual currents at both cells, with specified static hypercharges at the
boundaries. An integer q denotes `(q_left,q_right)=(0,q)`. A pair specifies
both charges. These boundary registers are explicit compensators; they
are not dynamical particle identifications.

For each boundary and `(L,2Sz)` the basis code returns a rational matrix B
with full column rank, its positive Gram matrix `G=B^dagger B`, and a
coordinate-row map E satisfying `E B=I`. Exact RREF and independent
SU(2) character counts certify the dimension. The complete physical
counts are 1053 for zero boundary, 718 for right charge minus one,
783 for right charge plus one, and 718 for left charge plus one.

For efficient real arithmetic introduce the diagonal unitary
`D|mask,link,hleft,hright> = i^(N_left)|mask,link,hleft,hright>`.
It commutes with every Gauss current and both conserved labels. In
`H_real=D^dagger H D`, hopping coefficients are real, left Yukawa
coefficients reverse sign, and every remaining coefficient is real.
The actual physical embedding is **D B**, including those phases.
Dropping D in an observable would change the calculation.

The code calculates

`C=E H_real B`, `K=G C`,

and checks exactly

`H_real B=B C`, `K=K^T`, `C^T G=G C`.

Therefore `(K,G)` is the complete generalized Hermitian eigenproblem of
this physical block. The first identity checks that no outgoing component
has been discarded. Enumerating every conserved block covers the full
physical Gauss domain. An empty-state trajectory cannot replace this
enumeration, since H conserves L and the physical ground need not be in
its L=0 subspace.

## 3. Rigorous complete-spectrum numerical enclosure

A floating generalized eigensolver proposes a full eigenbasis and energies
only in the explicit `--generate-witnesses` command. The compressed
`eigenbasis_witnesses.npz` artifact stores integer numerators with denominator
`2^44`, without pickle. Its manifest binds every proposal to the exact
rational K and G matrices, normalized boundary charges, conserved labels
and dimension. The receipt includes the artifact's SHA256 digest.
Both `--write` and `--check` replay these stored proposals without running
the eigensolver. A missing or mismatched witness fails replay. This fixes
the proposal inside each degenerate eigenspace across numerical platforms.

Its rank, orthogonality and residual are not accepted on floating-point
trust. The stored basis and diagonal energy matrix define dyadic rational
matrices Q and Lambda with denominator `2^44`. Replay computes the
following using exact integer arithmetic:

`F=Q^T G Q`, `R=K Q-G Q Lambda`,

`delta >= ||F-I||_2`, `r >= ||Q||_F ||R||_F`,

`M=||Lambda||_2`.

The symmetric row-sum norm gives delta. Integer squared Frobenius norms
and an upward-rounded integer square root give r. The implementation
checks `delta<1`. Hence F is positive definite and Q is invertible: its
columns span the entire block, including any states a low-energy Krylov
calculation could miss.

Set `A=Q^T K Q`. Then

`A-Lambda = Q^T R + (F-I)Lambda`,

so `||A-Lambda|| <= r+delta*M`. The exact physical eigenvalues are those
of the Hermitian matrix `F^(-1/2) A F^(-1/2)`. If
`a=(1-delta)^(-1/2)-1`, then

`||F^(-1/2)-I|| <= a`, `2a+a^2=delta/(1-delta)`.

Consequently

`||F^(-1/2) A F^(-1/2)-Lambda||
 <= (r+2 M delta)/(1-delta) = epsilon`.

Hermitian minmax, equivalently the ordered Weyl bound, encloses **every**
ordered physical eigenvalue within epsilon of its ordered proposed
energy. This is a full-spectrum certificate, not a residual claim about
one guessed eigenvector.

The exported intervals have dyadic centers on the `2^-32` grid and fixed
radius `2^-28`, about `3.73e-9` in the declared clock's energy units. The
code verifies the exact rational inequality

`|center-Lambda_j|+epsilon <= 2^-28`

for each eigenvalue before exporting it. The fixed outward intervals
avoid reporting solver last bits as physical precision. All occupation
blocks and all eigenvalues are included in `sector_certificate`.

## 4. The physical reference and charged-sector boundaries

For zero boundary charge the unique full physical ground is in
`(L,2Sz)=(2,0)`. Its certified energy interval is

`[-4890707863/1073741824, -4890707855/1073741824]`,

approximately `-4.55482663484`. Every other physical eigenvalue is
separated from this ground by at least

`211959005/536870912`, approximately `0.394804412`.

The ground is defined by the spectral projector of this same H in the
entire zero-boundary Gauss space. Its simplicity follows from the complete
interval separation, not from a chosen charge-zero or vacuum-cyclic block.

For right boundary charge minus one the lowest certified window has two
states, in `(L,2Sz)=(1,-1)` and `(1,1)`, near `-3.50278278709`. The receipt
calls this a ground-window rank; numerical overlap alone does not prove
exact degeneracy. An exact symmetry may separately identify the two
energies. For right boundary charge plus one the unique full ground lies
in `(L,2Sz)=(4,0)`, near `-3.71236223433`.

A single e^c addition to the neutral reference changes L from two to one.
A single removal changes it from two to three. Both change source fermion
parity from even to odd. Therefore the right-positive-boundary L=4 ground
is dark to this removal instrument. The lowest L=3 states are instead near
`-3.64502800307`. A charged-sector ground difference and a specified probe's
first nonzero spectral line are different quantities. The complete readout
must calculate the actual matrix elements and retain zero responses.

The reflected left-positive boundary has the same certified spectral
centers as the right-negative boundary, with the corresponding ground
window in L=3. Its relation to charged probes requires the actual spatial
and particle-hole symmetry, including CAR signs and scalar/link actions.

## 5. Preparation and interpretation

The unique neutral reference can be prepared by a complete energy PVM
from the physical maximally mixed state `P/1053`. Its original success
probability is `1/1053`; all other energy outcomes are retained. This is
a specified preparation controller. It is not the evolution of an empty
L=0 state into an L=2 ground under the conserving H.

If the apparatus must conserve the total L label, purify `P/1053` with a
retained reference whose label is minus the source label and whose parity
matches it. Source and reference then have total L zero and even parity.
A Gauss-invariant preparation unitary exists inside that joint sector.
The successful reference retains source L two and apparatus L minus two.
This refinement preserves the complete energy-PVM probabilities; it does
not erase the preparation reference or claim idle vacuum selection.

These results supply a fully calculated finite interacting reference and
controlled energies on one source. They do not establish a relativistic
vacuum, a momentum dispersion relation, measured species, a rest-mass
interpretation or a GeV calibration. The energy zero is the actual
specified physical reference. An interacting spatial refinement and its
error control must use the same quantum link, scalar state and clock.
