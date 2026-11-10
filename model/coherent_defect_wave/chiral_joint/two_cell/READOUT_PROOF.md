# Complete charged spectral records

The source, Hamiltonian and clock are those in `SPECTRUM_PROOF.md`.
All scalar and link registers remain quantum. The input to the measurement
sequence is the maximally mixed state on the complete 1,053-dimensional
neutral physical subspace, together with the retained preparation reference
specified in `INSTRUMENT_PROOF.md`. The binary ground filter succeeds with
probability 1/1053; its complementary outcome has probability 1052/1053.
Conditioned on success, the source is the unique ground of this same
Hamiltonian, in lepton-number sector L=2 and twice-spin-projection sector 0.
This preparation is a declared instrument, with no assertion of efficient
native idle cooling.

## Physical probe and its complete outcomes

The source e^c mode at either cell is swapped with a fermionic reservoir.
The reservoir and compensator start equally occupied, either both full or
both empty with probability 1/2. Every ancillary outcome is retained.
On the even source ground the four Kraus operators, in order of addition,
failed addition, removal and failed removal, are

    a†/sqrt(2), -n/sqrt(2), -a/sqrt(2), (I-n)/sqrt(2).

The minus signs follow from the full fermionic permutation. They are not
extensions of this formula to arbitrary source parity. The boundary charge
is the final reservoir occupancy minus its initial occupancy, placed at the
probed cell. The actual matrices satisfy all local Gauss constraints on
source plus retained apparatus. Addition reaches L=1 and removal L=3;
their twice-spin projections follow the probe mode. Both failed swaps
remain in the neutral L=2, twice-spin-projection 0 block.

In `readout.py`, each probe C is applied to the complete physical embedding
E_s=D_s B_s. The returned matrix T obeys the exact identity

    C E_s = E_d T.

No nonphysical image is dropped. The coordinate metrics G_s and G_d obey
sum T† G_d T/2 = G_s. Thus the four outcomes form a complete instrument
on the entire prepared conserved block, independently of any approximate
ground vector.

## Certified probabilities

The rationalized ground vector w has norm squared N=w†G_s w. The residual
and certified separation from every other physical energy bound its trace
distance from the exact pure ground by

    eta <= ||(H-d_0)w|| / (sqrt(N) separation).

Every destination block has a complete certified spectrum. Overlapping
energy enclosures are grouped into one detector window, so degeneracy is
never resolved by a numerical eigenbasis. Let Q contain the rationalized
vectors for such a window, delta bound ||Q†GQ-I||, and R their column
residuals. Separation s from all energies outside the window gives

    ||P_window - Q Q†|| <= min(1, ||R||_F/(s sqrt(1-delta))) + delta

in orthonormal physical coordinates. The first term bounds the angle
between equal-dimensional subspaces; the second bounds the replacement
of the approximate frame by its orthogonal projector. It must be added
after capping the subspace-angle term. The certified complete spectrum
establishes the required rank and separation.

For x=T w, the rational approximation to the conditional probability is
sum_j |q_j†G_d x|²/(2N). Its error is at most (eta+projector_error)/2:
each unscaled probe is a contraction and the original reservoir setting has
probability 1/2. Probabilities are rounded outward on a 2^-27 grid. The
lower probabilities sum to at most one and upper probabilities sum to at
least one. No normalization of approximate output is performed. Other
conserved destination blocks have exactly zero effect on this prepared
input by the exact probe intertwiner. The receipt certifies their common
zero-effect rule; it does not enumerate those redundant zero effects as
additional spectral records.

For the right-cell spin-up probe there are 460 conditional energy-window
records. The first addition record has energy transfer in
[1.052043840289116, 1.052043855190277] and probability in
[0.21838924288749695, 0.21838925033807755]. The first removal record has
transfer in [0.9097986244596541, 0.9097986393608153] and probability in
[0.23089293390512466, 0.23089294135570526]. Decimal endpoints here summarize
the outward rational bounds in the executable receipt. Positivity of each
lower probability proves visibility of that threshold. An interval
containing zero does not establish a dark state.

The right-positive boundary sector's absolute ground lies in L=4. It is
exactly dark to removal from the neutral L=2 ground, because a single
removal reaches L=3. An energy difference to that L=4 ground therefore
is not this detector's threshold.

## Original attempt law and conjugation

The unconditional instrument retains one failed-preparation outcome and
all 460 successful-preparation spectral windows. Each spectral probability
is multiplied by 1/1053; failed preparation has probability 1052/1053.
This gives 461 retained outcomes on the original attempt law. A failure
to prepare is recorded without claiming an energy measurement occurred.
The complete Gauss-constrained source plus reservoir/compensator has
boundary values 0 and +/-1. Their certified global minima are all at or
above the neutral E_0, so H_source-E_0 is positive on that full apparatus
domain when the apparatus Hamiltonian is declared zero. The swap is a
driven intervention; its work cost is not calculated.

The exact antiunitary in `SYMMETRIES.md` reflects the cells and conjugates
charge. It pairs right addition with left removal, and right removal
with left addition. Spin reversal gives the equivalent opposite-spin
experiments. All corresponding probability and energy enclosures overlap.
Comparing addition with removal at the same boundary does not use that
symmetry. Their distinct thresholds do not establish a particle versus
antiparticle mass difference.

These are finite, interacting, probe-weighted energy transfers. The open
two-cell source has no momentum dispersion relation or relativistic pole
mass. No conversion to GeV and no electron, muon or tau identification
follows from these records.
