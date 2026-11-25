# Three signed probes on one polarizable ring

Three identically calibrated material probes on a closed, frustrated ring can each prefer the opposite sign at either of the other two probes. On a ring of three guarded Interval Echo dimers a spectrometer starting from a completely mixed material state measured **+0.24609 ± 0.03126** exchange units for all three signed pair energies. On a longer ring crossing five owners, all three chosen pairs again gave positive signed energies, but a fourth port reversed the sign. Three equal charges together also cost substantially more than the pair measurements predict. This is a measurable collective polarization effect with a finite material range experiment; a universal electric potential requires different records. [Protocol](PROTOCOL.md) · [Mixed-source clock protocol](CLOCK_PROTOCOL.md) · [Native calculations and tests](test_ring.py)

## Preparing and reading the signs

A dimer has two ports of equal federation degree on one native seam. Its two orientations place an amount-two load on rail A or rail B. A guard preserves that two-state sector. Join the B rail of each dimer to the A rail of its successor, including the closing bond. Pin each external test load next to an A rail; the sign is its deviation from one. The sign convention and contact type are identical at every tested terminal. Each terminal has its own unit reference on a distinct adjacent port. An independent calibration writes an echo, reads the two outgoing sharp material loads and adds them. Subtracting two from that sum recovers the original sign even though the write has redistributed the loads. The private collar label is unused. All 16,384 calibration launches at each source setting gave the correct sign. The reference ports are used in separate calibrations. Other seams are shuttered during spectroscopy. Ranging uses an independent blank route and a travelling double; a held swap eliminates arrival. [Physical layouts](geometry.py) · [Calibration bins and raw phase counts](records.json)

The three-dimer ring fits inside one icosahedral cell. The eleven-dimer ring crosses five owners of the 20-cell federation. Its three designated terminals sit at ring sites 0, 3 and 6, giving shorter arcs of three, three and five dimer bonds. A fourth terminal at site 2 tests an even arc from site 0. Unused terminals carry the neutral load. Intact-board ranges between the designated physical pins are three, five and eight links on the longer ring. The fixed geometry and station menu are apparatus choices.

## Why closing the loop changes the sign test

Write Z=+1 when the amount-two load occupies A and Z=−1 when it occupies B. A connector favors equal physical loads, which requires *opposite* Z at its two dimers. An odd ring cannot satisfy all its connectors. Reversing a local dimer label flips both incident preferred orientations, so the product of preferred signs round the circuit is invariant and equals −1. Reversing any subset of local dimer names cannot remove this frustration. The two routes between terminals favor opposite relative signs.

Projecting the registered blank seam columns onto the guarded material words gives the following generator up to a sum of known constants and single-terminal terms:

\[
H=-h\sum_i X_i + J\sum_{i\sim i+1}Z_iZ_{i+1}
       -J\sum_a q_a Z_{s_a},\qquad h=1/2,\quad J=16.25 .
\]

These coefficients follow from the measured native operator with contrast strength eight. Every interior connector compares B to A, and every pin compares its signed load to A. This accounts for the opposite signs of the two diagonal terms. Degree completion on a dimer is a scalar because its two ports have equal degree. The executable restores every scalar, including the one-port degree costs, before taking absolute phase readings. All eight small-ring material words are also computed directly from the full 99-state seam matrix on the enrolled physical addresses. On the larger ring the ground energy is certified by componentwise bounds: for a strictly positive trial vector v the minimum and maximum of (Hv)ᵢ/vᵢ enclose the ground energy, since the hopping has nonpositive off-diagonal elements and the flip graph is connected. The executable verifies a width below 5×10⁻⁸ and that the enclosure lies inside one 2⁻²⁰ energy bin. The verifier represents every component of the proposed positive vector by its exact dyadic value and checks both row inequalities with integer arithmetic. [Native phase experiment](measure.py) · [Independent microscopic word matrix](mixed_clock.py)

Removing exchange leaves h=0 and J=16. For any two individually pinned signs, one can satisfy both terminal preferences while leaving just one unsatisfied bond. Of the two connecting arcs exactly one has the wrong parity; place the unsatisfied bond on it. One unsatisfied bond is a lower bound on the whole odd ring. All four sign settings therefore have exactly the same interaction contribution at their minimum. The two-source signed energy vanishes identically in this control at any distinct terminal pair. The nonzero observed pair preference requires mobile polarization. With exchange active and the bond on the short arc cut, the eleven-ring 0,3 pair gave **−15.77502**. Cutting the bond on the other arc gave **+15.84675**. The two signs were read at the original two physical terminals in every bank.

## Measured pair and triple energies

For each pair let C=(E₊₊+E₋₋−E₊₋−E₋₊)/4, with the other terminals held neutral. For three designated probes the three-ring gives C=0.245487, 0.245484 and 0.245486. The longer ring gives **0.116377** at sites (0,3), **0.116375** at (3,6), and **0.017568** at (0,6). Between the first and the fourth probe, sites (0,2), the same C is **−0.224410**. All contrast errors are below 0.000022 exchange units simultaneously at 99%. A strength calibrated at the first longer-ring pair predicts 0.06983 for the equally strong (3,6) pair at its five-link range, and 0.04364 at (0,6) at eight links. These predictions fail within the counting bounds.

C measures the part odd under reversal of each probe. In the small ring the neutral-subtracted pair energy itself is about +0.4015 for equal signs and −0.0895 for opposite signs, including an even contribution of +0.1560. At the five-bond arc of the longer ring the corresponding pair energies are +0.0660 and +0.0309, both positive despite the positive C. A positive C is thus a relative preference. It does not by itself establish an attractive force, a charge-product law or a distance-dependent potential of a single charge.

For a further test use all 27 source settings and form the third mixed difference

\[
T(q)=E(q_1,q_2,q_3)-\sum_{i<j}E(q_i,q_j;0)
           +\sum_i E(q_i;0,0)-E(0,0,0).
\]

Every pairwise-additive energy has T=0. For three equal signs T is **+31.36719** on the small ring and **+32.01094** on the larger ring; each has a conservative error below 0.000175. This failure of superposition is visible even though all three isolated pair contrasts are positive. A large mixed difference for strongly perturbed, closely spaced material probes does not decide weak-field electromagnetism in another medium.

The stationary-arm measurements use a supplied lowest-state selector, a phase reference and a geometric sequence of 18 calibrated exposures per setting. All absolute-energy quadrature outcomes enter the receipt, including neutral, cut and changed-law panels. Each measured phase starts within the unwrap window and its uncertainty during doubling is less than π/3. Hoeffding's bound over 30,000 quadrature components gives each component an error below 0.001930; converting to phase and including the certified energy mesh gives the stated absolute-energy and pair errors. There are 156 charge panels across the two rings. [Complete counts](records.json)

## Reading the floor without selecting a ground state

Independently prepare each of the three dimers in a fair mixture of its two sharp rails. The resulting density matrix is the identity divided by eight. A supplied local clock controls powers of native free evolution and reads its 262,144-way Fourier output. Conditional on the spectral line, successive clock measurements commute with the generator and have independent pointer noise. Tracing out the initial material gives the equal-weight mixture of the eight line laws, including their multiplicities. The reader uses only its clock outputs. For each original launch it takes the median of 25 clock readings, and for each source setting the minimum of 256 launch medians.

The Fourier line has bin width 1/256. At a circular offset d in bins its probability per bin is bounded by 1/(4d²) for 0<|d|≤131,072. The total probability beyond eight bins is at most 1/14. The binomial upper tail for at least 13 bad readings therefore bounds each median error. The probability of missing every ground input is at most (7/8)²⁵⁶ per setting. The sum of both bounds over all 27 settings is below 0.000020. A fixed interior spectral window removes the wrap ambiguity. The readout allowance includes half the 2⁻²⁵ energy grid. Each recovered ground energy and each of the three resulting C measurements has error below 0.03126 with this joint coverage. The measured C was +0.24609375 for all three pairs. The independently recovered equal-sign three-source residual was +31.37109 with uncertainty below 0.25001. The large-ring spectroscopy uses the stationary-arm source. [All 172,800 clock exits](clock.npz) · [Finite clock and mixed input](mixed_clock.py)

The complete OpenWave criterion asks for a far-field inverse-distance potential at elementary charge and an independently produced strength for any claim about α. The ring uses material excesses and vacancies on a guarded finite medium; it has neither an elementary-electric-charge calibration nor a universal far-field response. The full criterion is Open. [Official row](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#simplest-test-per-criterion)

Run from the repository root:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_20_electric.frustrated_ring
nice -n 15 .venv/bin/python -B -O -W error -m "$pkg.measure"
nice -n 15 .venv/bin/python -B -O -W error -m "$pkg.mixed_clock"
nice -n 15 .venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_20_electric/frustrated_ring
```
