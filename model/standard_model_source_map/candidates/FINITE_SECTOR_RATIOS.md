# Finite source threshold ratios

The six conserved occupation sectors have rigorous lower and upper
bounds on their lowest energies. Their three charge-conjugate pairs give
three distinct threshold intervals on the same member-zero Wilson source,
state, reference and clock. They are gauge-neutral composites, not three
charged-lepton generations. This calculation compares their finite source
energies without introducing a physical mass assignment.

Use the [full-space sector bounds](../../particle_primary_sector_threshold/PROOF.md),
the [admitted two-link upper bound](../../particle_primary_branch_certificate/PROOF.md)
for the weak trion, and [exact charge conjugation](../../particle_primary_charge_conjugation/PROOF.md).
For the common ground `E0 in [L0,U0]` and sector bottoms
`Ei in [Li,Ui]`, put `Di=Ei-E0>0`. The interval calculation preserves the
same unknown ground in numerator and denominator:

    min_{e=L0,U0} (Li-e)/(Uj-e) <= Di/Dj
        <= max_{e=L0,U0} (Ui-e)/(Lj-e).

Each bounding fractional-linear function is monotone or constant on the
allowed ground interval, so checking its two endpoints is sufficient.
The ratio within a sector, or between its exact conjugate thresholds, is
one. A positive common clock or unit conversion cancels. Independent gap
interval division would be valid but less informative because it discards
the common energy-zero correlation.

The outward-rounded finite excitation intervals are:

| Existing source sector ID, with its opposite increment | Native increment in charge order (-4,-3,1,2,6) | CAR parity | Threshold above full ground, source units |
| --- | --- | --- | --- |
| `weak_pair_charge_six` | (0,2,0,0,1) | Odd | [1.5467,12.3981] |
| `color_triple` | (1,0,0,2,0) | Odd | [0.5664,12.7244] |
| `color_determinant_charge_six_hole` | (0,0,0,3,-1) | Even | [3.5668,15.7156] |

Every ordered ratio among these three families is less than21.16; among
the two odd families it is less than15.31. The exact bounds are retained
in `finite_sector_ratios.certificate()`. Consequently these particular
member-zero thresholds cannot directly be the charged-lepton hierarchy
under one common mass conversion: the observed muon/electron ratio is
206.7682827 with standard uncertainty0.0000046
([2022 CODATA](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)). The
reference value is used only for comparison, not for choosing the source
or optimizing a bound.

This is a scoped exclusion of that direct assignment. It does not exclude
other excitations, another capacity, or a controlled effective mass
construction. None of these intervals supplies a momentum-dependent
branch, a charged threshold with its compensator removed, a Lorentz spin
assignment or a Standard Model species. The source ratio bounds are not
point eigenvalue estimates or physical rest-mass ratios. Their numerical
midpoints have no special status.
