# What occupies the two neutral ports?

At full conversion the receiving apparatus recorded **2,493,659 neutral arrivals** across twelve banks and both charge signs. Every arrival has one excess and one vacancy at the two inspected ports. With the pulse off, all 2,494,290 arrivals were charged. The programme here reanalyzes the published ternary words from the same starts, including source refusals and missed arrivals. It uses no hidden winding or particle labels. [Measured words and apparatus](../../../../../experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge/neutral_ladder/README.md) · [Reconstruction](measure.py) · [Complete bins](results.json)

## From the charge meter to an observable spatial profile

The meter reports loads 0, 1 or 2. Calibration on the filled host assigns signed charges −1, 0 and +1. For the addressed hub and target let these be h and t. Their net Q = h+t and oriented difference D = t−h recover the signed readings exactly, since 2h = Q−D and 2t = Q+D. The donor is read separately. This catalogues all 27 possible hub, target and donor words. Omitting the ordering of the first two ports leaves 18 distinct signatures. Source refusal and no arrival are separate recorded outcomes.

Observers can compare spatial charge profiles even if they describe the edge in opposite directions. For a phase step θ across one edge, the power of the reconstructed charge density is |h+t exp(iθ)|². It is invariant under swapping the coordinate origin to the other end, with the direction reversed. The uniform reading (θ = 0) measures Q² and the alternating reading (θ = π) measures D². On this ternary carrier every neutral word is either locally filled, with both powers zero, or a resolved pair, with powers 0 and 4. Giving both dipole orientations equal weight cancels their mean signed difference but retains alternating power 4. These are functions of joint local charge measurements. A scattering experiment would also need a calibrated coupling to a probe. [Exact reconstruction and neutral-word theorem](ChargeRecords.lean)

Each full-pulse sign has six independent simulated acquisition banks of 2³⁸ original starts. The negative input gives 1,245,274 resolved pairs and the positive input gives 1,248,385. The source and receiver recorded 3,298,532,389,669 refusals or missed arrivals on those starts combined. The neighboring half-pulse banks resolve both charged and neutral outcomes. A depleted donor gives zero resolved pairs in both independently prepared controls.

This fixes the material composition at the stated spatial resolution and appointment. Identifying a neutrino at longer wavelengths needs the propagating response of this same product and a flavour measurement identified by its charged partner. [Observer particle catalogue](../../../../../../docs/registers/PARTICLE_RECORDS.md)

To reproduce from the repository root, run the following.

~~~sh
p=model.candidate_laws.interval_echo.openwave.row_09_neutral.charge_profile
python -m "$p.measure" --write
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_09_neutral/charge_profile
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_09_neutral/charge_profile/ChargeRecords.lean)
~~~
