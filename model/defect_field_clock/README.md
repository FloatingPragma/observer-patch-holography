# One co-driven source clock and energy calibration

This package appends an operational clock to the retained complete output of
`defect_field_bridge.dynamics_packet`. Two clock hands, two untouched purifiers,
and two final binary pointers remain in the output. Native C/E records, their
zero and abort branches, internal charge pointers, source addresses, and the
previous field outputs are retained through their bound parent arrays. No
clock or field state is normalized after a measurement.

The explicit design convention is **co-driven action time**. In one physical
interval with native current j, the entire field Hamiltonian and both clock
hands have the same positive factor lambda_j=1+j/4. The clock Hamiltonian is
H_c=Omega_ref Z/2, with the single declared calibration Omega_ref=1/4 GeV.
The hand phase is Omega_ref times tau_eff, where tau_eff is the sum of
lambda_j delta_tau_j. The clock advances once per complete interval, not once
per term in a product-formula implementation. Native repair attempts,
compilation pulse counts, and host runtime do not define this time.

The faithful reference populations (1/6,5/6) come from the existing native
clock calibration. Its centered modular generator is log(5) Z/2, so
-log(rho_ref)=constant I+[log(5)/Omega_ref] H_c. A source-compiled preparation
purifies each reference and applies a Hadamard. Under exp(-i tau_eff H_c),
the X and Y expectations are respectively (-2/3) cos(theta) and
(-2/3) sin(theta). Source-compiled final Pauli instruments retain all four
outcomes. The check acts on all four work-input columns; an independent
reference test verifies complete trace preservation. The source uses its
specified floating numerical semantics, with a declared 1e-8 complete
instrument/preparation admission; this is not a physical hardware guarantee.

`infer_time` receives only the reconstructed four-outcome probabilities, the
reference calibration, and their error bound. `infer_energy_gap` additionally
receives a field Ramsey probability and its error bound. Neither accepts the
programmed duration or the field energy spectrum. Exact rational alternating
Taylor bounds and 60 bisections enclose the clock phase and the Ramsey phase;
their ratio gives an energy-gap interval in the same GeV scale for every mode.
The retained packet rounds reconstructed probabilities to a declared 1e-9
grid and adds its worst-case rounding error before this exact inversion.
Both inverse problems explicitly assume their phase lies in [0,3/2]. They do
not resolve phase aliases or an absolute common vacuum-energy offset. The
packet compares its inferred interval to the supplied calibrated spectrum
only after inversion. Such comparison is a postdiction, not an independent
prediction of the calibration inputs.

These are ensemble-probability inverses, not single-shot time readings.
`sampling_envelope` supplies a finite-shot option for independent repetitions
of the same conditional native branch: Hoeffding and a union bound bound all
four clock-cell frequencies and one Ramsey frequency. The exponential tail
has a conservative exact rational upper bound. Add the corresponding sampling
errors to the separate source-service errors before inversion. The receipt
contains an example budget but no collected experimental shots.

The SI conversion is t_seconds=hbar_GeV_s tau_eff. The package uses the exact
SI values h=6.62607015e-34 J s and e=1.602176634e-19 C, with hbar=h/(2 pi),
and an outward rational Machin-series enclosure for pi. The constants are
specified by the [BIPM SI defining constants](https://www.bipm.org/en/measurement-units/si-defining-constants).
This is an explicit common unit calibration, not evidence that this designed
model clock has been matched to a laboratory clock.

The clock has no intermediate feedback. It acts on separate factors, so its
source execution can be contracted separately and tensor-composed with each
native branch, including piecewise noncommuting field dynamics. This permits
use with propagation and Higgs interaction protocols sharing the same
lambda_j convention; their full combined outputs require their own receipt.

Run `python3 -m pytest model/defect_field_clock -q` and
`python3 -m model.defect_field_clock --check`. Use `--write` to regenerate the
JSON and compact retained clock-array receipt after an intentional change.
