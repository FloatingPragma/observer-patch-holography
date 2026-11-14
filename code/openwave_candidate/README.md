# Reproduce the monitored observer construction

From the public RER checkout, install and execute the portable observer model:

```sh
python -m pip install './code/openwave_candidate/candidate[integration,test]'
python -m pytest code/openwave_candidate/candidate/tests
python code/openwave_candidate/integration/assess.py --output /tmp/oph-assessment.json
python code/openwave_candidate/integration/assess.py --verify /tmp/oph-assessment.json
python code/openwave_candidate/platform_packet/validate.py
```

The bounded one-load instrument preserves its integer Gauss sector and retains
all measured outcomes. Complete physical criteria require particle identification,
interaction and refinement results beyond these finite checks. All 31 official
OpenWave rows are untested. The classical native registers, bare carrier controls
and monitored joint state retain their declared domains.

## Independent source/tag/link controls


Run from the published RER root with its declared Python dependencies:

```sh
python3 -m model.native_charged_propagation --check
python3 -m model.native_criterion_packet --check
python3 -m model.native_electric_gauss --check
python3 -m model.native_gauss_clock --check
python3 -m model.native_gauss_walk --check
python3 -m model.native_moving_clock --check
python3 -m model.native_observer_program --check
python3 -m model.native_spatial_refinement --check
python3 -m model.native_transport_spectrum --check
python3 -m tools.build_joint_observer_claims --check
```

Scientific commands recompute outputs from source. Receipt hashes authenticate
their inputs; they do not substitute for execution. The literal obligations
are retained in [the frozen platform packet](platform_packet/packet.json).

## Scoped observer experiments

The bundle explicitly includes three experiment groups under `validation/observables` and the scoped source-to-field bridge under `validation/bridges`.
Their code, receipts, complete declared readouts, reviews and plots are copied
byte-for-byte and authenticated by `export_manifest.json`. Other validation
work and private plans are excluded. These experiments do not change the
assessment or award a physical criterion.

The modular experiment conditions a declared two-port observer corner from
an actual native output. It retains the full native density and both terminal
P/Q measurement branches in a compressed archive. A faithful conditioned
state generates finite modular correlations, with no added Gibbs state or
noise floor. The dimensionless modular parameter has no automatic physical
clock or subjective-phenomenology interpretation.

The stacked experiment adds a declared third direction to the original 2D
ring screen, prepares a positive-band packet with uniform modular records,
and evolves at zero commit rate. Its saved seven-offset position apparatus
includes ports and seam midpoints. The centroid is compared with a reference
ray; anisotropy, periodic wrap and the supplied flat metric are explicit.
No hidden particle path, gravity or transfer of a 2D force kernel is derived.

From the public RER root, after installing the candidate and the repository's
plotting dependencies, run:

```sh
python code/openwave_candidate/validation/observables/modular_time/run_modular_time.py --output /tmp/modular-time.json --plot-data /tmp/modular-correlations.csv --native-state /tmp/native-state.npz
python code/openwave_candidate/validation/observables/modular_time/review_modular_time.py
python code/openwave_candidate/validation/observables/spacetime/stacked_wavepacket.py --output-dir /tmp/stacked-wavepacket
python -m pytest -q code/openwave_candidate/validation/observables/modular_time/test_modular_time.py code/openwave_candidate/validation/observables/spacetime/test_independent_stacked.py code/openwave_candidate/validation/observables/spacetime/test_stacked_assumptions.py
```

The independent modular review checks the frozen adjacent receipt and plot
custody; the first command emits a separate fresh native result. Stacked
archive tests inspect the adjacent frozen complete position distributions.
See [modular method](validation/observables/modular_time/METHOD.md) and
[stacked method](validation/observables/spacetime/METHOD.md) for every supplied
piece, observable and numerical limitation. Rendering requires Matplotlib,
which is declared in the public repository's `requirements.txt`; the candidate
package itself provides the numerical experiment dependencies.


## Shared preparations, occupations and calibrated energy readouts

The integrated runner includes the pinned-charge H_m orbital source, regional
continuity current, shared nonnegative Bose occupations and explicit uniform
Gauss-fiber preparation. The added preparation uses ideal Fourier and reversible
modular arithmetic gates. Continuing record measurements retain the complete
postmeasurement state, including its component outside the uniform code.

The GeV experiment preserves every energy branch. Its supplied 511 keV half-gap
anchor is a calibration, not a prediction of electron mass. A separate complete
sequential energy/native-record/energy instrument retains source and battery
energy changes. Its apparatus and work reference are declared resources.

```sh
python -m oph_openwave.runner pinned_charge
python -m oph_openwave.runner native_currents
python -m oph_openwave.runner many_load_native
python -m oph_openwave.runner fiber_preparation
python code/openwave_candidate/validation/observables/energy_gev/run_energy_gev.py
python code/openwave_candidate/validation/observables/energy_gev/spectral_apparatus.py
python code/openwave_candidate/validation/observables/energy_gev/review_energy_gev.py
python code/openwave_candidate/validation/reviews/review_energy_gev.py
python code/openwave_candidate/validation/reviews/review_native_apparatus.py
python code/openwave_candidate/validation/review_many_load_native.py
python code/openwave_candidate/validation/review_fiber_preparation.py
```

These commands regenerate adjacent scientific receipts. Reviews are separate
implementations where their methods state that; producer checks do not count as
independent physical audits. The shared-occupation control uses fixed total
nonnegative Bose load, and does not import signed-pair or fermion-sea results.

## Current driving the existing field executor

The bridge reuses the packaged, byte-identical Maxwell provider with its license
and provenance. It maps actual coherent native currents to a supplied spatial
Fourier mode, checks continuity and Gauss, and accounts for field energy as
supplied source work. It retains source-hold discretization error. There is no
reciprocal source/field coupling or distant radiation claim.

```sh
python code/openwave_candidate/validation/bridges/native_maxwell_source.py
python code/openwave_candidate/validation/bridges/review_native_maxwell_source.py
python -m pytest -q code/openwave_candidate/validation/bridges/test_native_maxwell_source.py code/openwave_candidate/validation/observables/energy_gev/test_independent_energy_gev.py code/openwave_candidate/validation/observables/energy_gev/test_spectral_apparatus.py
```
