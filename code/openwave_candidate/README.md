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

Scientific commands recompute outputs from source. Receipt hashes authenticate their inputs; they do not substitute for execution. Complete finite domains and physical limitations are in [WORK_ITEMS.md](WORK_ITEMS.md), with all 31 literal obligations in [SCOREBOARD.md](SCOREBOARD.md).

The published RER revision, six PDFs, ledger and live Results/status assets are authenticated. The flagship has 59 pages.
