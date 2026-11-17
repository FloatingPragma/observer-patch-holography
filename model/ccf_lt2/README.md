# Matter, fields and energy response

This is the official CCF-LT-2 executable model, the common construction
specified by [CCF-2](../SPEC.md). Every graph vertex carries 96 native
fermionic modes. One Hamiltonian combines signed charged and neutral
transport, weak contacts, family shifts, electric and magnetic fields,
and independent local capacity responses to their total energy.

The [source proof](PROOF.md) explains exact charge conservation and the
complete tensor action. [RESULTS.json](RESULTS.json) records the executed
checks. These establish finite algebra and response; physical spectra,
continuum limits and observed constants require their own calculations.

## Run

From the repository root, with NumPy, SciPy and SymPy installed:

```sh
python3 -B -m model.ccf_lt2.verify
python3 -B -m model.ccf_lt2.verify --full
python3 -B -m pytest -q model/ccf_lt2/test_model.py
```

The first command executes small complete-carrier checks and compares every
result with `HEADLINE.json`. The full command also recompiles all 114
retained source factors, verifies the source commitment, checks both
192-mode edge and 288-mode triangle polynomials, all weak contacts and
complete native columns through the per-vertex response and capacity factors. It compares the
result with `RESULTS.json`. A discrepancy exits with failure. Use
`--output PATH` to regenerate a current result file.

The package reads only its local Python and JSON scientific inputs. No
adjacent workspace, network service or particular Git checkout is needed.
Exact supplied binary coefficients and the empirical inputs are retained
in `source_snapshot.json` and `calibration_inputs.json`. Recompilation must
recover the specified matter-polynomial SHA256 before constructing the
combined Hamiltonian. Calibration agreement is not an independent prediction.

## Using the source

```python
from model.ccf_lt2 import source, spatial

law = source.assemble(spatial.triangle())
# Native Fock occupation, then one integer flux per oriented edge.
column = source.apply_source_column(law, occupation=0, flux=(0, 0, 0), cap=6)
current = source.phase_derivative(law, edge=0, matter_only=True)
contact = source.phase_derivative(law, edge=0, other_edge=0, matter_only=True)
# One capacity index per vertex, Phi index, native basis label.
key = ((0, 0, 0), (0, (0, (0, 0, 0))))
joined = source.apply_total_column(law, key, cap=6)
total_current = source.apply_total_phase_column(law, key, edge=0, cap=6)
```

Each returned mapping includes every nonzero output of that literal basis
column. Supported chart faces are simple closed boundaries with at least three
distinct vertices; repeated or backtracking face words are rejected. A native bit mask includes every one of the `96V` occupations;
finite flux shifts vanish at the hard wall. The reference is an idle tensor
factor, so every column retains its input reference label. Full arbitrary
state propagation on `2^(96V)` dimensions is not performed. The executed
finite matrix witnesses declare their own carriers and preparations.

The configurable local-response API requires nonnegative lapse coupling; the
specified value is `1/4`, with strictly positive lapse on every response label.

The component arguments are in [signed transport](DIRAC_PROOF.md),
[weak interactions](INTERACTIONS_PROOF.md), [capacity register](GEOMETRY_PROOF.md),
[local response](LOCAL_GEOMETRY_PROOF.md), and the
[complete tensor action](JOINT_GEOMETRY_PROOF.md). The
[spacetime bridge](EINSTEIN_BRIDGE.md) states the premises for consuming
RER's conditional geometric theorems.
