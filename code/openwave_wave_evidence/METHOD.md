# Wave evidence and reproduction

For a connected weighted graph with total edge weight $W$, graph Laplacian
$L$, and projection $\Pi$ off the constant mode, the declared one-unit
process and normalized occupation covariance satisfy

$$
P=I-L/(2W),\qquad
C_n=\frac{\Pi+2\sum_{k=1}^{n-1}(1-k/n)\Pi P^k\Pi}{4W},\qquad
C_\infty=L^+-\Pi/(4W).
$$

The scalar slice energy uses
$q_a(p)=\sum_j4a^{-2}\sin^2(ap_j/2)$ and
$E_a=2a^{-1}\operatorname{asinh}(a\sqrt{\lambda+q_a}/2)$.
For $S=\lambda+|p|^2$, its stated error bound is
$0\leq S-E_a^2\leq a^2[\sum_jp_j^4+S^2]/12$.
The internal gaps are $0,5-\sqrt5,6,5+\sqrt5$, with multiplicities
$1,3,5,3$. They are native band labels without measured-species assignments.
The field-strength reflection form has a positive rank-two quotient at nonzero
momentum. Its finite contact terms and separated-support domain are retained.

## What the package contains

`evidence/` preserves three sealed packages, including failures and development
attempts. Their manifest hashes and entry counts are in `SEALS.json`.
The source bundle has 420 entries and six ancestor manifests; the review has
60 entries; the replay amendment has 124 entries. Sealed files retain their
bytes. Their historical instructions, timestamps, local paths and pin proposals
are provenance. The commands below are the portable reproduction interface.

The original source lineage is tracking commit
`60519f4d50cd4b16f7eeaad1a2d1e92544a037dd`. That identifier records ancestry;
execution uses the bundled bytes and requires neither that repository nor Git.
The complete candidate is specified in
`evidence/shared/unified_wave_submission_20261115/evidence/unified/PIN_SPEC.md`
and `READOUT_SPEC.md` beside it. Publication of this evidence is separate from
adoption of that candidate as the OpenWave model pin.

## Reproduction

Use Python 3.12 or 3.13 and NumPy. The sealed comparison receipts cover
Python 3.12.3 / NumPy 2.3.5 and Python 3.13.2 / NumPy 2.5.1. A clean environment
can install `numpy==2.3.5`. All output paths must be fresh and outside this
package. From the repository root:

```sh
python -B code/openwave_wave_evidence/reproduce.py custody
python -B code/openwave_wave_evidence/reproduce.py independent --output /path/to/independent.json
python -B code/openwave_wave_evidence/reproduce.py replay --output /path/to/wave-run
python -B code/openwave_wave_evidence/reproduce.py tests --run /path/to/wave-run --output /path/to/tests.json
python -B code/openwave_wave_evidence/reproduce.py verify --run /path/to/wave-run --output /path/to/verification.json
python -B code/openwave_wave_evidence/test_transport.py
```

The independent check reconstructs 37,908 branch records and 131,072 actual
current records, occupation and clock readouts, 15 analytic populations and
2,000 dispersion cases. Its own implementation is the review's sealed
`independent.py`; it is separate from the producer's computation. Re-execution
by a package maintainer supplies no additional reviewing party.

The full replay executes the free, orbital, unified and boost leaves. Its
1,009,734 native calls include overlapping ancestry. The analytic population
block length of 100,000,000 is not an actual tape length. The finite contact
eigenvalue, all clock outcomes and the biased outer-coin control are retained.
The source review's 100 consistency checks and amendment's Windows/Linux
receipts accompany the code.

## Comparison policy

The sealed amendment's `POLICY.json` requires exact source hashes, inventory,
JSON structure, scalar types and discrete values. Finite floating values obey
`abs(x-y) <= 1e-12 + 1e-12 * max(abs(x), abs(y))`. Each embedded raw-record hash
must match its own decompressed bytes. Compression bytes and floating
serialization may differ. Original scientific thresholds are unchanged.
Optimized Python is rejected; scientific child processes explicitly verify
that their assertions are enabled.

The orbital dependency uses internal reference port 0 and a positive first
spatial ground coordinate. This is a declared readout convention. The
amendment reconstructs its expected diagnostics from the original occupation
records and the polynomial spectral projector. The historical argmax-selected
reference and its failed byte comparisons remain in the sealed evidence.
Adoption of this numerical policy does not follow from a successful replay:
the local literal-byte G4 requirement requires outside adjudication.

`reproduce.py` supplies `--bundle` to the sealed replay and independent checker.
For the sealed 21-test suite and saved-result verifier only, it replaces their
`extract_git` function with a byte-identical local copy operation. All science,
comparators, assertion checks and mutations are unchanged. `test_transport.py`
checks exact input, altered bytes, extra files, manifest substitution, path
escape and optimized entry. The transport checks are developer tests; they
do not constitute an independent adversarial audit.

## Review boundary

The sealed `REVIEW.md` and `PROOF_REVIEW.md` record the separate GPT-6 Astra
reviewer's assessment and independently implemented finite checks. The four
replay findings have an implementation amendment and retained Windows/Linux
receipts. The publication transport and the amendment require independent
adjudication before a submission claims the platform's audit requirement.
Same-group work does not satisfy the local G6 requirement.

The author's complete candidate declaration, its relationship to the pending
static-source proposal, outside-group review and a prospective platform task
are prerequisites for a row decision. This package supplies no static/wave
join, species identification, moving charged apparatus or minimal electric
coupling.

**Scope**: finite records on the declared four-dimensional routing family and
conditional free Gaussian/continuum reconstruction; this package confers no
OpenWave row award or approved candidate pin.
