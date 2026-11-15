# One calibrated source clock does not select light bending

This package gives a **constructive negative result** on an explicit positive
local two-field source-action class. It supplies an
analytic classification, actual finite source histories, independent replay
and a retrospective comparison to published gravitational measurements.
The [derivation](../../docs/research/PAIRED_SOURCE_GRAVITY.md) states the theorem and
its scope; [CONTRACT.md](CONTRACT.md) fixes the objective and exit.

After normalizing its clock susceptibility, every positive two-field edge
coefficient has the form
`[[1+kappa*gamma^2,-kappa*gamma],[-kappa*gamma,kappa]]`, kappa>0.
The full static clock response, Gaussian clock marginal and clock dynamics
are independent of gamma. The result survives the nonlinear resting-matter
source coupling `j exp(-u)`, with identical proper-time coupling for sources
and clocks. Yet the same universal metric/photon dictionary gives bending
amplitude `(1+gamma)/2` relative to GR. The conformal gamma=-1 source bends
no light, exactly; gamma=1 gives the Einstein coefficient. Positivity and
symmetry do not choose between them. Linearized physical exterior `G_ij=0`
selects gamma=1 in this class; scalar vacuum and the linear Ward identity do not.

The scalar support and twelve-port Gram addresses are inherited from the
existing source implementation. The additional fields, massless source
action, supplied update program and physical interpretation are hypotheses,
not a claim of native integer-repair admission or A1-A3 source selection.
The actual finite detectors consume versioned local source records; their
numerical clock and optical-index reads are not raw solar observations.

## Reproduce

From the repository root, with the pinned `requirements.txt` installed:

```sh
PYTHONPATH=code python -m paired_gravity.build
PYTHONPATH=code python -m paired_gravity.verify
PYTHONPATH=code python -m pytest code/paired_gravity
```

In PowerShell, set `$env:PYTHONPATH='code'` before the Python commands.
No network, simulator checkout or raw-data download is needed for replay.
The receipt contains final fields and hashes of every event, not a large
routing tape. Sixteen executions retain both source strengths, schedules,
transverse initializations and source laws, and all four constitutive
members. Rebuilding regenerates every write/read before checking the result.

The verifier independently reconstructs the operator from primitive Gram
tables, rounds using integer bisection rather than the producer's floating
seed, encloses nonlinear updates with a different exponential series,
recomputes the complete source chains and uses a maximum-principle residual
bound. It reconstructs linearized curvature symbolically and integrates the
optical force and nonlinear rays in different coordinates from the producer.
Source digests supplement these semantic checks. They do not replace them.
All explicit rejection checks remain active with `python -O`.

## Comparison and uncertainties

[DATA.md](DATA.md) records the primary measurement provenance and exposure.
All four source members predict unit static redshift ratio at laboratory
weak-field order. Pound-Snider report 0.9990 with statistical deviation
0.0076 and systematic bound 0.010; the latter is not a Gaussian variance.
The designated VLBA reduction gives total bending amplitude
0.99988 +/- 0.00016 relative to GR. Its other three overlapping reductions
are retained without treating them as independent data. Only gamma=1 among
the four controls lies within that designated one-deviation interval.
That agreement uses the supplied metric/probe law and is not independent
evidence selecting OPH's physical source. The exact gamma=-1 zero-bending
alternative demonstrates the unresolved prediction before any data fitting.

Exact finite residuals, outward exponential intervals, a compact source
with a C^3 potential, a continuum consistency estimate and a nonlinear asymptotic ray
remainder are provided. The last is not relabelled as a finite-observer
VLBA error budget. Original resonance counts, visibilities and likelihoods
are not available in this packet. Its completed exit is the constructive
nonidentifiability theorem, not a positive end-to-end physical postdiction.

## Requirements and retained obligations

| Requirement | Artifact and result |
| --- | --- |
| Static source, action, boundary, local records | Derivation sections 1-3; both source laws executed on the 64-site parent support |
| Inputs and common calibration | Same Schur normalization, source, units and probe dictionary throughout; physical maps explicitly supplied |
| Paired reads and spatial curvature | Frequency ratio and Fermat index derived from the same metric; joint amplitude range is all real values at first order as gamma varies |
| Error and useful comparison | Finite residual/clock bounds, continuum limit, asymptotic ray bound; published compressed comparison with stated limitations |
| Alternative-law test | Full positive constitutive classification, identical entire clock marginal, distinct optical outcomes, nonlinear source control |
| Measurement ancestry and GR | DATA.md and immutable primary transcription; same-input GR reference and all four overlapping reductions |
| Independent verification | Full source/read replay, independent curvature/ray calculations and hostile-input tests |
| Negative exit | Constructed nonidentifiability, with the exact exterior tensor selector and obligations below |

Physical curvature/stress matching and the Einstein branch, selection of the
source's spatial response or a justified exterior vacuum/null-stress
condition, mass/clock/unit calibration, their common physical history and the
source/outcome interpretation are open obligations (issues #729,
#779, #736, #740 and #730). This result discharges none of them and
does not falsify their additional physical Einstein premises. It establishes that more clock-only sampling on
this source cannot resolve the remaining spatial-response coefficient.
