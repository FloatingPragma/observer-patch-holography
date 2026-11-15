# Filled source vacuum and energy

The declared source uses bounded observer-like self-reading patches: local
states, ports, readback, feedback or repair operations, and public record
banks. The receipt and independent checker form its public evidence bundle.

The [derivation](../../docs/research/SOURCE_SEA_ENERGY.md) connects the actual massive
walk to a positive excitation generator on a deterministically prepared
filled vacuum. It proves boundary and local continuum control, a finite-energy
resolved clock read, and fixed-strength protection with the longer preparation
included. The [contract](CONTRACT.md) fixes the deliverables and interpretation.

The finite analytic clock witness has probability error below 0.038 and
observable swing above 0.923. Its enormous three-dimensional apparatus is
specified by bounds, not represented as an executed numerical experiment.
The actual finite controls have deliberately separate roles:

| Control | What is executed |
| --- | --- |
| Source spectrum | Eight complete Bloch matrices; all 256 Fock states for three cases; principal versus alternative generator branches |
| Reflecting boundary | Complete compiled operators on 1, 8 and 27 cells, both masses, all central covariance columns |
| Vacuum limit | 24 full eight-component projections, both carrier/valley labels, including the rejected ordinary-Dirac high-band replacement |
| Preparation/read | Four-mode source graph, routed rotations, both even/odd preparations, every loop syndrome, both returned QND outcomes and all spectators |
| Packet read | 256, 512 and 1024-site periodic controls, four times, two read phases, all Fourier modes, explicit sea background and excitation energy; whole-pattern translation compared with the wrong laboratory-anchored probe |
| Analytic witness | Outward interval bound for a compact three-dimensional reflecting experiment, including normalization, tails and boundary margin |
| Noise/resources | Exact integer schedules and polynomial exponents for the long preparation and complete diagnostic lifetime |

The small source graph tests the general preparation algorithm; it is not
called a four-mode discretization of the three-dimensional massive vacuum.
The periodic packet controls have a normalized constant transverse mode;
they are not evidence for infinite-volume Gaussian tails. Those tails and
the continuum limit have separate analytic bounds.

## Reproduce

Install the repository's pinned `requirements.txt`. From the repository root:

```powershell
$env:PYTHONPATH='code'
$env:OPENBLAS_NUM_THREADS='1'
python -m m1_sea_energy.verify
python -m pytest -q code/m1_sea_energy
```

On a POSIX shell, prefix commands with `PYTHONPATH=code OPENBLAS_NUM_THREADS=1`.
`python -m m1_sea_energy.build` produces a candidate, independently replays
its semantics and the inherited receipts, and then writes source/claim pins.
Do not regenerate a receipt merely to bless changed hashes.

The checker imports no producer. It reconstructs the logarithm using two
Hermitian eigensystems instead of a complex Schur decomposition; constructs
Fock generators directly from occupation bits; replays actual edge-code
rotations and the native QND helper; propagates packet states by binary
powers instead of eigensystem powers; translates probes with a Fourier
multiplier instead of copying the position-space carrier formula; and checks
interval and integer bounds. The complete QND channel's excitation-energy
change is also bounded on the full Fock space, retaining both outcomes.
Tests block producer imports under `python -O`, reject structural and semantic
corruptions, and accept occupied-orbital gauge changes and isometry phases.
Every positive scalar parameter uses a relative comparison, including tiny
boundary errors. Matrix roundoff tolerances are not treated as noise thresholds.

Custody recursively verifies the massive-clock and fixed-strength
archive results. The receipt stores small matrices and covariance columns;
full finite operators and Fock matrices are recomputed instead of archived.

## Interpretation

This is a filled-vacuum experiment, with a different physical read from
the empty-vacuum click detector. The generator is exactly positive after
normal ordering and produces the same integer-time source channels. It is
not a claim of unique laboratory energy: two explicit Floquet branch changes
preserve those channels and positivity while changing energy assignments.
No interacting-sea stability, autonomous schedule, source selection from
A1--A3 or physical work cost of the supplied controls is inferred.
