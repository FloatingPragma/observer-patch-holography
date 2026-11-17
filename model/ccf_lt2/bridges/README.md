# Finite bridges for CCF-LT-2

These calculations use the local [official source](../../SPEC.md), candidate
`A-ccf2-local-capacity-v2`. They distinguish full native polynomial identities,
literal native basis columns, complete smaller matrix diagnostics, and
conditional physical theorems. No sibling checkout is required.

Run `python -m model.ccf_lt2.bridges.verify --output model/ccf_lt2/bridges/RESULTS.json`
after the full native runtime receipt is current. Omit `--output` for strict
replay. The receipt authenticates the runtime source and complete bridge
inputs before and after execution; assertions must be enabled.

| Calculation | Exact scope |
| --- | --- |
| [Kinematics](kinematics/bridge.py) | Complete native charge currents, hard-wall Ampere identity, compact quadratic symbol, CAR/spin and actual mass compatibility. |
| [Phase instrument](kinematics/phase.py) | Three-outcome complete PVM, actual full initial native column, slopes −9/32,+9/32,0; retained complement and dephasing/Dirac-deletion controls. |
| [Matter](matter/bridge.py) | Full H_S charge/baryon/lepton/parity certificate, disconnected source additivity, colour-singlet seed and complete moment columns. |
| [Reference](reference/check.py) | Exact full Gauss-zero support rank, beta=1/1024 positive Gibbs floor, regional entropy fixture and actual stationarity countervariation. |
| [Coframe](gravity/coframe.py) | Constant real background family of H_S, six strain and three shift derivatives. |
| [Refinement](gravity/refinement.py) | Angular and radial shell chain identities and geometric mesh bound. |
| [Energy Ward](gravity/ward.py) | Complete local-capacity matrix partition, currents and coarse cancellation; analytic full-source identity in [PROOF.md](PROOF.md). |

The [finite-source Lean module](kinematics/FiniteSource.lean) and
[sector module](matter/SectorBridge.lean) state their external matrix-extraction
hypotheses. Compile from `Lean/` with `lake env lean ../model/ccf_lt2/bridges/kinematics/FiniteSource.lean`
and the corresponding matter path. The combined `python -m model.ccf_lt2.bridges.verify_lean`
checks both against [LEAN_RESULTS.json](LEAN_RESULTS.json). They do not formalize Python or identify
physical particles. [SCALING.md](SCALING.md) is an analytic design, without a
numerical family receipt. [EINSTEIN_BRIDGE.md](../EINSTEIN_BRIDGE.md) gives the
remaining same-source physical hypotheses.
