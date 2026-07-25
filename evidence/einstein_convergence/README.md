# Evidence E1: Einstein-Cone Convergence Ladder

Compressed primary evidence for the Lorentzian-signature emergence
measurements cited by the papers as [E1]. Each rung is one deterministic run
of the fixed federated source capture; the stored artifact carries the event
chart, the causal and spacelike pair samples, the fitted quadratic form, and
the flux and normalization vectors needed to recheck every printed number
without rerunning the simulation. `manifest.json` binds every file by sha256.

## Provenance

- Simulator: [oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim),
  commits `07e2faca` (three-rung ladder) and `4f0169c` (density-corrected top rung).
- Producer: `scripts/einstein_convergence_ladder.py` (deterministic; seed
  20260751; canonical capture path with `observer_cross_reads`,
  `snapshot_coverage=spanning`, `geometry_transport=held_out_flow`).
- Full per-rung configuration is embedded in each `rung_*.json` summary.
- Reproduction: `.venv/bin/python scripts/einstein_convergence_ladder.py`
  in the simulator repository regenerates every artifact bit for bit.

## Measured ladder

| Rung | Observers | Support | Cross edges | Held-out inertia | Cone margin | Coupling spread |
| --- | --- | --- | --- | --- | --- | --- |
| 16,384 | 128 | 96 | 348 | (1,3) | -5.62 | 0.1875 |
| 65,536 | 256 | 96 | 312 | (1,3) | -3.22 | 0.1860 |
| 262,144 | 512 | 96 | 312 | (2,2) | -2.49 | 0.1766 |
| 262,144 | 512 | 384 | 1,062 | (1,3) | -1.41 | 0.1766 |

The support-adjusted path selects the first, second, and fourth rows. Its
held-out signature is Lorentzian (1,3), its cone margins are -5.62, -3.22,
and -1.41, and its coupling spread decreases from 0.1875 to 0.1766. The path
does not hold observer density, support width, or cross-edge density fixed.

The two 262k rows form the direct control. At fixed carrier count, observer
count, chain depth, and seed, increasing support width from 96 to 384 changes
the recorded cross-observer edges from 312 to 1,062 and the held-out inertia
from (2,2) to (1,3). This measures sensitivity to the support and cross-read
structure. It does not isolate a unique density variable, establish an
invariant-density convergence law, or determine an infinite-scale limit.

## Claim boundary

These are finite measurements of the declared repair dynamics under frozen
instruments with a same-size support-width control (see the simulator's
Einstein-branch guide). They license statements about the measured event
forms on the four archived configurations. They do not by themselves
construct the continuum Einstein equation, whose named receipts are tracked
in the open program ledger.
