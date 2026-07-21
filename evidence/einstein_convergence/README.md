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

At constant coupling density the held-out signature is Lorentzian (1,3) at
every rung and the cone margin halves per rung (-5.62, -3.22, -1.41):
geometric convergence toward the Einstein cone, with the coupling spread
decreasing monotonically beside it. The unscaled-support 256k row is the
density control: with cross-observer edges stagnant while observers double,
the signature degrades to (2,2), isolating coupling density as the mechanism
of cone merging. Extrapolating the measured halving, the margin crosses zero
in the low millions of carriers; that extrapolation is a stated projection,
not a measurement.

## Claim boundary

These are finite measurements of the declared repair dynamics under frozen
instruments with adversarial controls (see the simulator's Einstein-branch
guide). They license the statements the papers make about measured signature
emergence and monotone convergence; they do not by themselves construct the
continuum Einstein equation, whose remaining named receipts are tracked in
the open program ledger.
