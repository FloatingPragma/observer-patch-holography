# Start Here: The Skeptical Technical Entry Point

Observer Patch Holography asks whether public physics can be reconstructed
from finite observers that compare overlapping records and repair
disagreement.

This repository contains:

- a sorry-free Lean 4 library of more than 800 checked theorems and lemmas, each public
  theorem carrying a per-theorem axiom report;
- reproducible finite simulations with an adversarial negative control;
- machine-certified interval and uniqueness certificates for the declared
  numerical maps;
- an explicit list of open proof and physical-construction gaps.

Nothing here asks you to accept the full theory before checking an artifact.
Start with one theorem, one receipt, or one counterexample.

## The Hypothesis In Two Sentences

Objective reality is modeled as the fixed point of finite observer patches
that keep records, compare them on overlaps, and repair disagreement until the
same record is recoverable from either side. The claim under test is that
spacetime signature, gauge structure, and specific dimensionless constants are
forced by that consistency requirement, on stated premises, without a table of
adjustable constants.

## Three Things You Can Check

### 1. Build the Lean library

Requires [elan](https://github.com/leanprover/elan); the toolchain is pinned
in `Lean/lean-toolchain`.

```bash
cd Lean
lake exe cache get
lake build
```

The library build is the proof receipt. Coverage and premises are indexed in
[Lean/docs/PROOF_INDEX.md](../Lean/docs/PROOF_INDEX.md).
Continuum geometry, asymptotic tails, and physical identification are explicit
premises, stated as such in [Lean/README.md](../Lean/README.md).

### 2. Reproduce one numerical certificate

The pixel closure $P=\varphi+\sqrt\pi/A_T(P)$ is claimed to have exactly one
fixed point per declared map on the declared domain. The certificate is
interval arithmetic, so a failure is a hard failure.

```bash
python -m venv .venv && . .venv/bin/activate   # CPython 3.12+
pip install -r requirements.txt
cd code/P_derivation
python3 global_uniqueness_certificate.py --mp-dps 40 --iv-dps 40 \
    --su2-cutoff 120 --su3-cutoff 90 --initial-pieces 256
python3 -m pytest test_global_uniqueness_certificate.py -q
```

The gauge-width root sits $2.5\times10^{-6}$ from the measured
$\alpha^{-1}=137.035999177(21)$. The remaining gap has a stated address, the
open hadronic transport term; see
[code/P_derivation/README.md](../code/P_derivation/README.md) for the claim
boundary.

### 3. Run the finite core and one negative control

From the repository root, with the same virtual environment:

```bash
python tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
python3 code/particles/calibration/strict_one_loop_pole_map/run_all.py
```

The last command regenerates a conditional receipt, runs its adversarial
suite, validates both JSON Schemas, and checks the receipt without importing
the producer. The simulation negative control lives in
[evidence/einstein_convergence](../evidence/einstein_convergence/): at
constant coupling density the held-out event form carries Lorentzian signature
$(1,3)$ at 16k, 65k, and 262k carriers with the cone margin halving per rung,
and the density-starved control run degrades the signature to $(2,2)$ on cue.
Every stored number is sha256-bound and regenerates bit for bit from the
companion [oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim)
repository.

## What Is Proved, Conditional, Measured, And Open

| Status | Contents |
| --- | --- |
| Proved, machine-checked | More than 800 Lean theorems and lemmas: finite consensus core, gauge identifiability, Einstein-branch composition, and a negative result against the naive version of the program's own claim |
| Proved, certified numerics | Fixed-point existence and uniqueness for each declared $P$ map on the declared domain, via interval contraction and adaptive subdivision certificates |
| Measured | The Einstein-cone convergence ladder above, with its adversarial density control |
| Conditional | The strict one-loop W/Z pole map (proved and machine checked, fixture is a post-exposure regression), the $N$ extension, the $N$–Higgs bridge, and the $N_g=3$ selection; each conditional branch names its open producer |
| Open | Physical Thomson transport for $P$, the physical $N$ packet, physical family attachment, the common-domain gravity tower, and the full list in the [README](../README.md#open-proof-obligations-and-falsification-boundary) |

The [claim registry](../claims/claim_registry.yaml) links prose claims to
artifacts and is checked by `tools/check_claim_registry.py`. The
[selection ledger](SELECTION_LEDGER.md) lists exactly what is proved and
exactly what is left.

## Five Scoped Contribution Opportunities

1. **Model counting for $N$.** The reversible public-checkpoint packet
   reduces record capacity to $M_0=|X_{\rm reach}|$, computable by exact CSP
   or model counting. An independent recomputation, or a faster counter, is a
   self-contained project in [code/capacity_readback](../code/capacity_readback/).
2. **Break a certificate.** The uniqueness certificate covers
   $\alpha\in[0.005,0.01]$ with forward-mode interval AD and tail majorants.
   Find a declared map and a subdomain where the Lipschitz verdict or a
   majorant is wrong. One verified counterexample is a publishable result
   against the program.
3. **Shrink a Lean premise.** Pick any public theorem, read its axiom report,
   and remove or weaken a premise; or formalize one of the named open lemmas
   in the [proof index](../Lean/docs/PROOF_INDEX.md).
4. **Extend the convergence ladder.** Reproduce the $(2,2)$ density control,
   or run a rung above 262k carriers in
   [oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim). The
   extrapolated zero crossing of the cone margin sits in the low millions of
   carriers and is a stated projection, unmeasured.
5. **Bound the hadronic transport term.** The $2.5\times10^{-6}$ gap between
   the $P$ closure and the measured $\alpha^{-1}$ is attributed to an open
   Thomson transport contribution. An independent QED/QCD bound on that term,
   in either direction, moves the claim.

## Licensing

All software, including the Lean library, `code/`, and `tools/`, is licensed
under [Apache-2.0](../code/LICENSE). Papers, the book, documentation,
figures, and data are [CC BY-NC-SA 4.0](../LICENSE). Commercial use, forks,
and derived libraries of the software are welcome under the Apache terms.

## Find An Error

The fastest way to contribute is to falsify something. The
[falsification program](OPH_FALSIFICATION_PROGRAM.md) collects the mature
claims with their evaluation boundaries; the
[common objections document](COMMON_OBJECTIONS.md) records the objections
raised so far and their current state. If a build fails, a certificate does
not verify, a receipt does not regenerate, or a stated premise hides a
stronger assumption, open an issue on the
[tracker](https://github.com/FloatingPragma/observer-patch-holography/issues)
with the command you ran and its output. Counterexamples are treated as
contributions, with the same standing as proofs.
