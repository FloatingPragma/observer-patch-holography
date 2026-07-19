# Correctable Public-Record Capacity

Exact finite contracts, controls, and archived candidate audits for the OPH
capacity map. The canonical Pro5 producer is

```text
N = log M_0(U_N),
F_set,r,epsilon(D) = {M_epsilon(q): q in Omega_tilde(r,D)},
M_0(q) = alpha(G_q).
```

The first line is the official universe-level closure equation; the next two
lines are its typed finite implementation.

The operational resolution, electroweak/Higgs bridge, and measured
cosmological constant are independent downstream comparisons. They never
define the direct map. The physical N closure remains open.

## Canonical lane

- [`F_READBACK_SPEC.md`](F_READBACK_SPEC.md) is the Pro5 acceptance contract:
  complete terminal fiber, atom readouts, endogenous reachability, frozen
  publicness, global joint kernels, compound confusability graph, exact and
  approximate correctable capacity, carrier saturation, scalarization,
  refinement, finite-size slack, receipts, and controls.
- [`correctable_public_record_capacity.py`](correctable_public_record_capacity.py)
  evaluates finite public checkpoint packets. It computes global
  sections, exact maximum independent sets, receipt-scale worst-input
  approximate capacities, support semigroups, carrier bounds, terminal-fiber
  scalarization, no-new-confusability, greatest fixed points, and unique
  slack-zero certificates.
- [`test_correctable_public_record_capacity.py`](test_correctable_public_record_capacity.py)
  covers saturation, cyclic permutation, joint-coupling nonidentifiability,
  approximate capacity, ambiguous fibers, order countermodels, target taint,
  and carrier failures.
- [`reversible_public_checkpoint_packet.py`](reversible_public_checkpoint_packet.py)
  supplies the finite twelve-port icosahedral reference packet requested by
  the proof program. It verifies every checkpoint generator is a permutation,
  certifies `M_0(q)=|X_reach(q)|`, and emits the exact rank-one saturation
  receipt. It is a schema/control, not the missing physical source packet.
- [`test_reversible_public_checkpoint_packet.py`](test_reversible_public_checkpoint_packet.py)
  checks the 12 vertices, 30 interfaces, exact reversible capacity identity,
  noninjective failure, and target-taint failure.
- [`public_record_capacity.py`](public_record_capacity.py) and its tests retain
  the superseded Pro4 checkpoint-fixed projection branch as a control. A cyclic
  permutation proves that it is not the canonical capacity definition.

The identity family fixes every dimension; the erasure family fixes only the
bottom dimension. Monotonicity and deflation therefore do not select the
cosmic value. Physical closure requires a source-derived complete packet and
an exact finite-size slack law with one regulator-stable zero.

## Independent geometric estimator

- [`operational_readback_contract.py`](operational_readback_contract.py)
  evaluates frozen scale-discrimination errors, preserves pre/post-checkpoint
  accounting, requires complete-fiber agreement for `rho_op`, and compares
  `log M_0` with `pi/rho_op^2` only after direct capacity exists.
- [`test_operational_readback_contract.py`](test_operational_readback_contract.py)
  checks discrimination endpoints, all-coarser thresholding, capped error
  accounting, complete-fiber agreement, and independence controls.

The diagnostic residual is

```text
R_rho = log M_0 - pi/rho_op^2.
```

Defining `rho_op` from `M_0`, or using capacity to select its protocol,
invalidates the comparison.

## Downstream bridges

After robust direct closure:

- identifying the correctable record carrier with the de Sitter horizon may
  identify `log D_star` with
  `A/(4 ell_star^2)` and yield `Lambda ell_star^2=3*pi/N_star`;
- a positive refinement-natural carrier map may identify the source-normalized screen load with
  the four-copy weak load and test
  `N_bridge=pi*exp(6*pi/(P*alpha_U(P)))`.

Neither bridge constructs the direct map. The exterior package proves the weak
multiplicity four, but that integer alone does not identify a physical load
carrier.

## Archived lanes

The dated construction notes and `F_candidate_*.py` files preserve historical
count, affine, and Banach candidates. They have diagnostic value only. The
`CP*` and `G2_GAP_1` notes likewise do not supply the Pro5 packet or the exact
finite-size selector.

## Open gates

- source-derived complete capacity-indexed trial universes and terminal fibers;
- total atom readouts, endogenous reachability, and frozen publicness;
- complete global joint checkpoint kernels with marginal consistency;
- compound-graph/MIS and carrier-projection receipts;
- no-new-confusability extension and refinement embeddings;
- an exact finite-size slack law with one regulator-stable zero;
- independent horizon, EW/Higgs, and operational-resolution bridge tests.

## Usage

```bash
python3 -m pytest test_correctable_public_record_capacity.py -q
python3 -m pytest test_reversible_public_checkpoint_packet.py -q
python3 -m pytest test_operational_readback_contract.py -q
python3 -m pytest test_public_record_capacity.py -q
```
