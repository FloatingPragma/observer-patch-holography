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
define the direct map. The universe-level physical `N` closure remains open.

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
- [`public_record_csp.py`](public_record_csp.py) is the exact
  constraint-propagating global-section backend. It is extensionally
  equivalent to Cartesian enumeration, but it model-counts the connected
  twelve-observer, twenty-four-atom source packet without exploring `24^12`
  assignments.
- [`test_correctable_public_record_capacity.py`](test_correctable_public_record_capacity.py)
  covers saturation, cyclic permutation, joint-coupling nonidentifiability,
  approximate capacity, ambiguous fibers, order countermodels, target taint,
  and carrier failures.
- [`reversible_public_checkpoint_packet.py`](reversible_public_checkpoint_packet.py)
  retains the finite twelve-port icosahedral reference control. It verifies
  every checkpoint generator is a permutation, certifies
  `M_0(q)=|X_reach(q)|`, and emits exact rank-one saturation. It remains
  explicitly nonphysical.
- [`test_reversible_public_checkpoint_packet.py`](test_reversible_public_checkpoint_packet.py)
  checks the 12 vertices, 30 interfaces, exact reversible capacity identity,
  noninjective failure, and target-taint failure.
- [`source_derived_public_checkpoint_packet.py`](source_derived_public_checkpoint_packet.py)
  defines the first source-derived fixed-cutoff physical packet for issue #548.
  It freezes the carrier to the twelve edge-center ports with reversible
  write/check orientation, so `D=|P_12 x {write,check}|=24`. The producer emits:
  - a complete 67-world structural one-fault trial manifest with one terminal
    world, fully materialized candidates, SHA-256 completeness receipts, and an
    output-blind membership predicate;
  - total observer/interface atom readouts and 24 endogenous reachability
    histories;
  - a frozen universal publicness policy;
  - the complete 40-element `D5 x C2 x C2` joint checkpoint family, all 1,600
    support-relation compositions, and independently checked local marginals;
  - an empty compound graph, a 24-record independent set, inverse decoders,
    worst-input and total-variation receipts;
  - 24 orthogonal rank-one carrier projections, proving `M_epsilon <= 24` and
    exact saturation `M_0=24`;
  - empty/incomplete/ambiguous/singleton fiber controls; isomorphism, cyclic,
    alternative-coupling, tiny-noise, circular-definition, taint, identity,
    erasure, and finite-suffix controls; and separate extension/refinement
    no-new-confusability injections with negative controls.
- [`test_source_derived_public_checkpoint_packet.py`](test_source_derived_public_checkpoint_packet.py)
  checks the full issue #548 acceptance surface.
- [`ISSUE_548_SOLUTION.md`](ISSUE_548_SOLUTION.md) maps every acceptance item
  to the executable receipt.
- [`public_record_capacity.py`](public_record_capacity.py) and its tests retain
  the superseded Pro4 checkpoint-fixed projection branch as a control. A cyclic
  permutation proves that it is not the canonical capacity definition.

The identity family fixes every dimension; the erasure family fixes only the
bottom dimension. Monotonicity and deflation therefore do not select the
cosmic value. The `D=24` artifact is the first physical finite packet,
not an exact finite-size slack law or a universe-level selection theorem.

## Generate the issue #548 receipts

```bash
python3 source_derived_public_checkpoint_packet.py --output-dir runtime
```

This writes the complete terminal-fiber manifest, public checkpoint packet,
and certificate as canonical JSON.

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
`CP*` and `G2_GAP_1` notes likewise do not supply the exact finite-size
selector.

## Open gates after issue #548

- extend the source construction from the first frozen `D=24` packet to the
  declared capacity-indexed regulator family;
- construct an exact finite-size slack law with one regulator-stable physical
  zero;
- independently certify the horizon-record, EW/Higgs load-carrier, and
  operational-resolution bridges;
- supply public hardware-realization evidence if a carrier implementation is
  claimed.

## Usage

```bash
python3 -m pytest test_correctable_public_record_capacity.py -q
python3 -m pytest test_public_record_csp.py -q
python3 -m pytest test_reversible_public_checkpoint_packet.py -q
python3 -m pytest test_source_derived_public_checkpoint_packet.py -q
python3 -m pytest test_operational_readback_contract.py -q
python3 -m pytest test_public_record_capacity.py -q
```
