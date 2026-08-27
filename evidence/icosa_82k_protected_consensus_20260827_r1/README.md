# Protected-authority 81,920-row consensus archive

This directory contains the claim-bearing data for the OPH-FPE run
`icosa_82k_protected_consensus_20260827_r1`. The run was produced from clean
checkouts of
[`oph-physics-sim` at `ce17921`](https://github.com/muellerberndt/oph-physics-sim/tree/ce17921eb7504106fef1ba445e1349b5367aa676)
and the research repository at `ee137dea6fd48487e30863db53bb787d2e5f45e8`.
The simulator commit, research commit, configuration, seed streams, primitive
arrays, terminal arrays, reports, patch state, committed records, and observer
rows are bound by `archive_manifest.json`.

## Result

The finite regulator uses the level-six cell rung of the geodesic
icosahedral refinement tower. It has 81,920 patch rows and 122,880 undirected
cell adjacencies. Every row contains twelve local port slots. Three slots are
routed through the cell adjacencies; nine remain exposed or reserved.

The repair kernel receives one immutable, pairwise-distinct signed integer
authority for every patch row. The authority array is created from the named
`repair_authority` seed stream before repair begins and is included in the
protected source hash. On each inconsistent oriented seam, the endpoint with
the larger authority value is preserved. The other endpoint receives its
unique gauge-transported matching label. Gauge links remain fixed.

The archived source has 102,415 gauge-covariant mismatches. Every enabled
repair removes exactly one mismatch. All sixteen shuffled replays accept
102,415 moves, finish with zero mismatches, and produce the same
authority-bound terminal hash:

```text
sha256:824a17d3315cb1ab703480d70fbb1368e319d6851c4e94014456a33c2f33fafb
```

The replay reports no strict-descent, global-increase, local-diamond,
disjoint-commutation, frame-covariance, move-contract, repair-completeness, or
link-mutation violation. It checks 512 shared-node pairs, 512 disjoint pairs,
and 16 local-frame relabelings. The exact argument is stronger than the sample:
every routed edge slot has a single authority-selected normal form, and the
archived patch map confirms that no routed carrier-port slot is shared by two
seams.

The driven trace applies 81,920 repairs in cycle zero and the remaining 20,495
in cycle one. All 81,920 patch records are committed by cycle eight and remain
committed. The archive retains the full twelve-coordinate patch state, 2,048
observer neighborhoods of size 96, and two cap-control views.

## Independent replay

`verify_archive.py` imports no simulator code. It constructs the six elements
of (S_3) directly from permutations, recomputes the initial mismatch set,
derives the authority-selected terminal state from the primitive source
arrays, checks the final state and full patch-state custody, and recomputes the
source, quotient, authority, protected-source, and terminal hashes. It also
checks every archived byte and decompresses the observer stream.

With Python 3, NumPy, and the `zstd` executable installed:

```bash
python3 verify_archive.py
```

## Scope

This is an exact finite edge-slot consensus result conditional on the archived
authority-decorated source. The present axioms do not select that authority
order. A fixed authority realization does not remain pointwise invariant
under arbitrary patch relabeling; covariance holds when authority travels
with the patch. A refinement-natural authority source has not been
constructed.

The run uses the declared identification between patch rows and the cells of
the icosahedral support. Its physical carrier-to-support realization receipt
is false. The exact normal form does not cover additional patch-local
constraints that couple several port slots, continuum dynamics, or laboratory
physics. The observer rows certify record custody and neighborhood readout.
They cause no later port writes, and the run's observer-like causal-feedback
receipt is false. Literal record-conditioned feedback is established by the
separate
[`vertex12_signed_record_feedback_20260827`](../vertex12_signed_record_feedback_20260827/)
archive and has not yet been integrated into this large run.

## Reproduce

In a clean clone of the simulator:

```bash
git checkout ce17921eb7504106fef1ba445e1349b5367aa676
python3 -m pip install -e '.[dev]'
python3 -m oph_fpe.cli run-bw-array \
  --config configs/icosa_82k_protected_consensus.yml \
  --out-dir runs
```

The seed is `20260827`. Reproductions should compare the archived primitive
and terminal hashes and disclose platform or dependency changes.
