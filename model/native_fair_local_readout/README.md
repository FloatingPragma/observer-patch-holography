# Complete local fair-source readout

Native E records recover every nonconstant initial occupation pattern on a
connected local window. The decoder tracks every incident current, including
boundary repairs, and acquires initial contrasts when chosen spanning-tree
edges first appear. It retains terminal disturbance and explicit coverage
and empty/full erasure flags across all local number sectors.

The complete output with any untouched reference differs from an initial
occupation-label oracle followed by the same native instrument by exactly
`p_coverage + (1-p_coverage)*eta`, where `eta` is the input weight of the
empty/full local subspace. For a canonical global number sector, `eta` has
an exact hypergeometric formula. Uniform selector coverage has certified
finite bounds; the fixed intracell tree needs eleven attempts. No success
branch is normalized and no empty/full coherence is silently measured.
An optional source-bound global-total check resolves constant patterns when
one is infeasible; globally empty or full sources need no tree coverage.

The uniform alternative is admitted at growing levels within the existing
attempt ceiling. Finite routed windows, original expanded execution limits,
and bounded exact raw-word queries have distinct domains. The code does not
supply a coherent service embedding or a physical clock. The proof and
failure controls are in [PROOF.md](PROOF.md) and
[test_readout.py](test_readout.py).

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --strict-markers --strict-config model/native_fair_local_readout
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_fair_local_readout.packet --check
```
