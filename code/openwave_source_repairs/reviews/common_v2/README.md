# Independent common-source replay

Run with Python3 and NumPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B check.py ../selected_source_common_v2
```

The default verifies saved evidence without writing. The explicit `--record`
mode creates a result once. The checker imports no producer code and fetches
nothing. `REVIEW.md` states its scientific scope and exact execution boundary.
`PREREG.json`, `RESULT.json` and the final seal binding identify the reviewed
producer. All review is same-group Codex; no other-group G6 is claimed.
