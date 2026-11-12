# Meson emission and outgoing volume

The same native bound meson and weak current have an exact short-time pair
emission coefficient `3/512`, independent of outgoing-ring size. A six-state
parent memory equation retains every emitted state and its interaction. The
large two-site conversion peak falls as the outgoing volume increases, so it
does not establish an irreversible lifetime. A 32-site calculation has a
rigorous finite-volume vector error below `2.547e-12` at source time0.01;
floating propagation roundoff is separate.

[Proof](PROOF.md), [complete records](results.json), and
[paper result](PAPER_RESULT.tex) specify the source and physical scope.

Run `python3 -m model.openwave_meson_ring --check` and
`python3 -m pytest -q model/openwave_meson_ring/test_ring.py` from the repository
root with single-threaded numerical libraries.
