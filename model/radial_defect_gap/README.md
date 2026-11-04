# An isolated pair of localized lattice modes

The supplied bounded radial lattice field has exactly two zero modes
when its mass magnitude exceeds the full kinetic band norm. Every other
energy is separated from zero by a positive gap independent of box size.
The modes remain concentrated near the origin, and the finite-box zero
projector converges to its infinite-lattice value at an explicit geometric
rate. The source, proof and conditions are stated in [PROOF.md](PROOF.md).

Choose a rational hopping bound `k>=sqrt(3)*v/a`, with `m>k`, and set
`d=m-k`. The gap is at least d. Each normalized zero orbital has origin
weight at least `1/(1+3v^2/(2a^2 d^2))`. Its local contact-overlap weight
has a volume-independent positive lower bound.

At `a=v=1`, `k=2`, the exact examples are:

| Supplied mass | Gap lower bound | Origin weight lower bound | Contact-overlap lower bound |
| --- | --- | --- | --- |
| `4` | `2` | `8/11` | `64/121` |
| `64` | `62` | `7688/7691` | `(7688/7691)^2` |

These are bounds on the same [bounded radial source](../bounded_defect_index/README.md).
The strong mass condition excludes refinement to zero spacing at fixed
physical mass and velocity. A dynamical origin of the supplied texture
and a volume-uniform interacting particle spectrum require additional
construction.

```sh
python -W error -m pytest -q model/radial_defect_gap
python -W error -m model.radial_defect_gap --write
python -W error -m model.radial_defect_gap --check
```

The [cold receipt](results.json) binds the proof and complete source
manifest. Independent tests use every finite source row, the actual
Schur resolvent, local zero frames and walk expansions. No physical
species assignment, native repair protection or completion gate is
accepted.
