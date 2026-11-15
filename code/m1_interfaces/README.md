# Sparse reads with variational area and a uniform field action

Fixed planar area agreement leaves a real microscopic defect. Half-residue
partitions of the spaced-ball read graph have macroscopic volume and entropy
tending to zero on the torus. Inside a clipped cube they can beat the
continuum isoperimetric value by the factor 2^(-1/3). Thin connecting paths
do not remove the lack of compactness.

An explicit short-ball addition closes this defect. For

    q=2^(16t), m=2^(8t), K=2^(6t), R=mK, r=2^(7t), t>=1,

take all m-spaced radius-K displacements, all shorter dyadic axis steps,
and every integer vector in the radius-r ball. Every unordered crossing
contributes one nat. Bounded normalized crossing entropy gives
strong compactness of binary partitions. The complete perimeter Gamma limit
holds, with binary recovery and exact finite volume constraints. Finite
minima and near-minimizers converge; inside a cube with sufficient room,
their limits are balls with area (36 pi)^(1/3) V^(2/3).

The threshold also controls the uniform action's hidden modes. For a short
radius r, put s=r/K. When s tends to a positive finite gamma, the extra
alias labels survive with limiting eigenvalues

    4 pi^2 (|j|^2 + gamma^5 |ell|^2),  j,ell in Z^3.

The full positive-mass thermal partition, normal-ordered energy and Gibbs
entropy then approach those of a product of two three-tori. This is a
spectral property of the regulator, not a claim of six physical spacetime
dimensions. When s tends to zero the thermal energy diverges; when s tends
to infinity the extra sector disappears. The last condition is exactly the
one required for cut compactness in this short-ball repair class.

The explicit family above has s=2^t and proves both desired limits. It uses
one common positive action coefficient on every read; no specially weighted
nearest-neighbor action is needed. The full clipped-box free thermal and
smooth Gaussian quantum limits follow from its whole-spectrum bound.
Causal/count limits on interior, non-wrapping diamonds are retained because
every original charged route remains available. Boundary-truncated and
wrapping intervals need their own volumes, not the flat diamond formula.
Logical reads over a fixed horizon grow as q^(71/16), compared
with q^5 for the original dense construction, at the same q. The larger
radius converges more slowly; this is not an equal-accuracy or native cost
claim. Exact finite M1 tuples and physical source/entropy/horizon attachments
are not inferred.

## Reproduce

From the repository root, with the pinned requirements installed:

```sh
PYTHONPATH=code python -m m1_interfaces.build
PYTHONPATH=code python -m m1_interfaces.verify
python -m pytest -q code/m1_interfaces
```

PowerShell can set `$env:PYTHONPATH='code'` before the first two commands.
Builds write a candidate receipt; verification recomputes its scientific
fields independently. The verifier imports neither the producer nor NumPy.
It rejects duplicate JSON keys, nonfinite numbers, incomplete catalogs,
booleans substituted for integer evidence, forged cuts, doubled unordered
counts, false connectivity, missing alias regimes and changed source custody.
An actual optimized Python invocation rejects a semantic forgery.
The field boundary control compares actual shifted-array actions with the
independent unordered-edge energy and fixed-zero exterior potential. It
rejects substituting a clipped graph Laplacian, which has a constant zero
mode and does not support the stated Dirichlet gap.

## Evidence and proof scope

The finite reference has 12 independently reconstructed graph menus and 96
cut cases, on three small generic parameter sets with periodic and clipped
boundaries. These menus contain 78,592 site instances and 3,340,024 ordered
read incidences; all eight label patterns are checked on each graph. This
is a cut census, not a native process history. The generic small controls
do not pretend to be the enormous asymptotic populations.

Complete exact moment and cut certificates at q=65536 compare the raw,
critical and repaired stencils using vertical column sums. Their degrees
are 1,097,965, 2,195,839 and 9,880,701. No q^3 history is stored or claimed
executed. Separate whole yz-plane censuses test the alias eigenvalue
enclosures. Twelve scale levels check the exact parameter identities;
the universal results come from the analytic proof, not extrapolation.

`DERIVATION.md` proves the counterexamples, compactness threshold, perimeter
and exact-volume recovery, minimizer limits, ball-symbol bound and full
thermal conclusions. Eight Lean reductions prove the finite cut-change
bound and the normalization/scale algebra, with transitive standard-axiom
audits. The universal BV, isoperimetric and spectral proofs are analytic.
`CONTRACT.md` states the objective and exit.

The spaced-ball stencil is defined in full here; nothing is imported from
`code/m1_necessity`, which introduces the same stencil in its own comparison.
The fixed-section area and source-causal claims keep their complete
statements and physical premise classifications. This work proves
a stronger interface and field comparison without asserting that the
coordinate cut problem is a physical horizon law.
