# OPH-native RG representation frontier

This package is the first positive source-side work product for GitHub issue
#32 and work item P594-07. It does not use the external Standard Model packet
from issue #593.

The finite rank-fifteen matter receipt fixes the exact quadratic
representation indices for one matter copy:

```text
sum_Weyl T_SU3 = 2
sum_Weyl T_SU2 = 2
sum_Weyl dim(SU3) dim(SU2) Y^2 = 10/3
```

After importing the standard four-dimensional one-loop QFT functional, the
coefficient law in the convention
`d g_i/d ln(mu) = b_i g_i^3/(16 pi^2)` is

```text
b_Y = (20/9) N_g + (1/6) N_H
b_2 = -22/3 + (4/3) N_g + (1/6) N_H
b_3 = -11 + (4/3) N_g
```

The evaluation `(N_g,N_H)=(3,1)` gives
`(41/6,-19/6,-7)`. It is a declared conditional completion, not a physical
OPH selection. The registered finite reduct admits three, four, and five
family copies. Its scalar grammar admits zero, duplicate, and inert-doublet
completions. The rank-three screen band has no completed physical matter-pole
attachment.

A direct summand with zero color, weak, and hypercharge indices changes none
of the three one-loop gauge coefficients. This exact result only quotients the
gauge-index calculation. A complete W/Z census must still include every field
with a W/Z, Yukawa, scalar, or mass-mixing vertex. A summand may be omitted
only after its complete interaction and mixing vertices are proved zero.

## Status

```text
PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN
```

The following objects are deliberately represented as `not_emitted`, never
as empty lists or trusted booleans:

- ordered EFT intervals;
- threshold locations;
- decoupling maps;
- finite scheme maps;
- Jacobians;
- finite-order term masks;
- certified vector remainders.

The complete issue remains open. A source-complete action and light/heavy
census for every W/Z-coupled field and operator modulo proved zero-vertex
decoupling, physical family, scalar, and local carrier attachment, mass
spectrum, scheme selection, and remainder-producing RG engine are required
for closure. These positive inputs are owned by #634, #569, #630, #631, and
#632.

## Evidence

- `data/source_rg_policy_v1.json` is the closed source allowlist and import
  boundary.
- `build_rg_representation_frontier.py` is the deterministic producer.
- `check_rg_representation_frontier.py` independently resolves every source
  pin and redoes the arithmetic without importing producer code.
- `outputs/rg_representation_frontier.json` is the committed receipt.
- `Lean/Screen/RGRepresentationFrontier.lean` checks the rational identities
  and copy-shift laws without axioms.
- `tests/test_rg_representation_frontier.py` rejects source-pin drift, target
  ancestry, normalization and sign mutations, multiplicity promotion,
  matching self-attestation, scheme-witness collapse, and Lean drift.

Run:

```bash
python3 run_all.py
```
