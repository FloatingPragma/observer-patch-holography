# Finite Born frame (B11 / issue #687)

This directory gives a bounded exact answer to the finite Born-frame question.
It does **not** prove Gleason's or Busch's theorem and does not claim that the
twelve-port source produces a physical quantum measurement instrument.

## Two different twelve-port objects

The audit found two objects that must not be identified.

1. The source manifest declares twelve pairwise orthogonal central port atoms
   summing to one. They form a single classical twelve-outcome context.
   Normalized weights have affine dimension 11 and determine a unique state on
   the commutative algebra `C^12`. They do not determine a density matrix in an
   ambient full matrix algebra: `I/12` and
   `I/12 + (E_01 + E_10)/24` are distinct positive trace-one matrices with the
   same twelve atom weights.
2. The binary-icosahedral receipt constructs a qubit Pauli/projector adapter
   from source-derived spin geometry. Its own status is “mathematical
   construction from the source-derived spin lift”; physical promotion is
   false. The twelve rank-one projectors form six disjoint antipodal binary
   contexts.

The second object is the strongest available finite projective frame, so the
certificate computes its complete linear answer while retaining its declared
attachment boundary.

## Exact result

Write `w_i` for the weight of projector `i`, with antipode `11-i`. Context
normalization is exactly

```text
w_i + w_(11-i) = 1,  i = 0,...,5.
```

The six equations have rank six. Thus normalized context-additive weights form
a six-dimensional affine space, with nonnegative weights giving `[0,1]^6`.
The contexts do not interlock: every rank-one projector occurs in only one
binary context.

Put `c_i = 2 w_i - 1` and `phi = (1+sqrt(5))/2`. A weight lies in the
trace-one Hermitian/Born slice exactly when

```text
c0 + c1 = phi (c2 + c3)
c2 - c3 = phi (c4 + c5)
c4 - c5 = phi (c0 - c1).
```

These three independent equations leave affine dimension three. When they
hold, tomography is unique:

```text
x = (c4+c5)/2,  y = (c0-c1)/2,  z = (c2+c3)/2.
```

The corresponding trace-one Hermitian matrix is positive precisely when

```text
(phi+2) (x^2+y^2+z^2) <= 1.
```

Therefore:

- a density representation is unique whenever it exists;
- not every normalized context-additive weight has even a Hermitian trace
  representation (dimension six versus three);
- positivity on the twelve listed projectors does not imply matrix
  positivity. The exact certificate includes a represented unit-interval
  weight with `(x,y,z)=(3/5,0,0)` and density norm squared
  `9/10 + 9/50*sqrt(5) > 1`.

This finite family does not derive the Born rule. A positive continuation
would need a source-produced, physically attached, interlocking effect family
rich enough to control the full positive cone.

## Reproduction

From the repository root:

```bash
python3 code/born_frame/finite_born_frame_certificate.py --check
python3 code/born_frame/verify_finite_born_frame_independent.py
python3 -m pytest -q code/born_frame/test_finite_born_frame_certificate.py
cd Lean && lake env lean EventAlgebra/FiniteBornFrame.lean
```

The producer and independent verifier use separate implementations of exact
`Q(sqrt(5))` arithmetic. The frozen output is
`runtime/finite_born_frame_certificate.json`.
