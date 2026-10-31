# One common internal representation and actual two-particle exchange

The retained chiral field's internal space is the same
`Lambda^2(C5) + Lambda^4(C5)` module used by `defect_identity`. `bridge.py`
constructs the exact permutation between their actual exterior bases and
checks all twelve source currents with an independent call to the original
Fock exterior derivative. The retained Higgs bilinears transform by
`P.T K P`, not by a Hermitian similarity. The chosen Spin2 and family3
factors then give the same ninety one-particle components as the existing
chiral action. These are permitted design choices. The code does not infer
particle identities from the number ninety.

There are three different constructions here. They must remain distinct.

1. The fixed negative one-particle encoding from `defect_identity` extends
   to the ninety-component space. It matches the one-particle spin central
   sign, and all current operators intertwine exactly.
2. The optional parity-correlated code `E=[I;P]/sqrt(2)`, `P=(-1)^N`, uses
   both ordering-memory sectors. It satisfies `Q E=E P`, restricts to the
   negative encoding on one-particle states, and uses the positive sector
   on vacuum and two-particle states. Even field operations preserve it;
   `(Z_memory tensor c_a)E=E c_a` provides same-code odd CAR operators.
   This is a designed coupling to **global parity**, not a theorem equating
   that parity with exchange of two particles. The two-mode control checks
   that a true mode swap sends the doubly occupied state to minus itself,
   while global parity is positive there. A charged odd interpolator also
   still needs its physical charge reference for an observable read.
3. `colored.py` gives an actual CAR two-particle exchange construction with
   a separate internal state for each occupied port. This is the relevant
   construction when two identical field particles are exchanged.

## The colored native repair law

Let `V=C^d`, and label the source pair by `x<y`. The exact isometry into
the two-particle CAR space is

```text
|{x,y}; a,b>  ->  c†_(x,a) c†_(y,b) |vacuum>.
```

It covers the sector with at most one particle at each of the twelve native
ports. If a native move crosses the occupied spectator, ordering the two
creation operators introduces a minus sign and reverses the internal
slots. If it does not cross, neither operation occurs. Consequently the
actual complete native Kraus map is

```text
K_h = (1/sqrt(2)) sum_x |f_h(x)><x| tensor (-SWAP)^p_h(x).
```

Every branch is a source partial injection times a unitary internal map.
Its effect is the unchanged native domain projector divided by two, so
the complete instrument sums to identity for every seam. All zero and
abort slots remain. The actual three-move closed odd exchange is `-SWAP`;
two identical internal states acquire minus one. `rho(g) tensor rho(g)`
commutes with this swap for every common/global internal action, including the actual
current representation and chosen Spin and family actions. Independent
bitstring creation/annihilation tests check every complete C/E branch at
`d=2` and `d=3`, including relative signs and internal ordering.
Site-local gauge covariance additionally requires the corresponding link
transports; it is not implied by this global representation identity.

The connection to the actual ordered-cover source is also explicit. With
`J=-SWAP`, the isometry `E_col=[I;J]/sqrt(2)` obeys `X E_col=E_col J`.
Thus each source orientation flip becomes exactly the colored CAR exchange
on the encoded pair. The natural design uses a `d^2`-dimensional extra
fiber and retains the orthogonal complement of this code. Its negative
exchange eigenspace is `Sym^2(V)`, dimension `d(d+1)/2`, and its positive
eigenspace is `Lambda^2(V)`, dimension `d(d-1)/2`. At `d=90` these dimensions
are 4,095 and 4,005. Both native characters therefore occur for identical
fermions with internal states. The earlier uniform-negative fifteen-label
construction is a different restricted encoding, not the full two-particle
internal space or a universal fermionic memory bound.

At `d=90`, this sector has dimension `66*90^2=534600`. The full two-particle
CAR space on `12*90` modes has dimension `582660`; the `48060` same-port
collision states are its explicit complement. No law on that complement,
other particle numbers, or a full interacting continuum is claimed. In
particular the old chiral Higgs pairing changes particle number; it cannot
be inserted into this fixed-two-particle repair law without a further
specified extension. Per-particle number-preserving actions do compose.
For a static quadratic Higgs background, a separately verified canonical
Bogoliubov transformation can instead define conserved quasiparticle number.
Its two-quasiparticle sector may then use this colored law, provided all
Nambu partners, zero modes, transformed currents and the chosen reference
are retained. This is a different restriction from dropping the bare-CAR
pairing terms. Time-dependent backgrounds and general gauge/Higgs operations
that create quasiparticle pairs require a further extension.

## Complete instruments and evidence

For any unitary field control `U_h`, the conditional joint operation
`K_h tensor U_h` is complete because its effect is `K_h† K_h tensor I`.
For the parity code the field control must also be even. This construction
does not use the false assertion that exterior second quantization of an
arbitrary Kraus family remains trace preserving. A test shows the failure
already on the vacuum of a one-mode example.

`small_complete_instrument` gives a complete, exact native-record-controlled
field instrument on the full one/two-mode Fock space, with a real field
phase selected by the recorded native coin. It retains every native outcome
and proves branch equality on the correlated code. Exact linear-map
equalities extend to arbitrary untouched references and finite adaptive
histories. The execution module supplies the separately counted native
source-word realization; mathematical completeness alone is not a claim
that every larger field program was run.

From the repository root:

```sh
python3 -m model.defect_field_bridge.bridge --check
python3 -m model.defect_field_bridge.colored --check
python3 -m pytest -q model/defect_field_bridge/test_bridge.py model/defect_field_bridge/test_colored.py
```

The retained receipts are `representation_results.json` and
`colored_results.json`. Physical masses, charge calibration, local spatial
identification, field-particle counting beyond this hard-core sector and
the full Lorentz spin-statistics connection require their own common-law
tests. A calibrated designed action is allowed; its calibration data must
be named as inputs rather than reported as independent predictions.
