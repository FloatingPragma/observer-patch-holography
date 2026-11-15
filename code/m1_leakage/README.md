# Full-interface leakage and general-noise source refinement

This package removes code preservation from the quantum-noise assumptions of
the declared source. Its existing complete failure/reset services
provide noisy leakage reduction. The [proof](../../docs/research/LEAKAGE_SOURCE_REFINEMENT.md)
also permits a fixed nonzero elementary error below the credited threshold:
growing concatenation depth replaces the parent's decreasing per-location
error, with the additional local drive, inventory and lifetime costs charged.
It also composes with the noisy-record compiler: general Markov
quantum noise, native leakage, noisy runtime control and a live noisy public
archive coexist in one constructed process at fixed finite physical rates.
The [contract](CONTRACT.md) states the scope, evidence requirements and
interpretation of these self-reading source patches and their retained records.

Finite evidence includes:

- Both native proper codes on every full-M6 matrix unit, with public outcomes
  separated from private Kraus indices, and two complete two-carrier gates on
  all 1,296 input matrix units each.
- 26,880 exact integer decoder identities spanning every single-site map from
  a code qubit to its full six-level interface. Twenty-one complete channels
  check erasure, dissipative leakage and coherent leakage on every code site.
- Six continuous full-M6 noise channels, with noncommuting drives and return
  from leaked states; a common-bath memory example; an analytic witness that
  amplitude damping cannot contain a positive ideal identity component.
- A corrected full-carrier accounting effect, a two-site failure witness,
  a nonvanishing fractional-control-error example, and exact level/resource
  inequalities for both refinement regimes.
- A coherent fault-tail identity replacing the archive's probabilistic
  estimate, its exact rational bounds, and the combined six-level resource
  and convergence calculation; recursive replay of the computation code,
  adaptive instruments and 8,100 archive fault cases from that construction.
- A coherent counterexample to an additive bad-region norm bound, with the
  required intersection terms and disjoint location ownership checked. The
  joint-history proof derives its product bound explicitly and counts every
  physical interval once, including export-to-refresh handoffs.

```sh
PYTHONPATH=code python -m m1_leakage.build
PYTHONPATH=code python -m m1_leakage.verify
python -m pytest -q code/m1_leakage
```

PowerShell: set `$env:PYTHONPATH='code'` first. The verifier imports no local
producer. Full channel matrices are serialized sparsely and checked within
`3e-10`; source custody is byte-exact, and integer identities are exact.
It does not interpret a Bell-state trace distance as a diamond norm: the
diamond upper bound follows from the generator estimate in the proof.

The arbitrary-size leakage/general-noise accuracy results are credited
Aliferis--Terhal and AGP theorems; the channel-to-dilation estimate is credited
Kretschmann--Schlingemann--Werner. The finite controls do not execute an entire
universal concatenated machine or measure its threshold. Five fixed levels
suffice for the conservative full-Markov route, four for bounded coherent
couplings; a fixed subthreshold elementary strength uses O(log log q) levels.
Both yield joint decoded/named-history error O(q^-12 log^6 q) and additional
accounting error O(q^-8 log^6 q), within the declared source family.

The combined fixed-rate construction uses the noisy-record parent's larger
measurement-free computation code, six protection levels and its actual
noisy classical archive. A norm fault-path argument handles damping without
inventing an ideal-operation probability. Joint error is O(q^-20 log^16 q),
accounting error O(q^-16 log^16 q), and complete diagnostic storage-time
O(q^6 log^25 q). This is an analytic composition, not a finite simulation of
the entire six-level universal machine.

Reliable runtime classical control is retained only in the separate
fixed-per-operation-error and common-bath variants, not in that combined
fixed-rate result. The static schedule, source/CCG capability and accounting
reference are declared. Extra diagnostic histories and raw counts need
not match the ideal experiment. Unknown raw-input encoding, autonomous clocks
and physical energy or source selection are not claimed.
