# Higgs and Yukawa source frontier

This package is the bounded, nonpromoting work product for issue #630. It
collects the finite source information shared by the RG, Higgs/top, and W/Z
producer lanes. It does not emit a positive scalar or Yukawa action.

The source allowlist contains eight hash-bound inputs:

- the exact non-gating #763 source-history custody projection: authenticated
  finite source order at independently generated cutoffs with exact induced
  restrictions between adjacent cutoffs;
- the finite port-current algebra from #566;
- the rank-fifteen chiral matter module and conditional invariant channels
  from #314;
- the finite global-form packet from #567;
- the conditional rank-three family packet from #569;
- the scalar multiplicity boundary from #616 and #617;
- the limited scalar-chain boundary from #623;
- the exact current-reduct Higgs/top separation implementation from #521.

No particle reference store or #593 calculation packet is in the ancestry.
The independent checker resolves every pin and repeats the structured target
firewall.

The #763 projection is non-gating receipt ancestry. It contributes only
authenticated finite source order and independent-cutoff induced-restriction
custody. It supplies no physical causal attachment, faithful embedding,
3+1-dimensional manifoldlikeness or dimension, count--volume calibration, or
smooth continuum. Issue #634 supplies a bounded finite local field and operator
domain. It does not supply a continuum Lorentzian/Spin quantum EFT. That
transfer is owned by #635.

## Exact conditional classification

On the declared one-doublet, rank-three family branch, with four-dimensional
power counting through degree four and a positive local scalar kinetic form,
the scalar and Yukawa operator space is classified modulo field-independent
constants and total derivatives:

```text
kinetic normalization:       one positive-real coordinate
scalar potential:            two real coordinates
each invariant family channel: Mat_3(C)
all three Yukawa channels:   Mat_3(C)^3, complex dimension 27
```

This is an operator-basis theorem under named premises. No point of this
coefficient space is selected. The packet contains no kinetic coefficient,
potential coefficient, Yukawa entry, vacuum coordinate, FJ map, mass, or
physical prediction.

The family-basis and field-normalization quotients remain open because the
finite parents do not select the required bases or equivalence maps.

## Completion-fiber controls

The receipt retains four exact boundaries:

- an empty scalar extension and a declared one-doublet extension share the
  finite parent projection and differ in scalar multiplicity;
- two formal points of `Mat_3(C)^3` share the finite parent projection and
  differ under a family-matrix coordinate projection;
- two invertible local coordinate one-jets share the finite parent projection
  and differ under the `v_chart` to `v_F` coordinate projection;
- on the #521 current reduct and `0 < u < 1`, the linear and Born lifts have
  exact formal separations `2u(1-u)` and `u(1-u)/2`.

The first three witnesses are completion-space statements. The #521 row is a
dimensionless formal readout statement. None is a positive physical source
action.

## Status

```text
BOUNDED_NONPROMOTING_FRONTIER__POSITIVE_SOURCE_ACTION_OPEN
```

The open positive objects and their owners are:

- #636: physical scalar carrier and multiplicity, canonical scalar kinetic
  normalization, scalar potential coefficients, and a stable vacuum branch;
- #637: complete `Yu`, `Yd`, and `Ye` matrices;
- #638: the scheme- and scale-defined `v_chart` to `v_F` map.

Issue #630 integrates the three final packets. The physical family and matter
attachment remains in #569 and is consumed by #637. The finite #634 domain is
resolved context, not an open positive dependency.

## Replay

```bash
python3 code/particles/hierarchy/higgs_yukawa_source_frontier/build_higgs_yukawa_source_frontier.py
python3 code/particles/hierarchy/higgs_yukawa_source_frontier/check_higgs_yukawa_source_frontier.py
python3 -m pytest -q code/particles/hierarchy/higgs_yukawa_source_frontier/tests/test_higgs_yukawa_source_frontier.py
```

The mutation suite rejects parent hash drift, target-path injection, physical
promotion, coefficient-value fields, collapsed completion witnesses, changed
#521 separation identities, coefficient-space dimension drift, status drift,
policy escape hatches, and subject-digest drift.
