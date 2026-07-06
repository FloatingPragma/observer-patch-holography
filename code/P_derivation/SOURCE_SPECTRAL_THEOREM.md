# Ward-Projected Source Spectral Theorem

Status: source-spectral reduction theorem emitted. The source spectral measure
payload is absent from the current corpus.

## Statement

Fix the D10 source family and the realized electromagnetic charge

```text
Q = T3 + Y.
```

Let the Ward-projected current-current correlator be

```text
Pi_Q^{mu nu}(q;P) = (q^mu q^nu - q^2 eta^{mu nu}) Pi_Q(q^2;P).
```

If the OPH source branch emits a positive `U(1)_Q` hadronic spectral measure
`rho_Q(s;P)`, a same-subtraction scheme lock tying the measure to
`a0(P)=alpha_em^-1(m_Z^2;P)`, and a source electroweak remainder
`Delta_EW_src(P)`, then the source-only Thomson endpoint is

```text
A_T(P) =
  a0(P)
  + Delta_lep_src(P)
  + mZ(P)^2/(3*pi) * integral rho_Q(s;P)/(s*(s+mZ(P)^2)) ds
  + Delta_EW_src(P).
```

The fixed-point lane is

```text
P = phi + sqrt(pi)/A_T(P).
```

An interval certificate for this map gives a unique pixel value and the
corresponding fine-structure constant.

## Proof

Ward projection selects the transverse electromagnetic channel, so the current
two-point function has the scalar form above. Positivity follows from the
source Hilbert-space spectral representation. A once-subtracted dispersion
relation between the `m_Z` anchor and the zero-momentum endpoint gives the
kernel

```text
mZ(P)^2/(3*pi*s*(s+mZ(P)^2)).
```

The charged-lepton part is the same closed one-loop kernel used by the D10
transport code. The confined hadronic part is not a screened free-quark sum. It
is the Ward-projected spectral moment of the source-emitted color-singlet
measure. The electroweak remainder must be a zero theorem in the declared
scheme or a source-bounded interval.

With those objects supplied on the same D10 branch, `A_T(P)` is a source-only
function. The Banach or monotonicity certificate then gives uniqueness of the
pixel fixed point.

## Current Corpus Boundary

The local hadron files define the production export contract, but the finite
volume levels, Ward-projected residues, continuum pushforward, and systematics
are not populated. The artifact

```text
code/P_derivation/runtime/source_spectral_theorem_current.json
```

therefore records the reduction theorem and blocks exact-alpha promotion from
the current corpus.

Fitted scalars such as `c_Q`, `S_required`, or a residual inverse-alpha packet
are rejected as source inputs. They are target diagnostics, not spectral data.

## Constructive No-External-Input Closeout

The current corpus cannot determine exact `alpha` without a populated source
spectral payload. The limitation has an explicit countermodel.

Use the dimensionless spectral coordinate

```text
y = s/mZ(P)^2
```

and remove the common positive factor `1/(3*pi*mZ(P)^2)` from the Thomson
kernel. The remaining kernel is

```text
k(y) = 1/(y*(1+y)).
```

Consider two positive atomic Ward-projected measure completions:

```text
mu_A = delta(y-2)
mu_B = delta(y-3).
```

Both project to the same invariants emitted by the corpus: the D10
family, the `U(1)_Q` lane, the anchor `a0(P)`, `mZ(P)`, the lepton kernel, and
the fact that no finite-volume levels, current residues, continuum
pushforward, or same-scheme source remainder have been emitted. They also have
the same total positive weight.

Their Thomson moments differ:

```text
Integral k dmu_A = 1/6
Integral k dmu_B = 1/12.
```

Therefore the current emitted invariant packet does not determine the Thomson
endpoint functional. Exact fine-structure promotion requires a source-emitted
spectral measure and same-scheme remainder; a scalar fitted to the measured
endpoint cannot replace those objects.

## Why the Hadronic Payload Is Unavailable Here

The missing object is nonperturbative QCD data in the Ward-projected
electromagnetic channel. It is not produced by increasing decimal precision in
the D10 scalar map. The required payload must contain:

- nonempty finite-volume vector-channel levels on the source family
- Ward-projected current residues or weights with current normalization
- positivity, threshold support, and the pushforward to `rho_Q(s;P)`
- statistical, continuum, finite-volume, chiral, current-matching,
  quadrature, and endpoint-remainder bounds
- a same-subtraction scheme lock to `a0(P)`
- a directed-rounding interval certificate for the fixed point

The local workspace contains contracts and skeletons for this route. The
production finite-volume correlator dump and the derived spectral-measure
artifact are absent. The viable source route is a first-principles
unquenched QCD/HVP computation on the OPH-emitted family, followed by the
interval certificate above.

## Empirical Display Surface

For numeric tables and plots, the public endpoint surface is carried through

```text
code/P_derivation/runtime/measured_endpoint_calibration_current.json
```

It uses

```text
alpha^-1(0) = 137.035999177(21)
P = 1.63096820940395932487927984778...
```

from the outer equation. This row class is OPH plus empirical hadron closure,
with a separate `e+e- -> hadrons` payload class for the empirical hadronic
endpoint contribution through the policy in `../../HADRON.md`. The source-audit
payload remains recorded in the runtime ledgers and theorem contract files. It
cannot satisfy the source-spectral theorem gate.
