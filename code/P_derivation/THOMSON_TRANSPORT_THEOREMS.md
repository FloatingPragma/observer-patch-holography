# Thomson Transport Theorem Suite

This file states the theorem package needed to turn the `P/alpha`
fixed-point witness into a measured fine-structure derivation.

The theorem package is conditional. The lepton one-loop transport is
implemented as a numerical kernel. The endpoint package computes the residual
inverse-alpha packet that a source-only Ward-projected transport theorem must
emit. The source-spectral reduction theorem is emitted in
`SOURCE_SPECTRAL_THEOREM.md` and
`runtime/source_spectral_theorem_current.json`. The Ward-projected spectral
measure payload, the electroweak scheme remainder, and the interval certificate
remain theorem burdens.

## Objects

Let

```text
P = phi + alpha * sqrt(pi)
```

and let the D10 source solve emit

```text
Q_src(P) = (alpha_2(m_Z;P), alpha_Y(m_Z;P), alpha_3(m_Z;P), v(P), a0(P))
```

where

```text
a0(P) = alpha_em^-1(m_Z^2;P).
```

The desired Thomson endpoint is

```text
alpha_Th^-1(P) = a0(P) + Delta_Th(P).
```

The missing source-only theorem is a derivation of `Delta_Th(P)` from the same
Ward-projected `U(1)_Q` family as `a0(P)`.

## Theorem 0: No Detuning-Only Thomson Bypass

**Status:** closed no-go.

The outer law

```text
P = phi + sqrt(pi) / A
```

defines the scalar fixed-point coordinate. It does not identify `A` with the
Thomson endpoint unless the Ward-projected `U(1)_Q` endpoint map has been
emitted on the same source branch.

Proof sketch. The D10 source branch fixes `a0(P)`, the electromagnetic inverse
coupling at `m_Z^2`. The Thomson endpoint is

```text
A_T(P) = a0(P) + Delta_Th(P).
```

The hadronic part of `Delta_Th(P)` is a subtracted dispersion functional of a
positive Ward-projected spectral density `rho_had(s;P)`. Ward transversality,
positivity, and the outer detuning equation constrain the form of the current
correlator, but they do not determine that spectral density. Two admissible
spectral densities can share the same D10 anchor and outer detuning coordinate
while giving different zero-momentum transport integrals. Therefore the
detuning coordinate alone cannot be promoted to `alpha_Th^-1(P)`.

## Theorem 1: Ward-Projected Source Lock

**Status:** closed at criterion level in the paper, not a populated transport
calculation.

Let the D10 source family be fixed by the forward pixel solve, and let the
realized low-energy charge operator be

```text
Q = T_3 + Y.
```

Let the electroweak transport kernel be

```text
K_D10^EW(q^2;P) = {Pi_AA, Pi_AZ, Pi_ZZ, Pi_WW}
```

with a Ward projector `W_Q` onto the unbroken electromagnetic channel. If

```text
W_Q[Pi_AZ(0;P)] = 0
```

and the projected `U(1)_Q` edge-sector probabilities satisfy

```text
p_n(q^2;P) proportional to exp(-t_Q(q^2;P) n^2),
```

then the electromagnetic readout on this lane is

```text
alpha_em^-1(q^2;P) = 8*pi^2 / t_Q(q^2;P).
```

Therefore

```text
alpha_Th^-1(P) = lim_{q^2 -> 0} 8*pi^2 / t_Q(q^2;P).
```

This theorem identifies the correct lane. It does not by itself compute
`t_Q(0;P)`.

## Theorem 2: Leptonic One-Loop Source Transport

**Status:** implemented with the closed-form one-loop kernel, conditional on the
Stage-5 charged-lepton mass emitter.

Assume the source branch emits charged-lepton masses

```text
m_e(P), m_mu(P), m_tau(P)
```

and the same Ward-projected scheme is used at `m_Z^2` and at the Thomson
endpoint. Then the perturbative charged-lepton contribution is

```text
Delta_lep(P) =
  sum_f K_f(m_Z(P)^2; m_f(P), Q_f^2 = 1, N_c = 1)
```

where

```text
K_f(Q^2;m_f,Q_f^2,N_c) =
  (2 N_c Q_f^2 / pi) * integral_0^1 x(1-x)
  log(1 + Q^2 x(1-x) / m_f^2) dx.
```

This is the exact one-loop kernel used by `paper_math.py`.

What remains to promote this to theorem grade:

1. prove that the Stage-5 lepton mass emitter is the source-only emitter on the
   same branch used by `a0(P)`;
2. evaluate the closed form in a directed-rounding interval backend;
3. prove the scheme match between `a0(P)` and the zero-momentum readout.

## Theorem 3: Hadronic Spectral Transport

**Status:** source-spectral reduction theorem emitted; source measure payload
absent.

The quark part cannot be theorem-grade if it is only a free-quark sum multiplied
by a simple screening factor. A source-only theorem must emit a positive
Ward-projected hadronic spectral density

```text
rho_had(s;P) >= 0
```

for the electromagnetic current-current correlator on the same D10 branch.

The required theorem is:

If `rho_had(s;P)` is emitted by the source branch, satisfies Ward positivity,
has the correct threshold support, and matches the high-energy quark/OPE tail
of the same `U(1)_Q` current, then the hadronic transport contribution is the
subtracted dispersion transport

```text
Delta_had(P) = Integral[ W_had(s, m_Z(P)^2) * rho_had(s;P) ds ],
```

with the subtraction chosen so that the same `a0(P)` scheme is used at
`m_Z^2`.

The reduction theorem replaces the free-quark screened ansatz once a populated
source measure exists. The local corpus contains the schema and contract, but
does not contain the finite-volume levels, Ward-projected residues, current
normalization, continuum pushforward, or systematics needed for the numerical
spectral moment.

Constructive implementation target:

- `thomson_endpoint_contract.py`
- `thomson_endpoint_package.py`
- `source_spectral_theorem.py`
- `runtime/thomson_endpoint_contract_current.json`
- `runtime/thomson_endpoint_package_current.json`
- `runtime/source_spectral_theorem_current.json`
- `../particles/hadron/ward_projected_spectral_measure.schema.json`

Workers should not return obstruction-only text for this branch. If the
free-quark screened route fails, the required replacement is a populated
Ward-projected spectral-measure export or a code/schema patch that moves that
export toward the endpoint builder.

## Theorem 4: Electroweak Matching Remainder

**Status:** open as a bound.

A full endpoint theorem also needs a residual matching term

```text
Delta_EW(P)
```

covering the difference between the D10 electroweak anchor convention and the
zero-momentum electromagnetic Thomson convention. The theorem must show either

```text
Delta_EW(P) = 0
```

in the declared scheme, or provide a source-only formula and an explicit error
bound.

## Theorem 5: Full Thomson Endpoint

**Status:** endpoint package closed for issue #223; source-residual
non-identifiability boundary closed for issue #235. Exact-alpha promotion
requires the populated source spectral measure payload and interval
certificate.

If Theorems 2, 3, and 4 are closed on the same source family and scheme, define

```text
Delta_Th(P) = Delta_lep(P) + Delta_had(P) + Delta_EW(P).
```

Then

```text
alpha_Th^-1(P) = a0(P) + Delta_Th(P)
```

is the source-only Thomson endpoint on the Ward-projected `U(1)_Q` lane.

At the runtime report root,

```text
a0(P) = 128.308268045165213892552005990181778935...
```

The public endpoint value requires

```text
Delta_Th(P) = 8.727731131834786107447994009818221065...
```

The implemented exact one-loop continuation gives

```text
Delta_impl(P) = 8.686567119456435565397988595605414327...
```

so the source-only theorem must account for

```text
Delta_missing(P) = 0.041164012378350542050005414212806738...
```

without using the measured endpoint as an input.

The measured endpoint can be used as a temporary calibration insert for
numeric tables and plots. The calibration artifact is

```text
runtime/measured_endpoint_calibration_current.json
```

and it carries `external_input_used=true`, `promotion_allowed=false`, and
`exact_alpha_promoted=false`. It does not close this theorem. Its purpose is to
keep the displayed numerical surface coherent while the nonperturbative
Ward-projected source spectral payload remains absent.

The endpoint package also evaluates the same question at the pixel obtained by
the outer equation from the CODATA/NIST comparison value:

## Theorem 6: Corpus-Limited Non-Identifiability of `R_Q(P)`

**Status:** closed no-go for the current corpus; issue #235 is closable as
first missing lemma isolated, not as an alpha promotion.

Let the current source corpus consist of the OPH axioms, realized Standard
Model branch, D10 source packet, outer detuning equation, Ward-projected
`U(1)_Q` lane criterion, and implemented closed-form baseline `B(P)`. This
corpus determines `a0(P)` and the admissible electromagnetic lane. It does not
determine a unique source residual

```text
R_Q(P) = Delta_had_src(P) + Delta_EW_src(P)
         - Delta_q_screened_impl(P).
```

Reason. The hadronic endpoint is a subtracted dispersion functional of a
positive Ward-projected spectral measure:

```text
Delta_had_src(P) = integral K_P(s) dmu_had^Q(s;P).
```

The current-corpus no-go is constructive. In the dimensionless coordinate
`y=s/mZ(P)^2`, remove the common positive factor `1/(3*pi*mZ(P)^2)` from the
Thomson kernel and write

```text
k(y) = 1/(y*(1+y)).
```

The two positive atomic completions

```text
mu_A = delta(y-2)
mu_B = delta(y-3)
```

have identical projection to the emitted current-corpus invariants: the D10
source family, `U(1)_Q` lane, `a0(P)`, `mZ(P)`, lepton kernel, and the absence
of populated finite-volume levels, residues, continuum pushforward, and
same-scheme remainder. They even share the same total positive weight. But

```text
integral k dmu_A = 1/6,
integral k dmu_B = 1/12.
```

Thus exact `alpha` is not determined by the current corpus. Any exact endpoint
theorem must emit the missing measure/remainder data; it cannot infer them from
the emitted invariant packet.

Ward transversality fixes the tensor form of the current correlator, positivity
places `dmu_had^Q` in a positive cone, threshold support restricts its support,
and the high-energy tail gives asymptotic conditions. Any finite package of
such constraints can be written as finitely many linear conditions

```text
integral f_i(s;P) dmu(s;P) = c_i(P).
```

Choose more support points than there are constraints. There is a nonzero
signed atomic perturbation `nu` that leaves all listed constraints fixed. Unless
the endpoint kernel `K_P` lies in the span of the constraint kernels, which the
current corpus does not prove, `integral K_P dnu` is nonzero. For a sufficiently
small perturbation around a strictly positive baseline measure, both
`mu + epsilon nu` and `mu - epsilon nu` remain admissible but produce different
Thomson endpoint transports.

The runtime screening witness gives the same obstruction in finite form:

```text
S_lambda(P) = 1 - x(P) + lambda x(P)^2,
x(P) = N_c alpha_3(m_Z;P) / pi.
```

The source packet hash is independent of `lambda`, while different `lambda`
values change the endpoint transport at the scale of the required residual.
Therefore neither `S_required`, `c_Q`, nor a public-endpoint scalar residual can
stand alone as an OPH derivation.

The emitted reduction theorem is:

```text
WardProjectedHadronicSpectralEmission_Q:
  OPH axioms + realized SM branch + D10 source packet
  + source spectral measure payload
  -> Delta_had_src(P), Delta_EW_src(P), scheme lock,
     quadrature/tail bound, derivative bound.
```

Only after the measure payload is populated can the conditional interval
certificate be upgraded to a theorem-grade proof for

```text
G(P) = phi + sqrt(pi) / (B(P) + R_Q(P)).
```

This closes the #235 boundary as a blocker-isolation theorem. It does not close
the exact fine-structure derivation.

```text
P_C = 1.630968209403959324879279847782648941...
a0(P_C) = 128.307965473286248209961108741756716187...
Delta_required(P_C) = 8.728033703713751790038891258243283813...
Delta_impl_exact(P_C) = 8.686567842708528400985442542885969768...
Delta_source_residual(P_C) = 0.041465861005223389053448715357314044...
```

Equivalently, if \(x=N_c\alpha_3(m_Z;P_C)/\pi\) and the implemented quark
screen is \(1-x\), the endpoint package requires

```text
S_required(P_C) = 0.895400132647658797805800283181670641...
c_Q(P_C) = 0.658025759927155435638230170232360050...
```

where \(c_Q=(S_{\rm required}-(1-x))/x^2\). This is a compact scalar form of
the missing Ward-projected QCD screening and endpoint-remainder map. It is not a
source-only derivation of that map.

## Theorem 6: Fixed-Point Closure With Certified Transport

**Status:** open.

Once `Delta_Th(P)` is emitted, define

```text
G(alpha) =
  1 / (a0(phi + alpha * sqrt(pi)) + Delta_Th(phi + alpha * sqrt(pi))).
```

A theorem-grade fine-structure derivation requires an interval `I` such that

```text
G(I) subset I
sup_{alpha in I} |G'(alpha)| < 1
```

or another certified uniqueness argument. Then the fixed point

```text
alpha_* = G(alpha_*)
```

exists uniquely in `I`, and

```text
P_* = phi + alpha_* * sqrt(pi)
```

is the source-only OPH pixel closure.

`fixed_point_certificate.py` is a local numerical certificate for
the implemented map only. It is not this final theorem.

## Promotion Rule

No theorem in this file may be promoted to a fine-structure endpoint derivation unless
all of the following are true:

1. `rho_had(s;P)` is emitted from the OPH source branch, not imported from a
   measured endpoint.
2. The same renormalization and matching scheme is used by `a0(P)` and
   `Delta_Th(P)`.
3. The transport integral has a certified numerical error bound.
4. The final fixed-point map has an interval-level existence and uniqueness
   certificate.
5. The public reference value is kept out of the source solver.
