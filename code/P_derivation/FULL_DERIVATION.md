# Full P/Alpha Derivation Contract

This note records the complete derivation shape for the OPH `P <-> alpha`
closure. It also records the part that is not closed by the code.

The implementation is a reproducible numerical witness. It is not a measured
fine-structure derivation. The open term is the same-family transport
from the electroweak source anchor at `m_Z^2` to the Thomson limit.

## Target Value

For comparison only, the 2022 CODATA/NIST value is

```text
alpha(0)    = 7.297 352 5643(11) x 10^-3
alpha^-1(0) = 137.035 999 177(21)
```

NIST source: https://physics.nist.gov/cuu/pdf/wall_2022.pdf

This value must not enter the solver. It may enter only after the OPH map has
emitted a candidate.

## Outer Closure

The pixel ratio is the screen-cell area in Planck units:

```text
P = a_cell / l_P^2
```

The outside reading of the cell is

```text
P = phi + alpha_in(P) * sqrt(pi)
```

where `alpha_in(P)` is the electromagnetic observation strength emitted by the
inside readout of the same cell.

For a proposed inverse fine-structure value `A = alpha^-1`, the outer equation
alone gives

```text
P(A) = phi + sqrt(pi) / A
```

Using the CODATA/NIST central value as compare-only metadata gives

```text
P_compare = 1.630968209403959324879279847782648941335982851627925...
```

## D10 Source Map

For a trial `P`, the code implements the D10 forward map:

```text
M_U(P)    = E_P * exp(-2*pi) * P^(1/6)
E_cell(P) = E_P / sqrt(P)
```

Then it solves the one-dimensional D10 pixel-closure equation for `alpha_U(P)`:

```text
ellbar_SU(2)(t2) + ellbar_SU(3)(t3) - P/4 = 0
```

with

```text
t2 = 4*pi^2 * alpha_2(m_Z; P)
t3 = 4*pi^2 * alpha_3(m_Z; P)
```

The same source branch gives the electroweak source anchor:

```text
a0(P) = alpha_em^-1(m_Z^2; P)
```

At the p80 closure solution, the tracked report has

```text
P                         = 1.630972095856551047687367577695292870658496913060878...
a0(P)                    = 128.308268057804920945587248022104657547038365445844...
alpha_U(P)               = 0.04112424744187026746915764925464448456993064904099...
alpha_1(m_Z; P)          = 0.01688566757069415835339363997836739395680047865877...
alpha_2(m_Z; P)          = 0.03377781410924628652300163613093328591835071186996...
alpha_3(m_Z; P)          = 0.11833586195787320830210183862144293744230714961799...
```

## Thomson Transport

A measured-alpha derivation needs the source-locked Thomson endpoint:

```text
alpha_Th^-1(P) = a0(P) + Delta_Th(P)
```

The implementation uses an internal Stage-5 charged-spectrum
continuation and an exact one-loop fermion transport kernel:

```text
K_f(Q^2;m_f) =
  (2 N_c Q_f^2 / pi) * integral_0^1 x(1-x)
  log(1 + Q^2 x(1-x) / m_f^2) dx
```

with a simple quark screening factor:

```text
1 - N_c * alpha_s(P) / pi
```

That gives

```text
Delta_impl(P) = 8.686567144452067333731552158679540632253260231830378...
alpha_impl^-1 = a0(P) + Delta_impl(P)
              = 136.994835202256988279318800180784198179291625677674792...
```

The closure residual is small:

```text
alpha_fixed_point_residual = -1.1689e-29
```

So the fixed-point algebra has converged.

## Compare-Only Gap

To hit the CODATA/NIST central value from the same source anchor, the transport
term would need to be

```text
Delta_required(P) = 137.035999177 - a0(P)
                  = 8.727731119195079054412751977895342452961634554155586...
```

The transport term is short by

```text
Delta_missing(P) = 0.041163974743011720681199819215801820708374322325208...
```

Equivalently,

```text
alpha^-1 gap = 300.388036648990592569245139235606485151919638673... ppm
P gap        = 0.000003886452591722808087729912643929322514061433...
```

This is too large to call a precision match.

## Checks That Pass

The comparison does not indicate a failure of the outer equation:

```text
P = phi + alpha * sqrt(pi)
```

It also does not indicate an obvious failure of numerical convergence. The
the p80 report contains 80 bisection steps and a tiny fixed-point residual.

The missing piece is the low-energy transport/readout term
`Delta_Th(P)`.

## What A Full Derivation Must Add

`THOMSON_TRANSPORT_THEOREMS.md` states the theorem suite for the missing layer.
The summary is below.

A full measured-alpha derivation must replace the structured-running
ansatz with a source-only transport theorem. The required object is a map

```text
Delta_Th(P)
```

derived from the same Ward-projected `U(1)_Q` source family as `a0(P)`, with no
CODATA endpoint and no imported compare bundle.

The missing theorem has to specify:

1. The renormalization scheme and matching convention that connects the D10
   source anchor `a0(P)` to a zero-momentum electromagnetic readout.
2. The charged lepton contribution from the same source branch.
3. The quark and hadronic vacuum-polarization contribution, including threshold
   and confinement effects, without replacing it by a fitted screen factor.
4. Any electroweak matching terms that belong to the chosen scheme.
5. An interval or certified numerical bound showing that the resulting
   `alpha -> alpha` map has the selected root.

Only after that object is closed should the final comparison be made:

```text
F(alpha) =
  1 / (a0(phi + alpha * sqrt(pi)) + Delta_Th(phi + alpha * sqrt(pi)))
  - alpha

F(alpha_*) = 0
```

Then

```text
alpha_*^-1
```

can be compared with `137.035999177(21)`.

## Audit Command

After generating a full report, run:

```bash
python3 alpha_gap_audit.py --report runtime/full_p_alpha_report_p80.json
python3 thomson_endpoint_package.py --report runtime/full_p_alpha_report_p80.json
python3 transport_theorem_manifest.py --report runtime/full_p_alpha_report_p80.json
```

The command prints the implemented transport term, the required compare-only
transport term, the missing inverse-alpha contribution, and the theorem-status
manifest. This makes any future replacement of `Delta_Th(P)` easy to check
without letting the measured constant feed the solver.

`thomson_endpoint_package.py` adds the conditional endpoint packet. At the
CODATA-mapped pixel point it reports

```text
P_C = 1.630968209403959324879279847782648941...
a0(P_C) = 128.307965473286248209948959819190019918...
Delta_required(P_C) = 8.728033703713751790051040180809980082...
Delta_impl_exact(P_C) = 8.686567867734823108913580310939963101...
Delta_source_residual(P_C) = 0.041465835978928681137459869870016982...
S_required(P_C) = 0.895400127551185647132725678585532880...
c_Q(P_C) = 0.658025360816792342502465198049036592...
```

The scalar \(c_Q\) is defined by \(S_{\rm required}=1-x+c_Qx^2\), with
\(x=N_c\alpha_3(m_Z;P_C)/\pi\). It is the compact endpoint target for a
source-only Ward-projected QCD screening and endpoint-remainder map.

## Status

```text
closed:   D10 source map P -> a0(P)
closed:   outer/inner numerical fixed-point witness for the implemented map
closed:   endpoint-package blocker isolation for issue #223
open:     source-only same-family residual map for Delta_Th(P), tracked by issue #235
open:     interval-wide proof for the final full transport map, tracked by issue #235
```
