# Full P/Alpha Derivation Contract

This note records the complete derivation shape for the OPH `P <-> alpha`
closure and the source-audit packet attached to the executable code.

The public fine-structure surface is the local fixed point of the same screen
cell read from the outside and from the electromagnetic inside channel. The
runtime code also keeps the finite source-audit trunk visible, including the
same-family transport remainder from the electroweak source anchor at `m_Z^2`
to the Thomson limit.

The public numeric row has output class `oph_plus_empirical_hadron_closure`.
The source-only audit row has output class `source_only_oph`. The two rows are
kept separate.

## Fine-Structure Readout

The OPH fixed-point readout is

```text
alpha(0)    = 7.297 352 5643(11) x 10^-3
alpha^-1(0) = 137.035 999 177(21)
```

Reference for the empirical endpoint value:
https://physics.nist.gov/cuu/pdf/wall_2022.pdf

## First-Principles Computation

The computation has one local unknown, the pixel ratio `P`. The source chain is:

```text
phi = (1 + sqrt(5)) / 2
P = a_cell / l_P^2
outer detuning = (P - phi) / sqrt(pi)
P -> M_U(P) -> alpha_U(P) -> alpha_i(m_Z;P) -> a0(P)
a0(P) -> A_T(P) by Ward-projected U(1)_Q Thomson transport
P = phi + sqrt(pi) / A_T(P)
alpha(0) = 1 / A_T(P)
```

The public root is

```text
P          = 1.630968209403959324879279847782648941...
alpha(0)  = 0.007297352564331425030245795264691683...
alpha^-1  = 137.035999177(21)
```

The source-only audit row emitted by the available transport package is

```text
P_source_audit          = 1.63097209585889737696451390350695562985390
alpha_source_audit^-1   = 136.994835177412937295289429464436887879658
```

The OPH plus empirical hadron closure row uses a separate `e+e- -> hadrons`
payload class for the hadronic endpoint contribution. The source-only audit row
is recorded separately.

The electroweak-hierarchy bundle in `../particles/hierarchy` uses these same
two branch surfaces. On the public endpoint branch it records

```text
alpha_U(P) = 0.041124336195630495
v/E_star   = 2.0199803239725553e-17
```

On the source-audit branch it records

```text
alpha_U(P_source_audit) = 0.04112424744557487
v(P_source_audit)/E_star = 2.0198114150099223e-17
```

The bundle records an earlier unconverged forward run. Converged precision-100
reruns supersede that source-audit branch beyond digit ~10 (CL-6, closed):
the certified forward values are
P_fwd = 1.630972095858897..., alpha_U(P_fwd) = 0.041124247441816685..., and
v(P_fwd)/E_star = 2.0198114078576331e-17. The public endpoint branch is
unchanged.

The `R_U` certificate proves a unique local source zero in the supplied
interval by Krawczyk inclusion and a derivative enclosure excluding zero. The
public-endpoint burden is the source-side Thomson transport object described in
this directory.

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

Using the fixed-point readout gives

```text
P = 1.630968209403959324879279847782648941335982851627925...
```

The public display surface uses

```text
alpha^-1(0) = 137.035999177(21)
alpha(0)   = 0.007297352564331425030245795264691683...
P          = 1.630968209403959324879279847782648941...
```

The finite code audit reports which source-side transport packet remains
outside the implemented D10 trunk.

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

At the tracked runtime candidate, the report has

```text
P                         = 1.63097209585889737696451390350695562985390...
a0(P)                    = 128.308268057987597347904057614084743934330...
alpha_U(P)               = 0.0411242474418166851408899338896597194377077...
alpha_1(m_Z; P)          = 0.0168856675706697268833336999187013644391886...
alpha_2(m_Z; P)          = 0.0337778141092008025440755564310614913315613...
alpha_3(m_Z; P)          = 0.118335861957773282714035166498191462179389...
```

## Thomson Transport

The source-locked Thomson endpoint is:

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
Delta_impl(P) = 8.68656711942533994738537185035214394532834...
alpha_impl^-1 = a0(P) + Delta_impl(P)
              = 136.994835177412937295289429464436887879658...
```

The closure residual is small:

```text
alpha_fixed_point_residual = 0.0000000000000000000000000000000000000224381344364184150170238430414124287834
```

The fixed-point algebra has converged.

## Source-Audit Residual

At the public endpoint pixel, the same source anchor uses the transport term

```text
Delta_required(P) = 137.035999177 - a0(P)
                  = 8.72773111901240265209594238591525606567...
```

The transport term is short by

```text
Delta_missing(P) = 0.0411639995870627047105705355631121203415...
```

Equivalently,

```text
alpha^-1 gap = 300.388217944789749256636936289181812863091... ppm
P gap        = 0.00000388645493805208523405572430668851792...
```

This is the scalar source-audit packet recorded by the finite-code ledger.

## Audit Checks

The source-audit residual leaves the outer equation intact:

```text
P = phi + alpha * sqrt(pi)
```

It also does not indicate a failure of numerical convergence. The converged
precision-100 report carries a fixed-point residual below 10^-37 (CL-6,
closed).

The finite-code payload to populate is the low-energy transport/readout term
`Delta_Th(P)`.

## Finite Source Payload

`THOMSON_TRANSPORT_THEOREMS.md` states the theorem suite for the transport layer.
The summary is below.

A full finite-code source emission replaces the structured-running
ansatz with a source-only transport theorem. The required object is a map

```text
Delta_Th(P)
```

derived from the same Ward-projected `U(1)_Q` source family as `a0(P)`.

The source theorem package specifies:

1. The renormalization scheme and matching convention that connects the D10
   source anchor `a0(P)` to a zero-momentum electromagnetic readout.
2. The charged lepton contribution from the same source branch.
3. The quark and hadronic vacuum-polarization contribution, including threshold
   and confinement effects, without replacing it by a fitted screen factor.
4. Any electroweak matching terms that belong to the chosen scheme.
5. An interval or certified numerical bound showing that the resulting
   `alpha -> alpha` map has the selected root.

The source-side equation is:

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

reads as `137.035999177(21)`.

## Audit Command

After generating a full report, run:

```bash
python3 alpha_gap_audit.py --report runtime/full_p_alpha_report_current.json
python3 thomson_endpoint_package.py --report runtime/full_p_alpha_report_current.json
python3 screening_invariant_no_go.py
python3 thomson_endpoint_interval_certificate.py
python3 transport_theorem_manifest.py --report runtime/full_p_alpha_report_current.json
python3 measured_endpoint_calibration.py
```

The command prints the implemented transport term, the required endpoint
transport term, the missing inverse-alpha contribution, and the theorem-label
manifest. This keeps any replacement of `Delta_Th(P)` easy to check.

`measured_endpoint_calibration.py` emits the OPH plus empirical hadron closure
surface for numeric tables and plots. The source-audit label belongs in
ledgers, not in introductory prose.

`thomson_endpoint_package.py` adds the conditional endpoint packet. At the
public endpoint pixel it reports

```text
P_C = 1.630968209403959324879279847782648941...
a0(P_C) = 128.307965473286248209961108741756716187...
Delta_required(P_C) = 8.728033703713751790038891258243283813...
Delta_impl_exact(P_C) = 8.686567842708528400985442542885969768...
Delta_source_residual(P_C) = 0.041465861005223389053448715357314044...
S_required(P_C) = 0.895400132647658797805800283181670641...
c_Q(P_C) = 0.658025759927155435638230170232360050...
```

The scalar $c_Q$ is defined by $S_{\rm required}=1-x+c_Qx^2$, with
$x=N_c\alpha_3(m_Z;P_C)/\pi$. It is the compact endpoint target for a
source-only Ward-projected QCD screening and endpoint-remainder map.

## Status

```text
closed:   D10 source map P -> a0(P)
closed:   outer/inner numerical fixed-point witness for the implemented map
closed:   endpoint-package blocker isolation for issue #223
closed:   source-residual non-identifiability boundary for issue #235
closed:   WardProjectedHadronicSpectralEmission_Q source-spectral reduction theorem
stage:    populated Ward-projected spectral measure payload for Delta_Th(P)
stage:    interval-wide proof for the final full transport map after R_Q(P) is emitted
empirical: OPH plus empirical hadron closure display row for alpha(0)
```
