# Publication-Grade Theorem Package

## Theorem A: Conditional OPH electroweak hierarchy

Assume:

1. `R_P`: pixel closure branch, either public endpoint or source-audit branch with status label.
2. `Pi_U`: frozen D10 source packet.
3. `DAG_U`: no forbidden measurement path into `Phi_U`, `alpha_U`, or `v/E_star`.
4. `R_U#`: interval/Krawczyk proof record.

Then the source-side weak hierarchy is:

```math
\frac{v}{E_\star}
=
P_\star^{-1/2}
\exp\left[-\frac{2\pi}{4\alpha_U(P_\star)}\right].
```

On the public endpoint branch:

```math
v/E_\star = 2.0199803239725553\times10^{-17}.
```

On the source-audit branch:

```math
v(P_{cand})/E_\star = 2.0198114150099223\times10^{-17}.
```

## Theorem B: Higgs D10/D11 declared surface

Assume the D10/D11 input interval box `I_HT`, interval extension `F_HT#`, D11 Jacobian interval map, and `DAG_HT`.

Then:

```math
m_H=125.1995304097179\,\mathrm{GeV},\qquad
m_t^{D11}=172.3523553288312\,\mathrm{GeV}.
```

This bundle records the formula/output surface but does not include the raw interval input box.

## Corollary: OPH Higgs naturalness

The selected exact source-to-Higgs comparison map has zero RG/coarse-graining
defect:

```math
\epsilon_H =
\sup_x |H_r(\rho_{sr}n_s(x))-H_r(n_r\rho_{sr}(x))|
=0.
```

Equivalently,

```math
\epsilon_H=\max\{\epsilon^n_{sH},\epsilon^h_{sH}\}=0,\qquad
\epsilon_H\in[0,0].
```

The certificate is `issue_332_rg_naturality_certificate.json`. It allows only
the OPH source branch, the D10 source `alpha_U(P_star)` interval, the upstream
repair-tick/tick-projection/joint-stability records, and the declared D10/D11
maps. It forbids measured weak-scale, Higgs, W/Z, gravity, Planck-area, Lambda,
and tuned bare-Higgs/cutoff-counterterm inputs. On the selected exact branch,
the bare/counterterm split is a regulator-coordinate split. Observer-facing
fine tuning does not appear as a proof datum.

## Theorem C: Joint `(P,N_CRC)` product fixed point

Let `x=log N`. The product-separated source branch has joint space
`X=I_P x I_x` and source map

```math
\mathcal J(P,x)=(\Gamma(P),\widehat C(x)).
```

Assume the local and global component contraction certificates

```math
|\Gamma(P)-\Gamma(Q)|\le q_P|P-Q|,\qquad
|\widehat C(x)-\widehat C(y)|\le q_N|x-y|,
\qquad q_P,q_N<1.
```

In the product metric,

```math
d((P,x),(Q,y))=\max\{|P-Q|,|x-y|\},
```

the joint map is a contraction with constant `q=max(qP,qN)<1`. Banach's theorem
therefore gives a unique stable joint fixed point
`(P_star,log N_CRC)`. The theorem does not source `N_CRC` from the weak scale.
The report `certificates/R_PN_joint_fixed_point_certificate_report.json` records
the rounded capacity display and flags the weak-scale backsolve as a circular
diagnostic only.

For a genuinely coupled source map

```math
\mathcal J(P,x)=(\Gamma(P,x),\widehat C(P,x)),
```

with derivative envelope `a,b,c,d`, the coupled branch is promoted only when
there is an `r>0` such that

```math
\max\{a+b/r,d+rc\}<1.
```

Without that bound, the coupled map has residual branch freedom and no
theorem-grade uniqueness claim.

## Theorem D: Electroweak tick-projection bridge

The local D10 transmutation law gives

```math
\frac{v}{E_{\rm cell}}
=
\exp\left[-\frac{2\pi}{4\alpha_U(P)}\right].
```

The global repair-tick lemma gives

```math
|g_*'|=(N/\pi)^{-1/48}.
```

Define the electroweak projection exponent by

```math
\Pi_{\rm EW}(P,N)
=
\frac{-\log(v/E_{\rm cell})}{-\log |g_*'|}
=
\frac{24\pi}{\alpha_U(P)\log(N/\pi)}.
```

Since \(0<|g_*'|<1\), this exponent is unique. The target projection
`Pi_EW(P_star,N_CRC)=4P_star` is equivalent to the bridge residual

```math
B_{\rm EW}(P,N)
=
\alpha_U(P)\log(N/\pi)-\frac{6\pi}{P}
=0.
```

Equivalently, the exact bridge target is

```math
N_{\rm EW}(P)=\pi\exp\left[\frac{6\pi}{P\alpha_U(P)}\right].
```

For the public endpoint branch this gives

```math
N_{\rm EW}
=
3.5323546226929906511187512962330547600462\times10^{122}.
```

The certificate is `certificates/R_EW_tick_projection_certificate.json`. It records the exact
bridge target and the rounded-capacity diagnostic. The rounded display `3.31e122` is a
capacity-scale label and fails the exact bridge residual.

The factor `4` is the electroweak transmutation multiplicity
`\beta_EW=N_c+1=4`. The factor `48` is `2*24` from the global repair tick. Hence
`12=48/4`, and the sampling exponent is `4P_star`.

## Theorem E: Icosahedral screen-sieve theorem

On the declared triangulated `S^2` screen branch, define the curvature charge

```math
q_v=6-\deg(v).
```

Euler and triangular incidence give

```math
V-E+F=2,\qquad 3F=2E,\qquad \sum_v q_v=12.
```

Convex positive defect cost selects twelve unit fivefold defects. Edge-center
collars expose those defects as central ports. The no-marked-point finite
isotropy rule places them on the `A5/C5` orbit, with size `60/5=12`.
Therefore a scalar screen load `X=log(N/pi)` is read locally as `X/12`.
With cell entropy `P/4` and `beta_EW=4`, this gives

```math
\Gamma_{\rm EW}
=
4\cdot\frac{P}{4}\cdot\frac{1}{12}\log(N/\pi)
=
\frac{P}{12}\log(N/\pi).
```

The certificate is `certificates/R_screen_sieve_icosahedral_certificate.json`.

## Theorem F: EW-refined exact capacity fixed point

Define the log-capacity map

```math
\mathcal C_{\rm EW}(P,x)
=
(1-\lambda)x+\lambda\frac{6\pi}{P\alpha_U(P)},
\qquad 0<\lambda\le 1.
```

For fixed `P=P_star` and `lambda=1/2`, this is a contraction with constant
`1/2`. Its unique fixed point is

```math
x_\star=\frac{6\pi}{P_\star\alpha_U(P_\star)},
\qquad
N_{\rm CRC}^{\rm EW}=\pi e^{x_\star}.
```

Therefore

```math
B_{\rm EW}(P_\star,N_{\rm CRC}^{\rm EW})=0.
```

The certificate is `certificates/R_EW_global_capacity_certificate.json`. It
records the fixed-point residual, the bridge residual, the residual contraction
check, and the rounded-capacity diagnostic. The rounded display `3.31e122`
remains a capacity-scale display and fails the exact EW bridge residual.

## Theorem G: Finite readback-resolution certificate

At fixed cutoff, define the finite readback pipeline

```math
F_r(N)=
\operatorname{Cap}_{\rm read}
\left(
\operatorname{Obs}
\left(
\operatorname{nf}_{r,N}(\mathfrak U_{r,N})
\right)
\right).
```

The finite delivery resolution is the positive D6 root

```math
\rho_{\rm read}(r,N)=\sqrt{\pi/F_r(N)}.
```

On the selected finite branch, the observer-facing capacity register is a
central finite record observable and the selected self-reading observer sector
has support on exactly one positive capacity atom. Therefore the finite branch
has one effective delivery resolution. If `F_r -> F` cofinally and
`F(N_CRC)=N_CRC>0`, then

```math
\rho_{\rm read}(r,N_{\rm CRC})
=
\sqrt{\pi/F_r(N_{\rm CRC})}
\to
\sqrt{\pi/N_{\rm CRC}}
=
(N_{\rm CRC}/\pi)^{-1/2}.
```

The certificate is `certificates/R_readback_resolution_certificate.json`. It
records the unique finite normal form, the central observer sector, the single
positive capacity atom, the positive-root extractor, the refinement bound, and
the forbidden-input ledger. It uses no measured weak/Higgs/gravity/cosmology
calibration.

## Theorem H: Representation-to-spectrum round count

On the OPH realized product branch, the observer-visible compact gauge repair
support is the product adjoint

```math
\mathfrak{su}(3)\oplus\mathfrak{su}(2)\oplus\mathfrak u(1).
```

Its unoriented dimension is

```math
8+3+1=12.
```

Reversible record repair carries two orientations for each adjoint generator,
the write/verify or action/coaction halves of the same public repair channel.
Therefore

```math
m_{\rm rep}
=2\dim(\mathfrak{su}(3)\oplus\mathfrak{su}(2)\oplus\mathfrak u(1))
=2(8+3+1)
=24.
```

The cyclic repair scheduler on this oriented support has period 24. The SU(5)
adjoint has the same single-orientation integer for a different support: it
contains X/Y mixed gauge bosons, which are excluded by the OPH product branch.
The certificate is `certificates/R_m_rep_24_certificate.json`.

## Lemma H: Global repair-tick lemma

Setting: the cosmic record-capacity fixed point `N_CRC = F(N_CRC)` with observed-branch
readout `N_CRC = S_dS` and D6 normalization `N_CRC = pi * (r_CRC/ell_star)^2`. Coordinate:
`rho = length / r_CRC`, so the horizon sits at `rho = 1` and the local cell at
`rho_star = ell_star / r_CRC = (N_CRC/pi)^(-1/2)`.

Definition (global repair cycle): `G_N` is the scale-free (homogeneous) readback transport
attached to the fixed point. Fixed-point closure, readback without deficit or slack, means
`G_N(1) = rho_star` exactly.

Realization (closure derived from the readback interface plus a declared counting model): the
corpus defines the readback map as `F(N) = Cap_read(Obs(nf_{r,N}(U_{r,N})))` in the refinement
limit, the capacity reconstructed by the stable self-reading observer sector from the horizon
record (`paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex`). Model `Cap_read` by the D6
area law evaluated at the effective delivery resolution `rho_read`, the canonical
scale-covariant extension of `N = pi (r/ell)^2`. The certificate records this as a declared
modeling identification. The corpus defines `Cap_read` as reconstructed capacity and does not
state this counting rule as one of its properties:

```math
F(N) = \pi / \rho_{read}^2 .
```

The fixed-point equation `N_CRC = F(N_CRC) = pi/rho_star^2` therefore forces
`rho_read = rho_star` (positive roots): deficit `rho_read > rho_star` gives `F(N) < N`, slack
`rho_read < rho_star` gives `F(N) > N`. So the closure condition `G_N(1) = rho_star` is
equivalent to the corpus fixed-point equation, and homogeneity gives

```math
G_N(\rho) = (N_{CRC}/\pi)^{-1/2}\rho.
```

Using the representation-to-spectrum theorem above, the selected repair branch
has `m_rep=24`. With `G_N = g_N^{24}` and the positive root, the one-tick normal
form is

```math
g_N(\rho) = (N_{CRC}/\pi)^{-1/48}\rho,
\qquad
|g_*'| = (N_{CRC}/\pi)^{-1/48}.
```

For a general `m`-tick decomposition the per-tick derivative is
`(N_CRC/pi)^(-1/(2m))`; the derived `m_rep = 24` gives the stated `-1/48`
exponent.

Boundary: the closure transport follows at the interface level of `F`, conditional on the
declared counting model above. The finite readback-resolution object is supplied by
`R_readback_resolution_certificate.json`. The representation-to-spectrum round count is
supplied by `R_m_rep_24_certificate.json`. The local/global hierarchy-resonance theorem
closes on the selected branch. The promoted hierarchy row does not use the local/global
resonance as an input.

## Boundary theorem: SI gravity/clock hierarchy

The electroweak hierarchy theorem supplies `R_U`, not the full clock bridge. SI gravity requires:

```math
R_\gamma =
R_U + R_\alpha + R_e^{abs} + R_{QCD/nuc}^{133Cs} + R_{atom}^{133Cs}.
```

Only after those components are source-side and no-G DAG-certified may one use:

```math
G_{SI} =
\frac{c^5}{4\pi^2\hbar\nu_{Cs}^2}\epsilon_{Cs}^2.
```
