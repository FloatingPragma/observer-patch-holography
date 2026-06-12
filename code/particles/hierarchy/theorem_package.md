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

If the RG/coarse-graining Higgs readout defect

```math
\epsilon_H =
\sup_x |H_r(\rho_{sr}n_s(x))-H_r(n_r\rho_{sr}(x))|
```

is zero, or interval-bounded within the declared Higgs output interval, then the bare/counterterm split is a regulator-coordinate split, not an observer-facing fine-tuning datum.

This bundle defines the defect but does not prove its numerical bound.

## Lemma C: Global repair-tick lemma

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
scale-covariant extension of `N = pi (r/ell)^2`, supplied here as a declared modeling
identification rather than a corpus-stated property of `Cap_read`:

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

On the declared 24-tick resonance branch (`G_N = g_N^{24}`, positive root) the one-tick
normal form is

```math
g_N(\rho) = (N_{CRC}/\pi)^{-1/48}\rho,
\qquad
|g_*'| = (N_{CRC}/\pi)^{-1/48}.
```

For a general `m`-tick decomposition the per-tick derivative is `(N_CRC/pi)^(-1/(2m))`;
the declared `m = 24` reproduces the stated `-1/48` exponent without correction.

Boundary: the closure transport is derived at the interface level of `F`, conditional on the
declared counting model above (the corpus marks the finite readback map `F_r` as schematic and
its refinement limit as conditional on existence, so the concrete finite-machinery
verification that `nf_{r,N}` delivers a single well-defined effective resolution, and that
`Cap_read` counts by the area law at that resolution, remains open). The round count `m = 24`
is the declared branch architecture, not a derived integer. That
finite-machinery verification, a representation-to-spectrum theorem for the round count, the
electroweak tick projection (assigned to later issues), the joint `(P,N)` stability theorem,
and the RG/coarse-graining naturality statement remain outside the promoted hierarchy theorem.

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
