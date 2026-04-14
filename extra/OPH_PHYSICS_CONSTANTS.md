# OPH Physics Constants Map

This file is a working map of the main constants and fixed parameters of modern physics, and how the OPH corpus derives them.

I am treating "important constants" broadly enough to include:

- fundamental couplings and scales
- particle masses
- flavor and mixing parameters
- the main gravity and cosmology constants
- the main discrete Standard Model numbers

I am not expanding every algebraically redundant quantity. For example:

- I list `alpha`-type couplings instead of duplicating them as both `alpha_i` and `g_i`.
- I list fermion masses rather than every Yukawa separately when the latter are just `y_f = sqrt(2) m_f / v`.
- I stop at the fundamental-physics layer plus flagship hadron outputs.

Tracker note:

- issue links in the `Tracker` column point to the GitHub tracker at `FloatingPragma/observer-patch-holography`

## OPH Fundamental Constants

These are OPH's own starting constants, so they are intentionally not listed as "not yet derived" rows in the tables below.

They are calibrated OPH implementation inputs, not theorem targets.

- `P = a_cell / l_P^2`
- `N_scr = log dim H_tot`

## Status Rule

- `fully derived`: fixed by the OPH structural theorem package once the realized branch is chosen, with no extra continuous implementation input beyond that branch.
- `conditional`: derived once OPH is supplied with one or more declared implementation inputs, usually `P` and/or `N_scr`, and sometimes a specific calibration or continuation branch.
- `unit convention`: a downstream familiar-unit conversion constant used to display OPH outputs in SI-facing units. Its numerical value is fixed by the declared readout convention rather than emitted today as a standalone OPH artifact.
- `not yet derived`: the corpus does not close the physical theorem lane. Compare-only adapters, same-family witnesses, and frozen simulation scaffolds do not count as closed derivations.

## Theorem-Chain Guide

This ledger is easiest to read as a small number of theorem chains. Each row says what the chain emits, what proof step is still missing, and which open GitHub issues carry that burden.

| Surface | Theorem chain | Remaining proof front | Open issues |
| --- | --- | --- | --- |
| Lorentz / BW geometry | overlap consistency -> fixed-cutoff collar control -> geometric subnet -> scaling-limit cap pair -> ordered cut-pair rigidity -> BW/Lorentz branch | canonical geometric cap-pair realization and ordered cut-pair rigidity on the realized limit | [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106), [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22) |
| Electroweak / gravity bridge | `P ->` edge entropy -> electroweak chain and Newton matching | identify the shared edge-entropy object, derive the transmutation and RG data non-circularly, and finish the local familiar-unit bridge | [#158](https://github.com/FloatingPragma/observer-patch-holography/issues/158), [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| Higgs / top | `P -> alpha_U -> (t_U, t_tr) -> v ->` one-scalar Higgs/top seed -> `(M_H, m_t)` | make the forward Higgs/top surface fully rigorous | [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34) |
| Charged leptons | exact same-family witness -> exact operator obstruction certificate -> post-promotion trace lift -> affine anchor -> `(m_e, m_mu, m_tau)` | close the post-promotion trace lift, affine anchor, and charged mass-emission route above the operator obstruction surface | [#151](https://github.com/FloatingPragma/observer-patch-holography/issues/151), [#150](https://github.com/FloatingPragma/observer-patch-holography/issues/150), [#152](https://github.com/FloatingPragma/observer-patch-holography/issues/152) |
| Quarks / CKM | selected public quark frame class `f_P` -> target-free D12 mass bridge -> public sigma descent on `f_P` -> affine mean readout -> exact running quark sextet and explicit exact forward Yukawas | class-uniform public quark-frame descent and quotient-intrinsic sigma classification | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| Neutrinos | weighted-cycle bridge rigidity + absolute attachment -> PMNS pattern + splittings + absolute masses | expose the remaining physical Majorana-phase surface | [#154](https://github.com/FloatingPragma/observer-patch-holography/issues/154) |
| Hadrons | quark descendants + `Lambda_MSbar` + nonperturbative execution -> hadron masses | unquench the hadron branch and publish systematics | [#153](https://github.com/FloatingPragma/observer-patch-holography/issues/153), [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157) |
| Familiar units / thermodynamic constants | Lorentz output + local unit bridge + global capacity bridge -> SI-facing `c`, `G`, `Lambda`, `H_dS`, temperature displays, with `hbar` / `k_B` kept in the readout-convention lane | close the local unit bridge and propagate the temperature readout honestly across the local and de Sitter surfaces | [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159), [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22) |

## Core Unit-Setting Constants

| Constant / parameter | Identified - summary | OPH constants required | How OPH derives it | Status | Tracker |
| --- | --- | --- | --- | --- | --- |
| `c` | Universal causal speed. | structural Lorentz branch only | Once the Lorentzian branch is recovered, OPH gets one invariant causal speed for all patches. The translated SI number then matches exactly by unit convention rather than by a second physical fit. | fully derived | [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `hbar` | Quantum of action. | familiar-unit bridge only | Not a standalone OPH theorem output on the declared public surface. In SI-facing readouts it is treated as a downstream exact unit convention that translates action/frequency normalizations into familiar energy-time units once the local bridge is fixed. | unit convention | [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `k_B` | Entropy-temperature conversion factor. | familiar-unit bridge only | Not a standalone OPH theorem output on the declared public surface. In SI-facing readouts it is treated as a downstream exact unit convention that translates entropy/energy normalizations into Kelvin temperature display once the local bridge is fixed. | unit convention | [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |

### Public Status Of `hbar` And `k_B`

On the OPH public constants surface, `hbar` and `k_B` are not carried as open theorem targets. They remain downstream familiar-unit readout conventions.

That means:

- OPH may recover the structural quantum and thermodynamic relations in which these symbols appear.
- the numerical SI-facing values of `hbar` and `k_B` are not claimed as standalone OPH-emitted artifacts.
- formulas such as `l_P = sqrt(hbar G / c^3)` and `T_dS = hbar H_dS / (2*pi k_B)` are familiar-unit readout formulas on the public surface, not extra OPH derivations.

## Gravity And Cosmology

| Constant / parameter | Identified - summary | OPH constants required | How OPH derives it | Status | Tracker |
| --- | --- | --- | --- | --- | --- |
| `G` | Newton's constant. Sets the strength of gravity. | UV edge entropy plus geometric area matching | OPH matches edge entropy to geometric area and gets `G = a_cell / (4 * lbar(t))`, where `lbar(t)` is the single-cell edge entropy from the heat-kernel edge distribution. | conditional | [#158](https://github.com/FloatingPragma/observer-patch-holography/issues/158), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `l_P` | Planck length. | `G`, `c`, plus the conventional `hbar` readout factor | Downstream familiar-unit definition `l_P = sqrt(hbar G / c^3)`, with `hbar` taken in the unit-convention sense above rather than as a standalone OPH-emitted constant. | conditional | — |
| `M_P` | Planck mass / Planck energy scale. | `G`, `c`, plus the conventional `hbar` readout factor | Downstream familiar-unit algebraic combination of `G`, `c`, and `hbar`, with `hbar` treated as a readout convention rather than as an independent OPH artifact. | conditional | — |
| `Lambda` | Cosmological constant. Benchmark `1.09e-52 m^-2`. | `G + N_scr` | OPH identifies `Lambda` as a global capacity parameter, not a local vacuum-energy output: `Lambda = 3*pi / (G N_scr)`. | conditional | — |
| `H_dS` | Asymptotic de Sitter Hubble scale. Benchmark `1.80706e-18 s^-1 = 55.76 km s^-1 Mpc^-1`. | `G + N_scr` | Downstream from `Lambda` via `H = c * sqrt(Lambda / 3)`. This is the OPH de Sitter `H`, not the time-dependent late-time `H_0`. | conditional | — |
| `r_dS` | de Sitter horizon radius. Benchmark `1.65900e26 m`. | `G + N_scr` | Downstream from `Lambda` via `r_dS = sqrt(3 / Lambda)`. | conditional | — |
| `A_dS` | de Sitter horizon area. Benchmark `3.45863e53 m^2`. | `G + N_scr` | Downstream from `r_dS` via `A_dS = 4*pi*r_dS^2`. | conditional | — |
| `S_dS` | de Sitter horizon entropy. Benchmark `3.30998e122` nats, or `4.77529e122` bits. | `N_scr` | In the static-patch reading this is the screen-capacity variable itself. Numerically it also matches `A_dS / (4 l_P^2)` once units are fixed. | conditional | — |
| `T_dS` | de Sitter temperature. Benchmark `2.19678e-30 K`. | `G + N_scr` | Downstream Kelvin readout from `H_dS`; on the public surface one may write `T_dS = hbar H / (2*pi k_B)`, with `hbar` and `k_B` functioning only as conventional familiar-unit readout factors rather than standalone OPH outputs. | conditional | — |
| `rho_Lambda` | Dark-energy density. Benchmark `5.84013e-27 kg m^-3`. | `G + N_scr` | Downstream algebraic quantity from `Lambda` and `G`. | conditional | — |
| `a_0^(OPH)` | MOND-like IR acceleration benchmark. Benchmark `1.02919e-10 m s^-2`. | `G + N_scr + modular-anomaly branch` | If the deep-IR response is controlled by the de Sitter scale and the fixed anomaly prefactor, OPH gives `a_0 = (15 / (8*pi^2)) c^2 sqrt(Lambda / 3)`. | conditional | — |

## Structural Standard Model Constants

| Constant / parameter | Identified - summary | OPH constants required | How OPH derives it | Status | Tracker |
| --- | --- | --- | --- | --- | --- |
| `G_phys = [SU(3) x SU(2) x U(1)] / Z_6` | Realized low-energy gauge group. | structural branch only | Recovered from compact gauge reconstruction plus MAR. | fully derived | — |
| `N_g = 3` | Number of generations. | structural branch only | CP-capability gives the lower bound, weak-sector UV-completability gives the upper window, and minimality picks `3`. | fully derived | — |
| `N_c = 3` | Number of colors. | structural branch only | Witten anomaly forces `N_c` odd once `N_g = 3`, and MAR picks the minimal genuinely nonabelian option `N_c = 3`. | fully derived | — |
| `Y_Q = 1/6` | Quark-doublet hypercharge. | structural branch only | Fixed exactly by anomaly freedom plus Yukawa invariance. | fully derived | — |
| `Y_L = -1/2` | Lepton-doublet hypercharge. | structural branch only | Fixed exactly by anomaly freedom plus Yukawa invariance. | fully derived | — |
| `Y_u = -2/3` | Up-type singlet hypercharge. | structural branch only | Fixed exactly by anomaly freedom plus Yukawa invariance. | fully derived | — |
| `Y_d = 1/3` | Down-type singlet hypercharge. | structural branch only | Fixed exactly by anomaly freedom plus Yukawa invariance. | fully derived | — |
| `Y_e = 1` | Charged-lepton singlet hypercharge. | structural branch only | Fixed exactly by anomaly freedom plus Yukawa invariance. | fully derived | — |
| `Y_H = 1/2` | Higgs hypercharge. | structural branch only | Fixed exactly by anomaly freedom plus Yukawa invariance. | fully derived | — |
| `m_gamma = 0` | Photon mass. | structural branch only | Unbroken `U(1)_em` overlap redundancy forbids a hard photon mass term. | fully derived | — |
| `m_gluon = 0` | Gluon mass. | structural branch only | Unbroken color gauge redundancy forbids a hard gluon mass term. | fully derived | — |
| `m_graviton = 0` | Graviton mass. | structural branch only on the dynamical-metric branch | Diffeomorphism redundancy forbids a hard graviton mass term. | fully derived | — |

## Gauge Couplings, Electroweak Scale, Higgs, And Top

| Constant / parameter | Identified - summary | OPH constants required | How OPH derives it | Status | Tracker |
| --- | --- | --- | --- | --- | --- |
| `alpha_U` | Unified gauge coupling. Reference value `0.04112498`. | `P` | The forward transmutation certificate solves `alpha_U` from the pixel-closure equation. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32) |
| `t_U` | Unified diffusion parameter. Reference value `1.62354915`. | `P` | Once `alpha_U` is fixed, `t_U = 4*pi^2*alpha_U`. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32) |
| `t_tr` | Electroweak transmutation exponent. Reference value `38.19567355`. | `P` | Once `alpha_U` is fixed, `t_tr = 2*pi / ((N_c + 1) * alpha_U)`. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32) |
| `alpha_1(m_Z)` | GUT-normalized hypercharge coupling at the electroweak reference scale. Reference value `0.01688600`. | `P` | Emitted on the target-free electroweak running tree. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `alpha_2(m_Z)` | Weak coupling at the electroweak reference scale. Reference value `0.03377844`. | `P` | Emitted on the target-free electroweak running tree and source basis. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `alpha_s(m_Z)` | Strong coupling at the electroweak reference scale. Reference value `0.11833723`. | `P` | Emitted on the target-free electroweak running tree and then used by the `Lambda_MSbar` descendant. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `alpha_em(m_Z)` | Electroweak-scale fine-structure constant. Reference value `alpha_em^-1 = 128.30576920`. | `P` | Emitted on the source-locked electroweak running family. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `alpha_em(0)` | Low-energy fine-structure constant. Reference value `alpha_em^-1(0) = 137.035999177`. | `P +` Ward-projected `U(1)_Q` transport theorem | Read as the Thomson endpoint of the Ward-projected electromagnetic transport family, anchored at `alpha_em^-1(m_Z^2) = 128.30576920`. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `sin^2 theta_W(m_Z)` | Weak mixing angle at the electroweak reference scale. Reference value `0.2307354235`. | `P` | Emitted on the source-locked electroweak running family. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `v` | Higgs vacuum expectation value. Reference value `246.76711733 GeV`. | `P` | The transmutation branch gives `v = E_cell(P) * exp(-t_tr)`. | conditional | — |
| `G_F` | Fermi constant. Reference value from the OPH `v` is `1.16121e-5 GeV^-2`. | `P` | Once `v` is fixed, `G_F = 1 / (sqrt(2) * v^2)` follows algebraically. | conditional | — |
| `M_W` | W-boson mass. Current public row `80.377000015 GeV`. | `P` | The active public theorem is the target-free source-only electroweak repair law. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `M_Z` | Z-boson mass. Current public row `91.187978078 GeV`. | `P` | Same target-free electroweak theorem as for `M_W`. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `M_H` | Higgs mass. Current public row `125.218922060 GeV`. | `P +` one-scalar Higgs/top branch | The forward Higgs/top seed closes the live predictive branch above the electroweak closure. | conditional | [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `m_t` | Top-quark pole mass. Current public row `172.388645595 GeV`. | `P +` one-scalar Higgs/top branch | Same branch as `M_H`. | conditional | [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) |
| `lambda_H` | Higgs quartic coupling. | `P +` one-scalar Higgs/top branch | Once `M_H` and `v` are fixed, this is a downstream algebraic descendant. | conditional | [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34) |
| `y_t` | Top Yukawa coupling. | `P +` one-scalar Higgs/top branch | Once `m_t` and `v` are fixed, this is a downstream algebraic descendant. | conditional | [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34) |

## Flavor And Neutrinos

| Constant / parameter | Identified - summary | OPH constants required | How OPH derives it | Status | Tracker |
| --- | --- | --- | --- | --- | --- |
| `m_e` | Electron mass. | `P +` charged-lepton absolute branch + missing final closure theorems | Exact same-family witnesses and an exact operator obstruction certificate exist, but the theorem lane lacks the post-promotion trace lift, the affine anchor, and the charged mass-emission step. | not yet derived | [#151](https://github.com/FloatingPragma/observer-patch-holography/issues/151), [#150](https://github.com/FloatingPragma/observer-patch-holography/issues/150), [#152](https://github.com/FloatingPragma/observer-patch-holography/issues/152) |
| `m_mu` | Muon mass. | `P +` charged-lepton absolute branch + missing final closure theorems | Same theorem boundary as `m_e`. | not yet derived | [#151](https://github.com/FloatingPragma/observer-patch-holography/issues/151), [#150](https://github.com/FloatingPragma/observer-patch-holography/issues/150), [#152](https://github.com/FloatingPragma/observer-patch-holography/issues/152) |
| `m_tau` | Tau mass. | `P +` charged-lepton absolute branch + missing final closure theorems | Same theorem boundary as `m_e`. | not yet derived | [#151](https://github.com/FloatingPragma/observer-patch-holography/issues/151), [#150](https://github.com/FloatingPragma/observer-patch-holography/issues/150), [#152](https://github.com/FloatingPragma/observer-patch-holography/issues/152) |
| `m_u` | Up-quark mass. | `P` on the selected public quark frame class `f_P` | The target-free D12 mass bridge, selected-class public sigma descent, and affine mean readout emit the exact PDG 2025 API running-quark value on `f_P`. Stronger class-uniform quark-frame classification is a separate theorem surface. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `m_d` | Down-quark mass. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `m_u`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `m_s` | Strange-quark mass. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `m_u`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `m_c` | Charm-quark mass. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `m_u`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `m_b` | Bottom-quark mass. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `m_u`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `theta_12^CKM` | First CKM mixing angle. | `P` on the selected public quark frame class `f_P` | The selected-class public sigma descent and explicit exact forward Yukawa package emit the physical CKM matrix on `f_P`; CKM angles are algebraic descendants on that scope. Stronger class-uniform quark-frame classification is a separate theorem surface. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `theta_13^CKM` | Second CKM mixing angle. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `theta_12^CKM`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `theta_23^CKM` | Third CKM mixing angle. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `theta_12^CKM`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `delta_CKM` | CKM CP phase. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `theta_12^CKM`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `J_CKM` | Jarlskog invariant. | `P` on the selected public quark frame class `f_P` | Same theorem surface as `theta_12^CKM`. | conditional | [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199) |
| `theta_12^PMNS` | Solar PMNS angle. Theorem-grade value `34.2259 deg`. | weighted-cycle neutrino theorem branch | The weighted-cycle branch closes the physical PMNS pattern on the active theorem lane. | conditional | — |
| `theta_13^PMNS` | Reactor PMNS angle. Theorem-grade value `8.68636 deg`. | weighted-cycle neutrino theorem branch | Same theorem lane as `theta_12^PMNS`. | conditional | — |
| `theta_23^PMNS` | Atmospheric PMNS angle. Theorem-grade value `49.7228 deg`. | weighted-cycle neutrino theorem branch | Same theorem lane as `theta_12^PMNS`. | conditional | — |
| `delta_PMNS` | PMNS Dirac CP phase. Theorem-grade value `305.581 deg`. | weighted-cycle neutrino theorem branch | Same theorem lane as `theta_12^PMNS`. | conditional | — |
| `J_PMNS` | PMNS Jarlskog invariant. Theorem-grade value `-0.0275312`. | weighted-cycle neutrino theorem branch | Same theorem lane as `theta_12^PMNS`. | conditional | — |
| `Delta m^2_21` | Solar neutrino splitting. Theorem-grade value `7.48806e-5 eV^2`. | `P + weighted-cycle bridge rigidity + absolute attachment` | The scale-free weighted-cycle branch closes first, then the absolute-attachment theorem emits the physical splitting. | conditional | — |
| `Delta m^2_31` | Larger ordered neutrino splitting. Theorem-grade value `2.51231e-3 eV^2`. | `P + weighted-cycle bridge rigidity + absolute attachment` | Same theorem lane as `Delta m^2_21`. | conditional | — |
| `Delta m^2_32` | Atmospheric neutrino splitting. Theorem-grade value `2.43743e-3 eV^2`. | `P + weighted-cycle bridge rigidity + absolute attachment` | Same theorem lane as `Delta m^2_21`. | conditional | — |
| `m_1` | Lightest neutrino mass. Theorem-grade value `0.0174547 eV`. | `P + weighted-cycle bridge rigidity + absolute attachment` | The absolute-attachment theorem emits `lambda_nu` and therefore the absolute family. | conditional | — |
| `m_2` | Second neutrino mass. Theorem-grade value `0.0194820 eV`. | `P + weighted-cycle bridge rigidity + absolute attachment` | Same theorem lane as `m_1`. | conditional | — |
| `m_3` | Heaviest neutrino mass. Theorem-grade value `0.0530752 eV`. | `P + weighted-cycle bridge rigidity + absolute attachment` | Same theorem lane as `m_1`. | conditional | — |
| `alpha_21^(Maj)` | First physical Majorana phase, if the neutrinos are Majorana. | open neutrino continuation | The corpus contains Majorana-holonomy machinery, but it does not yet expose final physical Majorana phases as public constants. | not yet derived | [#154](https://github.com/FloatingPragma/observer-patch-holography/issues/154) |
| `alpha_31^(Maj)` | Second physical Majorana phase, if the neutrinos are Majorana. | open neutrino continuation | Same theorem boundary as `alpha_21^(Maj)`. | not yet derived | [#154](https://github.com/FloatingPragma/observer-patch-holography/issues/154) |

## QCD And Hadrons

| Constant / parameter | Identified - summary | OPH constants required | How OPH derives it | Status | Tracker |
| --- | --- | --- | --- | --- | --- |
| `Lambda_MSbar` | QCD infrared scale in the `MSbar` scheme. Quoted `n_f = 3` value `0.334401707 GeV`. | `P + alpha_s(m_Z) + threshold matching` | Artifact `oph_qcd_lambda_msbar3` is a secondary descendant of the target-free electroweak strong-coupling branch plus the reference-free quark thresholds. | conditional | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) |
| `theta_QCD` | Strong-CP angle. | deferred continuation only | The SM/GR derivation paper explicitly treats strong-CP proposals as a later continuation branch, not part of the declared theorem package. | not yet derived | [#155](https://github.com/FloatingPragma/observer-patch-holography/issues/155) |
| `m_pi` | Pion mass. | `P + quark descendants + Lambda_MSbar + nonperturbative hadron closure` | The hadron lane is still execution-dependent and frozen behind theorem and production/systematics blockers. | not yet derived | [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157), [#153](https://github.com/FloatingPragma/observer-patch-holography/issues/153) |
| `m_K` | Kaon mass. | `P + quark descendants + Lambda_MSbar + nonperturbative hadron closure` | Same theorem boundary as `m_pi`. | not yet derived | [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157), [#153](https://github.com/FloatingPragma/observer-patch-holography/issues/153) |
| `m_p` | Proton mass. | `P + quark descendants + Lambda_MSbar + nonperturbative hadron closure` | Same theorem boundary as `m_pi`. | not yet derived | [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157), [#153](https://github.com/FloatingPragma/observer-patch-holography/issues/153) |
| `m_n` | Neutron mass. | `P + quark descendants + Lambda_MSbar + nonperturbative hadron closure` | Same theorem boundary as `m_pi`. | not yet derived | [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157), [#153](https://github.com/FloatingPragma/observer-patch-holography/issues/153) |
| `m_rho` | Rho-meson mass. | `P + quark descendants + Lambda_MSbar + nonperturbative hadron closure` | Same theorem boundary as `m_pi`. | not yet derived | [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157), [#153](https://github.com/FloatingPragma/observer-patch-holography/issues/153) |

## Bridge Program: `W/Z -> Edge Entropy -> G`

This is the unification statement worth pushing.

The clean claim is not:

- "`G` is numerically derived from `P`"

because `P = a_cell / l_P^2`, and `l_P^2` contains `G`.

The clean claim is:

- the same OPH edge-entropy law that calibrates the electroweak `W/Z` sector should also set Newton matching on the gravity side

### Electroweak Side

On the target-free electroweak surface, the forward transmutation certificate makes the electroweak chain explicit:

```text
P
-> alpha_U
-> (t_U, t_tr)
-> (alpha_2, alpha_Y, v)
-> (W, Z)
```

The explicit pixel-closure law on that surface is

```text
ellbar_SU2(t2_mz_run) + ellbar_SU3(t3_mz_run) = P / 4
```

So the `W/Z` calibration is controlled by heat-kernel edge entropy data, not by an unrelated fitted sector.

### Current Gravity Side

On the gravity side, the SM/GR derivation paper gives the Newton matching formula

```text
G = a_cell / (4 * ellbar(t))
```

where `ellbar(t)` is the single-cell edge entropy density from the heat-kernel edge distribution.

### Intended Unification Statement

What OPH wants to show is:

1. the `ellbar` controlling the electroweak pixel-closure law is the same OPH edge-entropy object that appears in Newton matching
2. the electroweak calibration therefore fixes the same edge microphysics that determines `G`
3. `W/Z` and `G` are then tied together by one common edge-entropy law rather than by two unrelated sectors

That would be the real Standard Model / GR bridge.

### What Is Still Missing

This is not yet a closed theorem-grade derivation. The remaining burden is to show, internally and non-circularly, that:

- the gravity-side `ellbar(t)` and the electroweak heat-kernel entropy entering `P / 4` are literally the same OPH object
- the fixed-cutoff edge heat-kernel / Casimir theorem is lifted and identified strongly enough to carry both branches
- the electroweak closure is promoted from a calibration surface to a theorem-grade predictive bridge

### Active Blockers

- [#158](https://github.com/FloatingPragma/observer-patch-holography/issues/158) identify the shared edge entropy object across electroweak closure and Newton matching
- [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31) fix the transmutation parameter from `P` without circular calibration
- [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32) derive RG matching and threshold structure from OPH
- [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) upgrade gauge/electroweak closure from consistency check to prediction
- [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) close the local unit bridge for familiar in-universe units
- [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22) unify local Einstein recovery with the global cosmological-capacity branch
- [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106) resolve the BW branch-selection and continuum algebra-type gap

Supporting upstream note:

- [#30](https://github.com/FloatingPragma/observer-patch-holography/issues/30) is still open, and its pipeline status is marked `completed` / `under audit` on the screen-microphysics paper. The remaining gap is not just "invent the heat-kernel law"; it is the lift into the SM/GR paper and the cross-branch identification.

### Best Present-Day Wording

The strongest accurate wording today is:

- "OPH has the skeleton of a `W/Z -> edge entropy -> G` bridge."
- "Electroweak calibration and Newton matching are controlled by the same edge-microphysics law, pending closure of the shared heat-kernel entropy object."

### Smallest Exact-Release Package

The local exact-release frontier consists of five explicit pieces:

1. a local familiar-unit readout package
2. a shared edge-entropy bridge
3. a strict classical-branch clause
4. a target-free electroweak identity surface
5. a closed one-scalar Higgs/top promotion surface

On that declared extension surface, the local release package is:

- `c = 299792458 m/s`
- `G = 6.674299995910528e-11 m^3 kg^-1 s^-2`
- `M_W = 80.377 GeV`
- `M_Z = 91.18797809193725 GeV`
- `M_H = 125.1995304097179 GeV`

with two important caveats:

- `c` is structural from the Lorentz branch rather than a particle-style `P` fit
- the `G` row is an exact emitted branch value on the local extension surface, not a literal zero-difference identity against the rounded display benchmark `6.6743e-11`

## Theorem Roadmap For Familiar In-Universe Units

The clean way to organize the familiar-unit problem is:

- dimensionless OPH outputs first
- then the local UV unit bridge from `a_cell`
- then the global de Sitter / capacity bridge from `N_scr`

Checked against the open `FloatingPragma/observer-patch-holography` tracker on `2026-04-04`, the practical roadmap is:

| Target family | OPH output | What still has to be proved or calculated | Varies with | Tracker | Conditional target surface |
| --- | --- | --- | --- | --- | --- |
| `c` | Structural Lorentz output: one invariant causal speed. | Keep the BW / Lorentz branch honest; for SI display, close the local unit bridge `c = ell_cell / tau_cell`. | invariant itself: neither; SI display: `a_cell` only | [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) | `c = 299792458 m/s`, with `ell_cell = 2.06409e-35 m` and `tau_cell = 6.88507e-44 s` |
| Massless bosons | Photon, gluon, and graviton mass zeros are structural OPH outputs. | No extra numerical solve is needed; only keep the dynamical-metric / Lorentz branch bookkeeping honest. | neither | [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106) | `m_gamma = 0`, `m_gluon = 0`, `m_graviton = 0` |
| `W`, `Z` | The target-free electroweak chain emits the public electroweak pair on the realized branch. | The remaining exact-live object is one target-free electroweak identity theorem on top of the forward `P -> t -> couplings` trunk and the explicit GeV bridge. | `a_cell` only for familiar units | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) | public live rows `M_W = 80.377000015 GeV`, `M_Z = 91.187978078 GeV`; exact codomain `80.377 GeV`, `91.18797809193725 GeV` |
| Higgs / top | The closed one-scalar Higgs/top chain emits a live Higgs/top branch above the electroweak closure. | The remaining exact-live object is one live-forward Higgs/top promotion theorem, with the local energy bridge kept explicit. | `a_cell` only for familiar units | [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) | public live rows `M_H = 125.218922060 GeV`, `m_t = 172.388645595 GeV`; exact codomain `125.1995304097179 GeV`, `172.3523553288312 GeV` |
| `G` | Gravity-side dictionary states `G = a_cell / (4 * ellbar(t))`. | The remaining exact-live objects are the shared electroweak/gravity edge-entropy bridge, the local readout package, and the strict classical-regime clause; the sharp primitive missing proof is the branch-preserving transport step behind the shared entropy object. | `a_cell` only | [#158](https://github.com/FloatingPragma/observer-patch-holography/issues/158), [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22), [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106), [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) | exact emitted branch value on the local extension surface `6.674299995910528e-11 m^3 kg^-1 s^-2` |
| Cell / Planck scales | Benchmark translation gives `a_cell`, `ell_cell`, `tau_cell`, and `E_cell`. | State one honest theorem-facing local unit bridge, with `hbar` and `k_B` kept explicitly in the downstream familiar-unit convention lane. | `a_cell` only | [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159) | `a_cell = 4.26047e-70 m^2`, `ell_cell = 2.06409e-35 m`, `tau_cell = 6.88507e-44 s`, `E_cell = 9.55999e18 GeV` |
| Fine-structure / weak couplings | The target-free electroweak chain emits dimensionless electroweak couplings on the realized branch. | Promote calibration to prediction and derive the matching structure from OPH rather than external RG packaging. | neither | [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33) | `alpha_U = 0.04112498`, `alpha_em(m_Z)^-1 = 128.30576920`, `sin^2 theta_W(m_Z) = 0.2307354235` |
| `Lambda`, `H_dS`, de Sitter observables | Once `G` and `N_scr` are supplied, OPH emits the de Sitter branch algebraically. | Unify the local gravity and global capacity stories into one theorem stack, then propagate the local unit bridge cleanly. | both `a_cell` and `N_scr`, except `S_dS` which is `N_scr` only | [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22), [#47](https://github.com/FloatingPragma/observer-patch-holography/issues/47), [#52](https://github.com/FloatingPragma/observer-patch-holography/issues/52) | `Lambda = 1.09e-52 m^-2`, `H_dS = 1.80706e-18 s^-1 = 55.76 km s^-1 Mpc^-1`, `r_dS = 1.65900e26 m` |
| Temperature-side quantities | Thermodynamic structure exists, and Kelvin displays use conventional `hbar` / `k_B` readout factors rather than standalone OPH constants. | Close the local familiar-unit bridge and propagate the public temperature readout honestly across the de Sitter branch. | local temperature scales: `a_cell`; de Sitter temperature: both `a_cell` and `N_scr` | [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159), [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22) | benchmark display: `T_cell = 1.10939e32 K`, `T_dS = 2.19678e-30 K` |

### Execution Order

If the goal is to maximize how many familiar constants OPH can quote numerically without overclaiming, the practical order is:

1. keep the Lorentz / BW branch honest enough that `c` and the local geometric language sit on a declared branch: [#106](https://github.com/FloatingPragma/observer-patch-holography/issues/106)
2. finish the shared edge-entropy bridge needed for non-circular `G`: [#158](https://github.com/FloatingPragma/observer-patch-holography/issues/158), [#22](https://github.com/FloatingPragma/observer-patch-holography/issues/22)
3. close the local familiar-unit bridge while keeping `hbar` / `k_B` in the downstream readout-convention lane: [#159](https://github.com/FloatingPragma/observer-patch-holography/issues/159)
4. prove the exact live electroweak chart identity so public `W/Z` rows land on the exact declared chart: [#31](https://github.com/FloatingPragma/observer-patch-holography/issues/31), [#32](https://github.com/FloatingPragma/observer-patch-holography/issues/32), [#33](https://github.com/FloatingPragma/observer-patch-holography/issues/33)
5. prove the exact live Higgs/top promotion on the forward branch: [#34](https://github.com/FloatingPragma/observer-patch-holography/issues/34)
6. once the local `c/G/W/Z/H` package is honest, propagate the same bridge into `Lambda`, `H_dS`, and the other de Sitter descendants: [#47](https://github.com/FloatingPragma/observer-patch-holography/issues/47), [#52](https://github.com/FloatingPragma/observer-patch-holography/issues/52)
7. only after that expand the same familiar-unit surface deeper into charged flavor, stronger ambient quark-frame classification, and hadrons: [#152](https://github.com/FloatingPragma/observer-patch-holography/issues/152), [#199](https://github.com/FloatingPragma/observer-patch-holography/issues/199), [#157](https://github.com/FloatingPragma/observer-patch-holography/issues/157)

### Conditional Familiar-Unit Targets If The Bridges Close

If the remaining theorem gaps close on the realized branch, the familiar-unit target surface to compare against is:

- local cell bridge: `a_cell = 4.2604724363543046e-70 m^2`, `ell_cell = 2.0640911889629063e-35 m`, `tau_cell = 6.885067098528897e-44 s`, `E_cell = 9.559993349864843e18 GeV`
- causal speed: `c = 299792458 m/s`
- electroweak bosons, public live rows: `M_W = 80.37700001539531 GeV`, `M_Z = 91.18797807794321 GeV`
- electroweak bosons, exact local extension surface: `M_W = 80.377 GeV`, `M_Z = 91.18797809193725 GeV`
- Higgs / top, public live rows: `M_H = 125.218922060 GeV`, `m_t = 172.388645595 GeV`
- Higgs / top, exact local extension surface: `M_H = 125.1995304097179 GeV`, `m_t = 172.3523553288312 GeV`
- dimensionless electroweak surface: `alpha_U = 0.04112498041477454`, `alpha_em(m_Z)^-1 = 128.30576920`, `sin^2 theta_W(m_Z) = 0.2307354235`
- Fermi constant: `G_F = 1.16120908925587e-5 GeV^-2`
- Newton bridge target: if `ellbar = P / 4`, then `G = 6.6743e-11 m^3 kg^-1 s^-2`
- cosmology benchmark surface: `Lambda = 1.09e-52 m^-2`, `H_dS = 1.80706e-18 s^-1 = 55.76 km s^-1 Mpc^-1`, `r_dS = 1.65900e26 m`, `rho_Lambda = 5.84013e-27 kg m^-3`

## Comparison To Official Reference Values

For the unit-setting constants, the comparison has to be separated into four categories:

- exact-by-definition SI matches
- conditional bridge matches
- downstream unit conventions rather than standalone OPH outputs
- not apples-to-apples quantities

### Exact-By-Definition SI Match

- `c`: the translated OPH familiar-unit value is `299792458 m/s`, which matches the official SI / NIST value exactly. This is not an independent particle-style postdiction; it is exact because the SI meter is defined from the speed of light once the invariant causal speed has been chosen as the unit bridge.

### Conditional Bridge Match

- `G`: on the sharpened local extension surface the emitted branch value is `6.674299995910528e-11 m^3 kg^-1 s^-2`, which matches the CODATA / NIST value `6.67430(15)e-11 m^3 kg^-1 s^-2` at the quoted digits. This is numerically strong, but it is not a theorem-grade independent prediction on the corpus surface because the shared `ellbar` identification, the strict classical-regime clause, and the local familiar-unit bridge remain open.

### Downstream Unit Conventions, Not Standalone OPH Outputs

- `hbar`: the SI-facing value is exact by SI unit definition, but on the OPH public surface it is used only as a downstream familiar-unit conversion factor, not as a standalone emitted OPH constant.
- `k_B`: likewise, the SI-facing value is exact by SI unit definition, but on the OPH public surface it is used only as a downstream entropy-to-Kelvin conversion factor, not as a standalone emitted OPH constant.

### Not Apples-To-Apples With A Present-Day Measurement

- `H_dS`: OPH’s `55.76 km s^-1 Mpc^-1` benchmark is the asymptotic de Sitter Hubble scale, not the observed late-time `H_0`, so it should not be compared directly to late-time observational `H_0` fits.

## Numerical Cosmology Outputs On This Surface

Using the benchmark `Lambda = 1.09e-52 m^-2`, OPH gives:

- `H_dS = c * sqrt(Lambda / 3) = 1.80706e-18 s^-1 = 55.76 km s^-1 Mpc^-1`
- `r_dS = sqrt(3 / Lambda) = 1.65900e26 m`
- `A_dS = 4*pi*r_dS^2 = 3.45863e53 m^2`
- `S_dS = A_dS / (4 l_P^2) = 3.30998e122` nats `= 4.77529e122` bits
- `T_dS = hbar H_dS / (2*pi k_B) = 2.19678e-30 K`, with `hbar` and `k_B` understood here as downstream familiar-unit readout conventions rather than standalone OPH outputs
- `rho_Lambda = Lambda c^2 / (8*pi G) = 5.84013e-27 kg m^-3`
- `a_0^(OPH) = (15 / (8*pi^2)) c^2 sqrt(Lambda / 3) = 1.02919e-10 m s^-2`

## Fast Reading Summary

Today, the OPH constant story looks like this:

- fully derived: the discrete Standard Model structure (`G_phys`, `N_g`, `N_c`, exact hypercharges) and the exact structural zeros (`m_gamma = 0`, `m_gluon = 0`, `m_graviton = 0`)
- conditional on `P`: the target-free electroweak and one-scalar Higgs/top stack, the QCD `Lambda_MSbar` descendant, and the theorem-grade neutrino branch
- conditional on `N_scr`: `Lambda` and its de Sitter descendants
- public familiar-unit conventions rather than standalone OPH outputs: `hbar` and `k_B`
- not yet derived: the physical charged-lepton theorem lane, stronger ambient quark-frame classification, physical Majorana phases, strong CP, and the hadron spectrum

## Primary Source Surface

- [SM/GR derivation paper](reverse-engineering-reality/paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex)
- [master paper fragment](reverse-engineering-reality/paper/tex_fragments/PAPER.tex)
- [particle paper](reverse-engineering-reality/paper/deriving_the_particle_zoo_from_observer_consistency.tex)
- [particle status table](reverse-engineering-reality/code/particles/RESULTS_STATUS.md)
- [exact non-hadron bundle](reverse-engineering-reality/code/particles/EXACT_NONHADRON_MASSES.md)
- [particle claim ledger](reverse-engineering-reality/code/particles/ledger.yaml)
