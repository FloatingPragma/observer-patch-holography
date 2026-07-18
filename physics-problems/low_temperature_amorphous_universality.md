# Universality of Low-Temperature Amorphous Solids

## Motivating Result

This note entered the queue because the low-temperature glass tables contained
a sharp anomaly: many chemically different amorphous solids showed
nearly the same acoustic and thermal attenuation band, while some amorphous
silicon and germanium films landed more than two orders of magnitude lower
([Pohl, Liu, and Thompson, Rev. Mod. Phys. 74, 991, 2002](https://link.aps.org/doi/10.1103/RevModPhys.74.991)).
Later elastic measurements in dense amorphous silicon even reported the
suppression of tunneling states
([Phys. Rev. Lett. 113, 025503, 2014](https://link.aps.org/doi/10.1103/PhysRevLett.113.025503)).
The OPH question is why a broad ordinary branch has a nearly universal
readout, while special films can escape it.

**Status:** conditional OPH phenomenological branch and simulator specification.
The standard TLS reduction is recovered, while the twelve-port/twenty-four-slot
acoustic bridge is an explicit model assumption awaiting independent material
receipts. The numerical output is therefore a within-model benchmark, not yet a
source-only material theorem.

Date: 2026-07-08

## Introduction

Low-temperature amorphous solids pose a microscopic-selection problem because many
different disordered materials show nearly the same acoustic and thermal
phonon attenuation. Composition, bond topology, density, sound speed, elastic
coupling, quench history, and defect history vary widely. The measured
dimensionless phonon mean-free-path ratio often lands in a narrow band. The
empirical problem is the coexistence of material-specific microscopic history
with near-universal low-energy readout.

In OPH the object is a self-reading rigid material patch. The patch has local
state, strain-visible ports, two-basin records, repair or tunneling moves,
readback by long-wavelength elastic carriers, and public acoustic or thermal
evidence receipts. Material details live in the source branch. The universal
claim belongs to the quotient that removes hidden labels and projects onto the
ordinary rigid-record normal form.

On the illustrative single-active-pair branch, the declared Bernoulli event
probability per wavelength-sized cell is

```math
p_{\rm cell}
=
\frac{e^{-P_{\rm src}/24}}{12\cdot24}
=
3.2440989\times10^{-3}.
```

If, in addition, this cell event is identified with exponential acoustic
scattering, the exact attenuation exponent and acoustic tunneling-strength
readout are

```math
\frac{\lambda}{\ell}
=
-\log(1-p_{\rm cell})
=
3.2493724\times10^{-3},
```

```math
C_\ast^{\rm OPH}
=
\frac{-\log(1-p_{\rm cell})}{\pi^2}
=
3.2923026\times10^{-4}.
```

Those values use the OPH source-side pixel, twelve exposed acoustic ports, a
twenty-four-state oriented repair register, and protected-reserve survival.
Those ingredients do not by themselves imply an acoustic cross section: the
event subset, product probability law, and material bridge are assumptions
stated below. Measured attenuation, internal friction, thermal conductivity,
and tunneling strength serve only as comparison data.

As a descriptive comparison, against the ordinary amorphous-solid table of
Pohl, Liu, and Thompson, the conditional value sits about $0.090\sigma$ from the
measured log-center. The same value sits $2.01\sigma$ above the separate
low-attenuation amorphous Si, Ge, and C film pool. The first number is the
location of the benchmark within the pooled historical spread; it is not a
prediction-significance or goodness-of-fit statistic. The second number is a
descriptive comparison with the known exception films, not a pre-registered
class-separation test.

## Standard-Theory Baseline and Unresolved Target

The standard tunneling-state model is empirically successful. It represents
low-energy defects in glasses as two-level systems coupled to strain and then
summarizes the acoustic response through

```math
C_i
=
\frac{\bar P\gamma_i^2}{\rho v_i^2}.
```

Here $\bar P$ is a density of tunneling states, $\gamma_i$ is an elastic
coupling, $\rho$ is mass density, and $v_i$ is sound speed. The difficulty is
that the combination is almost universal even though its factors are
material-dependent. Legacy physics can describe the readout once $C_i$ is
given. The source of the dimensionless band is underdetermined in that
parameterization.

Three ledgers need separation:

- the material source law that creates defects;
- the normal form of low-energy two-basin records;
- the quotient-visible acoustic probability seen by a long-wavelength phonon.

The standard tunneling model is not the whole conventional baseline.
Soft-potential, interacting-TLS, elastic-block, and collective-mode approaches
also seek a microscopic explanation of the narrow band. This note does not
prove those approaches impossible. Its narrower target is to make the source
law, TLS normal form, public attenuation readout, and promotion evidence
separate and auditable.

## Actual OPH Contribution

OPH separates the material source branch from the observer-visible readout.
The source branch keeps composition, bond topology, density, elastic tensor,
disorder, quench history, anneal history, and defect history. The quotient
removes hidden carrier labels, atom labels, basis choices, mesh labels, port
names, and repair schedules that leave the public strain response unchanged.

The low-temperature normal form is a rigid-record ensemble. Each active local
sector is a two-basin record with a repair amplitude and a strain dipole. A
long-wavelength phonon is insensitive to the hidden basin label. It records
whether one wavelength intersects an active quotient-visible rigid-record
event.

OPH's specific contribution is the self-reading-patch presentation: a bounded
material region has local state, strain-visible ports, durable TLS-like records,
feedback or repair moves, phonon readback, and public receipts. The
port/register construction is a provisional specialization of that
presentation, not a consequence of ordinary elasticity. On the illustrative
branch the visibility space and probability law are separately declared as

```math
\mathcal V_{\rm LTAS}
=
\mathcal P_{12}^{\rm ac}
\times
\mathcal R_{24}^{\rm or}.
```

```math
\mu_{\rm LTAS}
=
\mu_{12}^{\rm unif}\otimes\mu_{24}^{\rm unif}.
```

The Cartesian product specifies the classical outcome set. It does not license
multiplication of probabilities; the product measure is an additional
independence assumption. A declared event subset is also required before
$1/(12\cdot24)$ can be interpreted as an event probability.

## Abstract

This note gives a conditional OPH normal form for a saturated rigid-record
material quotient. Given a twelve-port acoustic presentation, a twenty-four-slot
repair register, a product probability law, a single active event pair, and the
protected-reserve survival assumption, it emits
$p_{\rm cell}=3.2440989\times10^{-3}$ and
$\lambda/\ell=3.2493724\times10^{-3}$. These are within-model consequences,
not a derivation that every glass realizes the branch. Low-attenuation films are
tests of independently declared membership receipts, not exceptions assigned
after their attenuation is inspected.

## Claim Boundary

The numerical value applies only to the declared uniform saturated LTAS model.
Its inputs are the candidate $P_{\rm src}$, a twelve-port acoustic
specialization of the screen sieve, a twenty-four-slot repair register, a
product law, a one-pair event subset, Poisson protected-reserve survival, and
the event-to-scattering bridge. Measured $\lambda/\ell$, internal friction,
thermal conductivity, and tunneling strength are excluded from the numerical
evaluation, but the choice of acoustic specialization requires an
independent no-target-leak audit.

The epistemic ledger is:

- **Assumed for this branch:** acoustic realization of the twelve ports,
  uniform product measure, one active port/register pair, independent reserve
  survival, scale-flat TLS action over declared cutoffs, and the material
  cross-section bridge.
- **Derived within the model:** the TLS two-state reduction, the event
  probability $p_{\rm cell}$, the exact Bernoulli attenuation exponent, and
  standard TLS readout conversions.
- **Empirical receipts:** pre-attenuation branch membership, TLS density and
  strain coupling, frequency/temperature plateau, polarization, processing
  history, attenuation, internal friction, and thermal conductivity.

Equal attenuation across every amorphous solid is outside the claim. The branch
claim is the quotient-visible leading scattering probability for

```math
\text{ordinary disordered rigid materials}
\quad
\text{with}
\quad
\text{saturated scale-flat strain-coupled rigid records}
```

on the saturated branch.

The branch excludes materials independently shown to lack scale-flat rigid records, with sparse
repair ensembles, with biased port or register activation, with non-uniform
protected-reserve behavior, or with suppressed elastic coupling to strain.

## Empirical Target

Pohl, Liu, and Thompson review the low-temperature glass data and report that
the phonon wavelength to mean-free-path ratio is typically in the band

```math
10^{-3}
\lesssim
\frac{\lambda}{\ell}
\lesssim
10^{-2},
```

The tunneling-strength band is

```math
10^{-4}
\lesssim
C
\lesssim
10^{-3}.
```

They also describe exceptions, including some amorphous Si, Ge, and C films
with much lower attenuation.

The acoustic readout convention is

```math
\frac{\lambda}{\ell}
=
\pi^2 C.
```

For the dominant transverse thermal-phonon convention,

```math
\frac{\lambda_{\rm dom}}{\ell}
\simeq
12.5 C_t.
```

The OPH target is this empirical band with no material-by-material fit of
$C$, $\bar P$, $\gamma$, $\rho$, or $v$.

## Source Branch

The source branch is

```math
\mathsf{OPH}_{\rm LTAS}
=
\mathsf{ObserverPatch}
+
\mathsf{ScreenMicrophysics}
+
\mathsf{ProtectedReserve}
+
\mathsf{RigidRecordTLS}
+
\mathsf{FiniteQuotientEnsemble}.
```

The OPH source-audit pixel candidate used here is

```math
P_{\rm src}
=
1.630972095858897\ldots.
```

The public comparison endpoint

```math
P_C
=
1.630968209403959\ldots
```

changes the result by about $1.6\times10^{-7}$ relatively. It is negligible at
four significant figures but visible in the long decimal displays, so the two
values must not be described as numerically identical. The local source ledger
classifies $P_{\rm src}$ as a source-audit witness rather than a
completed physical-endpoint proof. The measured Thomson endpoint remains a
comparison coordinate.

The bounded observer-like patch is a rigid material region with strain-visible
interfaces, local two-basin records, repair amplitudes, readback by acoustic or
thermal phonons, and public evidence receipts. Hidden atom labels and
microscopic carrier coordinates are silent unless they change the visible
record.

## Definition 1: LTAS Material Quotient

For a material $x$ and regulator $r$, define the LTAS material branch

```math
\mathcal B_{\rm LTAS}(x,r)
=
\left(
x,
G_x,
\{\mathcal A_i^x\},
H_x,
\Gamma_x,
\mathcal C_{\rm rr},
\mathcal C_{\rm el},
\mathcal P_{\rm proc}
\right).
```

Here $x$ contains composition, bond topology, density, elastic tensor,
coordination constraints, disorder, quench history, anneal history, and defect
history. $G_x$ is the finite structural graph. $\mathcal A_i^x$ are local
accessible algebras. $H_x$ is the source Hamiltonian or action. $\Gamma_x$
is the equivalence relation that removes hidden labels. $\mathcal C_{\rm rr}$
is the rigid-record class, $\mathcal C_{\rm el}$ is the elastic readout class,
and $\mathcal P_{\rm proc}$ records processing history.

The physical quotient is

```math
Q_{x,r}
=
\Sigma_{x,r}/\Gamma_{x,r}.
```

The source law is

```math
\mu_{x,r,T}(q)
=
Z_{x,r,T}^{-1}
m_{x,r}(q)e^{-S_{x,r,T}(q)}.
```

Material law selection belongs to the source branch. Normal-form projection
states the surviving ordinary low-energy rigid-record ensemble.

## Definition 2: Rigid-Record Pair

A rigid-record pair is a local two-basin quotient sector

```math
R_j
=
\{q_{j,+},q_{j,-}\}
\subset
Q_{x,r}.
```

Its asymmetry is

```math
\Delta_j
=
E(q_{j,+})-E(q_{j,-}).
```

Its repair or tunneling amplitude is

```math
\Delta_j^0.
```

For strain component $s$, its elastic dipole is

```math
\gamma_{j,s}
=
\frac12
\left(
\frac{\partial E(q_{j,+})}{\partial\varepsilon_s}
-
\frac{\partial E(q_{j,-})}{\partial\varepsilon_s}
\right).
```

The level splitting is

```math
E_j
=
\sqrt{\Delta_j^2+(\Delta_j^0)^2}.
```

The pair is a record because the two basins are locally stable, readably
different under strain, and connected by a permitted repair move.

## Definition 3: Rigid-Record Spectral Measure

Choose finite physical cutoffs

```math
|\Delta|\le \Delta_{\max},
\qquad
u_{\min}\le u\le u_{\max},
\qquad
u=-\log(\Delta^0/E_0).
```

The rigid-record spectral measure on this regulated domain is

```math
\nu_{x,r}
=
\frac1{V_r}
\sum_{j\in{\rm RR}(x,r)}
\delta_{(\Delta_j,u_j,\gamma_j;\rho_r,v_r)},
```

where

```math
u_j
=
-\log(\Delta_j^0/E_0).
```

The energy-resolved density of states is

```math
g_r(E)
=
\frac{d}{dE}
\mathbb E_{\mu_{x,r,0}}
\left[
\frac1{V_r}\#\{j:E_j\le E\}
\right].
```

It need not equal the joint TLS distribution coefficient. Denote that
coefficient by $\bar P_{0,r}$ through

```math
dN
=
\bar P_{0,r}\,dV\,d\Delta\,du
```

within the stated low-energy window. The cutoffs are essential: exact
translation invariance on an unbounded $u$ axis is not normalizable.

The elastic repair strength is

```math
C_{s,r}
=
\frac{\bar P_{0,r}\langle\gamma_s^2\rangle_r}
{\rho_r v_{s,r}^2}.
```

This is the standard material-facing TLS quantity. The conditional OPH branch
assigns its benchmark by first specifying a quotient-visible event probability;
the physical equality requires the material bridge below.

## Theorem 1: Rigid-Record TLS Reduction

**Statement.** In a localized-basin basis, after choosing the relative phase so
the leading tunneling amplitude is real and neglecting off-diagonal strain
coupling, an isolated two-basin rigid record reduces to

```math
H_j
=
\frac12
\left(
\Delta_j\sigma_j^z+\Delta_j^0\sigma_j^x
\right)
+
\sum_s
\gamma_{j,s}\varepsilon_s(x_j)\sigma_j^z.
```

**Proof.** A two-basin quotient sector has a two-dimensional low-energy
subspace. Any Hermitian operator on it is a scalar plus Pauli terms. A
strain-independent scalar is irrelevant to the level dynamics; a
strain-dependent common-mode term must instead be retained in an elastic-stress
calculation. In the declared localized basis the diagonal energy difference gives
$\Delta\sigma^z/2$. The off-diagonal repair move gives
$\Delta^0\sigma^x/2$. Long-wavelength strain perturbs the two basins through
opposite elastic dipoles, so it enters as $\gamma\varepsilon\sigma^z$.
Exchanging the two basin labels changes the corresponding parameter signs but
not the spectrum or readout. The displayed form is basis-fixed, not invariant
term by term under an arbitrary two-state rotation. $\square$

## Proposition 2: Scale-Flat Repair Action

**Statement.** If the regulated rigid-record density is smooth in $\Delta$
near zero and approximately translation-invariant in

```math
u
=
-\log(\Delta^0/E_0),
```

then

```math
dN
=
\bar P_0\,dV\,d\Delta\,|du|
=
\bar P_0\,dV\,d\Delta\,\frac{|d\Delta^0|}{\Delta^0},
```

or

```math
P(\Delta,\Delta^0)
=
\frac{\bar P_0}{\Delta^0}
```

over the declared finite window.

**Proof.** Smoothness near $\Delta=0$ gives a leading locally flat $d\Delta$
density. Approximate translation invariance in barrier action gives a leading
locally flat $|du|$ density on the finite window. Since
$|du|=|d\Delta^0|/\Delta^0$, the regulated tunneling distribution
follows. $\square$

## Assumption 3 and Proposition 3: Port-Slot Event Model

**Assumptions.** The LTAS visibility space factors as

```math
\mathcal V_{\rm LTAS}
=
\mathcal P_{12}^{\rm ac}
\times
\mathcal R_{24}^{\rm or}.
```

The probability law, which does not follow from the Cartesian-product sample space,
is the uniform product law

```math
\Pr(a,o)
=
\frac1{12}\frac1{24}.
```

A subset

```math
\mathcal A_{\rm scat}
\subseteq
\mathcal P_{12}^{\rm ac}\times\mathcal R_{24}^{\rm or}
```

must declare which pairs count as scattering events. The numerical benchmark
uses the **single-active-pair ansatz** $|\mathcal A_{\rm scat}|=1$. Reserve
survival is assumed independent of the selected pair and, on the uniform
Poisson product-thickening branch, has probability

```math
\Pr(\text{reserve survives})
=
e^{-P_{\rm src}/24}.
```

**Proposition.** Under those assumptions,

```math
p_{\rm cell}
=
\Pr[(a,o)\in\mathcal A_{\rm scat}]\,
\Pr(\text{reserve survives})
=
\frac{|\mathcal A_{\rm scat}|}{12\cdot24}
e^{-P_{\rm src}/24}.
```

For $|\mathcal A_{\rm scat}|=1$,

```math
p_{\rm cell}
=
3.244098917358505\times10^{-3}.
```

This is a conditional finite-event calculation. The port and slot counts alone
do not select $\mathcal A_{\rm scat}$ and do not establish that a real phonon
samples one independent cell per wavelength.

## Protected-Reserve Survival

On the uniform product-thickening branch with local Poisson reserve occupancy,
the protected-reserve factor is assumed to be

```math
\lambda_{\rm collar}
=
e^{-P/24}.
```

Using the source-audit pixel candidate,

```math
\lambda_{\rm collar}^{\rm src}
=
e^{-1.630972095858897/24}
=
0.9343004881928429.
```

The single-active-pair cell probability is therefore

```math
p_{\rm cell}
=
\frac{1}{12}
\frac{1}{24}
e^{-P_{\rm src}/24}
=
3.244098917358505\times10^{-3}.
```

This is a within-model cell probability, not yet a universal inverse relative
mean free path.

## Proposition 4: Event-to-Mean-Free-Path Bridge

**Assumption.** Let $n_{\rm TLS}$ be an active-defect number density and
$\sigma_{\rm TLS}(\omega,T)$ the total scattering cross section derived from
the TLS Hamiltonian. The material bridge is

```math
p_{\rm cell}(\omega,T)
=
1-\exp[-n_{\rm TLS}\sigma_{\rm TLS}(\omega,T)\lambda].
```

This equation is dimensionally coherent and makes explicit the object missing
from a purely combinatorial count. It must be derived or empirically certified
for a material before promotion.

**Proposition.** If successive wavelength-sized cells are independent and a
scattering event removes the carrier from the unscattered beam, then after $n$
cells the survival probability is

```math
(1-p_{\rm cell})^n
=
\exp[n\log(1-p_{\rm cell})].
```

Consequently the exponential attenuation convention gives

```math
\frac{\lambda}{\ell}
=
-\log(1-p_{\rm cell})
=
3.249372414504198\times10^{-3}.
```

Only at leading order for $p_{\rm cell}\ll1$ is
$\lambda/\ell\simeq p_{\rm cell}$. $\square$

## Proposition 5: Elastic TLS Readout

**Statement.** For the TLS ensemble of Theorem 1 and Proposition 2, in the weak-drive
relaxation plateau and the stated acoustic attenuation convention,

```math
Q_0^{-1}
=
\frac{\pi}{2}C,
```

and

```math
\frac{\lambda}{\ell}
=
\pi^2C
```

in the acoustic plateau convention. In the dominant transverse thermal-phonon
convention,

```math
\frac{\lambda_{\rm dom}}{\ell}
\simeq
12.5C_t.
```

The thermal coefficient additionally assumes transverse-dominant heat
transport, the dominant-phonon approximation, and the empirical
longitudinal/transverse velocity and coupling relations used in the cited
review. Under the event-to-attenuation bridge, the conditional OPH benchmark is

```math
C_\ast^{\rm OPH}
=
\frac{-\log(1-p_{\rm cell})}{\pi^2}.
```

**Proof.** The standard TLS attenuation calculation integrates the
strain-coupled TLS matrix element over the scale-flat distribution
$P(\Delta,\Delta_0)=\bar P_0/\Delta_0$. The leading material combination is

```math
C
=
\frac{\bar P_0\langle\gamma^2\rangle}{\rho v^2}.
```

The acoustic plateau relation maps $\lambda/\ell$ to $\pi^2C$. The dominant
thermal relation maps the same source branch to $12.5C_t$. Combining these
readout identities with Proposition 4 gives the conditional value of $C_\ast$.
This step does not derive $\bar P_0\langle\gamma^2\rangle/(\rho v^2)$ from the
port count; equality with the event model is the material bridge that requires
a receipt. $\square$

## Numerical Comparison

The conditional single-active-pair benchmark values are:

```math
\frac{\lambda}{\ell}
=
3.249372414504198\times10^{-3},
```

```math
C_\ast^{\rm OPH}
=
3.292302591323264\times10^{-4},
```

```math
Q_0^{-1}
=
\frac{\pi}{2}C_\ast
=
5.171536817147902\times10^{-4},
```

```math
\frac{\lambda_{\rm dom}}{\ell}
=
12.5C_\ast
=
4.115378239154080\times10^{-3}.
```

| quantity | OPH conditional benchmark | empirical target |
| --- | ---: | ---: |
| acoustic $\lambda/\ell$ | $3.2494\times10^{-3}$ | $10^{-3}$ to $10^{-2}$ |
| acoustic $C_\ast$ | $3.2923\times10^{-4}$ | $10^{-4}$ to $10^{-3}$ |
| acoustic plateau $Q_0^{-1}$ | $5.1715\times10^{-4}$ | typical glass plateau scale |
| dominant thermal $\lambda_{\rm dom}/\ell$ | $4.1154\times10^{-3}$ | $10^{-3}$ to $10^{-2}$ |

The benchmark lands inside the historical glass band without a
material-specific numerical fit of $\bar P_0$, $\gamma$, $\rho$, or $v$.
That coincidence is a comparison result, not proof of the acoustic bridge.

## Sigma Comparison With Measured Values

The sigma comparison uses the measured tunneling-strength entries in
Pohl, Liu, and Thompson, Tables I-IV. These tables list ordinary amorphous
solids and exclude the separate amorphous Si, Ge, and C exception-film class
of Table V. The entries used here are the numeric values in the $C_l$ and
$C_t$ columns, in units of $10^{-4}$.

Because the empirical spread is multiplicative, the comparison is done in
log space:

```math
\mu_{\log C}
=
\frac1N\sum_{j=1}^N\log_{10}C_j,
\qquad
\sigma_{\log C}
=
\sqrt{
\frac1{N-1}
\sum_{j=1}^N
\left(\log_{10}C_j-\mu_{\log C}\right)^2
}.
```

For the ordinary-amorphous pool,

```math
N=83,
\qquad
\mu_{\log C}=-3.512277129,
\qquad
\sigma_{\log C}=0.330270631\ {\rm dex}.
```

The corresponding geometric mean is

```math
C_{\rm meas,geo}
=
3.0741345\times10^{-4}.
```

The OPH value has

```math
\log_{10}C_\ast^{\rm OPH}
=
-3.482500256,
```

so its measured-table sigma displacement is

```math
z_C
=
\frac{\log_{10}C_\ast^{\rm OPH}-\mu_{\log C}}
{\sigma_{\log C}}
=
0.0902.
```

The acoustic $\lambda/\ell$ comparison uses the same sigma, since
$\lambda/\ell=\pi^2C$ shifts all log values by the same constant. The measured
ordinary-pool geometric mean maps to

```math
\left(\frac{\lambda}{\ell}\right)_{\rm meas,geo}
=
\pi^2C_{\rm meas,geo}
=
3.0340492\times10^{-3}.
```

The OPH value

```math
\left(\frac{\lambda}{\ell}\right)_{\rm OPH}
=
3.2493724\times10^{-3}
```

is again

```math
z_{\lambda/\ell}
=
0.0902.
```

Ordinary amorphous-solid pool, from Pohl/Liu/Thompson Tables I-IV:

- source columns: $C_l$ and $C_t$ numeric entries;
- entries: $N=83$;
- measured geometric mean: $C_{\rm meas,geo}=3.0741\times10^{-4}$;
- measured log spread: $\sigma_{\log C}=0.3303$ dex;
- OPH displacement: $0.0902\sigma$;
- acoustic readout mean: $\lambda/\ell=3.0340\times10^{-3}$;
- acoustic OPH displacement: $0.0902\sigma$.

The Table V exception films have a separate transverse-strength pool. For the
listed amorphous Si, Ge, and C films,

```math
C_{t,{\rm exc,geo}}
=
8.8169\times10^{-6},
\qquad
\sigma_{\log C,{\rm exc}}
=
0.7811\ {\rm dex}.
```

Relative to that exception pool, the conditional ordinary-branch benchmark
sits about $2.01\sigma$ above the exception-pool center.

The two standardized locations answer different descriptive questions. The
ordinary-pool value is

```math
z_{\rm ordinary}=0.0902.
```

The exception-pool statistic asks how far the conditional benchmark lies from
the historical low-attenuation film pool. The answer,

```math
z_{\rm exception}=2.01,
```

does not by itself establish class separation: the entries are pooled
measurements rather than independent material draws, include multiple
polarizations and processing states, and use the sample spread rather than a
prediction uncertainty. A reproducible comparison must publish the 83-row
extraction ledger, measurement qualifications, and a material-level
hierarchical or held-out analysis.

## Exception Boundary

The universal branch covers this class:

```math
\text{ordinary disordered rigid materials whose low-energy quotient has}
```

```math
\text{a saturated, scale-flat, strain-coupled rigid-record ensemble}
```

```math
\text{with reserve-unbiased port/register activation.}
```

The branch fails if any of the following fail:

```math
\text{scale-flat }P(\Delta,u),
```

```math
\mathcal V_{\rm LTAS}
=
\mathcal P_{12}^{\rm ac}
\times
\mathcal R_{24}^{\rm or},
```

```math
\lambda_{\rm collar}
=
e^{-P/24},
```

```math
S_{\rm rr}^{\rm sat}
\simeq
1,
```

```math
\gamma_s
\ne
0.
```

Here $S_{\rm rr}^{\rm sat}$ is the saturation score for the active
low-energy rigid-record ensemble. A material with sparse active records,
biased repair slots, strongly suppressed strain coupling, or non-scale-flat
barrier statistics lies outside the universal band.

Known low-attenuation samples such as certain amorphous Si, Ge, and C films may
belong outside the saturated/uniform branch only if receipts fixed independently
of the attenuation measurement fail the stated conditions. Otherwise they are
counterexamples to the branch prediction.

## Simulator Receipts

A future simulator should freeze the conditional cell probability

```math
p_{\rm cell}
=
\frac{e^{-P_{\rm src}/24}}{12\cdot24}
```

Experimental data enter after this benchmark and the branch-classification
rule are frozen. The following objects and receipts specify the required audit;
they supply no material receipts by themselves.

Required simulator objects:

- `LTASMaterialQuotient`
- `RigidRecordPair`
- `RigidRecordSpectralMeasure`
- `PortRegisterProjection`
- `ProtectedReserveSampler`
- `ElasticTLSReadout`
- `ExceptionClassifier`

Required receipts:

- `LTAS_SOURCE_BRANCH_RECEIPT`
- `QUOTIENT_INVARIANCE_RECEIPT`
- `SOURCE_ACTION_RECEIPT`
- `RIGID_RECORD_TLS_RECEIPT`
- `SCALE_FLAT_REPAIR_ACTION_RECEIPT`
- `PORT_12_UNIFORMITY_RECEIPT`
- `REGISTER_24_UNIFORMITY_RECEIPT`
- `PORT_REGISTER_TENSOR_RECEIPT`
- `PROTECTED_RESERVE_UNIFORMITY_RECEIPT`
- `MEAN_FREE_PATH_READOUT_RECEIPT`
- `EXCEPTION_CLASS_RECEIPT`

The simulator comparison must audit branch membership from independent inputs
and compare the measured attenuation against the frozen conditional benchmark.

## Falsifiers

The OPH LTAS branch fails if:

- a material passes all pre-registered branch receipts before attenuation is
  inspected but its attenuation disagrees with the benchmark outside the
  declared combined uncertainty;
- a material is classified as saturated even though the active rigid-record
  density is sparse or history-biased;
- port and register activation are empirically biased and the theorem uses the
  uniform $12\cdot24$ denominator anyway;
- branch membership is changed after a low-attenuation result is seen;
- the numerical value is re-fit after experimental data are loaded.

## Scope

OPH supplies a self-reading-patch organization and a conditional finite-event
benchmark; it does not yet close the microscopic universality problem. A
material-specific claim requires a declared source law, an identified
rigid-record ensemble, regulated scale-flat action, a derived acoustic
event/cross-section bridge, pre-registered port/register and reserve receipts,
and an uncertainty-aware comparison. Until then the $12\times24$ branch is a
useful visualization and falsifiable ansatz rather than a unique resolution.

## References

- R. O. Pohl, X. Liu, and E. Thompson, "Low-temperature thermal conductivity
  and acoustic attenuation in amorphous solids", Reviews of Modern Physics 74,
  991, 2002. https://link.aps.org/doi/10.1103/RevModPhys.74.991
- P. W. Anderson, B. I. Halperin, and C. M. Varma, "Anomalous low-temperature
  thermal properties of glasses and spin glasses", Philosophical Magazine 25,
  1, 1972. https://doi.org/10.1080/14786437208229210
- W. A. Phillips, "Two-level states in glasses", Reports on Progress in
  Physics 50, 1657, 1987. https://doi.org/10.1088/0034-4885/50/12/003
- X. Liu, D. R. Queen, T. H. Metcalf, J. E. Karel, and F. Hellman,
  "Hydrogen-Free Amorphous Silicon with No Tunneling States", Physical Review
  Letters 113, 025503, 2014. https://doi.org/10.1103/PhysRevLett.113.025503
