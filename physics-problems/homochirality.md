# Homochirality as Prebiotic Record-Branch Selection

## Motivating Result

This question entered the OPH problem queue as issue
[#492](https://github.com/FloatingPragma/observer-patch-holography/issues/492):
can biological handedness be stated precisely as branch selection in a
prebiotic replication-and-repair loop, while keeping autocatalytic, mineral,
circular-polarization, and parity-violation routes distinct?

The motivating laboratory result is Soai and collaborators' demonstration
that a small initial enantiomeric excess can be amplified by asymmetric
autocatalysis [Soai et al. (1995)](https://doi.org/10.1038/378767a0).
Independent experiments establish other pieces of a possible route:
face-selective amino-acid adsorption on calcite
[Hazen, Filley, and Goodfriend (2001)](https://doi.org/10.1073/pnas.101085998),
complete deracemization of conglomerate crystals under grinding and recycling
[Viedma (2005)](https://doi.org/10.1103/PhysRevLett.94.065504), and percent-level
enantiomeric excess after circularly polarized ultraviolet irradiation of
initially achiral interstellar-ice analogues
[de Marcellus et al. (2011)](https://doi.org/10.1088/2041-8205/727/2/L27).
These results establish viable seed and gain mechanisms. They do not yet form
one demonstrated path from geochemically plausible feedstock to the linked
amino-acid and sugar handedness of biology.

**Status:** conditional OPH audit and reduced-model derivation. The Frank-type
normal form, source/gain separation, route audit, and falsification conditions
are explicit under their stated assumptions. The note does not derive a unique
prebiotic pathway, the terrestrial sign, spatial planetary fixation, or a new
chiral force. Those remain empirical source, transport, and chemistry
boundaries.

**Date:** 2026-07-11.

## Abstract

OPH does not add a new chiral force and does not derive that terrestrial
proteins must use L-amino acids. Its useful contribution is a finite,
auditable decomposition of the homochirality problem. A prebiotic compartment
is a bounded observer-like patch when it has feed and transport ports, local
chemical state, a readback variable such as product enantiomeric excess,
persistent chiral templates or catalysts as records, and record-conditioned
copying, inhibition, recycling, or repair. In a symmetric Frank-type kinetic
normal form, the racemic state cannot choose a sign deterministically, but any
nonzero seed is amplified toward one of two stable chiral branches. In the
constant-coefficient, source-off error model below, non-racemic fixed branches
exist only when the effective amplification rate $\kappa$ exceeds
$2\mu$. This is not a universal prebiotic threshold. Mineral surfaces,
circularly polarized light, and electroweak parity violation are therefore
candidate source terms, not substitutes for amplification and repair.
Autocatalysis, mutual inhibition, and recycling are candidate gain terms, not
proofs of the original sign. The remaining scientific gap is an end-to-end,
prebiotically plausible, spatially resolved source-to-fixation receipt.

## Standard Baseline, Remaining Problem, and OPH Contribution

**Established chemistry.** Frank-type nonlinear kinetics, asymmetric
autocatalysis, mutual inhibition, crystallization/recycling, chiral surfaces,
circularly polarized photochemistry, stochastic fixation, and reaction--
diffusion models supply standard mechanisms for enrichment and
symmetry breaking. The cited experiments demonstrate particular seed or gain
mechanisms under particular conditions.

**Unresolved part.** No cited experiment supplies one geochemically plausible,
spatially resolved path from near-racemic feedstock through coupled amino-acid
and sugar selection to durable biological homochirality. The historical sign,
transport between compartments, survival under racemization and dilution, and
replicator architecture remain open.

**OPH contribution in this note.** OPH does not uniquely supply the kinetic
equations. It contributes the observer-like patch and audit decomposition:
boundary source, record-conditioned gain, readback, repair/error, transport,
checkpoint persistence, and public evidence are kept distinct. A specifically
OPH physical result would require the self-reading record architecture to make
a frozen prediction beyond the standard kinetic controls.

## 1. Claim boundary

Three claims must not be conflated:

1. **Molecular chirality:** a molecule has non-superposable mirror forms.
2. **Enantiomeric enrichment:** a process produces $L\ne R$ for a specified
   compound in a specified compartment.
3. **Biological homochirality:** a coupled replication and metabolism network
   persistently uses one amino-acid hand and the complementary sugar hand
   across copying, transport, and ecological succession.

The first is ordinary stereochemistry. The experiments cited here establish
instances of the second. The third requires a coupled, durable selection
process. OPH places a standard reduced normal form inside an observer-patch and
receipt audit for that selection process. It does not
identify the historically realized molecular network from observer
consistency alone.

The article therefore makes a conditional statement:

> If a bounded prebiotic network supplies a nonzero chiral seed, nonlinear
> record-conditioned gain, loss of the competing hand, and persistence across
> environmental checkpoints, then homochirality is fixation of one of two
> record branches. The selected sign and the rates are properties of the
> declared chemical source branch.

## 2. Source OPH branch

The source branch is the finite observer-patch and repair formalism in
[Screen Microphysics and Observer Synchronization](../paper/screen_microphysics_and_observer_synchronization.tex)
and the fixed-point/replicator formalism in
[Reality as Consensus Protocol](../paper/reality_as_consensus_protocol.tex).
The relevant OPH object is not a generic reaction vessel. It is a bounded,
self-reading chemical patch with the following operational quotient:

| OPH component | Prebiotic realization | Public readout |
|---|---|---|
| local state | concentrations of enantiomers, activated monomers, templates, catalysts, inhibitors, and waste | time-resolved concentration and speciation |
| boundary or ports | feedstock, mineral face, photon flux, wet-dry exchange, flow, or migration | calibrated flux, geometry, and residence time |
| readback | enantiomeric excess and copy yield produced by the current record | chiral chromatography or equivalent stereospecific assay |
| record | a persistent chiral template, catalyst, crystal population, or catalytic surface geometry | survival and re-read after transfer or dilution |
| repair move | autocatalytic copying, mutual inhibition, dissolution/recrystallization, degradation/resynthesis, or error correction | measured rate law and mass balance |
| checkpoint | continuation after cycling, dilution, compartment division, or transfer | preregistered pass/fail persistence test |
| evidence bundle | raw chromatograms, blanks, calibration, rates, face orientation, light helicity and dose, spatial samples, and analysis code | independently auditable receipt |

This is observer-like in the restricted operational sense: the patch's durable
chiral record changes its later boundary-visible products, those products read
back the record, and repair preferentially restores the same branch. No claim
of consciousness is involved.

The separation between **normal form** and **source selection** is essential.
OPH can classify the two fixed branches and the conditional model threshold. A
chemical rate law, mineral population, radiation field, or parity-violating
energy difference must supply the sign, magnitude, and physical realization.

## 3. Minimal branch-selection derivation

Let $L(t),R(t)\ge 0$ be the amounts of two enantiomeric record carriers in
one well-mixed compartment. The symmetric Frank normal form combines equal
net autocatalytic growth $g=a-d$ with mutual inhibition at rate $b>0$:

```math
\dot L=gL-bLR,\qquad
\dot R=gR-bLR.
\tag{1}
```

If $L$ and $R$ are amounts, $g$ has units of inverse time and $b$ has units of
inverse amount per unit time, so $\kappa=bx/2$ below has units of inverse time.

This is a deliberately minimal quotient of the chemistry, descended from
[Frank's asymmetric-synthesis model](https://doi.org/10.1016/0006-3002(53)90082-1).
It is not asserted to be a complete prebiotic reaction mechanism.

Define total record abundance and signed enantiomeric excess

```math
x=L+R,\qquad y=L-R,\qquad e=\frac{y}{x}\in[-1,1]
\quad (x>0).
\tag{2}
```

Since $LR=x^2(1-e^2)/4$, direct substitution gives

```math
\dot e=\frac{bx}{2}\,e(1-e^2)
       =\kappa(t)e(1-e^2),
\qquad \kappa(t):=\frac{b x(t)}{2}.
\tag{3}
```

### Proposition 1: symmetry does not choose a hand

For exactly racemic initial data, $e(0)=0$, equation (3) gives $e(t)=0$
for all $t$. The symmetric deterministic kinetics therefore cannot derive an
absolute sign. A fluctuation, asymmetric boundary condition, or explicit
source term is logically necessary.

The conclusion is stronger than the informal statement that one hand
"eventually wins": without a seed, the ideal symmetric model never leaves the
racemic branch. In a finite stochastic population, molecular shot noise can
provide a seed, but it selects either sign probabilistically unless an external
bias changes the fixation probabilities.

### Proposition 2: a nonzero seed has logarithmic amplification cost

Equation (3) does not generally have constant $\kappa$, because
$\kappa(t)=b x(t)/2$ changes as total abundance grows or depletes. Define the
cumulative dimensionless amplification exposure

```math
K(t):=\int_0^t\kappa(s)\,ds.
```

For $e_0\ne0$, equation (3) integrates to

```math
e(t)^2=
\left[1+\bigl(e_0^{-2}-1\bigr)e^{-2K(t)}\right]^{-1}.
\tag{4}
```

The sign of $e$ is conserved. The limit $|e|\to1$ follows only if
$K(t)\to\infty$; resource depletion can leave a finite final excess. The
exposure needed to reach a target $0<e_\star<1$ is

```math
K_\star=
\frac12\log\!\left[
\frac{e_\star^2}{1-e_\star^2}
\frac{1-e_0^2}{e_0^2}
\right].
\tag{5}
```

For constant $\kappa$, $K(t)=\kappa t$ is the special case used for a clock
estimate. Thus $e_0=10^{-2}\to e_\star=0.99$ requires $K_\star=6.554$,
while $e_0=10^{-6}\to0.99$ requires $K_\star=15.764$. A very small seed is not
automatically irrelevant: in this normal form, the required additional gain
grows only logarithmically. Whether a real system supplies that gain before
hydrolysis, racemization, depletion, or dispersal is an experimental question.

### Proposition 3: repair must beat chirality-erasing error

A useful biased, error-prone quotient is

```math
\dot e=(1-e^2)(\kappa e+h)-2\mu e.
\tag{6}
```

Here $h$ is a signed source bias and $\mu$ is an effective racemization,
wrong-handed copying, or record-erasure rate. Equation (6) is a phenomenological
normal form; each proposed chemistry must derive its own mapping to
$(\kappa,h,\mu)$. All three have units of inverse time in this normalization.
The factor $-2\mu e$ follows exactly from symmetric first-order
interconversion $L\rightleftarrows R$ with one-way rate $\mu$; using it
for copy errors, selective degradation, or template loss is an additional
effective-model assumption.

For constant $\kappa>0$ and $\mu\ge0$ with the source switched off
$(h=0)$, the non-racemic fixed points are

```math
e_\star=\pm\sqrt{1-\frac{2\mu}{\kappa}},
\tag{7}
```

and they exist only if

```math
\boxed{\kappa>2\mu.}
\tag{8}
```

Within this autonomous normal form, the linearization at $e=0$ is
$\dot e\simeq(\kappa-2\mu)e$. Below the threshold the racemic fixed point is
stable; above it the racemic fixed point is unstable and two stable chiral
fixed points exist. At equality, $\dot e=-\kappa e^3$, so the racemic
point is nonhyperbolic and attracting in this one-dimensional model.
This is a source-off pitchfork condition for the stated
model, not a universal repair threshold. If $h\ne0$, $e=0$ is not a fixed
point and the symmetric pitchfork is tilted. If $\kappa$, $h$, or $\mu$
varies with time, an instantaneous inequality $\kappa(t)>2\mu(t)$ does not by
itself prove fixation or persistence; the full exposure and checkpoint
dynamics must be integrated. The qualitative OPH audit statement is only that
measured gain must outpace measured chirality-erasing processes on the claimed
timescale.

## 4. Local fixation is not yet global homochirality

Let $e_i$ describe compartment $i$ in a graph of ponds, pores, films, or
droplets. A minimal spatial stochastic extension is

```math
de_i=\left[(1-e_i^2)(\kappa_i e_i+h_i)-2\mu_i e_i
+D\sum_j A_{ij}(e_j-e_i)\right]dt
+\sqrt{\frac{2\sigma_i^2}{N_i}(1-e_i^2)}\,dW_i.
\tag{9}
```

The adjacency matrix $A$ declares material exchange, $D$ its rate,
$N_i$ the effective molecular population, and $W_i$ independent noise
terms. Here $D$ and $\sigma_i^2$ have units of inverse time when $A$ is
dimensionless. The simple exchange drift $D\sum_jA_{ij}(e_j-e_i)$ assumes
symmetric transport and equal or externally fixed total record abundance in
connected compartments. If $x_i=L_i+R_i$ differs appreciably, transport must be
written for $L_i$ and $R_i$ and then reduced; abundance-weighted drift terms
result.
The factor $1-e_i^2$ is an illustrative bounded, Wright--Fisher-like
choice that makes the diffusion vanish at $e_i=\pm1$; the earlier additive
noise form would leave the physical interval with nonzero probability. Equation
(9) assumes symmetric effective abundance noise for the two hands,
independent compartments, and sufficiently large $N_i$. If $L$ and $R$ have
unequal copy/death noise, correlated feed fluctuations, or strongly varying
total abundance, the diffusion and Itô drift must be derived from the underlying
count process rather than assumed. Reflecting boundaries are an alternative but
must be stated.

Under those assumptions, equation (9) exposes a second limitation: local
amplification alone does not imply one planetary sign. Weakly coupled patches
can freeze into a mosaic of opposite chiral domains. A global claim therefore
needs at least one of the following, measured rather than assumed:

- a globally correlated signed seed;
- sufficiently rapid migration and competitive replacement to coarsen domains;
- a serial bottleneck in which one fixed patch seeds its descendants; or
- a shared environmental bias that consistently tilts fixation probabilities.

OPH calls agreement among patches overlap consensus, but the word
"consensus" is not a proof that equation (9) synchronizes. The source branch
must supply the graph, transport coefficients, noise, environmental cycling,
and a convergence or fixation receipt.

## 5. Audit of the four proposed source routes

| Route | Demonstrated role | OPH placement | What it does **not** yet establish |
|---|---|---|---|
| asymmetric autocatalysis and antagonism | Soai chemistry amplifies a small initial excess; Frank kinetics and Viedma-type recycling show how nonlinear gain plus removal of the rival hand can approach purity | primarily $\kappa$, with cross-inhibition or recycling as repair | prebiotic availability of the particular Soai reagents; the original sign; a coupled peptide/RNA/metabolic network |
| chiral mineral surfaces | mirror-related calcite faces adsorb D- and L-aspartate differently; the largest reported local excess was about 10% under the tested conditions | a boundary source $h_i$, a surface record, or a compartment-specific selector | opposite crystal faces favor opposite signs, so an unbiased mineral population need not yield a global excess; adsorption alone is not replication |
| circularly polarized light | irradiation of initially achiral interstellar-ice analogues produced alanine excess up to 1.34%, with sign and magnitude dependent on helicity and dose | an externally signed $h_i$ or initial $e_0$ | delivery, survival, global helicity history, amplification, or continuity to terrestrial biochemistry |
| electroweak parity violation | parity-violating energy differences between enantiomers are physically allowed; molecular calculations give extremely small and conformation-sensitive values | a possible microscopic contribution to $h$ | a robust universal sign for relevant aqueous networks or a macroscopic excess without nonequilibrium gain |

The mineral result is from
[Hazen, Filley, and Goodfriend (2001)](https://pmc.ncbi.nlm.nih.gov/articles/PMC33239/).
Its strongest lesson is local rather than planetary: the handed surface is a
boundary condition, and the abundance, orientation, weathering, and transport
of those boundaries must enter the ledger.

The circular-polarization result is similarly a credible seed receipt, not a
complete origin scenario. Extraterrestrial material can also contain much
larger excesses: L-isovaline excesses of $`18.5\pm2.6\,\%`$ in Murchison and
$`15.2\pm4.0\,\%`$ in Orgueil were reported by
[Glavin and Dworkin (2009)](https://doi.org/10.1073/pnas.0811618106).
Those authors associated the enrichment with parent-body aqueous alteration;
the observation does not identify circular polarization as the sole cause or
prove inheritance by terrestrial biology.

For scale, a representative calculated parity-violating energy difference for
amino-acid enantiomers is of order $10^{-14}\ \mathrm{J\,mol^{-1}}$, with
important conformational dependence
[Tranter (1985)](https://doi.org/10.1080/00268978500102741). At $300\ \mathrm K$,
a two-state equilibrium estimate gives

```math
|e_{\rm eq}|\simeq\frac{|\Delta E_{\rm PV}|}{2RT}
\sim2.0\times10^{-18}.
\tag{10}
```

The exact molecular value and even the preferred sign can depend on molecule,
conformation, solvation, and calculation. Equation (10) is therefore a scale
audit, not a universal sign prediction. Parity violation can at most be a tiny
seed in this account; a high-gain nonequilibrium mechanism has to pass
the threshold in equation (8).

## 6. Replication architecture is an assumption, not a free theorem

The common RNA-world argument is that opposite-handed monomers inhibit
template-directed extension, turning homochirality into a copying advantage.
That is an important candidate repair law, but it is not architecture
independent. A cross-chiral RNA polymerase ribozyme can catalyze the joining of
substrates of the opposite hand
[Sczepanski and Joyce (2014)](https://doi.org/10.1038/nature13900).
Consequently, OPH may not assume that every possible replicator rewards the
same-hand branch. The proposed prebiotic carrier must measure its own
stereochemical copying matrix:

```math
K=
\begin{pmatrix}
k_{L\leftarrow L} & k_{L\leftarrow R}\\
k_{R\leftarrow L} & k_{R\leftarrow R}
\end{pmatrix},
\tag{11}
```

together with degradation, inhibition, and resource-consumption rates. The
branch claim is earned only if the measured dynamics reduce to a positive
effective $\kappa$ and satisfy equation (8) over repeated checkpoints.

## 7. Evidence required for an end-to-end OPH certificate

A closing experiment should declare the source before inspecting the outcome
and publish enough information to reconstruct the branch ledger. At minimum:

1. **Input receipt:** achiral or measured near-racemic feedstock, contaminant
   bounds, isotopic provenance where relevant, and blind D/L controls.
2. **Seed receipt:** signed and calibrated mineral-face population, light
   Stokes parameters and dose, or a specified alternative source; a no-seed
   control must measure stochastic fixation.
3. **Kinetic receipt:** absolute $L,R$ concentrations, not enantiomeric excess
   alone; fitted autocatalysis, antagonism, racemization, and recycling rates
   with uncertainty and held-out predictive checks.
4. **Record/readback receipt:** demonstrate that a surviving chiral template,
   catalyst, or crystal population predicts the next cycle's products better
   than a shuffled-record control.
5. **Repair threshold:** for a source-off, approximately autonomous regime,
   independently estimate $\kappa$ and $\mu$, then test the model
   bifurcation near $\kappa=2\mu$ rather than merely fitting a final excess.
   For time-dependent chemistry, freeze and test the integrated exposure model
   instead of applying an instantaneous threshold.
6. **Spatial receipt:** sample multiple compartments and both mineral-face
   orientations; report domain size, exchange rate, sign concordance, and
   whether coarsening follows equation (9) or a declared alternative.
7. **Checkpoint receipt:** retain sign and copying advantage through dilution,
   wet-dry or freeze-thaw cycles, compartment division, and serial transfer.
8. **Coupling receipt:** show how the amino-acid branch becomes linked to the
   sugar/nucleotide branch instead of treating either sign in isolation.

Decisive falsifiers include: no amplification when $h$ or $e_0$ is blinded;
loss of the branch under realistic racemization; equal copying fitness after
resource accounting; opposite signs persisting in weakly coupled patches; or a
result explained by contamination, selective detection, or unrecorded
post-processing.

## 8. What the conditional model establishes and what remains

Under the declared reduced-model assumptions, the audit establishes:

- Homochirality can be represented as fixation of a persistent chiral record
  in a bounded replication-and-repair patch.
- Exact symmetry supplies no absolute hand.
- A nonzero seed is amplified by the Frank normal form according to cumulative
  exposure (4), with logarithmic target exposure (5); complete fixation also
  requires unbounded or sufficiently large exposure.
- The autonomous, constant-coefficient, source-off error model has a pitchfork
  threshold $\kappa>2\mu$; time-dependent or biased chemistry requires the
  full dynamics.
- Autocatalytic, mineral, circular-polarization, and parity-violation routes
  occupy different terms in the model and cannot be substituted for one
  another without evidence.
- Local fixation does not imply global homochirality; transport, bottlenecking,
  or correlated bias must close the spatial consensus step.

The audit specification stops without promoting an unknown historical pathway
to a theorem. Its empirical target is the following finite certificate:

> Produce one end-to-end source-to-fixation receipt: freeze a chemically and
> geochemically plausible network, measure
> $(\kappa,h,\mu)$ without target leakage, demonstrate record-conditioned
> copying and repair across serial checkpoints, and show spatial sign
> concordance plus amino-acid/sugar coupling against blinded controls.

This is an observational acceptance condition, not a completed physical
derivation. The OPH audit explains the proposed architecture by which a
microscopic or stochastic chiral difference becomes a durable biological
record and states what would test that account on a chosen chemical branch. It
deliberately does not claim that observer consistency fixes the
terrestrial L-amino-acid sign or reconstructs which combination of
astrophysical, mineral, and autocatalytic events occurred.

## References

- F. C. Frank, “On spontaneous asymmetric synthesis,” *Biochimica et
  Biophysica Acta* **11**, 459–463 (1953),
  [doi:10.1016/0006-3002(53)90082-1](https://doi.org/10.1016/0006-3002(53)90082-1).
- K. Soai, T. Shibata, H. Morioka, and K. Choji, “Asymmetric autocatalysis and
  amplification of enantiomeric excess of a chiral molecule,” *Nature* **378**,
  767–768 (1995), [doi:10.1038/378767a0](https://doi.org/10.1038/378767a0).
- R. M. Hazen, T. R. Filley, and G. A. Goodfriend, “Selective adsorption of
  L- and D-amino acids on calcite: Implications for biochemical
  homochirality,” *PNAS* **98**, 5487–5490 (2001),
  [doi:10.1073/pnas.101085998](https://doi.org/10.1073/pnas.101085998).
- C. Viedma, “Chiral symmetry breaking during crystallization: complete chiral
  purity induced by nonlinear autocatalysis and recycling,” *Physical Review
  Letters* **94**, 065504 (2005),
  [doi:10.1103/PhysRevLett.94.065504](https://doi.org/10.1103/PhysRevLett.94.065504).
- P. de Marcellus et al., “Non-racemic amino acid production by ultraviolet
  irradiation of achiral interstellar ice analogs with circularly polarized
  light,” *Astrophysical Journal Letters* **727**, L27 (2011),
  [doi:10.1088/2041-8205/727/2/L27](https://doi.org/10.1088/2041-8205/727/2/L27).
- D. P. Glavin and J. P. Dworkin, “Enrichment of the amino acid L-isovaline by
  aqueous alteration on CI and CM meteorite parent bodies,” *PNAS* **106**,
  5487–5492 (2009),
  [doi:10.1073/pnas.0811618106](https://doi.org/10.1073/pnas.0811618106).
- G. E. Tranter, “The parity violating energy differences between the
  enantiomers of alpha-amino acids,” *Molecular Physics* **56**, 825–838
  (1985),
  [doi:10.1080/00268978500102741](https://doi.org/10.1080/00268978500102741).
- J. T. Sczepanski and G. F. Joyce, “A cross-chiral RNA polymerase ribozyme,”
  *Nature* **515**, 440–442 (2014),
  [doi:10.1038/nature13900](https://doi.org/10.1038/nature13900).
