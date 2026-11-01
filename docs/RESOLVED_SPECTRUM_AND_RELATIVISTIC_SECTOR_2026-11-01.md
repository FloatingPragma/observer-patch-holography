# Resolved detector spectrum and a constructed massless relativistic sector

The [subsequent interacting-reference and mass extension](INTERACTING_VACUUM_AND_RELATIVISTIC_MASS_2026-11-01.md)
constructs a global cube ground and a one-dimensional massive relativistic
continuation with explicit coupling and reference choices. The results below
retain their earlier state domains and the original massless-limit conclusion.

Three new results answer different parts of the particle question. The
[resolved-spectrum calculation](../model/resolved_particle_spectrum/README.md)
recovers all eight weights in the protected cube's source-accessible response
from simulated pointer counts. The
[relativistic construction](../model/relativistic_sector_audit/README.md)
gives an explicit eight-valley massless Weyl limit of the original action,
with error bounds for the full interacting evolution on finite-particle
initial domains. A compiled momentum/occupation detector now identifies sixteen signed
model modes on a periodic lattice. These constructions retain distinct
carriers and approximation domains. A retrospective comparison finds no
direct match of the new massless flavors to the charged leptons, photon or
complete neutrino mass spectrum.

This extends the [five-issue campaign](THREE_DIMENSIONAL_PROGRESS_2026-11-01.md)
and [first internal-observer readout](INTERNAL_PARTICLE_MEASUREMENT_2026-11-01.md).
The new work advances the protected-sector and particle tasks
[#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139),
[#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143),
the field limit [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142),
and the contact-control prerequisite for
[#1145](https://github.com/FloatingPragma/observer-patch-holography/issues/1145).
The earlier adaptive-reference construction addresses the fifth selected
issue, [#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138).
No remote issue is closed or updated by this local work.

## All eight accessible response weights

The original cube Hamiltonian fixes the eight frequency intervals. On that
known support, both pointer quadratures at times `0,39,78,...,273` determine
the unknown weights. The inference takes counts, not the previously derived
weights, as input. A rational Lagrange inverse is independently certified
against the true frequency intervals. Its infinity norm is below two and
its residual below `10^-18`.

| Original H energy | Measured transition E-3 | Reconstructed normalized weight |
| ---: | ---: | ---: |
| 1.816498920498 | -1.183501079502 | 0.002457894673 |
| 1.860537501391 | -1.139462498609 | 0.005641397395 |
| 1.955137632057 | -1.044862367943 | 0.005890881572 |
| 2.108692846015 | -0.891307153985 | 0.007131698528 |
| 2.167639733928 | -0.832360266072 | 0.004528321839 |
| 3.015861345573 | +0.015861345573 | 0.421578515341 |
| 3.030769652593 | +0.030769652593 | 0.415796456292 |
| 3.044862367943 | +0.044862367943 | 0.136974834359 |

Each reconstructed weight has error below `0.000258394`, conditional on
the declared simultaneous sampling event. Every interval is strictly
positive and contains the independently derived source weight. The
sampling failure bound is below one percent for the specified independent
reset model. These normalized weights omit the separately known elastic
weight `1/16`; the remaining response is multiplied by `15/16` in the
actual pointer correlation.

The computation uses 16 seeded binomial samples with 4,831,838,208 trials
per setting, totaling 77,309,411,328 simulated reset trials. It samples a
binary64 proxy to rigorously enclosed ideal-word probabilities. It does
not individually execute the billions of native pulses or draw complete
individual native histories. The exact source word, every abort address
and all complete repair outcomes remain specified; a norm bound covers
their difference from the computed comparator. The maximum combined
source, splitting, certified arithmetic and probability-conversion budget
is below `0.000008569`.

This acceleration has a proof: consecutive spin-complete edges and the
contact layer intertwine through the eta lift. The ideal chronological
word therefore has an exact 36-dimensional singlet representation.
Every long-time power is evaluated with Gaussian integers and outward
error rounding; the final time uses 3,670,016 steps but only 24 matrix
multiplications. No equality between the rounded native word and its
compressed ideal comparator is assumed.

The result closes the earlier prospective *known-support weight*
reconstruction premise. It is not blind inference of arbitrary unknown
frequencies, nor spectral tomography of all 1,820 protected states. The
pointer measures differences relative to the excited reference energy
three. Those transition lines are not eight relativistic rest masses.

## Constructed relativistic energies and masses

Keep the original action and choose the scaling `c=v/a`, with `v>0` fixed:

```
H_a = dGamma[-(v/a) sigma.sin(k)] + U(a) sum_x n_up(x)n_down(x).
```

For each of the eight nodes `nu in {0,pi}^3`, put `k=nu+a p`. The explicit
isometry maps a continuum spinor to the corresponding lattice patch with
factor `a^(-3/2)`. Its exterior powers embed every fixed finite-particle
sector. All eight nodes and both energy signs are retained. The limit is

```
h_nu(p) = -v sum_j (-1)^(nu_j/pi) sigma_j p_j,
E_±(p) = ±v|p|,
E²-v²|p|² = 0.
```

The limiting spinors obey the massless Weyl equation. Exact SU(2) frame
checks give four valleys of each chirality; pairing opposite chiralities
satisfies the Lorentz Clifford algebra. This does not select four physical
Dirac species. The fixed regulator and bounded momentum domain are not
Lorentz invariant; the limiting differential expression has that symmetry.

At `v=1,a=1/64`, every valley has these energies in model units:

| Momentum p | Finite-lattice energies | Limiting relativistic energies | Limiting mass |
| --- | ---: | ---: | ---: |
| (0,0,0) | 0 | 0 | 0 |
| (1,0,0) | ±0.999959310393 | ±1 | 0 |
| (1,1,0) | ±1.414156018578 | ±sqrt(2) | 0 |
| (1,1,1) | ±1.731980331101 | ±sqrt(3) | 0 |

The receipt contains exact rational energy enclosures. The nonzero
finite-spacing mass-shell defect is a regulator error, not a positive
rest mass or a tachyon prediction. No joule, GeV or kilogram calibration
is inferred.

For `|p_j|<=Lambda`, `N<=Nmax`, nonnegative supplied `U(a)` and
`a Lambda<pi/2`, the comparison controls the full unprojected evolution:

```
||exp(-it H_a) W_a - W_a exp(-it H_Weyl)||
 <= min(2, |t| [Nmax v a² Lambda³/2
               + U(a) floor(Nmax/2) sqrt(8(a Lambda/pi)³)]).
```

A countable-lattice contact argument bounds all scattering out of the
initial band. The band is not artificially made invariant. For example,
`Nmax=2, Lambda=v=U=t=1, a=1/64` gives error below `0.001237`, including
all outgoing states. The error vanishes at fixed particle number and
momentum cutoff when `U(a)=o(a^(-3/2))`. For one particle, contact vanishes
exactly for every allowed coupling. The proved limiting interaction is
therefore free. The critically tuned moving-pair scaling `U~a^(-2)` is
outside this estimate.

## The detector now identifies internal particle flavors

The [compiled model-particle detector](../model/particle_species_detector/README.md)
measures momentum and spin occupation on the original periodic `8^3`
Weyl/contact lattice, retaining all 1,024 one-particle output modes. It
attaches sixteen exact energy eigenmodes to detector addresses: eight
emergent Weyl flavors and both energy signs. Their number is one, CAR
parity is odd and their exchange algebra is fermionic. At finite spacing,
the flavor is a valley-momentum label, not an extra microscopic field.

For the selected offset `k_z=nu_z+pi/4` and `c=1/8`, every flavor has
energies `E=±sqrt(2)/16=±0.0883883476483` in model units. The positive
branch has the following helicity in its derived spin frame:

| Valley flavor b | Chirality | Positive-energy helicity | Continuum mass |
| --- | ---: | ---: | ---: |
| 000 | -1 | -1/2 | 0 |
| 001 | +1 | +1/2 | 0 |
| 010 | +1 | +1/2 | 0 |
| 011 | -1 | -1/2 | 0 |
| 100 | +1 | +1/2 | 0 |
| 101 | -1 | -1/2 | 0 |
| 110 | -1 | -1/2 | 0 |
| 111 | +1 | +1/2 | 0 |

The negative-energy branch has the opposite helicity. It is retained as
a negative-energy branch, without assuming a physical antiparticle
interpretation. Exact source equations `H Pi=E Pi` make this a spectral
energy readout on the declared sixteen-mode family. The other 1,008
one-particle outcomes are retained and do not receive these energy labels.
The corresponding complement commutes with H.

The implementation uses a physical-mode Fourier network, initialized
single-occupation preparation and complete occupation pointers. Its
65,822,721 formal source pulses have a total complete-output error below
`8.31e-9`. Thus the specified native protocol reports the correct prepared
label with unconditional probability at least `1-8.31e-9`. This is proved
against the full native instrument, including off-number states and
aborts. The executed numerical array is the ideal Gaussian comparator on
all 1,024 one-particle rows, not a literal `2^1024`-row source simulation.
The original-H continuation is a declared model-law service and is not
included in that pulse count.

For the separate periodic refinement `a=2^-m, L=8*2^m, c=1/(8a)`, this
same detector family has
`E_±=±sin(a*pi/4)/(8a) -> ±pi/32≈±0.098174770425` at fixed model momentum.
This massless limit follows from its own exact dispersion. It does not
transfer the infinite-volume wave-packet theorem to a finite plane wave.

## Comparison with observed masses

The comparison uses a fixed **PDG 2025** baseline and is retrospective.
No measured mass enters the Hamiltonian, detector or continuum calculation.
A particle's energy depends on its momentum; the comparable invariant is
its rest mass, not its nonzero finite-momentum detector energy.

| Candidate physical identification | Observed rest mass | This sector's continuum rest mass | Result |
| --- | ---: | ---: | --- |
| Electron | 0.51099895000 MeV/c² | 0 under any positive unit conversion | Does not match |
| Muon | 105.6583755 MeV/c² | 0 under any positive unit conversion | Does not match |
| Tau | 1776.93 MeV/c² | 0 under any positive unit conversion | Does not match |

The numbers and their retained uncertainties come from the PDG
[electron](https://pdg.lbl.gov/2025/listings/rpp2025-list-electron.pdf),
[muon](https://pdg.lbl.gov/2025/listings/rpp2025-list-muon.pdf) and
[tau](https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf) listings.
A common positive change of units keeps zero mass zero. The model's mass
ratios are undefined `0/0`; treating the equal finite-momentum energies
as masses would instead give ratios one, also failing the observed
muon/electron and tau/electron hierarchy.

Masslessness alone also does not identify a photon: these branches have
helicity magnitude one half, while the photon has spin one. This excludes
a direct spin-preserving identification, not a composite photon or a
photon in another sector. See the [PDG photon listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-photon.pdf).

Assigning three of these massless flavors to the complete neutrino mass
spectrum would give all mass-squared differences zero. The selected
normal-ordering fit in PDG Table 14.7 instead gives approximately
`7.49e-5` and `2.459e-3 eV²` in `c=1` units. The retained published
three-sigma ranges exclude zero; they are statistical fit ranges, not
exact mathematical bounds. This does not exclude one massless neutrino
or determine any absolute neutrino mass. See the
[PDG neutrino review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-neutrino-mixing.pdf).

These are scoped mismatches of the constructed sector, not a theorem that
OPH cannot produce other particles. The new diagnostic adds no successful
physical identification or eligible prospective comparison. Exact input
values, source URLs and the computed comparisons are retained in
[physical_targets.json](../model/particle_species_detector/physical_targets.json)
and the detector receipt.

## Physical identification still to construct

The empty Fock reference is not a ground state because the original
one-particle action has negative energies. A stable physical vacuum and
its particle/antiparticle representation remain unconstructed here.
The new branch now has its own momentum/occupation detector; it is not
yet attached to the protected cube's five-seam repair service or its
neutral-response detector. Electric charge, gauge representation,
physical spin and a common physical clock and units still need a common
identification. There is no new assignment to an electron, neutrino,
photon or other observed species.

The separate moving pair retains its derived nonrelativistic inertial
mass two and binding energy gamma squared; neither number is imported as
a relativistic rest mass. No massive relativistic sector, interacting
relativistic vacuum theory, finite-density limit or 3+1-dimensional
Einstein coarse-graining theorem is claimed. The supplied simulator,
observer, scheduler, repair and memory laws remain explicit inputs.

The resolved-response and continuum results are recorded as FA-62 and
FA-63; the subsequent model-particle detector is FA-64. The registry has
59 realized scoped attachments among 64 reviewed domains. All 346 claim
classifications and all 84 canonical premise dispositions are preserved;
no eligible prospective comparison or integrated completion gate is
added. The PDG mismatch comparison is explicitly retrospective.

## Why the masses fail and what needs to change

The comparison failure is a limitation of this constructed sector, not
an arithmetic error. The original contact operator annihilates every
one-particle state, for every allowed U. Changing U while retaining the
empty-reference one-particle calculation therefore cannot generate a
mass. Changing c changes the slope of the cone, while every node stays
gapless. The proved finite-particle continuum scaling also makes the
contact contribution vanish.

The next calculation must first construct a stable interacting reference
for the original action and determine its excitation spectrum. That may
reveal a massive phase already present in the law; the present result
does not decide this. It is premature to assert that the microscopic law
must be changed. If no suitable phase exists on a controlled domain, an
additional term must be justified by the simulator or repair/memory law,
with its selection and coefficient recorded as new obligations.

A precise target is an effective opposite-chirality block
`H_eff(p)=v alpha.p + beta Delta`, where beta anticommutes with the three
kinetic alpha matrices and Delta acts on flavor. If Delta is a constant
Hermitian flavor matrix commuting with these spin matrices, its eigenvalues
delta_i give `E_i²=v²|p|²+delta_i²` and rest masses `|delta_i|/v²` in these
model coordinates. The existing algebra certificate supplies the possible
Clifford structure; the source does **not** yet supply such a derived Delta.
Writing down or fitting this matrix would not establish mass generation.

Required work is therefore: derive the interacting reference and any
mass-generating order or self-energy, prove stable relativistic excitation
poles with their errors, identify their gauge charges and detector
couplings, and test their dimensionless mass ratios. A common physical
unit can be assigned only after that dictionary is justified. The
massless physical gauge sector must also be constructed; retaining a
massless spinor is not a photon construction. This sequence preserves
the useful massless theorem while addressing the actual missing physics.
