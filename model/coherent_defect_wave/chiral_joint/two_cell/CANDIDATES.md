# Candidate channels on the exact finite source

The calculation identifies operator channels and complete finite energy
records. A channel is a specified way to prepare and interrogate the source;
it need not create one elementary particle. The same channel can overlap
with an elementary excitation, a composite, several levels, or an excitation
bound to its apparatus boundary. These alternatives require discriminating
measurements on the same Hamiltonian family.

The source is the uniform open two-cell control of dimension 614400, with
both charged-lepton Yukawa couplings, scalar penalties and link penalty one,
and both quark Yukawa couplings zero. Its reference is the complete neutral
physical ground at L=2 and 2Sz=0, prepared with original probability 1/1053.
Time uses this same dimensionless Hamiltonian. The source parameters and
finite memories are declared architecture. No measured mass enters them.

## What is identified exactly

| Channel | Actual source operator and charge | Established classification | Physical identity owed |
|---|---|---|---|
| e^c addition | Native R15 row 9 creator, Y=+1, weak/color singlet; source delta L=-1 | Odd CAR channel, complete compensated energy response, exact odd-L Kramers pairing | A physical positron or charged lepton; an elementary isolated pole |
| e^c removal | Adjoint of the same source operator, Y=-1; source delta L=+1 | Conjugate operator channel, reflected opposite-boundary response relation | A physical electron or other named antiparticle |
| Neutral scalar occupation | Local scalar-memory projector N_H | Gauge-neutral, even, source-spin-scalar operator with complete binary PVM | A propagating Higgs or any relativistic scalar particle |
| Neutral Yukawa composite | Local matter-pair/scalar interaction density H_Y | Gauge-neutral even operator, source-spin-scalar, nontrivial same-H dynamics | An elementary scalar, bound state or multiparticle response |

`candidate_data.py` reads the actual native current matrices. Row 9 has
eleven zero generator eigenvalues and hypercharge one, and both its color
and weak Casimirs vanish. The algebraic combination T3+Y equals one on that
row. This does not identify an emergent electromagnetic current: the finite
quantum scalar reference has not supplied that physical field connection.

The active native family contains the e^c singlet and L doublet exactly:

| Native R15 row | Channel | Y | T3 | Algebraic T3+Y | Weak Casimir | Source L |
|---|---|---:|---:|---:|---:|---:|
| 9 | e^c | 1 | 0 | 1 | 0 | -1 |
| 10 | L component | -1/2 | 1/2 | 0 | 3/4 | 1 |
| 11 | L component | -1/2 | -1/2 | -1 | 3/4 | 1 |

All three are color singlets. The actual off-diagonal weak generators join
rows 10 and 11 into one doublet; these are not merely matching diagonal
dimension counts. The full parent carrier retains the quark fields, while
this calculation sets their Yukawa couplings to zero and keeps their empty
reducing sector. It therefore calculates no quark-particle spectrum.

The four labeled e^c probes are two cells times two source spin components
of one internal channel. They are not four species or three generations.
Their CAR relations are checked on every one of the 4096 twelve-mode
occupations, including other occupied modes as spectators. This is the
declared fermionic source algebra, not an independent proof of continuum
particle statistics.

The same-boundary antiunitary Theta reverses 2Sz and squares to source
fermion parity, equal to `(-1)^L`. Odd-L energy levels therefore have exact
Kramers partners. The distinct antiunitary A conjugates charge and reflects
the two cells: right addition matches left removal with the opposite
reflected boundary and transformed reference. Neither statement identifies
Wigner spin or proves physical CPT. Indeed, the open one-axis Hamiltonian
does not commute with the full internal spin SU(2); a spin projection is
not a relativistic spin classification.

## Apparatus and neutral measurements

Every charged probe uses the retained reservoir and opposite-charge
compensator of [the complete instrument](INSTRUMENT_PROOF.md). Its final
boundary is `q=n_r,out-initial_setting`, at the selected cell. Source L and
2Sz changes are recorded with the corresponding apparatus changes. Addition
and removal failures remain part of the complete law. Energy records use
the spectrum module's phased physical embedding D B and the same H.

At either cell, scalar occupation has complete effects `{N_H,I-N_H}`.
The local composite interaction has the complete effects
`{(I+H_Y/8)/2,(I-H_Y/8)/2}`. Positivity follows from the eight unit-norm
pair/scalar terms in H_Y, giving norm at most eight. Their square-root
measurement operators preserve local Gauss, L, 2Sz and parity. The exact
320-state local calculation checks all twelve actual current commutators
and source-spin-scalar identities. Both probes have nonzero commutators
with the same local interacting Hamiltonian. Their neutral spectral weights
are not calculated by this classification module.

Connected correlations must subtract the reference expectation of each
neutral operator. A nonzero scalar-density response alone can be a pair
threshold or a composite response; it does not establish an elementary
scalar. Comparing the occupation, interaction-density and charged/Higgs
conversion channels is a way to test their shared or different spectral
overlap, without choosing a particle name first.

## Measurements which discriminate the alternatives

1. **Boundary localization.** Record the site-resolved charged response and
   source-energy distribution while moving the compensator and preparation
   farther apart. Repeat with reflected opposite charges as the exact A
   control. A level that follows the compensator or grows with an attached
   flux configuration is evidence for a boundary-bound excitation. A
   two-cell energy difference alone cannot separate that contribution from
   an intrinsic excitation energy.
2. **Size and cutoff dependence.** Specify a common Hamiltonian family,
   scale conventions and state preparation, then increase volume and the
   scalar/link cutoffs separately. Carry controlled truncation and numerical
   errors. Merely observing similar values at two sizes is not an error
   bound or a proof of a limiting particle.
3. **Translation and dispersion.** The open source has no proved conserved
   lattice momentum. Its one-sided finite raising link obstructs the stated
   local exchange translation. The separately declared two-link periodic
   source has exact finite translation phases and complete phase effects,
   but its carrier and Hamiltonian differ. Its charged Gauss-compatible
   apparatus translates the retained compensator with the source, so the
   phase is a joint source/reference observable. It is not an identified
   momentum of an isolated source particle. The separate
   [complete periodic calculation](../periodic_spectrum/PROOF.md) supplies
   its own ground and joint energy/reference-phase records. It does not
   provide a three-dimensional relativistic momentum measurement. On the two-cell
   antiperiodic odd sector there is only one momentum magnitude, so a fit
   of both rest parameter and limiting speed is underdetermined. A
   ring in the specified even antiperiodic family needs at least four cells
   for two magnitudes and six for three; the third permits a nontrivial
   shell test. These menu counts are necessary and do not suffice for a
   particle dispersion. They supply neither the required same-H branches
   and gauge-compatible states nor a controlled relativistic continuum.
4. **Spectral response and particle criterion.** Compute connected responses
   for the actual charged, quantum Higgs and neutral/composite operators.
   The [spectral-edge attachment](../spectral_edges/PROOF.md) reuses the
   existing finite threshold theorem with original preparation weights.
   Stable gapped particles can have poles; charged infrared sectors can
   require dressed thresholds and inclusive records; unstable particles
   require resonances; confined quarks use specified renormalized parameters.
   Finite discrete levels alone prove none of these limiting identifications.
   Operator overlap by itself does not uniquely distinguish an elementary
   field from a composite.
5. **Physical species.** Establish the spatial rotation/little-group action,
   interactions and decay or scattering channels on that same family, then
   compare its complete predictions and error bars with observations. A
   physical energy convention is required for GeV. The current packet gives
   no empirical confidence score, lepton-generation assignment or calibrated
   claim of agreement with measured masses.

The [basis proof](BASIS_PROOF.md), [spectrum proof](SPECTRUM_PROOF.md),
[instrument](INSTRUMENT_PROOF.md) and [symmetry proof](SYMMETRIES.md) describe
which finite parts are executed and where the relativistic interpretation
requires further results.
