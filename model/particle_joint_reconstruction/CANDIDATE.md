# A conditional gauge and matter construction

The candidate joins the selected finite observer program to an explicitly
assumed chiral gauge and Higgs extension. Its executable result is the exact
quadratic spectrum of that extension at inherited calibrated coefficients,
checked against two independent action producers. The construction supplies
an interface for an interacting four-dimensional limit; it does not supply
that limit or an embedding into the selected finite program.

Run `python3 -m model.particle_joint_reconstruction.candidate --write` to
recompute `candidate_results.json`, and use `--check` to replay it. The
receipt retains the original theorem premises, the exact selected primary
member, each calibration packet, and their source hashes. Every assumption
has explicit references and a false `proved_on_primary` field. The
[dependency graph](DEPENDENCIES.md) gives the joint reconstruction domains.

## The common source boundary

The selected member has 60 CAR annihilation modes per cell: the native
15-component internal representation tensored with four Dirac components.
Its finite-range Wilson action, source-Gauss polynomial preparation, conjugate
idle reference and source Hamiltonian clock are fixed by
[the primary contract](../particle_common_contract/PROOF.md). The clock has
no supplied conversion to GeV or proper time. The selected polynomial density
is stationary and has controlled ground energy excess. It is not an exact
vacuum, and its purified joint state is not stationary under source-only
evolution.

The candidate extension uses three all-left families and two-component
Weyl spin, giving 45 internal Weyl components and 90 static CAR coordinates,
with a complex Higgs doublet. It consumes the exact current and Yukawa
intertwiners from [the conditional action](../conditional_standard_model/PROOF.md).
Shared internal matrices are an algebraic compatibility condition. They do
not identify the two spatial Hamiltonians, carriers, states or clocks.
In particular, the 60-to-90 comparison is not a basis change or a proof of
three physical families. A local regulator, chiral measure, mirror control,
compatible state preparation, physical time conversion and source map are
explicit extension premises. The primary preparation and locality theorems
apply to their own source until a transfer theorem establishes these maps.

The global gauge group and charge lattice also require a declaration.
Anomaly-free Lie-algebra representations do not choose them. Gauge theta
angles and the additive vacuum energy have no selected physical values in
this candidate. They are not hidden fitted constants in the receipt.

## Two spectrum arms

The first arm calls `conditional_standard_model.spectrum.build_catalogue()`
without changing its inputs. It retains the inherited common-scale
modified-minimal-subtraction coefficients, the actual symmetric 45-row
Yukawa mass kernel, the antisymmetric 90-coordinate creation pairing and
the calibrated Cabibbo–Kobayashi–Maskawa matrix. The neutral gauge Hessian
and bosonic coefficients agree exactly with the symbolic producer in
`particle_candidate_atlas.conditional_standard_model`. The broken-phase
quadratic boson count is 28. These checks test action consistency at supplied
coefficients.

The inherited boson tree coefficients are approximately 80.848 GeV,
92.2264 GeV and 131.280 GeV for the charged weak boson, neutral weak boson
and radial Higgs mode. Their input scheme and scale do not make them
measured pole masses. Directly assigning this catalogue as the physical
spectrum fails: the minimal renormalizable action also has three massless
neutrinos and cannot reproduce nonzero oscillation splittings. Quark
running parameters are not free asymptotic quark poles. No interacting
width, confinement spectrum or inclusive detector probability is computed
by these quadratic identities.

The second arm retains those charged and boson inputs and separately calls
`calibrated_matter.neutrino_report()`. Its dimension-five Weinberg operator
uses the same Higgs and invariant current tensor. The code verifies 12
current identities, the symmetric electrically neutral mass kernel,
unitarity of the supplied Pontecorvo–Maki–Nakagawa–Sakata rotations and
both calibrated mass-squared splittings. The minimal catalogue is preserved
verbatim; its zero neutrino entries are not overwritten and are not used
as a fit to oscillations. Mass eigenstates and flavor labels remain distinct.

This neutrino packet uses its declared oscillation dataset and normal
ordering. The lightest mass and Majorana phases are declared unmeasured
choices. The packet supplies no certified low-energy-to-common-scale
running or matching error. A finite effective-theory window and its
truncation error are therefore required; an all-energy claim additionally
requires a specified ultraviolet completion. Recovering the fitted
splittings is a calibration identity. No independent mass prediction is
counted.

Physical compatibility with all reference targets is conditional on
radiative, infrared, resonance and neutrino matching, with their errors,
and on the detector dictionary. No complete numerical matching witness is
supplied. The receipt does not adjust coefficients to force agreement or
claim that unspecified corrections must succeed. The available exact
feasibility result is narrower: the inherited charged and neutrino tensors
can occur in one electrically neutral gauge and Higgs action, and the
symbolic and calibrated quadratic producers agree. It supplies a coherent
comparison action for testing those physical matching hypotheses.

## Conditional interacting reconstruction

Let the extension be a regulated family on a declared local observable
domain, with the supplied chiral action and calibrated effective coefficients
in a specified energy window. Assume a compatible physical vacuum limit,
one reference and time normalization, a gauge-invariant chiral measure and
mirror control. Assume its renormalized physical correlations converge in
the joint spacing, cutoff and infinite-volume limit and satisfy the
positivity, regularity, Euclidean covariance, locality and clustering
hypotheses needed for relativistic reconstruction on four-dimensional
spacetime. Require a surviving interaction witness in a fixed renormalized
generating-field and current dictionary, incompatible with a free/Wick
hierarchy under the admitted identifications, or an appropriate nonfree
scattering criterion. Connected higher correlations of nonlinear composites
in a free theory do not meet this requirement. These assumptions yield a
nontrivial interacting reconstructed theory on that domain. Nonzero bare
coefficients alone do not establish the interaction premise.

Identification with physical particles additionally requires conserved
current and Lorentz-spin intertwiners and a detector interpretation, with
radiative and infrared control. If a source-sector embedding is claimed,
its isometry must control the entire rectangular generator residual
\(H_{\mathrm{source}}J-JH_{\mathrm{field}}\) in the same energy reference and
clock. The inherited executable counterexample has identical compressed
Hamiltonians while a nonzero off-image residual persists. Compression
alone cannot attach the spectrum. Completeness of the physical spectrum,
including composite hadrons and unstable resonances, is a separate
statement from the elementary quadratic channel census.

The ten explicit assumption rows in `candidate.py` are the interface to
these statements. None is a discharge of the corresponding primary
construction premise. The inherited finite theorems constrain their
own inputs and supply exact algebraic steps within this conditional join.

## The gap belongs to a separate pure gauge branch

The photon null vector is exact. Under the assumed continuum dispersion,
photon excitations have energies \(|p|\) for arbitrarily small nonzero
momentum. The receipt exercises this with a decreasing momentum sequence.
Thus the photon-containing continuum scenario has no positive global
vacuum gap. This is a conditional consequence of that scenario, not a
spectrum calculation for the selected primary.

For a separate pure compact-simple gauge source,
[the dense-correlation criterion](../particle_yang_mills_gap_contract/PROOF.md)
provides a stronger useful conditional conclusion. A reconstructed
nonnegative Hamiltonian and dense centered vacuum vectors satisfying
\(C_A(kh)\le B_Aq^k\), with one common physical \(h>0\) and \(0<q<1\),
imply spectral exclusion below \(-\log(q)/h\). A surviving variance
\(v_0>0\) and energy bound \(M<\infty\) give a finite upper threshold
\(M/v_0\). Correlation convergence at fixed times permits this transfer
without an assumed equality of finite and continuum spectra.

The selected action scales as \(H_n=a_n^{-1}\widehat H_n\) at fixed electric
and magnetic coefficients. A constant dimensionless finite gap therefore
gives a divergent physical threshold; it does not establish a finite
continuum mass. The selected low-energy polynomial state does not directly
supply the positive vacuum correlations required by the criterion.

Selecting a Lie-algebra factor from the coupled product-group matter
source does not construct a pure gauge invariant Hilbert sector. Freezing
other links at the group identity is nonnormalizable in Haar space, and
projecting their wavefunctions to constants need not preserve the
Hamiltonian. A gap for the strong sector of the full theory needs a proved
sector or decoupling map, with its clock. The pure gauge branch therefore
has its own source and continuum premises. The receipt establishes neither
the Yang–Mills continuum nor the Clay mass-gap claim.
