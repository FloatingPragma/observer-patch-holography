# Controlled refinement and joint translation readouts

The [open quantum source](SPECTRUM_PROOF.md) has certified complete physical
energy spectra. The [periodic completion](SYMMETRIES.md#4-a-periodic-completion-retaining-the-quantum-links)
has a different Hamiltonian and an exact translation. This proof gives an
actual comparison between those sources, supplies a charged reference for
physical translation readout, and states quantitative acceptance conditions
for a dispersion and refinement claim. No periodic ground energy is taken
from the open source.

## 1. Exact failure of the blank-link spectral transfer

Let J append a blank second link to every open-source basis vector:

`J|mask,link,H_left,H_right> = |mask,link,0,H_left,H_right>`.

This is an isometry from dimension 614400 to dimension 3686400. The blank
link is a gauge singlet, so J preserves the physical Gauss constraints.
Both local quantum Higgs operators and scalar memories are unchanged.

Prepare the normalized **actual** left pair/scalar vector

`psi = H_Y,left |empty>/2`.

It has four nonzero coefficients of magnitude one half. Every local Gauss
generator annihilates psi and J psi. The code checks all 24 generators on
each carrier. Thus this is a physical input, not a naked charged source
mode or a comparison using different gauge constraints.

Evaluate the complete residual

`R psi = (H_periodic J - J H_open) psi`.

Exact sparse execution gives

`J^dagger R psi = 0`, `||R psi||^2 = 3/16`.

Every residual component has the second link excited and therefore lies
outside the image of J. The open raising link cannot transport this left
pair in the corresponding direction from a blank link. The second
oriented quantum link in the periodic construction can. The algebra
retains its complete CAR signs, spin components and scalar alternatives.
For any constant clock scale alpha and energy shift b,

`(I-J J^dagger) [H_periodic J-J(alpha H_open+b I)] psi = R psi`.

A rescaled clock or changed energy zero cannot remove the omitted coupling.
This disproves exact spectral transfer under **this** blank-link embedding.
It does not exclude another explicitly constructed dressed state map.

There is also a complete observable distinction. Let F project onto an
excited second link. Embedded open evolution has probability zero for F
at every time. On the periodic source, `F J psi=0` and
`||F H_periodic J psi||=sqrt(3)/4`. The whole-source norm bound is 30:
the two oriented hopping sums contribute at most ten, the two local
Yukawa sums at most sixteen, and the scalar/link projectors at most four.
The first-order unitary Taylor remainder is at most `t^2*30^2/2`. Hence

`p_periodic(F,t) >= (|t| sqrt(3)/4 - t^2*30^2/2)^2 >0`

whenever the parenthesis is positive. The receipt executes this bound at
`t=1/1000000` and retains `p_no=1-p_yes`. Zero compressed residual therefore
coexists with a nonzero complete detector discrepancy on the actual source.

The full residual still has a useful uniform **upper** bound. The original
link hopping changes by minus one half of its open value, contributing at
most five. On the appended blank link only the ten forward quarter-unit
monomials act, contributing at most five halves. All local scalar, Yukawa
and original electric terms cancel; the new electric projector annihilates
the image. Thus

`||H_periodic J-J H_open|| <=15/2`.

For matching preparations and compressed detector effects the complete
readout error is at most `min(1,15|t|/2)`. This controls short-time comparison.
It does not justify equality of stationary spectra, a transferred ground,
or a vanishing refinement error.

## 2. A physical charged packet with a retained reference

The periodic source's charged Gauss domain requires its compensating charge
to translate too. Supply a two-state apparatus register `|left>,|right>`
for the retained position of one fixed hypercharge minus one. Its local
charge operators are

`Q_left=diag(-1,0)`, `Q_right=diag(0,-1)`.

The register is a specified boundary/reference system, not an identified
source particle. Its swap S obeys `S Q_left S^dagger=Q_right`. Take

`H_joint=H_periodic tensor I_reference`,

`T_joint=T_source tensor S`.

The existing source translation exchanges the source Gauss generators;
S exchanges the reference charges. Therefore T_joint preserves the
complete joint Gauss domain, while H_joint commutes with it and with
T_joint. This is an operator statement on the entire joint carrier,
dimension 7372800, using the actual periodic translation theorem and the
complete two-state reference matrices.

Prepare one right spin-up e^c source occupation, blank quantum links and
scalar memories, and the compensator at the right cell. Source charge
plus one and reference charge minus one cancel locally. All 24 Gauss
constraints hold. For `T_joint^4=I` use the complete projectors

`Q_k=(1/4) sum_(r=0)^3 exp(-i k r) T_joint^r`,

with `k=0,pi/2,pi,3pi/2`. Exact execution gives

`(p_0,p_pi/2,p_pi,p_3pi/2)=(0,1/2,0,1/2)`.

The code retains both zero outcomes, checks the projector products and
sum, and verifies every local Gauss constraint in every projected branch.
The reference position remains in each quantum state. A pointer may copy
the orthogonal phase record; no copying of an arbitrary state is required.

**The measured phase belongs to the joint source/reference translation.**
It is not an independently established momentum of an isolated source
particle. The input is nonstationary. No periodic ground spectrum or
energy-dispersion sample is supplied by this momentum experiment. Separating
an intrinsic charged excitation from its retained boundary reference needs
an additional localization and reference-independence argument.

## 3. Why two cells cannot test the desired mass shell

On the odd-source-parity sector, `T_joint^2=-I`. With the declared principal
branch, only momenta `p=+pi/(2a)` and `p=-pi/(2a)` occur. Consequently

`P^2 = pi^2/(4a^2) I`.

The code checks this identity on the physical compensated preparation; the
identity on the entire odd sector follows from the two Fourier projectors.
It gives only one squared momentum. For a shell

`E(p)^2 = M + C p^2`, `M=m^2 c^4`, `C=c^2`,

a single squared-momentum sample has design rank one. It cannot determine
both M and C. Assigning c separately permits a conditional value of M, but
one sample cannot test the dispersion law. No zero momentum is present.

For a declared even N-cell antiperiodic translation family, the odd-sector
menu is `k=(2j+1)pi/N` in the principal branch. The exact counts are:

| Cells | Distinct squared momenta | Rank for M and C | Independent shell constraints after fitting two parameters |
|---|---|---|---|
| 2 | 1 | 1 | 0 |
| 4 | 2 | 2 | 0 |
| 6 | 3 | 2 | 1 |

These are necessary menu conditions, not computed spectra of those larger
sources. Four cells permit a two-parameter fit; six give the first extra
axial shell constraint. A single particle band must be identified across
the momenta by its currents and actual nonzero readout residues. Arbitrary
levels from different sectors cannot be assembled into a dispersion curve.
Three-dimensional identification also requires three compatible spatial
translations and directional tests, not only an axial momentum menu.

`refinement.shell_feasibility` supplies an exact rational acceptance test.
For each supplied squared-energy interval `[l_j,u_j]` and squared momentum
`x_j`, it tests the simultaneous linear inequalities

`l_j <= M+C x_j <= u_j`, `M>=0`, `C>=0`.

It enumerates the exact feasible vertices. For any positive x_j the region
is bounded, so vertex emptiness is equivalent to inconsistency. With only
zero x_j, C is unconstrained and is not reported determined. Code controls
use explicitly algebraic test data, never values presented as outputs of
the quantum source. An incompatible third point is rejected. Consistency
is a necessary shell test, not a proof of relativistic covariance or a
physical particle identity.

A controlled finite-volume dispersion estimate is useful in its own right.
One may report a finite effective mass interval from the same H, state,
translation and clock, with its finite-resolution and boundary errors.
A full continuum construction or infinite-time scattering theorem is not
a prerequisite for that finite estimate. Calling it a measured particle
mass requires the additional physical identification evidence.

If an energy estimate E has error dE and a vector momentum estimate p has
Euclidean error dp, a fixed c gives the shell-error allowance

`2|E| dE+dE^2+c^2(2|p| dp+dp^2)`.

Numerical, finite-volume, boundary, cutoff and state-transfer errors belong
inside that allowance. Reassigning the clock independently for each species
or momentum would change the construction being tested.

## 4. A fixed primitive link is not a noiseless classical transporter

The actual e^c link channel is `U=|raised><blank|`. For every normalized
link state psi and complex u,

`||(U-u)psi||^2+||(U^dagger-conj(u))psi||^2 >= |u|^2`.

To prove it, rephase the two channel states so `u=r>=0`. After subtracting
`r^2 I`, the sum-of-squares operator has two eigenvalues `(r-1)^2` and
`(r+1)^2`; its four other active-link eigenvalues are r^2. All are
nonnegative. The code verifies this exact two-channel polynomial and its
embedding in the actual six-state link.

Thus no choice of state in this fixed primitive memory makes both a
nonzero transporter and its adjoint simultaneously classical with vanishing
relative residual. Multiplying the hopping by a scalar g preserves the
lower bound `|g u|^2`. Keeping a nonzero effective hop while rescaling the
coupling does not eliminate this product-background obstruction.

This is a boundary for a **product-state replacement by a c-number**, not
a no-go theorem for an interacting quantum gauge field. A directly
calculated quantum spectrum, an enlarged collective link memory, or a
proved entangled/dressed effective map has a different acceptance test.
The scalar memory's possible large-memory approximation does not by itself
supply a classical link state.

## 5. Reusing the finite theorem stack with quantitative premises

The comparison route uses
[EffectiveQuantumComparison.lean](../../../../Lean/QFT/EffectiveQuantumComparison.lean)
and the [owning analytic theorem](../../../../paper/tex_fragments/EFFECTIVE_QUANTUM_COMPARISON.tex).
For finite self-adjoint generators H and H', an isometry J and a reference
spectral window P commuting with H, supply the **full** residual bound

`||(H'J-JH)P|| <=epsilon_H`.

This is an operator bound on the invariant reference window, hence a bound
along the entire reference trajectory. A residual computed only on the
initial vector does not supply it. Section 1 separately proves the full-
domain bound 15/2 and uses a Taylor remainder for its initial-vector control.

Matrix exponentials make the interaction-picture derivative identity
unproblematic in finite dimension. The theorem gives a vector propagation
error `|t|epsilon_H/hbar`. Its effect theorem has constant one, and its
tail theorem retains the omitted vector. With initial-state error epsilon_0,
initial vector tail tau, detector error epsilon_D and numerical error
nu, the complete event-probability bound is

`min(1,epsilon_0+|t|epsilon_H/hbar+2tau+epsilon_D+nu)`.

If the analytic comparison uses a nonnegative reference generator, shift
both finite Hamiltonians by the same sufficiently large constant. Their
full residual and all energy differences are unchanged; the associated
common overall phase has no effect on the recorded probabilities.

For mixed states, tau may be replaced by the square root of the spectral
tail probability, using the existing analytic conditioning proof. For a
whole finite outcome distribution, epsilon_D must hold uniformly over
unions of its retained outcomes, or one must sum the individual detector
bounds. Failures, compensator positions and reference records belong to
that common outcome domain. The exact leakage in section 1 demonstrates
why a compressed generator is insufficient here.

For a refinement family, fix the source program at each N, its clock and
energy zero, the physical state, all apparatus references, and compatible
maps between common record meanings. These are declared construction
choices; they do not require a uniqueness proof. If adjacent complete
readout errors satisfy

`delta_N(T)=epsilon_0,N+T epsilon_H,N/hbar+2tau_N+epsilon_D,N+nu_N`,

then the distance between levels N and M is at most
`sum_(n=N)^(M-1) delta_n(T)`. Summability gives a Cauchy family with explicit
tail bound `sum_(n=N)^infinity delta_n(T)`. This proves convergence of those
specified finite-window records. It does not automatically construct a
Wightman field, Lorentz covariance or interacting scattering.

Longer spectroscopy times impose the additional joint condition
`T_N epsilon_H,N ->0`, together with the state, detector and tail bounds.
A fixed finite-time accuracy cannot resolve arbitrarily narrow spectral
lines. Boundaries must also recede or have separately controlled effects
on the selected localized source-and-reference observables.

The finite theorem attachments have distinct roles:

- [PeriodicMomentum.lean](../../../../Lean/QFT/PeriodicMomentum.lean)
  transports the supplied `HT=TH` identity to Fourier projectors and a
  declared momentum logarithm. The periodic joint source/reference
  construction supplies that premise. Its commuting-polynomials theorem
  also applies once three commuting physical translations are supplied.
- [MassShellKinematics.lean](../../../../Lean/Geometry/MassShellKinematics.lean)
  transports a declared Lorentz-module mass and momentum and characterizes
  rest at zero spatial momentum. It does not identify our finite excitation
  energies and translation phases with that module.
- The [chiral free regulator](../../chiral_regulator/PROOF.md#controlled-free-field-refinement)
  has exact windowed Weyl evolution and controlled smearing tails.
  Reusing them for the quantum-link/Higgs source requires the full residual,
  state and detector bounds above; a blank link is not its free vacuum.
- [source_mass_shell](../../../source_mass_shell/PROOF.md) proves
  quantified three-dimensional Dirac shells for its specified one-body
  mass family and exterior packets. Its own proof distinguishes those
  results from an interacting charged pole. It is an available comparison
  target, not an interchangeable Hamiltonian or reference for this source.
- [FiniteUnitaryScatteringNoGo.lean](../../../../Lean/QFT/FiniteUnitaryScatteringNoGo.lean)
  excludes taking a direct large-time limit of a nontrivial fixed finite
  unitary. Its relative-evolution result explicitly permits a comparison
  dynamics. An actual scattering claim needs that comparison or a controlled
  order of volume, cutoff and time limits; the theorem neither forbids nor
  supplies scattering for this source.

A nonzero finite coupling alone does not prove a nontrivial continuum.
For example, let the four common recorded functions X_i obey |X_i|<=1
and have mean exactly zero under both compared laws. Total-variation error
delta changes their fourth connected moment by at
most `14 delta`: the four-point moment contributes `2 delta`, and the
three products of two-point moments contribute at most `4 delta` each.
A connected signal larger than the certified refinement tail therefore
persists in the limiting records. To interpret it as an interacting field
signal additionally requires its actual linear-field and vacuum ancestry;
a composite observable or non-Gaussian input can be connected even in a
free theory.

## 6. Concrete next construction

The direct route keeps the quantum links and quantum scalar memories,
specifies a periodic finite source family and computes its own physical
reference, joint energy/translation records and charged-field residues.
The spectral regulator's path routing and boundary signs must translate
with the whole program. Choosing a different local regulator is allowed,
but its chirality, extra branches and continuum estimates need their own
proofs. The two-cell momentum algebra alone does not select that family.

After a sufficiently rich momentum menu is available, the interval shell
test, reference-independence estimates and same-source refinement bound
above supply concrete acceptance checks. A complementary comparison route
can use the existing free-field theorems through a proved dressed isometry
and controlled quantum residual. The current blank-link and noiseless-link
shortcuts fail by the explicit calculations above.

The certified finite energies remain valid outputs. A finite
same-source dispersion estimate can precede a full continuum proof; it
must retain its box, cutoff and error domain. The fixed-cutoff large-time
no-go does not forbid finite-time collision or spectroscopy experiments.

Neither route is completed by changing names of the open-source energy
lines. A physical mass requires a common energy/momentum convention, an
identified stable field branch and controlled reference, boundary and
cutoff errors. Converting that established energy to GeV is a subsequent
unit/reference attachment. The finite joint momentum apparatus and
certified source energy spectra provide separate inputs to this test;
they are not combined across Hamiltonians without the transfer proof.
