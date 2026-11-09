# Phase records and the charged-response premises

The configured coherent source has an exact finite action spectrum. Its
phase is readable from the program specification. A complete eight-record
experiment also recovers the unordered spectrum from an original outcome
probability, with a certified interval error. Neither operation selects a
phase from the measured primitive repair law. The phase-zero default has
two equal action eigenvalues and therefore cannot describe the three
unequal charged-lepton masses through one common proportional mass map.

A valid theorem whose hypotheses have not been attached to a specified
physical source remains a valid theorem. The distinctions below concern
which source, state, operator, clock and observations instantiate its
hypotheses. They do not require every admissible architecture to select the
same source, and they do not prohibit calculations for configured inputs.

## The specified finite source

Use the literal native C/E record partition and family three-cycle in
[`koide_source.py`](../standard_model_source_map/koide_source.py). The 128
fine outcomes include 104 zero-weight outcomes. Their five exclusive group
weights are `(10,5,5,2,50)/72`. The actual product of the C/E half-turns is
conjugate to the cyclic permutation P on the three-dimensional family band.
For a supplied unit phase z=exp(i theta), define

    a=sqrt(10/72), b=sqrt(5/72),
    C=a I+b(z P+conjugate(z) P^2), Y=C^2/3.           (1)

The supplied coherent preparation retains eight flags: the five groups and
three padding outcomes. A Fourier analyzer acts on the three active flags.
Its flag-zero Kraus operator is C/sqrt(3). The complete active probability
on the tracial family state is 5/18; no division by that probability occurs
in (1). In particular, summing all active effects yields `(5/18)I` and
would erase the nontrivial action.

The [coherent source proof](../standard_model_source_map/KOIDE_SOURCE.md)
checks the full native analyzer and finite weak-round columns. The
calculation here uses its exact matrices; it does not execute the native
pulse word again. All preparation and detector services stated below are
finite logical services. Their approximation and clock costs are separate
from the exact matrix identities.

## Which records identify a phase

Before coherent flag mixing, each fine map has the form

    K_r(z)=sqrt(w_r) z^(e_r) R_r,
    e_r in {0,1,-1}, R_r in {I,P,P^2}.               (2)

For every input, including arbitrary reference entanglement,
`K_r(z) rho K_r(z)^* = w_r R_r rho R_r^*`. The executable checks all 128
maps on all nine family matrix units. Therefore every finite history built
from these measured maps and fixed, classically record-controlled
interventions is identical for every z. This follows by induction on the
complete conditional maps, so no independence or tracial-input assumption
is needed. Zero records and the full target output are retained.

This identity has a precise access boundary. A coherent operation mixing
flag amplitudes is outside the measured-map interface of (2). Changing
source flag phases while holding that analyzer fixed changes an experiment.
Transporting both the dilation and the analyzer by the same flag-basis
unitary is a coordinate change and preserves its probabilities. These are
different comparisons. The existing fresh-record proof in
`oph-models/coherent_law_selection_20261025/reviews/FRESH_RECORD_SELECTION.md`
proves this distinction for actual native repair maps as well.

Even the complete single Fourier-analyzer law on I3/3 is phase-blind:

    (5/54,5/54,5/54,1/36,25/36,0,0,0).              (3)

The exact programs z=1 and z=(99+20i)/101 have the same law (3) and the same
fine maps (2), but different Tr(Y^2). Both phases satisfy positive-root balance and the generic same-current
Higgs coefficient identities. Thus the available data (2)-(3) do not
identify their action spectra. This is an exact counterfamily to that
specified inference, not to identification with other coherent experiments.
For example, the pure P-invariant family input gives flag-zero probability
`(a+2b cos(theta))^2/3`, which distinguishes these two programs. Choosing
a known program phase does not require inferring it from any record.

## Original-probability spectral readout

Prepare the tracial Nambu/family input I2/2 tensor I3/3 and apply one exact
weak round from the same coherent source. The half-angle tangent u>0 is a
known clock/control parameter, h=2 atan(u), and

    kappa=4 sin^2(h/2)=4u^2/(1+u^2).

Retain all eight reset flags. With P_U the analyzer's conjugated flag-zero
projection, the round is exp(-i h X tensor P_U), followed by complete flag
readout. Its nonzero-flag effect on the family is
`kappa(Y-Y^2)`. The original event probability p therefore satisfies

    Tr(Y)=5/18,
    Tr(Y^2)=25[17+4sqrt(2) cos(3theta)]/7776,
    p=kappa[1735-100sqrt(2) cos(3theta)]/23328.        (4)

The second identity is checked as a Laurent-polynomial identity for every
unit z. The program independently evaluates all eight effects and checks
that their sum is I and their probabilities sum to one. The controls z=1,
z=(99+20i)/101 and its conjugate retain every flag. The two nonconjugate
controls have different p; conjugation preserves the unordered spectrum.
No successful branch is normalized and no rare outcome is discarded.

Equation (4) has the exact inverse

    c3=cos(3theta)=[1735-23328 p/kappa]/[100sqrt(2)].  (5)

Thus a single original probability supplies the invariant needed for the
unordered balanced spectrum. A supplied probability interval [p-,p+] and
an additional certified event-probability error epsilon give the interval
obtained from (5) on
`[max(0,p--epsilon),min(1,p++epsilon)]`. The implementation uses outward
rational bounds for sqrt(2) and intersects with [-1,1]. An empty intersection
refuses the source model. The inverse sensitivity is

    |Delta c3| <= 23328 |Delta p|/[100sqrt(2) kappa].  (6)

The interval function does not certify sampling coverage, preparation
accuracy or native arithmetic error. Those errors must be supplied in the
same probability convention. In particular, a reduced output half-trace
error bounds an event error directly; an unhalved trace error contributes
at most half its value. An operator-norm isometry error epsilon contributes
at most epsilon for a complete event probability on normalized inputs.
There is no justification for substituting a probability-normalization
residual for an error bound on p.

For C to be nonnegative, the phase lies in
`|theta|<=pi/12 modulo 2pi/3`; equivalently c3 is in `[1/sqrt(2),1]`.
Write r for the eigenvalues of C/a. They are precisely the three real roots

    r^3-3r^2+(3/2)r+1/2-c3/sqrt(2)=0.               (7)

The action eigenvalues are `(5/108)r_i^2`. For an interval strictly inside
`1/sqrt(2)<c3<1`, the implementation isolates the three simple positive
roots by rational Sturm intervals and propagates them to the two ratios
relative to the smallest action eigenvalue. The lower and upper endpoints
use the derivative signs of the three cubic branches, rather than a
point-valued numerical eigenvalue fit. The root/ratio routine refuses a
zero-root boundary or unresolved degeneracy. Such refusal is a numerical
conditioning boundary, not a failure of (4)-(7).

At the known default z=1, the exact roots are
`1-1/sqrt(2),1-1/sqrt(2),1+sqrt(2)`. The action ratios are

    1 : 1 : (34+24sqrt(2)).                          (8)

The double eigenvalue cannot be removed by a family basis change or common
mass unit. The charged-lepton target has three unequal masses. Equation
(8) therefore excludes the default as that target under the stated common
proportional action-mass interpretation. It does not exclude the continuum
of configured phases in (1), and the certificate uses no empirical mass.
At theta and -theta the unordered spectrum agrees; a chirality/orientation
measurement is needed only to resolve that signed phase convention.

The clock in (4) is the supplied finite weak-kick angle. A sequence of such
rounds has the controlled logical generator X tensor Y in the clock t=Nh,
as proved in KOIDE_SOURCE. The fixed analyzer/uncompute pulse overhead does
not vanish with h. This calculation does not equate t with global native
attempt time, construct a physical vacuum or produce relativistic poles.

## Geometric connections and phase selection

The finite holonomy results supply useful positive information with an
explicit connection premise.

* `Lean/Screen/HolonomyInterference.lean`,
  `relativeCharacterPhase_eq_closedLoopPhase`, takes edge labels, reversal
  compatibility, a unitary character and paths as inputs. It identifies the
  relative phase with their loop holonomy; it does not choose those labels.
* `Lean/Screen/SeamU1HolonomyClassification.lean`,
  `gaugeEquiv_iff_chordHolonomy_eq` and `chordHolonomy_surjective`, prove
  that the actual 30-seam graph has a nineteen-circle space of gauge classes.
  Its unique tree-gauge representative is unique for each supplied class.
  Gauge fixing does not choose one of the nineteen holonomies.
* `common_source/native_mass_holonomy.geometry()` constructs a particular
  spin-half coherent-line connection from the twelve native directions.
  All twenty oriented triangular Bargmann products have phase pi/10. Its
  exact unimodular face-boundary certificate covers every integral closed
  edge chain, whose phase is an integer multiple of pi/10. The spin-one
  symmetric-square lift doubles that phase. This is a derived phase for
  this geometric connection, not an arbitrary fitted angle.

Attaching a chosen closed-loop phase to an equal-link family circulant adds
`3theta=loop phase mod 2pi`. With the separate balance and positive mass-map
premises, that connection gives unordered positive phase classes
`0,pi/30,pi/15`. The existing certificate evaluates their conditional action
ratios. It does not identify which loop or spin lift controls the charged
response, or identify its connection with the phase control in (1).
Choosing the spin-half single-face loop gives pi/30 as a valid configured
geometric alternative. It does not follow from raw duty probabilities.
The much larger space of continuously supplied seam connections is not
excluded by this discrete coherent-line example.

The local search covered the relevant Lean, model and charged-lepton code,
and the sibling `leptons`, `pro`, `closure`, `oph-models`, `proof` and
correspondence results. The
strongest source-law rigidity submission is represented by
`code/particles/leptons/derive_charged_source_law_rigidity_conditional.py`:
its finite rank, trace, clock and path assumptions fix eight affine
coefficients. Its baseline coordinate `zeta=theta-2/9` assumes that phase;
the coefficient theorem does not derive its baseline. The explicit phase
mechanism program under `pro/Particles/04/lanes/L2_phase_mechanism` requires
one radian per selected event to turn an event frequency into the number
2/9 as an angle. That angular unit and the event-to-coherent-control map
are separate inputs. None of these candidate attachments supplies an
identity selecting the phase of (1) from the specified raw repair process.
This is a scoped conclusion about the examined source interfaces and
results, not a theorem excluding every possible phase-selection dynamics.

## Premise crosswalk for the same objects

| Result | Exact data consumed | Attachment for the finite source | What that result does not imply |
|---|---|---|---|
| `KoideCirculant.koide_quotient_phase_independent` | Real a,b with a nonzero; three cosine coordinates at 2pi/3 | Equation (1), b/a=1/sqrt(2); positivity makes C the positive root of 3Y | A phase value, a vacuum, or pole masses |
| Finite tracial-GNS orientation isometry in `derive_charged_koide_orientation_isometry.py` | M6 trace state, two rank-two admitted sectors, canonical square-root map and declared orientation resolution | Its admission is 2/3; its normalized weights give the balanced modulus. The coherent source (1) has a separate complete construction with admission 5/18 | Identification of these two instruments, or a map of GNS roots to physical masses |
| `SourcePhaseSelection.generated_effect_values_exact` | Declared two-dimensional representation, complexification adapter and ordered commutators | Produces the unordered Pauli-Y pair under that adapter | A source-forced first orientation or the C3 phase theta |
| `OperationalBornBridge.operational_record_probabilities_unique_born` | Complete weights/effects, effect-equivalent probabilities and universal coexistence refinements | The finite matrix services here carry explicit complete Kraus effects and their specified Born law | Verification of its universal operational premises from a finite sample of native classical records |
| `ColimitSelectedStateGNS` | A consensus tower and a compatible selected observer/state family | A fixed finite matrix state has a GNS representation; the specified trace input above is explicit | Selection of a physical vacuum or compatible growing spatial family |
| Native P90/Higgs/current attachment | Same fifteen internal matrices, separate odd CAR modes, specified kinetic operator and complex Yukawa matrix | Y from (1) can be inserted into that verified coefficient map; all generation factors are retained | That the primitive repair dynamics supplies this Hamiltonian, or that the logical weak kick is its full-CAR implementation |
| Effective quantum comparison and weak feedback | A same-clock residual bound on the whole retained state, including leakage/reference/output error | The logical weak-kick theorem supplies its finite norm budget; (4) supplies a complete original-probability readout | A source-clock identification, a momentum shell, or a field limit without their own hypotheses |
| Free Weyl/Lorentz and gauge identities | Specified spacetime translations, kinetic scaling, state and currents; separate operator-valued link law for the interacting case | The current-preserving kinetic and finite Gauss identities are exact in their declared architectures | A joint interacting vacuum, mirror removal or physical charged-lepton pole assignment |

The tracial input used in (4) has a trivial modular flow: its modular
Hamiltonian is a scalar. It therefore cannot be identified with the
nontrivial X tensor Y weak-kick clock by naming both a quantum clock.
`Lean/QFT/ModularClockComparison.lean` requires an actual frequency
intertwiner and forces a nonresonant coefficient to vanish. A physical
attachment must preserve the selected dynamics rather than silently
substitute this tracial modular flow.

If a physical charged-lepton mass map supplies positive masses obeying
Koide, two independently supplied masses determine the heavy third root:
with x=sqrt(me), y=sqrt(mmu),

    sqrt(mtau)=2(x+y)+sqrt(3(x^2+4xy+y^2)).           (9)

No separate phase premise is needed for that conditional comparison.
A calibrated sanity test may instead supply one independently measured ratio
to choose the readable phase and compare the other ratio as held-out data.
That is a conditional test of balance and the chosen mass map; the input
ratio and phase calibration must be reported and cannot be labeled a
source-only prediction. The default program in (8) is unchanged.
Equations (4)-(8) instead compute the action spectrum of a known coherent
source without mass inputs. Identifying this spectrum with the same
physical lepton poles requires the actual spatial state, clock and
interaction/renormalization attachment; the valid finite Koide identity
alone does not supply that identification.
