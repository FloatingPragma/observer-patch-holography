# A faithful thermal state and the same finite observer clock

The declared 307-state Higgs source has a faithful Gibbs reference whose
modular flow is exactly a rescaling of its Hamiltonian evolution. Complete
scalar-number records recover the same interaction gap without requiring
ground-state preparation. The inverse temperature and preparation service
are architectural choices; the density, detector probabilities and clock
relation follow from those choices on the same source.

## The full source Gibbs state

Use the actual Hamiltonian and projectors from
[`higgs_records.py`](higgs_records.py), with \(Y=\operatorname{diag}(1,2,3)\)
and pair cost \(d=1\). On its entire \(Q_*=N_F-2N_H=0\) sector,

\[
 H=-7\Pi_-+8\Pi_++D,
 \qquad \operatorname{rank}(\Pi_-,\Pi_+,D)=(1,1,305).
\]

Choose \(\beta=\log2\) in the declared Hamiltonian units. The exact
partition function includes every dark state:

\[
 Z=2^7+2^{-8}+305\,2^{-1}=\frac{71809}{256},
 \qquad
 \rho_\beta=\frac{32768\Pi_-+\Pi_++128D}{71809}.
\]

These three coefficients are probabilities **per state** in their
respective eigenspaces. The dark block's total probability is
\(305\times128/71809\), rather than \(128/71809\). The density has
trace one and smallest eigenvalue \(1/71809>0\), so it is faithful on
all 307 states. It commutes with the same \(H\) used during readout.
No state outside this conserved sector is assigned a thermal weight by
this construction.

The [finite Gibbs theorem](../../Lean/QFT/StructuralInheritance.lean)
proves normalization, positivity and the algebraic KMS boundary identity
for a finite Hermitian Hamiltonian. In particular,
`gibbsState_one`, `gibbsState_star_mul_self_nonneg` and `gibbsState_kms`
apply to this finite matrix. This package supplies a concrete Hermitian
source matrix and a faithful density; it does not infer a thermal state
from a public classical record algebra.

## The modular parameter and the Hamiltonian parameter

Functional calculus on the three complete spectral blocks gives

\[
 K=-\log\rho_\beta=\beta H+\log Z\,I.
\]

Fix the conventions

\[
 \sigma_s(A)=\rho_\beta^{is}A\rho_\beta^{-is},\qquad
 \alpha_\tau(A)=e^{i\tau H}Ae^{-i\tau H}.
\]

The scalar partition factor cancels in conjugation, yielding
\(\sigma_s=\alpha_{-\beta s}\). Thus the exact clock conversion is
\(\tau=-\beta s\) with these conventions. The negative sign reflects
the stated definitions, rather than a physical reversal of an observer's
records. Replacing the modular convention by its inverse reverses the sign.

For a block \(\Pi_aA\Pi_b\), the modular frequency is
\(\log(p_a/p_b)=-\beta(E_a-E_b)\). The executable checks all nine
block pairs, accounting for all \(307^2=94249\) matrix units, including
the complete degenerate dark algebra. The comparison is therefore an
operator identity on the full finite field, rather than a fit to one
oscillation. It is consistent with the finite spectral comparison in
[`ModularClockComparison.lean`](../../Lean/QFT/ModularClockComparison.lean).

The thermal parameter does not select seconds, kelvins or GeV. Both the
Hamiltonian unit and the inverse temperature belong to the declared
architecture. A different faithful reference generally has a different
modular generator; the displayed identity uses the Gibbs reference of
this same Hamiltonian.

## A complete, gauge-invariant preparation with a retained reference

Let \(C=\sqrt{\rho_\beta}\) and prepare

\[
 |\Psi_\beta\rangle=\sum_{i,j=1}^{307} C_{ij}
 |i\rangle\otimes|\bar j\rangle.
\]

The second register carries the conjugate representation. The coefficient
matrix gives reduced states \(CC^\dagger=\rho_\beta\) and
\(C^{\mathsf T}\bar C=\rho_\beta^{\mathsf T}\); its Schmidt rank is
307. The reference is retained throughout the experiment. An apparatus
acting only on the field has the same probabilities as the density
\(\rho_\beta\), while the joint evolution retains the reference and
its correlations.

The source Ward identity makes \(H\), and hence \(C\), commute with
its gauge representation \(U\). The joint transformation
\(U\otimes\bar U\) acts on coefficients as \(C\mapsto UCU^\dagger\),
so the purification is a singlet. This distinction is essential: the
full 307-state Gibbs density is invariant under gauge transformations,
but it has support beyond the field's Gauss-singlet subspace. In this
preparation route the conjugate reference carries compensating charge.
It is an indispensable part of the physical Gauss interpretation, rather
than a discarded auxiliary label. Projecting the field onto singlets
would define a different density and different weights.

The preparation can be specified as a complete finite unitary. Let the
blank be the invariant pair \(|b\rangle=|\Omega\rangle|\bar\Omega\rangle\)
and \(a=\langle b|\Psi_\beta\rangle\), which is real and satisfies
\(0<a<1\). Define

\[
 |u\rangle=\frac{|\Psi_\beta\rangle-|b\rangle}{\sqrt{2(1-a)}},
 \qquad R=I-2|u\rangle\langle u|.
\]

Then \(R^\dagger R=I\) on all 94249 joint inputs and
\(R|b\rangle=|\Psi_\beta\rangle\). Both \(|b\rangle\) and
\(|\Psi_\beta\rangle\) are singlets, so the exact reflection is
gauge invariant. This is preparation from a specified blank, not a reset
of an arbitrary unknown input.

The exact reflection can be generated by the declared preparation
Hamiltonian \(\pi|u\rangle\langle u|\) for one preparation duration. This
generator commutes with the joint gauge action. The compiler route gives
an approximation to the complete logical operation; it does not assert
that every intermediate pulse separately preserves the field Gauss sector.

One allowed architecture declares \(R\) as a coherent primitive. Another
uses the inherited [finite instrument compiler](READOUT_CONTRACT.md) on
two nine-bit field/reference encodings, with identity on unused code inputs.
`preparation_plan` reuses its 18-spin complete-channel bound. It reports
the coherent service assumptions and retains every compiler failure label;
no large compiled tape is executed here. Preparations and detectors are
separate services with autonomous field evolution paused or refocused.
During observation the reference is untouched and the field evolves under
the same \(H\).

This construction does not claim that passive repair equilibrates to this
Gibbs state or prepares a Lorentz-invariant vacuum. No uniqueness theorem
for the architecture is required to calculate its chosen preparation.

## Read the gap from every original detector outcome

Use the complete scalar-number instrument
\(\{P_0=|\Omega\rangle\langle\Omega|,I-P_0\}\). For an arbitrary
input state, write
\(\alpha=\operatorname{Tr}(\rho P_0)\),
\(b=\operatorname{Tr}(\rho B)\) and
\(\delta=\operatorname{Tr}(\rho D)\), with
\(\alpha+b+\delta=1\). The two measurements, separated by the same
Hamiltonian duration \(t\), have the four unconditional joint probabilities

\[
 \begin{array}{c|cc}
 &\text{empty second}&\text{occupied second}\\\hline
 \text{empty first}&\alpha(1-T)&\alpha T\\
 \text{occupied first}&bT&b(1-T)+\delta
 \end{array},
 \qquad
 T(t)=\frac{4g_{\rm pair}^2}{\Delta^2}\sin^2\!\left(\frac{\Delta t}{2}\right).
\]

This is the actual full-sector instrument calculation from
[HIGGS_RECORDS.md](HIGGS_RECORDS.md), with \(g_{\rm pair}^2=56\) and
\(\Delta=15\). Initial coherences are retained in the quantum instrument;
they cancel from these four pointer probabilities. An equilibrium input
does not make the sequential measurements nondisturbing, and these joint
probabilities are distinct from an ordered quantum two-point correlator.

For the Gibbs input,

\[
 (\alpha,b,\delta)=\frac{(262151,229384,585600)}{1077135}.
\]

At \(t=\pi/45\), \(T(t)=56/225\), and at \(2t\),
\(T(2t)=56/75\). The code exports all four probabilities at both
durations, including the entire dark contribution. The original joint
event \(J_{01}(t)=\alpha T(t)\) therefore gives

\[
 \frac{J_{01}(2t)}{4J_{01}(t)}=
 \cos^2\!\left(\frac{\Delta t}{2}\right)=\frac34,
 \qquad
 \Delta=\frac{2}{t}\arccos\sqrt{\frac34}=15.
\]

The cancellation works for every reproducible input with \(\alpha>0\)
and nonzero coupling, on the declared unaliased window
\(0<\Delta t<\pi\). Ground-state preparation is unnecessary. The same
input preparation must be used for both durations; changing \(\alpha\)
changes the ratio. A dark-only input has \(\alpha=0\) and supplies no
such gap readout. No experiment is conditioned on seeing the first empty
outcome: both runs retain all four original histories.

The [common macro-clock contract](READOUT_CONTRACT.md) uses one fixed
implementation of \(U(\pi/45)\), repeated twice for the second run.
Preparation, tick and detector errors compose on the complete record
instrument, with arbitrary untouched references. If each service has
complete error at most \(10^{-6}\), the two-tick, two-detector run has
error at most \(5\times10^{-6}\), including its one preparation.
Sampling and implementation uncertainties enter both probability
intervals before calling the owning `gap_interval_from_records` inverse.
That inverse rejects a denominator interval touching zero or an unresolved
alias window.

## Reproduction and scope

Run `python3 -m model.particle_observer_outputs.thermal_clock` for the
exact receipt and compiler bounds. Run
`python3 -m pytest -q model/particle_observer_outputs/test_thermal_clock.py`
for the full-sector density, modular sign, complete purification, direct
307-state Kraus comparison and preparation-independence checks.

The result connects a faithful reference, modular parameter, Hamiltonian
evolution and complete detector records on one declared interacting finite
source. It does not identify its gap with a Standard Model particle mass,
establish an interacting spatial continuum or assign physical units.
