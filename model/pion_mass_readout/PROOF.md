# Pion operators and conditional energy calibration

The source is the archived quenched Wilson calculation in
`model/quenched_hadron_binding`. Its production lattice is 12³ × 24 at
β = 5.7, with 128 configurations in 32 streams. The package replays reduced
correlator tables; it does not regenerate the absent raw gauge configurations.

## Isospin triplet

Write a colour-singlet pseudoscalar quark-antiquark operator as a flavour
matrix M, with row a and column b labelling the ket content q_a anti-q_b.
The normalized triplet is M+ = |u><d|, M− = |d><u| and
M0 = (|u><u| − |d><d|)/√2. All three have unit Hilbert-Schmidt norm and zero
trace. For Q = diag(2/3, −1/3), [Q,M±] = ±M± and [Q,M0] = 0.
The adjoint action of the isospin generators closes on this triplet.
The matrices describe leading valence operators, not a complete physical
Fock-state decomposition of a pion.

With identical u and d Wilson propagators and no electromagnetic field, the
connected flavour contraction is Tr(M_a† M_b) C(t) = δ_ab C(t).
The neutral disconnected contraction cancels because Tr(M0) = 0. Thus the
same fitted pseudoscalar decay energy serves all three charge labels. A
nondegenerate light-quark operator or electromagnetism invalidates this
equality. The tests include a nondegenerate-propagator counterexample to the
disconnected cancellation. No charged-neutral splitting is computed here.

## Calibration ancestry

The retained fit chain sets the light-quark point by the imported ratio
c = mπ/mρ = 0.1780. Let r and n denote the fitted lattice rho and nucleon
energies at that point. For a declared rho anchor R = 0.77526 GeV,
the conversion is a⁻¹ = R/r, mπ = cR = 0.13799628 GeV and mN = R(n/r).
The pion and rho values in this conversion are calibration inputs. Their
conditional fit variance vanishes by exact shared-scale cancellation;
this statement supplies no total physical uncertainty or independent mass
prediction. The nucleon ratio and its jackknife error come from paired
delete-one-stream replicas, not independent numerator and denominator errors.

An alternative conventional string anchor S = 0.44 GeV gives
a⁻¹ = S/√(a²σ), mρ = S r/√(a²σ), mπ = c mρ and
mN = S n/√(a²σ). Both conversions are retained. Their disagreement is
exposed rather than absorbed into a selected scale. Fit sensitivities are
finite envelopes over retained fit choices; they do not bound quenching,
finite-volume or continuum effects. The conventional string anchor has no
assigned Gaussian error here.

The complete four-kappa mass table is retained, including the unresolved
lightest nucleon estimate. Experimental charged and neutral pion references
are stored separately and do not enter either conversion. A physical
identification requires quark loops, volume and continuum control,
electromagnetic and light-quark splitting, and a transfer to the native
observer source and clock. A correlator fit is not a complete native
observer instrument. No such instrument or transfer is asserted here.
