# Magnetic refinement and finite-error precession on one source

## Fixed-volume and independent controls

The source is the positive kernel and transverse field energy in
[the stationary construction](../openwave_intrinsic_magnetic/PROOF.md).
The mass is 8, scalar mass-well coefficient 0.01 and Wilson coefficient
1/2. The spacing sequence has n=5,9,13,17, a=6/(n+1) and field side
N=7(n+1)/2. Both the matter Dirichlet extent 6 and periodic field extent
21 are therefore fixed exactly. The self-field coupling is 0.1. Every
stored state/field witness is checked against the reconstructed full
Dirac operator, actual Peierls current and Maxwell equations. The
separate regeneration command solves those sparse eigenproblems again.
The witness is a finite numerical certificate, without a claim that its
small floating residual alone identifies a global minimizer.

The bound-energy static g values are 1.47516855, 1.78874943, 1.88532471
and 1.92737696. Independent fixed-spacing matter-volume and field-volume
sequences distinguish those effects. Symmetric probes at B=0.02,0.01,
0.005 compare the total-energy derivative with the same zero-probe
moment. The first three spacings' quadratic extrapolation predicts the
fourth within 0.00180. Zero-spacing intercepts range from about 2.00021
to 2.07456 with fit order. These are diagnostics, not a selected
continuum value or a physical uncertainty interval. Matter-volume
changes are material.

## A refined dynamical member with explicit accuracy control

For a controlled precession experiment use the same n=17, a=1/3,
N=63 matter/field construction with kappa=10^-14 and external B=1/100.
This small supplied coupling controls self-field backreaction errors;
it is not an empirical fine-structure constant. It changes neither the
native Clifford matrices nor the derived magnetic-current coefficient.
The mean-field dynamics and its arbitrary-entangled-input limitation
are those of [the Larmor construction](../openwave_larmor_dynamics/PROOF.md).

Two numerical eigenlines are computed from the full 19,652-row sparse
Dirac matrix at zero self-field and the stated external probe. The
initial state is their normalized equal superposition, with relative
phase fixed by a spin-x readout. The full kernel retains every other
spectral mode, both signs and all four identical native blocks.
Specifying an input superposition is not a truncation of that kernel.
The same weak-coupling member also has an explicit zero-probe stationary-
residual witness: select positive spin in its actual zero-field ground
doublet and compute that state's self-field. Both the field fixed-point
defect and Maxwell residual are divided by kappa before comparison, so
setting the self-field to zero fails even at this small coupling. These
are numerical source-equation residuals, without an interval proof of a
particular exact nonlinear stationary eigenvector.

The link phases have rational angles. Their degree-18 exponential Taylor
polynomials and the binary numerical vector components are treated as
exact rational numbers. For |theta|<1/2, the omitted phase is bounded by
2|theta|^19/19!. Exact rational sparse matrix multiplication, a coordinate
norm bound and the link tail give an upper bound r_j on
||D v_j-E_j v_j||. The scalar inequality
(|lambda|-|E_j|)^2 <= (lambda-E_j)^2 then gives the same bound for the
full positive kernel. No positive-spectrum-only operator replaces |D|.
The exact binary squared norms exceed 1-10^-12. The residual bound is
divided by 1-10^-12 before outward rounding, so it bounds both the raw
vector residual and the normalized residual required for a spectral
root enclosure.

The scalar Wilson mass operator at zero probe separates into three
tridiagonal axes. Exact rational Sturm inertia counts enclose its first
two axis eigenvalues. For every vector,

    Re <beta v,Dv> = <v,Mv>,   ||Dv|| >= <v,Mv>/||v||.

The singular-value min-max principle therefore bounds the fifth
absolute Dirac eigenvalue from below by the second spatial eigenvalue of M,
including a norm bound for the external Peierls perturbation. At the
executed probe this lower bound is 8.162073472479, above both computed
low energies. Their exact residual intervals are disjoint. Chiral
symmetry supplies their negative-energy partners, accounting for all
four levels below that threshold. At zero probe, the exact antiunitary
symmetry with square -1 and chiral symmetry give a fourfold absolute-
energy quartet. The single isolated zero-probe root thus identifies the
bottom energy, not just an unspecified eigenvalue. Four identical native
blocks preserve this statement on all sixteen internal rows.

## Uniform self-field and dynamical error bounds

Let M bound ||D|| and let L bound the norm of the derivative of |D| with
respect to the complete edge potential. The same positive Sylvester
argument gives L <= (M/m)sqrt(2(1+r^2)); the executable uses a rational
upper bound. With transverse Laplacian lower bound lambda_1,

    ||A[psi]|| <= kappa L/(a^3 lambda_1),
    ||K(A[psi])-K(0)|| <= delta = kappa L^2/(a^3 lambda_1).

These bounds hold for every normalized state. The rational Hessian
bound makes the field map a contraction at the declared small coupling,
so its minimizer is unique. Duhamel's formula compares the nonlinear
quasistatic trajectory with the fixed unitary without an exponential
Gronwall factor:

    ||psi_qs(t)-exp(-itK(0))psi(0)|| <= |t| delta.

The normalized sum of the two witness vectors has norm denominator
larger than one. Hermiticity gives an exact bound on their overlap,
(r_1+r_2)/|E_1-E_2|. Consequently the numerical two-frequency preparation
has total state error at time T bounded by

    epsilon(T) <= T(r_1+r_2+delta).

This includes every unprepared mode. The executable evaluates the bound
as an exact rational number at T=6000. A small residual is not promoted
to a proof that a two-dimensional effective law was inserted into the
source.

The same uniform perturbation bounds the true zero-probe one-particle
minimum E_star between E_ground(0)-delta and E_ground(0). The zero-field
root enclosure gives exact rational endpoints for that rest-energy
reference. Supplied mass 8 and this bound-state energy are separate
normalizations.

## Complete spin and dipole observations

All canonical spin measurements retain both outcomes. Time samples
describe separate preparations evolved to the stated time;
no intermediate measurement is omitted from the evolution of one copy.
The exact rational
matrix elements of S_x+i S_y in the witness pair bound the transverse
spin radius uniformly above 0.49. Counterrotation and constant terms are
bounded separately. The uniform state error gives an unwrapped spin
phase error, hence an exact finite-duration interval for

    g_T = 2 E_star |arg<S_+(T)>-arg<S_+(0)>|/(B T).

The rational endpoints and outward decimal endpoints are printed. This
is a finite-duration precession normalization on the declared mean-field
source. It differs from the static convention dividing moment by the
measured canonical spin, whose length is slightly below one half.
Neither convention is silently identified with a continuum spatial
rotation generator or with a physical electron.

The bond matrices and their moment are obtained from the full native
current derivative evaluated on the prepared eigenlines. The small
negative-sign leakage is bounded using r_j/(E_j+m); it is not deleted.
Current errors include state error, field Hessian response and this
signed-branch correction. Complete weak dipole instruments have effects
(I+/-M_j/L_mu)/2, with L_mu an explicit full-operator norm bound. Their
probability errors are bounded. The intrinsic current produces its own
instantaneous transverse Green field, which follows the dipole. There
is no claim of a propagating photon field or a linear completely
positive state-dependent evolution on arbitrary entangled inputs.

The finite g is approximately 1.92095. The mathematical enclosure is
for that finite construction, with its supplied mass, trap, regulator,
weak coupling and probe. It is not an enclosure of the continuum or
observed electron g, and does not erase the finite offset from two.
