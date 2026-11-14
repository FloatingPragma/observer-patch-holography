# Finite projections for the full-strength compact N16 ground energy

This theorem concerns the autonomous N16 source df34a04c at h=1/8,
m=1, a=1/2 and b=3 in its kappa-zero compact harmonic sector. All
190 coexact modes, the three harmonic momenta, both spin bands, the
internal Coulomb term and the occupation-dependent twist are retained.
The common coexact vacuum scalar is subtracted. This is a finite-energy
projection theorem, not an executed native ground preparation.

Write

    H_lambda = F + D + K + lambda W,   W=n_0-2,   |lambda|<=4.

Here F is the nonnegative normal-ordered free field Hamiltonian in
physical time units. D is the diagonal mass plus Coulomb term. The
source proofs give D>=-16, ||K||<=384 and ||W||<=2. The compact
fiber has compact resolvent. The interacting ground energy E_0<-50
follows from the explicit two-dimensional variational trial, without
using a free-reference compression.

Cell translations preserve this fiber. Averaging the ground projector
over translations gives mean W zero. Its variational energy under
H_lambda is E_0, hence E_lambda<=E_0<-50 for every |lambda|<=4. Ground
degeneracy is allowed. For any normalized ground vector psi_lambda,

    <F> <= E_lambda + 16 + 384 + 2|lambda| < 358.

Let P_Lambda be the exact spectral projection of F onto [0,Lambda],
restricted to the N16 Gauss/twist domain. F is diagonal in the exact
normal-mode/Fourier basis, and D and W commute with P_Lambda. There are
finitely many modes and all coexact frequencies are strictly positive;
P_Lambda has finite rank. The harmonic momentum labels retain their
occupation-dependent parity. No change of harmonic sector is made.

For Lambda>358 put delta=||(1-P_Lambda)psi_lambda||^2. Markov's inequality
gives delta<=358/Lambda<1. The eigenvalue equation and the commutation
of the diagonal terms with P_Lambda give the exact form identity

    <P psi,(H_lambda-E_lambda)P psi> = -<P psi,K(1-P)psi>.

Consequently the lowest Ritz value E_lambda,Lambda of the **compressed
physical operator** P_Lambda H_lambda P_Lambda satisfies

    0 <= E_lambda,Lambda-E_lambda
       <= 384 sqrt(delta/(1-delta))
       <= 384 sqrt(358/(Lambda-358)).

The lower bound is the variational principle; the upper bound uses a
single exact ground vector only to prove existence of a finite trial.
It does not assume a spectral gap, unique ground or overlap with the
free filled state. A sufficient rational cutoff for error epsilon is

    Lambda = 358 [1 + (384/epsilon)^2].

The projected Hamiltonian is the same source at every lambda; only the
diagonal external probe changes. The actual neutral external cluster
q_star=-BB^T delta_root has H(q_star)=H_4+12. Its bare self-energy 12 is
separate throughout.

## Size and computability

The pinned exact curvature certificate implies every physical coexact
frequency is greater than one. Thus each of the 190 occupation numbers
in P_Lambda is at most floor(Lambda). The physical harmonic energy is
|p_phi|^2/6 at h=1/8 and L=2, so each integer Fourier coordinate has
absolute value at most sqrt(6 Lambda). With Q=ceil(sqrt(6 Lambda)),

    rank(P_Lambda) <= binom(32,16) [floor(Lambda)+1]^190 (2Q+1)^3.

This deliberately loose finite bound does not use a volume-independent
gap. It counts extra states that fail the total-energy or parity test;
the exact physical basis is obtained by imposing both tests. This bound
is enormous even for a modest energy error and is not a feasibility claim.

All entries of the finite Hamiltonian are computable from the declared
source. The diagonal mass/Coulomb and harmonic energies are rational.
The exact coexact frequencies and frame are algebraic spectral data of
the finite physical curvature operator. The matrix of a Weyl displacement
in that basis has products of one-mode entries

    exp(-gamma^2/2) sqrt(min(m,n)!/max(m,n)!)
      (i gamma)^|m-n| L_min(m,n)^|m-n|(gamma^2).

The integer harmonic shift and exact CAR signs select the other entries.
Certified algebraic and exponential enclosures therefore provide a
finite matrix oracle to arbitrary error. An interval matrix with entry
error eta/rank has operator error at most eta. A finite eigensolve with
residual/Ritz bounds can certify the lowest energy to any requested
additive precision. Degenerate eigenvalues need not be separated.

A native preparation must implement an approximate low-energy vector
of this finite matrix and its energy reader, with coefficient, reset,
export and pointer errors controlled. The universal finite-code native
services make this a concrete synthesis problem; this packet does not
claim that its very large program or matrix has been generated or run.
The state-dependent short-time bridge 1f510cbd begins from a different
preparation and is not itself that preparation theorem.

## Projection and truncated displacement

The projected operator uses P exp(i J.q) P. In general this is **not**
exp(i P J.q P). The latter is unitary on the truncated coordinate code
and is useful in the bounded short-time evolution comparison,
but it has no automatic variational ordering relative to the physical
ground energy. Thus its lowest eigenvalue must not be called an upper
bound for the physical ground energy.

An exact two-dimensional control makes the distinction explicit. Put
F=diag(0,1), q=(pi/2)X, K=-(exp(iq)+exp(-iq))/2=0, and P=|0><0|.
The physical and Ritz ground energies are zero, whereas replacing the
compressed displacement by exp(iPqP)=1 gives energy -1. This is an
algebraic negative control, not a simulation of the 190-mode source.

## Consequence for a screening experiment

The external-screening packet supplies, on L2,

    E_4-E_0 <= -1252815/135895533121 < -9.2189e-6.

With both finite Ritz errors at most epsilon,

    E_4,Lambda-E_0,Lambda <= -1252815/135895533121 + epsilon.

Only one epsilon is needed: E_0,Lambda>=E_0. If a prepared state in each
finite subspace has energy at most eta above its respective Ritz minimum,
the same difference of prepared energies is at most the right side plus
eta. A complete native energy measurement adds its own certified error
interval and failure probability; neither may be discarded by conditioning
on a successful readout. This gives explicit tolerances for a physical
integer-charge screening experiment without assuming an interacting gap.

It does not establish a momentum-dependent coupling, scale onset, an
infinite-volume state, a charged particle dispersion or any OpenWave row.
Those claims require separate evidence on the same physical source.

The accompanying exact-arithmetic driver consumes the frozen sibling-workspace sources named in its receipt. It is not a self-contained public export of those dependencies.
