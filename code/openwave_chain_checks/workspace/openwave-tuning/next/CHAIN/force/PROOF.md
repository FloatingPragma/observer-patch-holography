# A finite force measurement from native current records

The source is the frozen reciprocal two-copy classical process on the
side-two cubic routing. It has 96 ports and 288 seams per copy, a fair
native coin, two mobile loads and the signed-charge readout
rho=delta_x-delta_y. The full source has 9,216 position states. This
consumer constructs no full position generator or semigroup. Its measured
observable is an oriented seam-current score on a complete first-event
instrument, with native holds and the attempt-cap flag retained.

The source manifest is `auditor/native_field_bridge/FINAL_MANIFEST.json`,
SHA256 `4597dc9c226d203d890bc3d3a929bbe1be80e2dc57ec03898c77da4fe5b954e1`.
Its finite local field controller defines the symmetric neutral Green
operator G_K with K=512 and the energy

    U_K(x,y) = (kappa/2) rho^T G_K rho.

The real arithmetic field and its controller energy are architectural
inputs. They are distinct from the native integer current tape and an
integer face-load minimizer. The exact twelve-port control in
`coordination/INTEGER_FIELD_CONTROL.json` has integer minimum energy 1/2
and real minimum 11/60. No equality between those domains is used here.

## Energy, field and native transition rates

Write B_e=delta_a-delta_b for an oriented seam a to b. A move of charge
sigma along that seam changes rho by -sigma B_e. Expanding the same
quadratic energy gives

    Delta U_K = kappa[-sigma B_e^T phi
                     + (B_e^T G_K B_e)/2],   phi=G_K rho.       (1)

The self-response term is part of a finite hop. Substituting only the
present field B_e^T phi gives an incorrect energy difference. Equivalently,
Delta U_K is minus the charge times the midpoint field of the two charge
states, multiplied by kappa. The executable checks (1) on the measured
source and retains the bare-field error as a rejecting control.

Each seam of either copy has scheduled attempt rate

    r = nu exp(-beta Delta U_K/2).

The literal native rule moves the load on exactly one of its two fair
coin branches. Therefore the moving transition rate is

    w(s to s') = (nu/2) exp(-beta Delta U_K/2).                  (2)

The opposite move changes the sign of Delta U_K. Equations (1)-(2) give
the exact finite detailed-balance ratio

    w(s to s') / w(s' to s) = exp(-beta Delta U_K).              (3)

The rate nu is one common scalar multiplying every seam rate of both
copies. It has units of inverse declared source time. Holds cost an
attempt, including the many seams not incident to either load. At zero
coupling the joint attempt rate is 576 nu, and the total load-transfer
rate is 6 nu. Source time, event count and controller arithmetic work
are different observables. No relation to seconds is supplied.

## Complete finite instrument and drift

Fix the positive-copy seam (0,1). Prepare x=0 or x=1 with probability
one half each, retaining the preparation label, and prepare the other
copy at y=9. The source evolves with its actual joint scheduler until
time T=1/1000 or the first scheduled native attempt. The latter outcome
has `cap_exhausted` status and retains the native record and terminal
positions. This is the cap-one instance of the frozen complete capped
instrument. It does not report a capped position as a horizon position.

For a preparation s, let Lambda_s be the total rate of all 576 scheduled
seams and put

    g_s = [1-exp(-Lambda_s T)]/Lambda_s.

The no-event outcome has probability exp(-Lambda_s T). Every copy/seam/
coin mark has probability g_s r_e(s)/2. There are 1,153 outcomes for
each preparation: one no-event outcome and 1,152 marked cap outcomes.
The literal native endpoint routine is called for every mark. Their
probabilities sum to one, and the archive retains both complete laws.

Let j_e be the native current if the attempted mark belongs to the
selected seam of the positive copy, and zero otherwise. The score
j_e/g_s is zero on the no-event branch. Its magnitude is bounded by
1/min_s g_s. It uses every cap output. Its unconditional expectation is
the instantaneous oriented seam drift of the declared endpoint mixture:

    J = [w(0,y to 1,y)-w(1,y to 0,y)]/2
      = (nu/2) sinh(beta F/2),
    F = U_K(0,y)-U_K(1,y).                                    (4)

Thus F is a finite lattice energy-drop force, with a declared unit seam
length. This measurement determines a current response, rather than an
acceleration or inertial mass. A physical force would require an energy
unit and a physical seam length. The expression is specific to this
instantaneous local preparation; it is not a persistent drift of a
uniform or equilibrium source on the whole graph.

At zero coupling the second-current rate is nu/2. Define the local seam
diffusion coefficient by half that rate, D_0=nu/4. The discrete Einstein
relation is then

    J = beta D_0 F + R,
    |R| <= nu |beta F/2|^3 cosh(|beta F/2|)/12.                 (5)

This follows from the Taylor remainder of sinh and holds at finite F.
The executable measures both moments from the native outcomes and checks
the bound at kappa=0,0.01,0.1,1,2. No temperature, kappa, nu or microscopic
electromagnetic coupling is inferred from the relation.

At beta=4, kappa=nu=1, the measured energy drop is
0.034922889610389585 and the drift is 0.03495129136610145 per source time.
The source retains a nonzero unit-field readout at kappa=0 while its
measured drift is zero. Thus the native Green field alone does not supply
this response; the installed scheduler is doing observable work.

## Interventions and retained failures

Reflecting the first coordinate of the internal port frame induces the
permutation (0 1)(2 3)(8 10)(9 11) in every cell. On the side-two routing
it is an exact graph automorphism and preserves the pair-energy table.
It exchanges the measured seam endpoints and maps y=9 to y=11. The
measured force and drift reverse. This is a reflection of the source
preparation on the actual graph.

A second intervention executes an actual nonzero native move of the
other copy. The reciprocal controller changes the positive-copy drift.
The held-field comparison keeps that copy's rates at the field of the
other copy's initial port. The held-field drift has zero cross-response.
Each comparison uses its own total-rate exposure, preventing a change
in waiting time from masquerading as a force response. The held-field
comparison has a different controller; it is not labelled as the frozen
reciprocal law. Doubling nu doubles the measured drift.

Two short executions call the frozen source executor at caps one and
two, retain their cap failures and replay all waits and native currents.
The exact force instrument itself retains all first-event probabilities.
Deleting cap outcomes, deleting holds, conditioning on conclusive outputs,
reversing currents, dropping the fair-coin factor, using the wrong nu,
or relabelling held-field output as reciprocal output is rejected by a
separate implementation. The same author wrote that verification, so it
is not an external scientific audit.

## Quantitative transfer from a controlled arithmetic field

The field-limit manifest is pinned at SHA256
`612aedd40f13faea2b80b5c31d4afbb19cb2332ff22db5b459ac2ae936879008`.
Its direct-pair zero-start 512-round 40-bit solver has field-norm and
unit-energy error bounded by

    eta = 2 theta^512 + (408/2^40)(1-theta^512), theta=23/24.

Define an explicit comparison scheduler that solves every candidate
neutral pair separately with that arithmetic and uses its energy in
(2). It is named `direct-pair-zero-start512-fixedpoint40-electric-scheduler`.
It is different from the frozen floating 512-round cache. No source file
or frozen scheduler is changed. Only the 22 required candidate pairs of
the first-event experiment are compiled; the following analytic bound
is uniform on all pair states.

The ideal frozen controller's unit-energy error from the exact Green
energy is at most epsilon=2 theta^512. Let delta=eta+epsilon. Every
state energy differs by at most kappa delta between the two controllers.
An edge energy difference has error at most 2 kappa delta, so

    |log(r_variant/r_frozen)| <= ell=beta kappa delta.          (6)

The frozen cap certificate bounds each incident seam rate by 8 nu on
its beta<=4, kappa<=2 domain. There are twelve incident seams across
the two copies; all other rates agree. Coupling equal marks and coins
until their first rate disagreement bounds the total disagreement
intensity by 96 nu (exp(ell)-1). Stopping both paths at the same finite
attempt cap preserves the coupling. Consequently their complete marked
capped path laws, including waits and all cap outcomes, obey

    TV <= min(1,96 nu T [exp(ell)-1])
       <= min(1,96 nu T ell/(1-ell))   when ell<1.              (7)

This argument retains holds and requires no full pair semigroup. It
also bounds the expectation difference of any bounded detector by its
range times TV. For the frozen exposure score in (4), the error is at
most 2 TV/min_s g_s. The executable compares both complete first-event
laws directly and checks these bounds. It evaluates the identical
frozen score on both laws, avoiding an estimator change in the bound.
The uniform attempt-rate envelope Lambda<=660 also gives the rational
lower bound g_s>=T-660 T^2/2=67/100000. Using the outward rational TV
bound, the separate verifier obtains a score-error bound below 2.003e-6
without a floating exponential in that certificate.

For beta=4, kappa=nu=1 and T=1/1000 the exact rational upper expression
in (7) is below 6.70801e-10. Controller fixed-point arithmetic is exact.
Equations (6)-(7) compare ideal exact-real exponentials and the ideal
frozen stencil. Floating-library evaluation error has no general
certificate here and is stated separately from the arithmetic bound.

## Source and physical scope

The controlled field theorem supplies a static smoothed Coulomb limit
on a unit three-torus, with explicit Fourier and port-readout bounds.
It supplies no dynamic trajectory or common-clock refinement across
carrier sizes. This force measurement uses its side-two error theorem
and does not transfer the static continuum normalization into a
continuum force law. The native mean accumulated-current Green identity
requires uniform scheduling; it is not used for this state-dependent
electric scheduler.

This is the literal classical reciprocal branch. It is not the C10
384-state quantum source, its Hamiltonian, its reservoir preparation or
its marked quantum clock. The construction supplies a controlled finite
dissipative force response with explicit field-to-rate error transfer.
It supplies no coherent quantum-source transfer, Maxwell waves, Lorentz
limit, mechanical inertia, isolated-system energy conservation or
derived fine-structure constant. The signs, kappa, beta, nu, energy
functional, field register and centralized controller are supplied
architecture. No acceptance configuration or official row is changed.

Run from this directory with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 force_measurement.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify.py
    python3 finalize.py

`RESULT.json` and `FIRST_EVENT_OUTCOMES.npz` contain the executed laws.
`FORCE_CONTRACT.json` records the typed handoff and `MANIFEST.json` binds
the complete package.
