# A native first-hit controller creates an attractive internal-band well

This construction replaces a supplied numerical atomic weight by a native
record-dependent scheduler. It is an explicit new N+ architecture, not a
claim that OPH uniquely forces this controller or that a charged nucleus
has already been identified. No primitive oscillator or potential phase is
inserted.

## 1. The inner probe is a complete native record experiment

The inner source is one native unit on a three-dimensional box of the same
twelve-port graphs. The target is a disclosed central cell; escape means
first hitting the outer boundary. An observer knows the initial position
and reconstructs every later position from the selected seam and its current.
It selects an incident seam uniformly and applies the pinned fair endpoint
event, stopping on target or boundary. The entire-cell target makes its
hit probability h_B(x) independent of the internal port at starting cell x.

Internal moves and holds do not change cell. Conditional on a spatial move,
the six neighboring cells are equally likely. Therefore h_B is exactly the
simple cubic discrete harmonic function with target value1 and boundary0.
This is derived from native current probabilities. No numerical Green or
Coulomb table is supplied to the controller.

A probe can be prepared at x using only native events. Keep its known final
position from the previous probe, choose a declared coordinate path to x,
and repeatedly attempt each path edge until its current records the unit's
move. Each such relocation edge succeeds with probability1/2, so the mean
cost is two native calls per edge, with every failed attempt retained.
The internal port can remain the current known one. There is no hidden-state
read or externally reset location. The preparation terminates almost surely;
a bounded implementation must retain its timeout and resulting error.

## 2. Scheduler feedback on an outer source

The outer source is the common periodic four-axis native graph, with fixed
global seam proposals: internal weights1 and external weights w=a^(-2).
Let W be the total proposal weight, and disclose the outer unit's initial
position. Each completed outer macro-attempt does the following.

* On an external proposal, execute the proposed native endpoint event.
* On an internal proposal in cell x, prepare and run the inner probe from
  its spatial coordinate. If the probe hits the target, execute a native
  event on a disclosed zero-load auxiliary pair instead of the outer seam.
  If it escapes, execute the proposed outer event. Every inner event and
  the actual final selected seam are retained.

The zero-load pair stays zero, so repeated use needs no reset. If an outer
proposal has no occupied endpoint, the observer may skip its inner probe
and execute that zero-current outer event; continuity identifies this case
from records. This shortcut changes no outer transition. Variable inner
preparation/probe costs are charged as actual native calls, not erased by
naming one macro-attempt a unit of physical time.

After each exact preparation the inner target probability is h_B(x),
independent of earlier inner paths. Both endpoints of an internal seam are
in the same cell, so the acceptance factor1-h_B(x) is the same in both
directions. The outer marginal is thus a symmetric Markov transfer

    P_out = I - L_eff/(2W),
    L_eff = w L_external + sum_x (1-h_B(x)) L_internal,x.     (1)

Equation(1) follows from actual controller outcome probabilities, not from
choosing a desired potential and asking a solver to reproduce it. Since
0<=h_B<=1, P_out is a convex mixture of the ordinary native one-seam
projections and the identity. It is positive semidefinite and preserves the
uniform reference. Connectivity holds when some internal repair remains
available and the external routing connects the cells.

## 3. Exact native band consequence

The same internal scalar factor multiplies all thirty internal seams at x.
Consequently the native internal eigenspaces remain invariant. For an
internal eigenvalue lambda, the outer operator is exactly

    L_lambda = -Delta_time,a + S_lambda,a,
    S_lambda,a = -Delta_space,a + lambda(1-h_B(x)).            (2)

Its attractive term relative to the asymptotic threshold lambda is
-lambda h_B(x). Both the strength lambda and h_B come from native records:
lambda from the internal transition and h_B from the auxiliary first-hit
law. There is no independently selected gamma or numerical potential table.
The controller architecture and the target preparation remain disclosed.

Stability needs no added inequality: S_lambda,a is a sum of positive
kinetic and nonnegative multiplication forms. The lambda=0 band is
unchanged by the feedback. Thus the constant-sector free Maxwell source
retains its previous law, while the nonconstant scalar sectors respond.
This statement is about source sectors; it is not yet an assignment of
electric charge or a universal electromagnetic coupling.

The existing Gaussian occupation proof applies with the scale W in(1),
even though the effective edge weights sum to less than W: the remaining
proposal probability is the explicit identity/auxiliary-hold event.
Its limiting covariance is L_eff^+-Pi/(4W). In each internal band, the
time-slice transfer of a spatial eigenmode S psi=s psi is

    exp[-2 asinh(a sqrt(s)/2)],

with the same first-hit clock construction and positive scalar reflection.
The reconstructed standing mode has spatial profile psi and energy
E_a=2asinh(a sqrt(s)/2)/a. This is an effective field description of native
controller records; a primitive wave phase is not an extra source operation.

## 4. A small exact source-derived well

For the inner box[-2,2]^3, the three nonterminal cell orbits are100,110,111.
Writing their hit probabilities h1,h2,h3, native neighbor counts give

    6h1=1+4h2, 6h2=2h1+2h3, 6h3=3h2.

Thus h1=5/22, h2=1/11, h3=1/22. The conditional probability for an outer
unit at one endpoint of a proposed internal seam to move to the other is
(1-h_i)/2. At100 this is17/44, while feedback-disabled repair gives1/2.
At the target h=1 and the outer internal seam is always replaced by the
zero-load auxiliary event. These are direct record-controller predictions.

The transient orbit matrix per incident native attempt is

    Q=(1/22) [[16,4,0],[2,16,2],[0,3,16]].

Its maximum row sum is10/11. Hence a probe capped at T has unresolved
probability at most(10/11)^T. Its exact mean uncapped attempt counts are
31/3,10,26/3. A capped controller that holds the outer source on both target
hit and timeout has effective acceptance equal to the recorded escape
probability by T. Equation(1) then holds exactly with that finite acceptance;
its discrepancy from1-h_B is at most the retained timeout probability.

On a spatial3^3 outer torus, these27 target/neighbor values have mean14/99.
For every positive native lambda, a uniform spatial trial vector has
Rayleigh value85lambda/99. Therefore the lowest S_lambda eigenvalue is
strictly below lambda. It is strictly positive: zero kinetic energy would
force a constant vector, and the nonnegative potential vanishes only at
the target. This proves a nontrivial stable bound-well component before
any orbital eigensolver or hydrogen fit is used. The torus alone does not
establish a Rydberg ladder.

## 5. Long-distance target and scaling, with the remaining proof named

The infinite cubic first-hit limit is h_infinity(z)=G(z,0)/G(0,0), with
G the Green function of the source cubic Laplacian. Its far tail is
c_star/|z|+O(|z|^-3), c_star=1/(4pi G(0,0)). Finite probes take the boundary
to infinity; an unobservable event “escape to infinity” is never executed.
The finite-boundary error is bounded by the maximum boundary value of the
infinite h through the discrete maximum principle.

A single shrinking target does not give nonzero binding at a fixed physical
scale under a->0. The explicit weak-binding scale is y=a² z. On that scale
the proposed limiting operator is

    (S_lambda,a-lambda)/a² -> -Delta_y - lambda c_star/|y|.

Take the probe radius and outer box radius much larger than a^(-2), for
example order a^(-3), so that the boundary error divided by a² vanishes.
The candidate hydrogenic squared-level shifts are then
-a²(lambda c_star)²/(4n²), with shell multiplicity n², and the physical
energy shifts from sqrt(lambda) are
-a²lambda^(3/2)c_star²/(8n²) to leading order.

This paragraph is a next theorem target, not a sealed spectral conclusion.
A controlled discrete Coulomb form/min-max limit and finite ordinal/shell
checks are still required. The small controller experiment proves the source
interaction and well; it must not be reported as already proving hydrogenic
levels. Charge, nuclear species and the physical strength dictionary are
also unassigned. The source construction makes those remaining questions
concrete while retaining the observer-first standing-wave route.
