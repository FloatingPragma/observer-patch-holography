# Reciprocal response of two mobile native copies

`reciprocal_feedback.py` replaces the fixed compensator by a second mobile
native one-load carrier. Each copy has the same 96 ports, 288 seams, fair
coin and literal C repair. Their signs are declared readout labels. The
charge is rho=delta_x-delta_y, reconstructed from the two known initial
preparations and their actual current records. The source positions are
passed only to native repair, never to the scheduling controller.

Both coordinates use the same energy

    U_K(x,y) = (kappa/2)(delta_x-delta_y)*G_K(delta_x-delta_y).

The finite local relaxation circuit and exact spectral error bound are
those in `PROOF.md`. If Phi has columns G_K(delta_x-delta_0), the energy
table follows from

    U_K(x,y)/kappa = [Phi_xx+Phi_yy-Phi_xy-Phi_yx]/2.

This cancels the computational reference exactly. No pinned physical charge
remains. The code symmetrizes only floating roundoff after checking the
theoretical table symmetry and zero diagonal to 1e-14. Every source state
and both sign labels remain in the declared 9,216-state joint domain.
At coincidence both carrier loads survive; no annihilation rule is added.

## One joint clock and two response channels

A seam of either copy has attempt rate
nu exp[-beta Delta U_K/2], where Delta U_K changes that copy's coordinate
and holds the other coordinate fixed. Nonincident seams have rate nu.
The joint total rate is the sum over both copies. A scheduled mark specifies
the copy, seam and fresh fair coin. Both repairs therefore act in the same
clock and both rates depend on the actual recorded state of the other copy.

The joint generator obeys detailed balance with exp[-beta U_K(x,y)]. Its
complete uncapped position law has 9,216 outcomes. Exchange of the two
copies leaves U_K and the generator invariant. The receipt verifies the
full generator identity and equality of the two response contributions
under a uniform joint preparation. At beta=4, nu=1 and kappa=1, each copy
contributes -0.005318327866635085 to the initial rate of change of the
unit-kappa energy. At kappa=0 both contributions vanish. At kappa=2 each
contribution is -0.01106023698768859.

For every actual jump the work readout is the change of the shared energy
computed from the observed current. The two work sums satisfy, pathwise,

    W_plus + W_minus = U_K(final) - U_K(initial).

Holding events contribute zero and remain in the tape. The signed flow
F_plus-F_minus has divergence rho_initial-rho_final on every record.
These are energy and charge accounting identities in the declared
controller. The stochastic thermal scheduler can exchange energy with its
environment; the result does not assert isolated Hamiltonian conservation
or equal and opposite instantaneous mechanical momentum.

The executable runs 24 complete development paths containing 6,379 native
events, 6,321 zero currents and 31/27 moves of the two respective copies.
Every event is independently replayed with endpoint floor/ceiling repair,
including its waiting-time density, mark probability, energy increment and
signed divergence. The largest tested change in one copy's seam rate when
only the other copy's configuration changes is 0.3606004729302097.

An injected omission of the first nonzero current stops the control at that
event with a distinct flag. The observer's counter then disagrees with the
source, and the next scheduler settings differ by up to
0.17327752554994058. This verifies that actual current information is doing
work in the coupling; discarding it does not preserve the executed law.

The same finite-field approximation gives a bound below 9.014e-7 on total
variation of complete exact-arithmetic path laws at time1/2 and kappa=1.
The factor doubles because up to twelve incident seam rates can differ.
Floating rounding remains separately scoped. `CAPPED_INSTRUMENT.md`
defines the complete finite-cap instrument and its nonzero failure tail;
the uncapped position semigroup alone is not that instrument.

The source architecture now has two mobile native copies and a common
record feedback controller. It removes fixed-source dynamics from this
response experiment while retaining the supplied signs, kappa, beta,
energy functional and centralized rate-selection service. It is a
reciprocal dissipative electric response, without a transverse field,
Lorentz limit, physical coupling calibration or orbital spectrum.

The latest coordination read for this extension is Codex C029/Fable F009.
The prior review's empty N31 directory and absent WQ-11 receipt were
rechecked. A later producer receipt needs a fresh source-specific review.
