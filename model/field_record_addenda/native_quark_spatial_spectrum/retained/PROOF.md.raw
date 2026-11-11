# The spectral content of one native quark hop

The Hamiltonian is the original eight-cell source `h=H0+lambda T`,
`lambda=2^-34`, with thirty CAR modes, five scalar states per cell and
PW119 links. The preparation, partner source, full static reference
`4 N_ref` and every preparation record are those of
`interacting_ground_preparation` and `native_quark_response`.
The recorded atomic localization and compensated removal pump are those
of `native_quark_transport`. No hopping or string term is added.

There are two results. The selected first-hop spectral measure is exactly
the local quark removal measure shifted upward by one link-energy unit.
The exact visible initial and first-hop energies are nonresonant. A
complete ordered arrival-then-energy experiment has a strictly positive
full-H probability after one source-time unit. The arrival projection
changes the state before the energy instrument. Its recorded energy is
therefore not an autonomous energy gain or a simultaneous energy/position
PVM.

## Scalar spectral compression on the whole doublet

Write `p` for the actual local neutral ground doublet and `e` for its
energy. For family r let epsilon_j and w_j be all the native local removal
energies and one-component weights in `native_quark_response`, including
the exact zero residues. They satisfy `sum_j w_j=1/2`.

The local Hamiltonian commutes with the separate quark and lepton spin
actions. Its maximal quark Spin(5) module has multiplicity one, hence its
commuting quark spin is trivial. Gauge covariance and quark-spin
covariance give the internal and spin deltas below. Lepton spin commutes
with the quark probe and the Hamiltonian; Schur's lemma on the ground
doublet fixes the remaining operator to a scalar. The actual one-component
spectral measure determines that scalar:

`p c_(s,a)^dagger pi_j c_(t,b) p = delta_st delta_ab w_j p`. (1)

Here pi_j is the exact *local* Hamiltonian projector, rather than a
projector onto a prescribed quark species. Summing (1) gives the
`delta_st delta_ab p/2` compression. There is no mixed-state-only
restriction: (1) holds on every ground-spin input and its retained
environment. No quark state or weight is imported from a lepton channel.
Only the common neutral root e and full hopping norm are shared.

## Literal spectral route matrix

Let P0 put all eight cells in their true local ground doublets and all
links blank. Let B be the full compensated pump from the static-reference
vacuum, with coefficient `1/sqrt(2d_r)`. Let Pe be any of the six actual
incident routes: only that link is excited, every nonarrival cell is in
its local ground, and arrival matter and reference states are unrestricted.
These are the source projectors of `native_quark_transport`.
They obey `Pe B P0=0` and `Pe exp(-it H0) B P0=0` exactly. Thus every
selected arrival is absent before the wait and in the zero-hopping
control. A positive delayed arrival is a change in the source/flux
record, rather than a relabeling of its fixed boundary reference.

Let Qj be the initial atomic spectral projector after B and Rk the arrival
atomic spectral projector within Pe. On these supported vectors,

`(H0+4N_ref) Qj B P0 = (7e+epsilon_j+4) Qj B P0`,
`(H0+4N_ref) Rk Pe T B P0 = (7e+epsilon_k+5) Rk Pe T B P0`. (2)

The extra one is the literal `I-|blank><blank|` link penalty. The static
reference costs four before and after the hop and stays at the origin.
Every nonarrival cell returns to energy e. The matter at the destination
has the same local removal support as the original pump.

The full route calculation retains all reference components, the actual
Pauli hopping matrix, its coefficient of magnitude 1/4 and the PW blank
Gram matrix `I/d_r`. Applying (1) at the origin twice and the destination
once gives the exact cross-input effect

`P0 B^dagger Ql T Pe Rk Pe T Qj B P0`
`                  = w_l w_j w_k P0/16`.             (3)

The spin sum is two, the two internal sums give `d_r^2`, and the pump,
hopping and link factors are `1/(2d_r)`, `1/16` and `1/d_r`. Thus their
product is `1/16` after the spin sum. `spatial.py` executes every spin
contraction for all three Pauli matrices and all three family dimensions,
together with the actual source's six route and link-Gram controls.
The CAR parity factor is common to a selected route and cancels in (3).

In particular an initial energy record j followed by an arrival energy
record k has first-hop weight `w_j^2 w_k/16`. Initial energy indices must
be summed coherently when no initial energy measurement is performed.
Summing j and l in (3), using `sum w=1/2`, gives

`P0 B^dagger T Pe Rk Pe T B P0 = w_k P0/64`,
`P0 B^dagger T Pe T B P0 = P0/128`.                  (4)

Equivalently, for the local Stieltjes transform
`R_r(z)=sum_j w_j/(z-epsilon_j)`, the selected-route transform, with
`7e+4` removed from its energy variable, is `R_r(z-1)/64`. Its normalized
energy law is the initial pump's normalized law shifted by one. Its
conditional mean therefore increases by exactly one at the first-hop
coefficient level. Equation (3) also shows why the route is not diagonal
in the incoming local spectral lines.

## Nonresonance and the actual first-order clock

The two thirteen-root visible support polynomials obey
`gcd(q_r(z),q_r(z-1))=1`. Exact isolating intervals for every visible root
give, for all j,k in the probe support,

`|epsilon_k+1-epsilon_j| > 0.0749111786` for uc and dc,
`|epsilon_k+1-epsilon_j| > 0.0511420752` for Q.        (5)

The twelve dark Q roots are retained in the final complete energy
instrument; they have zero coefficient in the atomic first-hop map.
Equation (5) concerns these particular initial and arrival supports.
Other charged matter and scalar sectors are present in the full source.
Their energies are not bounded below by the neutral local energy by this
argument.

Let the actual state conditional on physical preparation be z, including
its environment, and let `a=||P0 z||^2`. For the actual pump,
`||B P0 z||^2=a/2`. Expand the complete source propagator once in lambda.
On a selected route, its first-order probability is exactly

`lambda^2 a/16 sum_k w_k |sum_j w_j F_kj(t)|^2`,
`F_kj(t)=[exp(-it epsilon_j)-exp(-it(epsilon_k+1))]`
`                                      /(epsilon_k+1-epsilon_j)`. (6)

The fixed phase from `7e+4` and the retained partner unitary cancel from
this probability. They are not removed from the physical clock. Formula
(3) proves (6) for an arbitrary spin-entangled P0 input, without a
comparison to a selected cube spin vector. As t tends to zero, (6)
reduces to `a(lambda t)^2/128`.

With Delta the certified separation (5),

`||first-order route vector||/sqrt(a)`
`       <=lambda sum_k sqrt(w_k)/(4 Delta)`.          (7)

The first-order contribution is bounded in time and has oscillatory
detunings, rather than a resonant term linear in time at fixed energy.
This identifies off-resonant one-flux dressing in the first-order channel.
It is not a proof that every full-H motion is virtual, nor a confinement
or string-tension theorem.

The full interaction-picture Dyson identity supplies a remainder using
unitarity of the retained full propagator:

`||R_2(t) B P0 z|| <= v^2 t^2 sqrt(a)/(2 sqrt(2))`,
`v=lambda ||T|| <=360 lambda`.                       (8)

No exponential norm loss, discarded complement or finite charged-sector
truncation is used. Equations (7) and (8) give a finite-time full-H upper
bound and separate higher-order processes from the first-order statement.

## Positive arrival followed by a full-H energy record

Choose the source wait `t=1`. This is a different declared observation
time from the short-time setting in `native_quark_transport`; the
Hamiltonian, preparation and pump are unchanged. For any outgoing line k,
the j=k term in (6) has detuning exactly one. The reverse triangle
inequality and `|F_kj(t)|<=t` give

`|sum_j w_j F_kj(t)|`
` >= max(0,2 w_k sin(t/2)-(1/2-w_k)t)`.              (9)

An alternating rational Taylor enclosure supplies sin(1/2). The strongest
line has weight greater than .495, so (9) is positive with a wide margin.
All other lines, including exact zeros, remain explicit records.

After the arrival projection, measure the complete spectral bins of the
actual `h_active+4N_ref`, centered on `7e+epsilon_k+5`, with radius `128v`.
The rational bins include the atomic root uncertainty and are disjoint.
Let dk be the distance from the exact atomic energy interval to its bin
complement, and dkl the distance from another atomic interval to this bin.
For a vector in an atomic energy interval,

`||(I-Pi_k) chi_k|| <= v ||chi_k||/dk`,
`||Pi_k chi_l|| <= v ||chi_l||/dkl`, l not equal k.   (10)

The atomic first-order components in (6) have lower norm
`sqrt(a) sqrt(w_k) L_k(t)/4` and upper norm
`sqrt(a) sqrt(w_k) t/8`, where Lk is (9). Define

`A_k = sqrt(w_k) L_k(t) sqrt(1-(v/dk)^2)/4`
`      -sum_(l != k) v sqrt(w_l) t/(8 dkl)`.

Combining (8)-(10), retaining the actual physical preparation success q,
gives the original-attempt lower bound

`q_lower a_lower max(0,lambda A_k-v^2 t^2/(2sqrt(2)))^2`. (11)

The implementation rounds every root, trigonometric series, weight and
spectral leakage outward. It also retains an upper endpoint for every
bin. Its a lower bound is the actual recorded localization bound from
the interacting preparation, not a replacement by an atomic trial.
The full doublet identity eliminates the need to subtract a trial-vector
error from this small transport amplitude.
The outside-all-bins amplitude is bounded by
`lambda sum_k v sqrt(w_k) t/(8 dk)+v^2t^2/(2sqrt(2))`;
its squared upper bound and the complete complementary outcome are retained.

This is the ordered Kraus map `Pi_k Pe exp(-it h) B P0`, with the source
partner and static reference retained. The operator Pe does not commute
with h. The energy record is of the state after the arrival instrument;
it is not a joint spectral measure of position and the premeasurement
energy. In particular, the shift in (4) does not show that autonomous
evolution created one unit of energy. The unitary wait conserves the
actual source energy; the declared intervention can exchange energy.

There is an exact same-H ordering control. The full-H bins centered on
the visible *unshifted* pump energies and the one-link bins have separation
at least `Delta-256v>0`. Their spectral projectors satisfy
`Pi_shifted exp(-it H_total) Pi_unshifted=0` at every time. An initial
full-H energy admission cannot acquire a disjoint energy record under the
autonomous wait alone. Inserting the arrival projector before the final
energy readout changes this composition. This operator identity is an
ordering control, not an additional executed preparation setting or a
claim that the localized pump is an exact full-H energy eigenstate.

## Complete law, apparatus and remaining domain

After the original physical preparation, a duration-one complete service
records localization failure, pump failure or success. A full autonomous
wait of one source unit follows. One duration-one detector implements
the complete sequential instrument with Kraus maps `P_no_arrival` and
`Pi_k Pe` for each route and every energy bin, including its outside
projector. Their adjoint products sum exactly to I. Failure branches run
matching no-op controls and retain their source states.

With n the number of atomic roots, the final detector has
`1+6(n+1)` outcomes. The complete original law has `4+6(n+1)` non-padded
coarse labels: physical preparation failure, localization failure, pump
failure, pump-success nonarrival, and each route/energy/outside label.
All unused pointer inputs and zero-probability records are retained.
Every fine preparation flag and every source, partner and reference is
present. Static reference charge stays at the origin; flux cancels it
there and cancels the conjugate removal charge at the destination.

The two specified services use Hermitian logarithms of finite unitary
completions. Constant controls `i Log U-H_total` explicitly cancel the
background during their duration-one slots. The raw two-source plus
thirty-mode-reference norm is at most `1080+720lambda`. Each control is
bounded by that value plus pi; gross switching work is at most four
times that bound. Extra recorded capacity is
`ln[3(1+6(n+1))]` nats beyond the retained inherited systems. Controller
synthesis workspace and efficient realization are separate premises.
The complete added elapsed source time is three on every branch.

The certificate gives a source-calculated spatial/energy record and an
exact off-resonant spectral route matrix. It supplies no lower bound
across all separated charged sectors, asymptotic colored state,
ballistic dispersion, relativistic speed, particle identity or GeV mass.
