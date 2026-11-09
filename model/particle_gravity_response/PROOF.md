# Same-program energy transport and action response

The [gauge-coordinate field program](../particle_observer_outputs/FIELD_ROUTE.md)
and [finite scalar source](../particle_observer_outputs/HIGGS_RECORDS.md)
admit energy partitions, energy-transfer detectors and fixed-state action
responses on their stated Hamiltonian clock. This construction supplies the
finite conservation input of the coarse-energy theorem. Its coefficients
are declared program inputs; no particle mass or gravitational measurement
is used. A source coefficient is not assigned a geometric meaning solely
because its derivative is observable.

## Energy and work from one Hamiltonian

Let the complete retained finite field carrier have a Hermitian Hamiltonian

\[
H=\sum_x h_x,\qquad h_x=h_x^\dagger,
\qquad U(t)=e^{-itH}.
\tag{1}
\]

The index can label spatial cells or energy components. A spatial
assignment must be specified independently of the continuity identity.
Define the action's energy current

\[
J_{xy}=i[h_x,h_y].\tag{2}
\]

Then `J_xy` is Hermitian, `J_yx=-J_xy`, and

\[
\frac{d}{dt}U(t)^\dagger h_xU(t)
 =-\sum_y U(t)^\dagger J_{xy}U(t).
\tag{3}
\]

This is an operator identity, so it holds in every density and with any
untouched reference. Summing over a region cancels its internal pairs and
retains every crossing pair. Integrating over a common time interval gives
the exact endpoint energy increment plus integrated outward flux equal
to zero. No approximation replaces the integral by an endpoint current.

The declared weighted energy-source family is `H_f(t)=sum_x f_x(t) h_x`,
with real control functions. Setting `e_x(t)=f_x(t)h_x` gives

\[
\frac{d}{dt}\langle e_x\rangle
 +\sum_y f_xf_y\langle J_{xy}\rangle
 =\dot f_x\langle h_x\rangle,\qquad
\frac{d}{dt}\langle H_f\rangle
 =\sum_x\dot f_x\langle h_x\rangle.
\tag{4}
\]

The right side accounts for controlled work. At fixed state,
`partial H_f/partial f_x=h_x`; this is the exact response to an energy-source
coefficient. If an independently defined metric dictionary identifies
`f_x` as a lapse, this becomes its finite lapse response. Equation (4)
alone defines neither proper time nor a symmetric spacetime stress tensor.
Detector and preparation work also belongs in the Hamiltonian whenever
those apparatus degrees of freedom participate in an energy experiment.

If every `h_x` commutes with a gauge generator, Jacobi's identity shows that
every `J_xy` does too. The current is formed from energy operators, including
the interaction energy. A conserved particle-number current is insufficient
to replace (2).

The source-to-theorem attachment is direct. For a state `rho`, set
`delta_x=Tr rho(U(T)^dagger h_x U(T)-h_x)` and
`Phi_xy=integral_0^T Tr rho U(t)^dagger J_xy U(t) dt`. Equation (3) supplies
the premise of
[`OPH.QFT.CoarseEnergyWard.coarse_continuity`](../../Lean/QFT/CoarseEnergyWard.lean).
`coarse_compose` preserves it under consecutive block maps;
`time_increments_telescope` preserves exact consecutive time intervals.
These finite sum theorems are reused. The source matrices are checked by
Python exact arithmetic; this package adds no Lean source encoding.

## Assignment on the declared spatial field

For a finite cubic instance of the coordinate-field program, divide each
Hermitian covariant hopping or overlap-polynomial term among its declared
support sites, divide each electric link term between its endpoints,
divide each plaquette term among its corner sites, and assign each scalar
number term and gauge-invariant pair vertex to its site. Any fixed real
weights summing to one per term define a valid Hermitian partition (1).
Its coefficients and full Hamiltonian are those of the chosen program.
No new interaction is added by this energy assignment.

On the uncompressed finite-range field, disjoint even CAR terms with
disjoint link supports commute. On the finite total-representation-degree
compression used for computation, `P_K h_x P_K` need not preserve this
support commutation. Consequently the compressed construction retains
every pair in (2). Inferring a bounded current range from the cutoff matrix
alone is invalid. The [coordinate-field locality proof](../particle_field_causality/PROOF.md)
supplies an exact cutoff counterexample and the parent field estimates. The separate
complete-history cutoff errors are needed to transport a locality claim.
The finite conservation identities hold for the compressed action without
such a claim.
The all-state operator statements in (1)-(4) concern finite retained
carriers. On the uncut coordinate Hilbert space the electric energies are
unbounded; an energy-transport statement there needs a common operator
core and finite-energy/domain control. The bounded-observable locality
theorem does not remove those energy-domain requirements.

The [Weyl/contact energy package](../field_stress/README.md) verifies the
same general identity on a different specified action. Its numerical
current matrix is not substituted for this program's current. The
executable controls below replay the actual coordinate-field and scalar
matrices. They check component exchange within a spatial cell, not the
full cubic spatial tensor or its continuum limit.

## The complete quantum-flux control

Use the literal `field_route.exact_rotor(6)`. Its carrier has sixteen CAR
occupations and all thirteen electric fluxes `-6,...,6`, hence dimension
208. Its charge is the native singlet value `q6=6`, and its Hamiltonian is

\[
H=H_E+H_M+H_B,\quad H_E=E^2/2,\quad
H_M=d\Gamma(\beta),\quad
H_B=d\Gamma(b)\otimes S_6+d\Gamma(b)^\dagger\otimes S_6^\dagger.
\tag{5}
\]

Here `b=(-beta+beta gamma)/2` is the owning source's exact one-cell
coefficient; the flux shift has its full cutoff endpoints. The prepared
state is the owning Gauss-zero state
`(|mask3,flux0>+i|mask3,flux6>)/sqrt(2)`. Every energy term preserves the
same Gauss generator. Exact operator arithmetic gives the initial values

| Component | Mean energy | Mean time derivative |
|---|---:|---:|
| Electric | 9 | 18 |
| Dirac onsite | 2 | 0 |
| Covariant hopping | 0 | -18 |

The total energy is 11. The nonzero current operators are `J_EB` and
`J_MB`; the latter has zero expectation in this particular input.
`J_EM=0`. No operator is dropped because its initial mean vanishes.
The same source has `d<E>/dt=6`. Under the declared perturbation
`H(lambda)=H+lambda H_B`, the exact initial slope satisfies
`partial_lambda partial_t <E(t)>|_(0,0)=6`. This is an observable response
to the hopping coefficient. No metric interpretation of that coefficient
is inserted into the proof.

## Scalar energy exchange with a complete detector

Use the actual three-family quantum-Higgs vertex with `Y=diag(1,2,3)` and
pair cost `d=1` on its full 307-dimensional `N_F-2N_H=0` sector. The owning
source supplies `g^2=4 Tr(Y^dagger Y)=56`, `Delta=sqrt(d^2+4g^2)=15`,
the empty vector `Omega`, normalized bright vector `b=W Omega/g`, and the
305-dimensional dark projector `D`. Its energy decomposition is

\[
H=S+V,\quad S=dN_H,\quad
V=g(|\Omega\rangle\langle b|+|b\rangle\langle\Omega|).
\tag{6}
\]

Choose the reproducible state
`psi=(Omega+i b)/sqrt(2)`, independent of `d`. The empty state and
`W Omega` are gauge singlets; their superposition is a singlet. The full
307-state carrier also contains nonsinglet dark states. They are retained
as operator inputs and detector outcomes, without asserting that every
state is a Gauss-invariant preparation.

The energy current is

\[
J_{SV}=idg(|b\rangle\langle\Omega|-|\Omega\rangle\langle b|).
\tag{7}
\]

Its complete spectral menu has values `(dg,-dg,0)` and projectors
`(P_Omega+P_b+J/(dg))/2`, `(P_Omega+P_b-J/(dg))/2`, and `D`, with ranks
`(1,1,305)`. The initial probabilities are `(1,0,0)`. The exact source
propagator gives two alternative complete detector menus:

\[
\begin{aligned}
p_{\rm empty}(t)&=\tfrac12+\frac g\Delta\sin(\Delta t),&
p_{\rm bright}(t)&=\tfrac12-\frac g\Delta\sin(\Delta t),&p_D(t)&=0,\\
p_{J,+}(t)&=\tfrac12(1+\cos\Delta t),&
p_{J,-}(t)&=\tfrac12(1-\cos\Delta t),&p_{J,0}(t)&=0.
\end{aligned}\tag{8}
\]

The zero-probability dark outcome is retained. It has probability one on
a normalized dark input. The two menus are alternative experiments;
equation (8) does not assign their noncommuting outcomes a joint classical
value. The full source propagator, including the dark phase, is defined
in the owning package. The bright reduction in this prepared experiment
follows from its invariant subspace, without deleting dark inputs.
The executable readout vector removes the common phase `exp(-i d t/2)`
of this prepared bright state. Its density and all record probabilities
are unchanged; this does not define a replacement propagator on dark inputs.

Their energy and integrated flux are

\[
\langle S(t)\rangle=\frac d2-\frac{dg}\Delta\sin(\Delta t),\quad
\langle V(t)\rangle=\frac{dg}\Delta\sin(\Delta t),\quad
\int_0^t\langle J_{SV}(u)\rangle du=\frac{dg}\Delta\sin(\Delta t).
\tag{9}
\]

The initial component rates are `(-2sqrt(14),2sqrt(14))`, and total energy
is `1/2`. Equations (8)-(9) supply explicit endpoint increments and
integrated fluxes for the finite coarse-energy theorem.

A combined control declares
`H_total=H_rotor tensor I+I tensor H_scalar` on dimension 63856. These
commuting terms evolve on the same chosen program clock. The scalar vertex
is interacting; no rotor-scalar interaction is declared in this product
control. Its initial complete flux/occupation detector has 26 outcomes.
Four have probability `1/4` and twenty have probability zero. All outcomes
and any reference are retained. This product control does not claim a
fully coupled spatial gauge-Higgs dynamics.

## A fixed-state action response and its error bound

Change the same source pair cost to `d+lambda`, keep `Y`, the preparation
`psi`, apparatus and time `t` fixed, and retain both occupation outcomes.
The Hamiltonian change is exactly `lambda N_H`. Differentiating its actual
occupation probability in (8) gives

\[
\frac{\partial p_{\rm occupied}(t)}{\partial d}
 =\frac{gd}{\Delta^3}\bigl(\sin\Delta t-\Delta t\cos\Delta t\bigr).
\tag{10}
\]

At `Y=diag(1,2,3)`, `d=1`, `t=pi/30`, it is `2sqrt(14)/3375`. This uses
the same input density for every perturbation, rather than silently
repreparing the perturbed ground state. The pair cost is an action
coefficient. Its relationship to a metric or strain is not defined by
(10).

More generally, Duhamel differentiation of `H_lambda=H+lambda V` gives

\[
\partial_\lambda\langle A(t)\rangle\big|_0
 =i\int_0^t\langle[V(s),A(t)]\rangle ds.
\tag{11}
\]

The formula uses the fixed initial state and the unperturbed Heisenberg
operators. Each nested commutator contributes at most `2 norm(V)` to the
operator norm. Time-ordered integration gives
`|d^3 <A(t)>/d lambda^3| <= 8 |t|^3 norm(V)^3 norm(A)` for real lambda.
The central difference at `+h,-h`, `h>0`, therefore has derivative error
at most

\[
\frac43 h^2|t|^3\|V\|^3\|A\|.\tag{12}
\]

For (10), `V=A=N_H` have norm one. If each complete detector mean has
absolute error at most `epsilon`, the central derivative also has error
at most `epsilon/h` from those records. This includes preparation,
evolution, cutoff, apparatus and sampling errors when their stated bounds
cover the complete unconditioned law. For labels bounded by `L`, a full
record `l1` error `delta` contributes at most `L delta` to a mean. Failures
must have explicit labels and errors, rather than being renormalized away.

The [finite readout contract](../particle_observer_outputs/READOUT_CONTRACT.md)
implements every specified finite spectral menu using a full unitary
pointer coupling, with padding and retained references. A primitive
apparatus is also an allowed design route. Repeating one fixed compiled
step defines the common program clock; operator error per step accumulates
under that contract. No equality with an unexamined repair-step count or
physical proper time is asserted.

## The Einstein consumer and the geometric gap

The [matched regional theorem](../matched_einstein_limit/PROOF.md) consumes
one source family with actual curvature, action-derived Lorentz stress,
regional entropy and area in a common proper-time normalization. For each
of ten future timelike observers it requires, with `C=4pi ell^4/15`,

\[
|s_r-2\pi C T_{uu}|\le e_b,\quad
|a_r+C(G_{uu}-\Lambda)|\le e_a,\quad
|s_r+\eta a_r|\le e_s.
\tag{13}
\]

The ten probes are `e0`, `(5e0+4ei)/3`, `(5e0-4ei)/3` for
`i=1,2,3`, and `3e0+2ei+2ej` for `i<j`, in signature `(-,+,+,+)`.
They are future unit timelike vectors. Their exact decoder has sharp
maximum-norm bound `17/8`. For
`R=G+Lambda g-(2pi/eta)T`, the same regional contrasts satisfy

\[
\eta C R_{uu}
 =(s_r-2\pi C T_{uu})
   +\eta[a_r+C(G_{uu}-\Lambda)]-(s_r+\eta a_r).
\tag{14}
\]

Thus the decoder bounds every component of R by
`255(e_b+eta e_a+e_s)/(32 pi eta ell^4)`. The exact inverse,
sharp constant and residual identity are reused from the
[owning finite proof](../matched_einstein_limit/PROOF.md#ten-observers-and-an-exact-inverse)
and [`small_ball_tensor_bound`](../../Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean).
No separate all-null balance premise is needed for this route. This code does not enter
invented values into that consumer. The realized finite energy identity
is only the energy part of a prospective stress conservation law.
An exact tomography null direction makes the distinction explicit:
changing only `T_11` leaves the resting observer's energy contraction
`T_00` unchanged while altering boosted contractions. The executable
reuses the owning ten-observer decoder to verify this direction. Energy
conservation supplies no missing nine tensor components.

The hypotheses required to invoke (13) are:

1. A same-source metric and proper clocks, volumes and boosted regional
   measurements, with actual curvature rather than a tensor defined by
   the desired field equation.
2. One symmetric Lorentz stress from that source's matter action, with
   independent covariant conservation and the measured energy currents
   attached to its appropriate components.
3. The regional modular diamond kernel for that same stress, state and
   observer, including the retained reference and edge contributions.
4. Fixed-proper-volume area response and generalized-entropy stationarity.
5. A common shrinking-ball family with all three errors `o(ell^4)` and
   one positive spacetime-independent limiting area coefficient.
6. A smooth metric limit, its contracted Bianchi identity, and matched
   state and detector convergence.

With `ell_n` tending to zero, the three errors divided by `ell_n^4`
tending to zero, common tensor limits and `eta_n` tending to a positive
spacetime-independent `eta_*`, equation (14) yields
`G+Lambda g=8pi G_N T`, where `G_N=1/(4eta_*)`. Independent matter
conservation and the contracted Bianchi identity make Lambda constant
on each connected component. These are source hypotheses on matched
geometric readouts. First variations about one reference give only the
linearized equation; the nonlinear conclusion requires the corresponding
finite regional contrasts for the admitted states.

### Reusable matter and entropy premises

The source-specific entropy construction supplies part of this interface
beyond the energy identity. The
[actual internal-family marginal](../particle_regional_geometry/PROOF.md)
of the same 307-state scalar Gibbs preparation has a completely positive
restriction, faithful positive support and a gauge-invariant tangent with
`delta S=Tr(delta rho K)`. Its induced pair term and quadratic entropy
remainder are calculated on that state. The explicit tangent lifts to a
positive global preparation with a retained conjugate reference. This
realizes a finite matter entropy-pairing premise; the selected family
is not a spatial ball and its modular generator is not the geometric
kernel in (13). Its tracial regional expectation also fails to preserve
the interacting Gibbs state, so state-preserving modular restriction
cannot be inferred from trace preservation.

| Existing construction | Reusable result | Same-source boundary |
| --- | --- | --- |
| [Canonical fixed-number repair state](../native_gravity_identification/PROOF.md#a-realized-finite-bulk-and-edge-entropy-split) | Exact direct-sum bulk/edge entropy convention and finite relative-entropy first law. | Its regional modular generator is central; nontrivial native attempt evolution is not this modular flow. No geometric area follows from the edge multiplicity. |
| [Free eight-valley continuum](../weyl_geometric_modular/README.md#attachment-to-the-existing-einstein-stack) and [shrinking lattice entropy](../shrinking_modular_entropy/README.md) | Same-vacuum Lorentz stress, geometrically normalized diamond modular charge, matter first law and controlled matter small-ball coefficient on their stated free-source domains. | The free source, vacuum and scaling differ from the interacting scalar and gauge-coordinate programs. Their modular, stress and lattice error bounds cannot be transferred without a state/action/readout map. |
| [Reciprocal area dynamics](../reciprocal_area_dynamics/DERIVATION.md#one-metric-energy-and-constrained-entropy-family) | A responsive two-geometry process, positive thermal entropy/area coefficient and complete energy/entropy balance. | Its full-state balance retains excess matter energy. The thermal coefficient is not a universal Newton constant, and two geometry labels do not supply a spatial shrinking-ball family. |

These results remove the need to reconstruct finite entropy algebra or
tensor tomography. They do not combine their different preparations into
one instance of (13). The remaining attachment is spatial modular/stress
response, edge/area response and fixed-volume entropy balance on the same
interacting state, with the common refinement errors and proper clocks.

The [native gravity identification proof](../native_gravity_identification/PROOF.md)
shows why static cuts and stationary repair entropy do not supply these
geometric identifications. The faithful finite Gibbs construction gives a
state-specific modular clock on its own algebra; it does not establish the
regional diamond kernel in (13). Neither entropy stationarity nor the
Einstein equations are used as definitions of the outputs here.

These requirements concern the observable consequences of a declared
architecture. No uniqueness theorem for that architecture is required.
The finite energy, work and detector construction is valid independently
of whether a particular continuation meets the geometric conditions.
