# Propagating native registers and the common interaction clock

This packet gives two declared completions of the
[finite native-routing join](../native_routing_join/PROOF.md). Both move
native occupations between spatial cells while retaining the actual link
register on which the native repair acts. The first retains all 7,200
primitive native modes. The second adds explicit four-component routing
memory and gives native occupations a controlled free Dirac window on the
same clock as the paired flux. Their finite interaction is calculated;
the interacting relativistic limit is not inferred from their local gates.

## 1. Native transport on the unchanged occupation space

The logical code at each flux lane remains `C11 tensor F_native`, with
the occupation-dependent Gauss encoding

\[
 E|r,b;n\rangle=
 (|r,b;n\rangle+\chi(n)|r\mathbin{\mathrm{xor}}STAR,b;n\rangle)/\sqrt2.
\]

Two independent native registers are retained per cell. On a neighboring
pair of registers, set `A=(1-i)/2`, `B=(1+i)/2` and take the full exterior
lift of the one-body gate `A I+B SWAP`. This transports native modes and
fixes flux and routing labels. It includes the CAR permutation sign and
Pauli exclusion in every occupation sector. Its source realization is
Gauss dressed through E; it does not reset a link or copy a field.
`literal_pair_edge` applies the dressed operation directly to both literal
link words and verifies its intertwiner with the logical native exchange.

Even/odd nearest-neighbor matchings along the three axes give six layers.
A six-state scheduler applies one matching `T_j` per tick, then advances
`j` modulo six. The common tick is

\[
 J_{a,\lambda,j}=e^{-ia\lambda K}T_jR_a,
 \qquad K=\sum_\ell(I-\operatorname{Re}U_{F,\ell}).
\]

Here `R_a` is the complete paired-flux route, native U_F reads the same
retained link word, and all ticks have duration a. The scheduler is explicit
source memory. The flat/native-empty data reference tensored with the
uniform scheduler state is stationary. Total native number is conserved;
the individual cell occupations are not. The K bound `0<=K<=2N` and the
controlled Taylor/error theorem transfer by their actual algebraic identities.

The executed side-four experiment starts the flux at `(1,1,1)` and a native
wave at `(3,2,2)`, lane zero. The first native matching sends the wave to
`(2,2,2)` with probability one half. Its inverse is checked. After routing
and matching, the state has 32 amplitudes; K produces 2,880 amplitudes.
The joint complete readout asks whether the native wave arrived and whether
the flux changed from its prepared internal face. With `epsilon=a lambda`,

\[
 \|K\psi\|^2=283/72,\qquad
 \|P_{\rm arrived,changed}K\psi\|^2=1/38400.
\]

For `a=pi/12`, `lambda=3/(2500pi)`, use `epsilon=1/10000` and Taylor error
`delta=1/50000000`. Since `195<sqrt(38400)<196`, the exact changed-face
probability p lies in

\[
 [(\epsilon/196-\delta)^2,(\epsilon/195+\delta)^2].
\]

All four exact records, ordered as `(arrived,changed)=00,01,10,11`, are
`1/2,0,1/2-p,p`. Native number at the receiver is conserved during this
weak gate, which proves these exact marginal relations. Approximate Taylor
weights are recorded separately and are never renormalized into exact records.

## 2. Why the primitive exchange is not a relativistic native kinetic law

For one axis, the two matching layers have a two-cell Bloch symbol

\[
 T(q)=(A I+B\begin{pmatrix}0&e^{-iq}\\e^{iq}&0\end{pmatrix})
       (A I+B\begin{pmatrix}0&1\\1&0\end{pmatrix}).
\]

Its trace is `-2i sin^2(q/2)`, determinant minus one and eigenvalues
`exp(-i theta(q)),-exp(i theta(q))`, where
`theta(q)=arcsin(sin^2(q/2))`. The low branch over a six-tick cubic cycle
has frequency `sum_i theta(2a p_i)/(6a)` because the clock duration is 6a.
For `|a p_i|<=1/2`, `theta(2a p_i)<=2a^2 p_i^2`, so this frequency is at
most `a |p|^2/3`. It tends to zero on every fixed physical momentum window.
The other branches have carrier phases and the same quadratic deviations.
Relabelling this quadratic branch cannot make it a native relativistic mode
on the flux route's a-clock. This is an exact kinetic restriction of this
specified matching law, without an exclusion of other architectures.

## 3. An explicit spinful native extension with the same link

The second completion adds a four-dimensional native routing memory.
There are 28,800 native one-body modes per link lane. Primitive port,
internal charge, repair memory and seam address retain their original
meaning. They are not inferred to be physical spin.

The shared-link native wave operation is the exterior lift of `W_g tensor I4`.
Each native controlled link transport uses endpoint-number parity summed
over all four routing components. This is one common link register, rather
than four independent link backgrounds or sequential unrelated source ticks.
The code in `spinful.py` defines this operation on every CAR occupation.
Its inverse reverses the controls and wave lift. It agrees with the actual
native one-wave U in each spin component, which is checked by full source
re-encoding. The Gauss sign counts all native port-five doublet occupations.
Native spin rotations and primitive `U(3) x U(2)` currents commute with
the wave and the summed-parity controls. These internal covariance statements
do not alone prove Lorentz covariance of the interaction.

The execution reuses the proved mobile-source one-object isometry when a
flux configuration contains one object. General flux configurations use
the complete tensor-lane routing rule. This avoids repeatedly applying
identity coins to empty cells without discarding collision configurations.

Use the complete exterior lift `Gamma(D_a^0)` of the three directional
Clifford shifts, with native mass angle zero. It acts on the additional
memory and transports the actual native occupations. The local spin coins
and conditional translations have number-preserving exterior lifts on the
full finite Fock space. No two-particle collision sector is discarded.
The same occupation-dependent E lifts them into the source Gauss sector.
Thus a second explicit clock is

\[
 J^K_a=e^{-ia\lambda K}\,\Gamma(D_a^0)R_a.
\]

The native and flux kinetic factors commute, acting on different logical
registers. Their source implementations retain the same link custody via E.
With native interaction switched off, the one-flux/one-native sector has
free frequencies `+/-sqrt(|p_F|^2+4)` and `+/-|p_N|`, on one spacing/time
convention. For momentum windows `|p_F,i|<=P_F`, `|p_N,i|<=P_N`, the
combined free approximation error is at most

\[
 ta(3P_F^2+6P_F+3P_N^2).
\]

This follows from the existing Clifford product estimate and the tensor
product norm inequality. An N-native exterior lift multiplies its one-body
error by at most N on the corresponding fixed-number momentum window.
The native-empty sector recovers the complete flux route exactly.
The empty native Fock vacuum is stationary, without being proved the
positive-energy vacuum of the signed Dirac spectrum.

## 4. Bulk native energy and an explicit flat-background subtraction

On a flat flux orbit, the native positive generator is a nonzero operator
`K0`. Extend this actual flat-orbit native operator identically across the
flux labels, retaining the same Gauss encoding. On one native occupation
it is a uniform primitive-mode operator kappa tensor spin identity.
Consequently the asymptotic native branches of `J^K` are
`lambda*kappa_j +/- |p|` in the continuum window. Nonzero additive offsets
are not relativistic mass terms. Kinetic spin memory alone does not fix this
bulk consequence of the literal positive-clock choice.

One may instead declare the explicitly different clock

\[
 V=K-K_0,\qquad J^V_a=e^{-ia\lambda V}\,\Gamma(D_a^0)R_a.
\]

This is a flat-background subtraction, not a change of scalar energy zero,
an identity between the two clocks, or Wick normal ordering. It retains the
literal link-dependent native response. It vanishes when no flux is present
or the native register is empty. Its free comparison clock is exactly the
two kinetic factors above. On fixed native number N, `||V||<=4N` is a valid
bound. A positive many-body reference is not established by this subtraction.

The executed native incoming packet is the exact inverse native Dirac step
of a localized receiver state. It has eight amplitudes and zero receiver
occupation initially. The forward native step places it at the receiver
with probability one. Simultaneously the flux executes its own route.
The exact joint K action then gives changed-face squared norm `1/19200`.
K0 cannot change a face, so the V action gives the same changed-face norm.
The exact complete arrival/face records and controlled weak-gate bounds
are retained. Turning off native propagation leaves its incoming support
disjoint from the routed flux, giving zero changed-face response. The
calculation therefore involves two propagated source objects, not a pinned
native register. This sharply localized finite preparation does not itself
define a fixed-energy continuum state.

## 5. A bounded-contact survival estimate on the free reference

Take one flux object and one native occupation. Let C project onto their
coincidence in the same spatial link lane. Q and native number are preserved
by each local repair, so

\[
 V=CVC,\qquad \|V\|\le4.
\]

For `J^K`, choose the free reference gate with K replaced by K0; for `J^V`,
choose the gate with V replaced by zero. Duhamel's exact integral identity
and telescoping the complete finite ticks give

\[
 \|(J_n\cdots J_1-J_{0,n}\cdots J_{0,1})\psi\|
 \le4a|\lambda|\sum_{j=1}^n\|C\psi_{0,j}\|.
\]

Here psi_0,j is the freely evolved reference just before the weak gate.
The unitary future factors are those of the interacting law; only their
norm is used. This estimate makes no assumption that the interacting state
has the reference's coincidence probabilities. K0 commutes with C inside
each weak step, which permits the displayed contact bound for `J^K` too.

On a fixed physical torus of side ell, let the native state lie in M Bloch
momenta, with all its finite internal components retained. The reference
evolution preserves that window. For the one-cell Dirac route, the exact
compressed position-projector bound gives
`||C psi_0,j||^2 <= M a^3/ell^3`. This also holds for entangled initial
flux/native states in that native momentum subspace. For the two-cell
matching route a conservative bound has an additional factor eight.
Therefore, at physical time `t=na`,

\[
 \|\psi(t)-\psi_0(t)\|
 \le4|\lambda(a)|t\sqrt{M/\ell^3}\,a^{3/2}
\]

for the Dirac route. The same bound controls the total variation of any
common complete final detector. In particular, bounded lambda produces
the free records in this smooth, fixed-particle, fixed-volume window.
A nonvanishing interaction must escape these hypotheses, for example
through coupling scaling, populated density/reference or an extended
interaction range. The necessary failure of
`|lambda(a)| a^(3/2)->0` within this bound is not a sufficient prescription
for a renormalized interacting theory. Dense seas, cutoff-scale packets
and unbounded native number are outside this few-particle estimate.

## 6. A transformation test for the retained vertex

The added native memory supplies a physical-spin candidate for the free
Clifford action. The actual link-dependent vertex is independent of that
spin and preserves the two flux routing lanes. In one owner-lane channel
its spin factor is `P_lane tensor I4`. Under the standard Dirac spinor
boost `S(eta)=exp(eta alpha_i/2)`, a Lorentz-scalar bilinear must use a
matrix M satisfying \(S^\dagger M S=M\). The density matrix I fails this test;
beta satisfies it. For the two-spinor vertex, the infinitesimal test yields

\[
 G_i=\tfrac12\{\alpha_i,P_{lane}\}\otimes I
       +P_{lane}\otimes\alpha_i.
\]

The executed Clifford matrices give \(\operatorname{Tr}(G_i^\dagger G_i)=16\), whereas the scalar
`beta tensor beta` has zero infinitesimal residual. Thus this nonzero
density-type vertex is not a scalar interaction under that field dictionary.
A covariant current contraction, an explicitly scalar field coupling or a
different proved field dictionary is additional architecture. Adding beta
alone does not supply a covariant dynamics for the coefficient field.

The packet consequently supplies complete finite native motion, a constructive
common-speed free kinetic alternative and controlled finite interaction.
It also gives exact bulk, transformation and weak-coupling boundaries. It
does not identify Standard Model particles, establish an interacting
relativistic vacuum, convert a source frequency to GeV, or equate added
spin memory with the primitive internal representation.
The flux factor retains the complete 121-state tensor-lane rule. A separate
graded flux field and its filled negative-band reference require an explicit
transfer; neither is supplied by adding CAR propagation to native registers.
