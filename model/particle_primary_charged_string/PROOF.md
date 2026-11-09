# Physical charged endpoints with original energy weight

This construction uses member zero of `oph.open_wilson_gauss_observer.v1`:
the eight-cell open cube, its complete physical carrier, polynomially
prepared state, idle conjugate reference and source Hamiltonian clock.
All products are parent products before word compression. The finite
controls below check the operator algebra; they do not execute the full
primary-state histogram.

## 1. Gauge-constrained creation, without changing the preparation

The native one-dimensional weak/color singlet has integer hypercharge
label `Q6=6` and trace-normalized Casimir `3/10`. On an interior link
oriented from `y` to `x`, let `chi(U)` be its character, with the orientation
chosen to compensate the endpoint matter transformation. For either of
the two matched beta-spin components define

\[
 B_{xy\alpha}=c^\dagger(p_{x\alpha})c(m_{y\alpha})\chi(U)P_{0,xy}.
\]

`P0` is the Haar constant-link projection. For the reverse orientation use
the conjugate character. There are 12 interior links, two directions and
two components, hence 48 such probes. Each preserves Gauss support and
every global native occupation. Each has even CAR parity. Their matter
charge selection rules are

\[
 [Q_6(z),B_{xy\alpha}]=6(\delta_{zx}-\delta_{zy})B_{xy\alpha}.
\]

These are opposite endpoint increments, not an identification with an
electron or with measured electromagnetic charge. In particular, an even
pair-creating operator is not itself a single fermionic particle.

The complete instrument has 48 successes `K_b=B_b/sqrt(48)` and one
failure `K_f=(I-sum_b K_b^*K_b)^(1/2)`. Indeed,

\[
 B_b^*B_b=P_{0,xy}n_{m_y}(1-n_{p_x})\le I.
\]

All failures and the reference are retained. The link projection is part
of this intervention; the original state is not replaced by a state with
constant links. Character degree one and the constant character are in
the member-zero word carrier, so `B` commutes with its word projection.
Consequently its compressed Hamiltonian commutator is the compression of
the parent commutator, without an omitted cutoff boundary term.

A complete local-charge experiment can measure the initial endpoint
charges, apply this instrument and measure the final endpoint charges.
The displayed selection rule holds for every successful history. The
initial charge measurement generally disturbs the selected state. The
energy bounds in the following sections apply to the original state,
not to its charge-dephased replacement. They are separate apparatuses.

## 2. Full-space original success probability

Retain fraction `r=7/10` of the native-six onsite penalty and fraction
`s=8/45` of the electric penalty. These are allocations in an operator
inequality, not changes to the source Hamiltonian. For the five native
factors with dimensions `d_j`, Casimirs `C_j` and charges `q_j`, put

\[
 m_j=1-r\,1_{q_j=6},\qquad
 \Delta=4\sum_jd_j\bigl(\sqrt{m_j^2+1/4}-m_j\bigr),\qquad
 w_j=\bigl(105(1-s)C_j/2-\Delta\bigr)^{-1}.
\]

The smallest nonconstant-link margin is strictly positive. The
[full-Fock Schur bound](../particle_primary_occupation_bound/PROOF.md),
with the [exact native-six block bound](../particle_primary_refinement_response/PROOF.md),
gives

\[
 H\ge-720-\mathcal L+\frac{21}{10}D_6+
                  \frac{28}{3}\sum_\ell C_\ell,
 \qquad \mathcal L=12\sum_jc_j(m_j,w_j)+\frac67<10.612.
\]

Here `D6` counts native-six particles above, and holes below, the beta
reference. The `c_j` are bounds for the entire endpoint Fock block, not
just its vacuum vector. Electric terms on ghost links are nonnegative
and vanish on their physical Gauss support.

For each probe, the three commuting missing-success projections satisfy

\[
 1-B_b^*B_b\le n_{p_x}+(1-n_{m_y})+(1-P_{0,xy}).
\]

Every endpoint component occurs three times in the 48-term average, and
each link occurs four times. Since the first nonconstant Casimir is at
least `3/10`,

\[
 1-\frac1{48}\sum_b\langle B_b^*B_b\rangle
 \le\left\langle\frac{D_6}{16}+\frac5{18}\sum_\ell C_\ell\right\rangle.
\]

Both terms are controlled with the same multiplier
`k=48r=189s=168/5`. The unchanged selected state has
`<H> <= U0+1/16`, where the independently replayed physical variational
bound is `U0 < -723.7678199426`. Thus the original total success weight is

\[
 w\ge1-\frac{U_0+1/16+720+\mathcal L}{168/5}
       >0.79446654>\frac{79}{100}.
\]

No symmetry of the prepared state is needed for this uniform average.

## 3. Centered energy moment on the same physical source

The exact onsite and electric commutator is

\[
 [H_{\rm onsite}+H_E,B_b]=\left(6+\frac{63}{4}\right)B_b
                         =\frac{87}{4}B_b.
\]

The electric identity uses `C chi=(3/10)chi` and `C P0=0`.
It is a reference transfer frequency, not a measured rest mass.

For the selected bond let `T` contain all five native hopping factors.
Haar orthogonality gives `P0 T P0=0` and
`P0 T^2 P0=sum_j Q_j`. The full-Fock quartic bounds
`||Q_j||<=2(d_j+1)` sum to 40. Hence `||T P0||<=sqrt(40)`.
Scalar multiplication by `chi` commutes with the parent `T`; with
`Pchi=chi P0 chi^*` the same bound holds for `T Pchi`, and
`Pchi T Pchi=0`. In `[T,B]=TB-BT`, the first summand has input in `P0`
and output orthogonal to `Pchi`, while the second has input orthogonal
to `P0` and output in `Pchi`. Both input and output supports are
orthogonal, so the norm of their sum is at most `sqrt(40)`, without a
factor of two. Compression cannot increase this bound.

Among the other bonds only the four incident native-six bonds fail to
commute with the matter pair. Write their one-body action on its two
modes as `u=h_rest p_x` and `v=h_rest m_y`. On the eight-cell cube these
vectors each have norm one, have disjoint support and are orthogonal to
both original modes. The CAR commutator is the sum of two disjoint
normalized raising pairs, with exact Fock norm `sqrt(2)`. This argument
uses the actual Wilson rank-two hopping matrices, whose column norms
are `1/sqrt(2)` on each of the two remaining incident bonds. It is not
a bound obtained from an extensive global Hamiltonian norm.

Only two magnetic plaquettes meet the link. Each magnetic term has norm
at most `1/7`, so their commutator norm is at most `4/7`. Therefore

\[
 \|[H,B_b]-(87/4)B_b\|
 \le\sqrt{40}+\sqrt2+\frac47<\sqrt{70}.
\]

For the stationary selected density matrix, the original two-energy
transfer measure averaged over the 48 probes consequently obeys

\[
 \int(\omega-87/4)^2\,d\mu(\omega)\le70.
\]

The projection and selected-link estimates hold on the entire physical
carrier. They do not substitute the variational seed for the selected
state.

## 4. Finite complete records and their meaning

Chebyshev's inequality with radius 12 gives

\[
 \mu([39/4,135/4])\ge w-70/144>3/10.
\]

The existing finite ground-record tail theorem, with initial excess cut
one and selected mean excess at most `1/16`, puts original weight at
least `19/80` in final excitation energies `[39/4,139/4]` above the same
source ground. This is an inference about a final excitation, not an
absolute ground-energy measurement.

The complete source-clock plan uses resolution `1/16`, tail allowance
`1/128` per clock and implementation allowance `1/128` per clock.
Initial energy, probe result, final energy and reference are all retained:
there are `49 J^2` fine outcomes. The total event loss is at most `1/32`,
and the original recorded transfer in `[77/8,271/8]` has weight greater
than `1/4`. No unsuccessful outcome is conditioned away.

The construction supplies a physical oppositely charged endpoint pair
and a populated energy band from the configured source. It supplies no
unique spectral line, momentum dispersion, isolated charged sector,
spin/statistics identification, continuum particle or GeV rest mass.
The link-supported pair response can be used as a seed in a separate
translation-resolved apparatus; its single-member weight is not a
uniform-capacity lower bound for that different apparatus.
