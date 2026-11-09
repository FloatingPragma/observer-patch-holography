# Cubic and parity channels with original spectral weight

The actual primary member zero has a populated scalar-odd pair channel
and nine other populated cubic/parity isotypic channels. All thirty-two
coherent spatial Walsh and internal pair records also have positive
original weight in an explicit positive transfer band. These statements
resolve operator channels of a finite interacting response. They do not
isolate a pole, identify a continuum spin or supply a particle rest mass.

## 1. Reused source results and their domains

Fix [the primary contract](../particle_common_contract/PROOF.md), its
member \(n=0\), spacing \(a=1\), eight cells, complete 60-mode matter
factors, quantum links, Gauss restriction, polynomial state \(\rho\),
idle conjugate reference and original source clock. The interaction
coefficient is one. [The occupation estimate](../particle_primary_occupation_bound/PROOF.md)
proves \(\langle D\rangle_\rho<5\), where \(D\) counts all occupied
positive-beta modes and holes in negative-beta modes. Write

\[
 w=1-\frac{\langle D\rangle_{\rm upper}}{16}
 >\frac{11}{16}.
 \tag{1}
\]

[The local pair response](../particle_primary_response_moments/PROOF.md)
supplies the full-Fock residual inequality
\(\|[H,c^\dagger(p)c(m)]-6c^\dagger(p)c(m)\|^2\le3\)
for any normalized positive- and negative-beta vectors at a corner.
Its proof retains nonorthogonal neighbor images, the product link
projection and every matter occupation. Its complete energy apparatus
is reused here through [joint spectroscopy](../particle_primary_joint_spectroscopy/PROOF.md).

[The physical cube action](../particle_primary_channel_labels/PROOF.md)
acts on the interior graph after the one-link ghost Gauss constraints
force dangling links to constants. A geometric action on the asymmetric
ambient owned-link allocation is not assumed. Only the exact binary cube
enumeration and characters from `cubic_particle_channels/group.py` are
reused; that package's comparison Hamiltonian and spectra are absent.

`Lean/QFT/PeriodicMomentum.lean`, `polynomial_conserved`, propagates a
supplied relation \([H,T]=0\); it supplies no translation of this open
source. `Lean/QFT/GroundRecordTail.lean`, `downwardWeight_le_excess_div`,
and [the mass readout](../particle_primary_mass_readout/PROOF.md) retain
original weights when converting source transfers to final excitation.
These finite results do not assume a unique vacuum or a continuum field.

## 2. Four exact local pair components

Let \(p_1,p_2\) and \(m_1,m_2\) be the actual beta eigenvectors in the
response proof, assembled into isometries \(P,M:\mathbb C^2\to\mathbb C^4\).
The source Clifford chart fixes
\(\sigma_j=-iP^\dagger\gamma_{j+1}\gamma_{j+2}P\), with cyclic indices,
and \(\sigma_0=I_2\). These are an orthonormal Pauli basis after division
by \(\sqrt2\). At cell \(x\), in the one-dimensional native \(Q_6=6\)
block, define

\[
 C_{xa}=\frac1{\sqrt2}\sum_{u,v=1}^2
       (\sigma_a)_{uv}c^\dagger(p_{xu})c(m_{xv}),
 \quad a=0,1,2,3.
 \tag{2}
\]

Every \(C_{xa}\) is even, gauge invariant and preserves all five native
occupation counts. It is a contraction on the full four-mode Fock factor,
including the one-, two- and three-particle sectors. Its effect eigenvalues
are \(0,1/2,1\). If \(D_x^{(6)}\) counts the four defects in this native
block, direct CAR algebra gives

\[
 [C_{xa}^\dagger,C_{xa}]=I-\frac12D_x^{(6)},\qquad
 C_{xa}^\dagger C_{xa}\ge I-\frac12D_x^{(6)}.      \tag{3}
\]

The second inequality follows by adding the positive operator
\(C_{xa}C_{xa}^\dagger\). Thus each separate Pauli component, averaged
uniformly over the eight cells, has original success weight at least
\(w\). No equality of individual site probabilities is required.

The coefficient \(\sigma_a/\sqrt2\) has two singular values
\(1/\sqrt2\), so each residual is a sum of two rank-one pair residuals.
The inherited all-occupation bound and the triangle inequality give

\[
 R_{xa}=[H,C_{xa}]-6C_{xa},\qquad \|R_{xa}\|^2\le6.
 \tag{4}
\]

For each spatially averaged Pauli measure, the centered second moment
at frequency six is therefore at most six. Its original transfer weight
in \([2,10]\) is at least \(w-6/16>5/16\). This spatially averaged
instrument is distinct from a coherent zero-wavevector component.

## 3. Physical rotations and inversion

The source spin lift obeys
\(S(q)\beta_D S(q)^\dagger=\beta_D\) and rotates \(\gamma_j\) by the
proper cube matrix \(R(q)\). The literal Clifford matrices give

\[
 U_q C_{x0}U_q^\dagger=C_{qx,0},\qquad
 U_q C_{xa}U_q^\dagger=\sum_{b=1}^3R(q)_{ba}C_{qx,b}.
 \tag{5}
\]

Spatial inversion maps each interior cell to \(1-x\), reverses link
orientations and acts on matter by \(\beta_D\). Haar measure, the word
alphabet and the tensor cutoff are preserved under link inversion. The
literal hopping matrix satisfies \(\beta_D S_j\beta_D=S_j^\dagger\).
The electric Casimir and real complete plaquettes are invariant, and
the Gauss constraints are permuted. This defines an exact physical
involution commuting with \(H\) and \(\rho\). It sends

\[
                       C_{xa}\longmapsto-C_{1-x,a}.          \tag{6}
\]

The internal pair components consequently transform as \(A_{1u}+T_{1u}\).
The label \(u\) denotes odd inversion parity. The reference is idle;
state covariance concerns the source density. The canonical purification
has the corresponding covariance under the source unitary tensored with
its conjugate, without a source-only invariance claim.

The eight vertex coefficients decompose as
\(A_{1g}+T_{1u}+T_{2g}+A_{2u}\). Their tensor product with (5)-(6) gives
the complete thirty-two-dimensional pair coefficient space:

\[
 \bigoplus_{\pi=g,u}
       \left(A_{1\pi}+A_{2\pi}+E_\pi+2T_{1\pi}+2T_{2\pi}\right).
 \tag{7}
\]

The ranks of its ten isotypic projectors are therefore
\(1,1,2,6,6\) for each parity. `model.py` constructs each projector by
the exact finite character sum, checks orthogonality and their sum to the
identity, and constructs an orthonormal adapted coefficient basis.

These coefficients define independent physical operators. Applied to the
filled-negative-beta, constant-link physical anchor, the thirty-two
operators in (2) have orthonormal particle-hole images. The anchor is only
an injectivity witness; it does not replace the selected state.

## 4. Coherent records and their normalization

Index all local components by \(i=(x,a)\), so \(N=4V=32\). For each
contraction declare the complete two-column source service

\[
 C_{0i}=C_i,\qquad C_{1i}=\sqrt{I-C_i^\dagger C_i}.
\]

It is the blank column of a finite contraction dilation, implementable
under the existing exact algebraic apparatus service. For any specified
unitary coefficient matrix \(W\), a uniform neutral address followed by
the controlled columns and \(W\) gives

\[
 K_{bj}=\frac1{\sqrt{32}}\sum_iW_{ji}C_{bi},\qquad
 \sum_{b,j}K_{bj}^\dagger K_{bj}=I.                \tag{8}
\]

All thirty-two success and thirty-two failure records are retained,
together with arbitrary source-reference entanglement. Failure maps are
the ones in (8); a different coarse failure instrument is not identified
with them. Exact root service specification does not assert efficient
compilation or execution of the full spatial histogram.
The isotypic classification below concerns the success operators. The
complementary square roots need not transform in that coefficient
representation, so their retained labels do not assert source irreps.

The spatial Walsh choice is
\(W_{(k,a),(x,b)}=8^{-1/2}(-1)^{k\cdot x}\delta_{ab}\),
\(k\in\{0,1\}^3\). Its success Kraus operator is
\(K_{0,k,a}=\frac1{16}\sum_x(-1)^{k\cdot x}C_{xa}\).
For every one of these thirty-two records, (3) gives

\[
 \operatorname{Tr}\rho[K_{0,k,a}^\dagger,K_{0,k,a}]
     =\frac1{32}\left(1-\frac{\langle D^{(6)}\rangle}{16}\right)
     \ge\frac w{32}>\frac{11}{512}.              \tag{9}
\]

This lower-bounds its original success weight. An apparatus selecting one
Pauli component in advance would instead have factor \(w/8\). Equation
(9) includes all four components in the same complete instrument and
does not divide out their relative probability.

The alternate choice of \(W\) is the exact cubic/parity adapted basis.
Let \(P_\lambda\) be its isotypic coefficient projector of rank \(d\).
For each cell the exact projector satisfies

\[
 \frac12\sum_{a,b=0}^3(P_\lambda)_{xa,xb}\sigma_a^\dagger\sigma_b
                         =\frac d{16}I_2.         \tag{10}
\]

It follows either by the finite character sum or by vertex transitivity,
the real symmetric projector and its separate scalar/vector blocks.
The code checks (10) at every cell with exact rational entries. Contracting
the CAR before taking any state expectation yields the original
isotypic commutator lower bound

\[
 \chi_\lambda:=\sum_{j\in\lambda}
  \operatorname{Tr}\rho[K_{0j}^\dagger,K_{0j}]
       \ge\frac{d w}{32}>\frac{11d}{512}.         \tag{11}
\]

The corresponding weight is an aggregate over the indicated isotypic
subspace. In particular the rank-six vector channels retain two
multiplicity spaces. No lower bound on each multiplicity vector is
inferred. The rank-one \(A_{1u}\) projector is exactly the coherent
zero-Walsh scalar mode. The three zero-Walsh vector components form a
proper subspace of the rank-six \(T_{1u}\) isotypic component.

## 5. A positive-frequency band in every channel

On the physical source spectrum the selected density is

\[
 p(E)=Z^{-1}T_m(x(E))^{2s},\quad
 x(E)=\frac{\beta(U-E)}{2s}\ge0,\quad
 T_m(x)=\sum_{j=0}^m\frac{x^j}{j!}.
\]

Its logarithmic slope obeys the exact polynomial identity

\[
 0\le-\frac{d}{dE}\log p(E)
       =\beta\frac{T_{m-1}(x(E))}{T_m(x(E))}\le\beta.          \tag{12}
\]

This is a bound for the actual polynomial state, without an exact Gibbs
substitution. For any success subinstrument define its original positive
transfer measure
\(\mu=\sum_{j,i,k}p_i|(K_j)_{ki}|^2\delta_{E_k-E_i}\).
The corresponding commutator is

\[
 \chi=\sum_{j,i,k}(p_i-p_k)|(K_j)_{ki}|^2.
 \tag{13}
\]

Negative transfers contribute at most zero to (13), and zero transfers
contribute zero. For \(0<\omega<\delta\), (12) implies
\(1-p(E+\omega)/p(E)\le\beta\delta\). If \(\mu(\mathbb R)\le M_0\)
and \(\int(\omega-c)^2d\mu\le M_2\), then for \(R>\max(c,\delta)>0\),

\[
 \mu([\delta,R])\ge
     \max\left(0,\chi-\beta\delta M_0-
                            \frac{M_2}{(R-c)^2}\right).     \tag{14}
\]

The inequality discards the negative contribution in the upper bound on
\(\chi\), in the direction required to infer a positive-band lower bound.
It retains the original measure throughout. `passive_window` implements
(14) for reuse by authenticated source consumers; its arithmetic alone
does not prove their commutator or moment premises.

For either complete coefficient transform in (8), Parseval's identity and
(4) give the total success centered moment

\[
 \sum_j\operatorname{Tr}\rho([H,K_{0j}]-6K_{0j})^\dagger
                       ([H,K_{0j}]-6K_{0j})\le6.              \tag{15}
\]

Every subinstrument has total mass at most one and centered moment at
most six. Choose the common explicit window

\[
 \delta=\frac{11}{2048\beta},\qquad R=40.
 \tag{16}
\]

The low-frequency allowance is \(11/2048\); the high-frequency allowance
is \(6/34^2<11/2048\). Equations (11) and (14) prove

\[
 \boxed{\mu_\lambda([\delta,40])\ge\frac{11d}{1024}}.
 \tag{17}
\]

The stronger arithmetic value is retained in the receipt. Equation (9)
gives the same rank-one bound \(11/1024\) for each coherent Walsh/pair
record, including zero wavevector. These are substantial original
probability floors independent of the exponentially small faithful-state
floor. Their lower endpoint depends on the actual cooling parameter.

The selected excess is \(\epsilon=1/16\). Reusing original-weight
excitation transport with initial-excitation cutoff eight loses at most
\(1/128\). Each isotypic channel therefore has final ground-relative
excitation in \([\delta,48]\) with original weight at least
\((11d-8)/1024>0\). This inference requires no numerical ground-energy
subtraction and supplies no rest-frame momentum.

## 6. Complete finite clocks and physical meaning

Use the unchanged source-only clock plan at resolution \(\eta=\delta/4\),
with each complete tail and each complete implementation allowance at
most \(1/4096\). The two clocks together lose at most \(1/1024\), including
both tails. The complete instrument has \(64J^2\) fine records. Initial
energy readout preserves the stationary source marginal and changes its
joint purification; the initial pointer and reference are retained.

Its reported transfer window is
\([\delta/2,40+\delta/2]\). The original isotypic reported probability
is at least \((11d-1)/1024\). Any additional probe-compiler or sampling
event error must be deducted separately. The receipt specifies a finite
clock plan and proves its probability floor; no full spatial matrix,
compiled tape or primary histogram is executed.

The spatial phases in (8) are finite-aperture Walsh labels. The open
Hamiltonian has no proved commuting lattice translations, and the label
zero does not establish a physical rest frame. Cubic irreps and parity
classify operators. A scalar operator in a stationary mixed state can
connect equal non-scalar state irreps and different multiplicity spaces.
The result does not assign continuum spin zero or one, an isolated
energy pole, a unique threshold, a dispersion law or a named species.

A first finite-capacity GeV mass comparison needs a resolved same-state
branch or experimentally defined threshold/resonance, an operational
momentum and initial-momentum bound, a controlled energy-momentum relation
on its stated domain, and complete prediction/readout errors. The owning
mass-readout and joint-spectroscopy proofs supply conditional interval
transport, without those missing spectral premises. An exact continuum
claim additionally needs fixed physical band survival and the joint field
limit. Here \(\delta\) scales as \(1/\beta\), and the occupation input is
authenticated only at member zero.

Unit conversion reuses `readout.declared_SI_comparison`. The receipt gives
coefficients for transfer energy in GeV times the declared seconds per
source tick; the tick itself is unspecified. This is an energy conversion,
without an empirical mass input or an implicit rest-mass assignment. The
primary physical-mass table in GeV is empty.
