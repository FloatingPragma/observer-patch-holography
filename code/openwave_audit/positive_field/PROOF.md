# A neutral three-quark composite coupled to a positive Maxwell field

This construction gives an energy bound on breakup without a spectral ground-state preparation. The matter source is the finite three-quark colour chain of `color_chain_native`, with an actual electromagnetic covariant hopping term. Its field is the common LOCAL12 Maxwell field. The construction is an explicit amendment; no equivalence to the fixed-charge Dirac source or the autonomous filled CAR source is asserted. All quantities below are source quantities. Native realization requires the preparation, evolution and observation contract in the sibling implementation package.

## Source and clock

Take the complete three-quark SU(3) singlet with total spin projection one half and weak-up number one on an open path of length L in {2,3,4,5}. Each weak-up constituent has electric charge 2/3 and each weak-down constituent has charge -1/3 in the same elementary unit. The colour epsilon isometry gives the symmetric spin/weak/orbital factor. The SU(3) tree-Gauss image and local and prefix Casimirs have the definitions of the finite chain source. Choose t=1/8, J=1/4, U=1 and kappa in {1/2,1,2}.

Embed the path as distinct vertices of the actual finite LOCAL12 complex. The concrete five-vertex path is [0,1,5,4,2], using oriented edges [(0,+1),(6,+1),(22,-1),(12,-1)] of the inherited 96-vertex,288-edge complex; shorter members use its prefixes. Let B be its oriented incidence matrix and let P(q)=B* (B B*)^+ q be the longitudinal edge field for a neutral occupation charge vector q. Orthogonal coexact and harmonic coordinates carry the transverse field. With a=1/2, b=3 and the inherited canonical field normalization, use the vacuum-subtracted coexact oscillator Hamiltonian and the noncompact harmonic kinetic energy:

`F = sum_j omega_j N_j + (1/4) sum_(r=1)^3 p_r^2 >= 0`.

The commutators are [q_j,p_k]=i delta_jk. Normal ordering subtracts the same finite coexact vacuum constant in every matter sector; it changes no source evolution or energy difference. There is no normalizable zero-momentum harmonic ground vector.

Pin h_star=1/2. The interaction uses the calibrated coefficient `g=sqrt(12 pi alpha_input)`, with `alpha_input=72973525643/10^13`. The stability proof also holds for any real g. The longitudinal term is

`V_C(q) = (a g^2 / (2 h_star)) ||P(q)||^2 >= 0`.

Every allowed constituent hop has the actual unitary Weyl factor of its charge times g and its LOCAL12 route. Fix B_e=delta_old-delta_new. A right hop of charge q_const changes density by -q_const B_e; under B E=g q its transverse Weyl factor is `W=exp(-i g q_const J_e dot A)` and the longitudinal particular field changes by `-g q_const P_div J_e`. The negative elementary charge q_const=-1 therefore has the inherited +g J route sign. These factors implement the charge change and preserve the full Gauss constraint under the inherited particular/transverse field dictionary. Their actual coefficients and generator order are implementation obligations; an arbitrary hopping phase does not establish this dictionary.

Write K_g for the positive covariant path Laplacian on the three constituents. Each edge contribution has the quadratic form `sum_s (a_y,s - W_s a_x,s)* (a_y,s - W_s a_x,s)` in the full matter/field space, restricted to the fixed three-particle domain. The W_s are unitary and preserve weak species. Colour dressing uses the same exact tree-Gauss isometry. Thus K_g is positive even if different W_s fail to commute. The complete source is

`H = 3 + K_g/8 + (1-kappa/2) sum_x C_x^2 + kappa sum_e E_e^2 + (1/4) sum_x S_x^2 + F/h_star + V_C`.

The finite matter terms and Weyl hopping are bounded. F is self-adjoint on its oscillator/free-particle domain. The displayed sum is self-adjoint on that domain by bounded perturbation and is bounded below. Source time is tau in exp(-i tau H); the native reference clock reads relative phase tau/8. The field evolves for tau/h_star cell-time units. Pulse counts do not define this clock. Under the calibrated unit map, tau0=2T_g and one source-energy unit is h_star U_g; one field-cell energy unit is U_g. The chosen colour coefficients are architectural inputs. An empirical baryon energy is not used to choose them, and the endpoint mean energy is not identified with an empirical baryon rest mass.

## Block bounds unaffected by the added field

Let P project onto all three constituents at one path site, including all allowed spin/weak states and every field vector. Let Q=1-P. Neutrality makes the charge vector q identically zero on P, so V_C P=0. Both field energy and the Coulomb term commute with P. The colour and spin terms also preserve P. Only covariant hopping connects P to Q.

On P, the smallest local spin is one half. The path degree is at least one, giving

`PHP >= a0 P`, with `a0=3+3/16+3/8=57/16`.

For a separated configuration of span s>=1, the local colour cost is at least 8/3, the spin cost at least 3/16 and the tree electric term is 4 kappa s/3. The combination of local and tree colour terms is at least `8/3+(4 kappa/3)(s-1)`. Covariant K_g, F and V_C are positive, so

`QHQ >= (281/48) Q`, and `281/48-a0=55/24`.

The cross-block norm remains bounded by `b0^2=3/32`. To see this without treating field coefficients as scalars, start with a colocated occupation vector at x. A single hop gives a two-at-x/one-at-y pattern. That pattern identifies its source x, destination y and the moved spin/weak label; moving the particle back identifies the original occupation. Distinct such outputs are orthogonal. On returning along the same hop, W_s* W_s=1 exactly, and the squared occupation coefficient is n_s. Thus, on the whole colocated subspace with arbitrary field/reference entanglement,

`P H Q H P = (3/64) deg(x) P`

blockwise in the original colocated occupation labels. Since deg(x)<=2, `||QHP||^2<=3/32`. This identity requires three constituents: the two-at-source pattern distinguishes the direction. It is not obtained by replacing the Weyl factor with its vacuum expectation. Adding code padding at a diagonal energy at least281/48, with no cross terms, preserves these bounds.

For a neutral external probe let xi(x) be any real function with adjacent path increments at most one in absolute value. Put `D=sum_i q_i xi(x_i)`. Neutrality gives D P=0 and `|D|<=2s/3`. If |f|<=1/4, then

`Q(H+fD)Q >= [281/48-1/6]Q = (91/16)Q`.

Indeed its span-dependent lower bound is `281/48+(4 kappa/3)(s-1)-(2|f|/3)s`, minimized at s=1. The same a0 and b0 bounds hold, with separated excess `d=91/16-a0=17/8`. A linear spatial potential is included when its restriction to the embedded path has the stated increment bound.

## Explicit finite-energy initial state

At one endpoint prepare the normalized local spin-half vector. In the neutral orbital basis it is `(|0,1,3> - sqrt(2)|1,1,2>)/sqrt(3)`, with normalized symmetric occupation vectors. Prepare the coexact oscillators in their vacuum and each of the three harmonic pairs in a pure centered Gaussian with

`Var(p)=1/1000`, `Var(q)=250`, `Cov(q,p)=0`.

Their uncertainty product is exactly1/4. These are normalized finite-energy states on the noncompact harmonic coordinate, with no assumed momentum eigenvector. The initial matter state is colocated, so its Coulomb and external-probe energies vanish. Off-diagonal hopping has zero expectation; its endpoint diagonal contribution is3/8. The total initial energy is

`E_initial = a0 + delta`, where `delta = 3/(4000 h_star)=3/2000`.

No spectral projection, cooling process, reflected mixture or ground-state fidelity estimate is needed for this initial state. Its charge neutrality is essential to the displayed initial Coulomb energy.

## Uniform breakup theorem

Let rho(tau) evolve under H+fD from that state, and let `q_sep=Tr(Q rho(tau))`. The same argument applies to mixtures or external references. Block Cauchy-Schwarz and energy conservation give

`delta >= d q_sep - 2 b0 sqrt(q_sep(1-q_sep))`.

For no external probe, d=55/24. If q_sep>=7/100, then

`2 b0 sqrt((1-q_sep)/q_sep) <= sqrt(279/56) < 9/4`.

The right-hand side of the energy inequality would exceed `q_sep/24>=7/2400>3/2000`, a contradiction. Therefore **q_sep<7/100 for every source time** under H.

For every |f|<=1/4, d=17/8. If q_sep>=2/25, then

`2 b0 sqrt((1-q_sep)/q_sep) <= sqrt(69)/4 < 133/64`.

The last inequality is certified by133^2>69*256. Consequently the energy right-hand side exceeds `3 q_sep/64>=3/800>3/2000`, another contradiction. Therefore **q_sep<2/25 for every source time** under H+fD. This leaves a probability margin of1/50 below the breakup gate1/10.

For comparison, zero initial excess gives the exact bound `q_sep<=4 b0^2/(d^2+4 b0^2)`:216/3241 without a probe and24/313 with a probe. This zero-excess statement is an abstract energy case; it does not assert a normalized zero-energy state of the noncompact harmonic sector.

## Observation transfer and boundary

The breakup effect is the complete occupation projector Q, with failed preparation, reset, export, invalid code and unreadable outputs assigned to the bad outcome. If a native finite-time experiment approximates the displayed complete preparation/evolution/occupation instrument within total variation epsilon, its bad-outcome probability is at most2/25+epsilon for the weak-probe domain. A bound epsilon<=1/1000 leaves it below81/1000<1/10. This comparison uses norm/channel control of the prepared Gaussian and the evolution; it does not infer small energy error from trace distance for an unbounded field Hamiltonian. A uniform source-time theorem and a separately priced program for each finite time window are distinct claims.

The theorem depends on positive field and Coulomb energies, unitary routed hopping, the exact colour Casimirs, and the neutral endpoint preparation. A resonantly lowered separated sector invalidates the bound and permits complete transfer in a two-state control. Charged proton preparation, derived electromagnetic strength, empirical masses, spontaneous formation, deuteron binding, a relativistic interacting continuum and a common source satisfying all31 criteria do not follow from this energy estimate. The explicit noncompact Gaussian/native route implementation supplies the remaining observation-transfer obligation for this construction.
