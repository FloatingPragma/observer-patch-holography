# A local strange source and an open energy shell

The complete two-cell strange action supplies calibrated K and Lambda
energies and a nonzero native nonleptonic vertex. Its global quark–antiquark
count/spin terms prevent cluster additivity: a separated baryon and meson
cannot be assigned the sum of their independent masses under that action.
The source here is an explicit finite-range extension. Its bulk bands have
an energy-conserving proton/pion shell. No outgoing weak-decay theorem or
lifetime is inferred from that kinematic result.

## Full native finite-range action

For the complete two-cell Hamiltonian h_xy, let h_x and h_y be its exact
vacuum restrictions to one cell. Extend its matter terms as
sum_x h_x+sum_edges(h_xy−h_x−h_y). The colour-electric term is the parent
Gauss action (U−kappa/2)sum C_x²+kappa sum E_e². On a two-cell singlet its
extra electric contribution cancels identically, preserving the full
calibrated operator. The degree coefficients are the inherited exact
coefficients and vanish on two cells.

The construction localizes every global pair/spin contribution from the
two-cell action. The quark–antiquark spin term between distinct cells is
h[3nq_x na_y/2+2S_q,x dot S_a,y], and its exchanged copy. Spin addition
bounds each by h[nq_x na_y−min(nq_x,na_y)]>=0. Strange-number and strange-pair
terms retain their positive coefficients. The charge-conjugate completion is

    M=sum_x min(nq_x,na_x)
      +sum_edges[min(nq_x+nq_y,na_x+na_y)
                 −min(nq_x,na_x)−min(nq_y,na_y)].

Every bond summand is nonnegative. On a path,
M=sum_edges min(nq_x+nq_y,na_x+na_y)−sum_x(deg(x)−1)min(nq_x,na_x)<=N.
Thus its negative contribution is at worst−cN. The inherited full-Fock
count proof reserves−cN/2; paying another−cN/2 leaves the strictly positive
number floor reported in the receipt. Local positive pair interactions,
Gauss terms and the degree correction satisfy the same count inequalities.
This proves stability at arbitrary family capacity and every path length.

The source has no global min(Nq,Na). Interactions of separated colour-singlet
clusters vanish once their supports are more than the interaction range
apart. The intervening Gauss flux is zero. The spectrum of isolated bulk
clusters must be recalculated; a preserved two-cell fit is not automatically
a bulk pole.

## Native Lambda and meson fibres

For uds there is one quark of each flavour. The complete position/spin
carrier on L cells has dimension(2L)^3. The explicit epsilon isometry is
inherited from the native symmetric cube. `lambda_source.py` checks the
entire64-coordinate two-cell uds matrix against the literal native parent.
Spin permutations and the light-isospin exchange generate its exact
operators. The compact span-at-most-one frame has56 coordinates before
spin/flavour reduction, and seven at S=1/2,I=0,Sz=1/2. Exact rational
nullspaces and the complete matrix give the bulk rest interval.

For three separated quarks, the same colour/number estimates as the native
baryon theorem give the two rational Q floors stored in `lambda_source.py`.
The local light-isospin subtraction is bounded by2max(0,eta−Jw) in the
split case. For three occupied cells it is exactly3(Jw−eta)/2. The full
hopping norm is at most6t. Hence the full Lambda ground lies in
[e_P−(6t)^2/(C−e_P,upper),e_P]. Every relative configuration is retained.
The interval is approximately[1.115684628607,1.115684878081] GeV.

For spin-zero charged pions and kaons, all relative separations r in Z
remain. The compact coordinates are r=0,+1,−1. Their complete two-cell
four-position native spin-zero operators are checked against the256-state
native meson source. The two tails have linearly increasing electric
energy. In the bulk, the common degree law adds2dn to every pair and
8dU/3 to split and tail configurations. The returned parameters from
`constants(*degree_coefficients())` already include these corrections;
the independent two-cell control uses the unchanged parameters. The
receipt and regression test retain this ancestry explicitly. Their lower
bound discards no kinetic term without paying its
norm4t. Each tail endpoint couples to P with norm at most2t, so the full
squared off-block norm is at most4t². The ground interval follows from
[e_P−4t²/(C−e_P,upper),e_P]. The pion bulk rest interval is approximately
[0.139596209155,0.139596220909] GeV. The retained shift from its two-cell fit
is a consequence of the common extension, rather than a new fit.

## Positive propagation and the shell test

Apply the same inverse-rest-operator rule as `openwave_baryon_dispersion`
to all compact triquark channels and all compact quark–antiquark channels.
For a quark–antiquark pair the annihilator is an even native CAR polynomial
with its Wilson transport. The sum DeltaC* J DeltaC is positive on full
Fock space. It vanishes on the complete one-cluster k=0 fibre. Meson
annihilators vanish in a pure three-quark sector; triquark annihilators
vanish in a one-pair sector. Both rules therefore coexist on one declared
strong action without changing the quoted isolated rest fibres.

Far-separated cluster actions are additive. The outgoing ground bands have
opposite momenta k and−k and energy E_p(k)+E_pi(−k). Both full Hamiltonian
fibres depend continuously in bounded norm on k, so their ground energies
are continuous. All electric tails and native internal complements remain.
At tan(k/2)=1/6, the upper pair-energy enclosure is strictly below the
lower Lambda-rest enclosure. At tan(k/2)=1/5, the lower pair-energy enclosure
is strictly above the upper Lambda enclosure. The intermediate-value
argument establishes an energy-conserving shell between these momenta.
Alternating rational arctangent sums give the physical-momentum bracket.
This removes both the artificial separated-cluster energy and insufficient
bandwidth obstructions in this declared member.

The parent native weak current is not replaced by a hand-selected decay
projector. Its nonzero two-cell vertex and closed-channel control retain
their original domains. Transfer to interacting incoming and outgoing
scattering states, with the full weak vacuum and all other channels, is
not proved here. A shell crossing therefore receives no claim of an
executed strange lifetime or physical weak-decay validation.
