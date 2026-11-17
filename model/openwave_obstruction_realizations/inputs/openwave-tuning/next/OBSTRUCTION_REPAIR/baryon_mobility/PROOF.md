# A complete composite-motion term from the unchanged quark hopping

The frozen shared-auxiliary source has all24 CAR modes on two cells,
both25-state auxiliaries, the total colour Gauss constraint, and inherited
t=g=1/8. Define H0 as its electric and auxiliary-gap terms. P denotes both
auxiliary vacua and local quark colour singlets on both cells. Its dimension
is92 squared=8464. This calculation derives a term on that complete P
space; it does not identify P states with the full interacting vacuum.

The consumed local frames are independently checked against all eight
colour generators. Their Gram matrices are exactly invertible, and the
local dimensions1,20,50,20,1 at N=0,3,6,9,12 exhaust the local singlet
representations. No spin, isospin or local occupation is selected out.

## The coefficient and its actual denominators

Write K=sum_mu(c_Lmu^dagger c_Rmu+h.c.) and H_hop=-t K. A single hop from
P has colour-electric energy Delta=4kappa/3. This is an actual eigenvalue
on the complete one-hop image, as verified in the shared-auxiliary
denominator packet. Consequently the off-diagonal third-order term is

    P H_hop Q H0^-1 Q H_hop Q H0^-1 Q H_hop P
      = -t^3/Delta^2 P K^3 P.

Both outside resolvents act on one-hop images. Possible other colour
representations in a two-hop expansion cannot alter this identity: only
their component that can return through one hop to P contributes.
There is no colour-singlet one-hop P-to-P matrix element. The subtracted
third-order normalization terms involving P V P cannot transfer three
quarks, since P V P here is the diagonal mass operator. Mixed three-factor
terms containing auxiliary couplings or the mass do not transfer three
quarks between local singlets. There can be additional diagonal terms at
the same perturbative order; they are not omitted from a claimed full
H_eff3, because no such claim is made.

## All triple operators and the fermion sign

Let T=sum_mu c_Lmu^dagger c_Rmu. Three hops ending in local singlets must
all go in the same direction. Thus P K^3 P=P(T^3+(T^dagger)^3)P. The CAR give

    T^3=6 sum_(i<j<k) c_Li^dagger c_Lj^dagger c_Lk^dagger c_Rk c_Rj c_Ri.

All220 mode triples are constructed directly. Exactly64 have nonzero
local singlet projections; they collapse, including their signs, into
20 symmetric flavour triples. The multiplicities are1,3 or6 and sum to64.
For local projected creation B_I^dagger and annihilation B_I, and local
fermion parity Pi=(-1)^N, the ordinary product-basis expression is

    P K^3 P=6 sum_I [B_I^dagger Pi tensor B_I + Pi B_I tensor B_I^dagger].

The parity is essential. These composites are odd in the quark CAR;
replacing them by independent commuting site operators changes the source.
An independent raw 24-mode CAR control starts with three quarks of one
flavour on the left and three of another on the right. The matrix element
of T cubed into the six-quark left state is +6. Dropping the parity in the
fixed local product frame gives -6. This is a rejecting control for those
fixed preparation and probe maps, not a claim that a simultaneous change
of every state and operator convention would change the physics.
The exact local matrices and every triple sign are retained in
`LOCAL_TRIPLES.json`. The full8464-dimensional rational transfer matrix
has denominator36 and is saved with its product Gram form. Exact tests
verify physical self-adjointness, total-number conservation, a change of
three in the left quark number, and covariance under spin and isospin
raising and diagonal generators.

The complete number-sector dimensions are

    N:         0,  3,   6,    9,   12,   15,  18, 21, 24
    dimension: 1, 40, 500, 2040, 3302, 2040, 500, 40,  1.

The exact left-baryon-current operator is i[H_transfer,N_L/3]. It is
nonzero and retained through the same full matrix. Coordinate row-sum
bounds in `RESULT.json` bound the physical spectral norm because the
operator is self-adjoint in the positive Gram metric. They are conservative
bounds, not claimed dispersion curves or physical velocities.

## Exact complete hopping control

At N3 on the empty reference, P K^3 P has square36 I, so the derived
one-baryon transfer amplitude is

    t_B=6t^3/Delta^2=27/(4096 kappa^2).

The N21 hole control gives the same exact square. These identities do not
make the empty reference the actual N12 vacuum.

A separate g=0,m=0 control can be solved exactly before eliminating its
colour-flux states. Starting with three quarks of one spin/flavour at the
right cell, raw24-mode CAR hopping generates four invariant colour-singlet
vectors with norms1,3,12,36. In their normalized basis the diagonal entries
are0,Delta,Delta,0 and the three off-diagonal entries are
-sqrt(3)t,-2t,-sqrt(3)t. The even and odd low eigenvalues are

    E_even=(Delta-2t-sqrt((Delta-2t)^2+12t^2))/2,
    E_odd =(Delta+2t-sqrt((Delta+2t)^2+12t^2))/2.

The exact half splitting tends to t_B at weak hopping. At the three
inherited kappas1/2,1,2 the exact-to-third-order ratios are approximately
0.82614,0.95537,0.98895. Rational square-root enclosures retain the exact
control intervals. The leading approximation is measurably imperfect at
the smallest coupling; it is not treated as the exact transfer rate.

`CONTROL_DYNAMICS.json` records the full unitary control, with all four
left-quark-count outcomes at each source time. The intermediate one- and
two-quark placements carry colour flux and are retained as outcomes.
Independent even/odd amplitudes check the same probabilities. This is a
coherent transfer experiment, not exponential decay or a nuclear lifetime,
and no literal native program execution is claimed for it.

## Remaining physical obligations

The actual shared source retains nonzero g and mass, its half-filled
vacuum, all auxiliary complements, and the previously stated effective
remainder limitations. A physical baryon dispersion requires that source's
states and spectral error bounds. A halo additionally needs an identified
core, separation thresholds, spatial tails, and the correct motion-to-binding
ratio. Its lifetime requires the weak current and all outgoing channels
on the same clock. This packet supplies an exact source-generated motion
operator and complete finite controls; it does not supply those missing
nuclear predictions or a free coefficient chosen to fit them.
