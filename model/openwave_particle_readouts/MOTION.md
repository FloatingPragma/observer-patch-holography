# Finite charged excursions after energy and position preparation

## The same source and a complete position intervention

Use the full source, physical preparation, reference and energy bin of
[PROOF.md](PROOF.md). Its Hamiltonian is
`H=H0+lambda T+4 N_ref`, with `lambda=2^-34`, `||T||<=360`.
Every link has the declared unit excitation penalty. The full static
CAR reference stays at cell zero. Write `chi=P_I(H) B_r Psi` for the
unnormalized prepared line, conditional only on physical preparation
success. The retained-environment comparison is

`||chi-phi||<=epsilon`, `||phi||^2=w`,

where phi is the lowest actual atomic charged component of `B_r u`.
Its total atomic energy is `a=7E0+epsilon_r+4`. The vector u is the
same native-coordinate cube trial used in the interacting-ground
comparison. Its spin factors may be entangled.

Define S0 as the tensor product of blank-link projectors and the true
local ground-doublet projector at every cell except zero. It is identity
on the complete origin CAR/scalar space and the retained reference.
The two-outcome instrument `{S0,I-S0}` preserves every joint Gauss
constraint and parity. It follows the energy-bin measurement. Its failure
is retained with subsequent arrival controls off. Set `psi=S0 chi` and
`s=||psi||^2`; intermediate branches are unnormalized.

Since `S0 phi=phi`,

`||psi-phi||<=epsilon`, `p_line-s<=epsilon^2`,

and the same line bounds give

`s_-=(sqrt(w_-)-epsilon)^2 <= s <= (sqrt(w_+)+epsilon)^2=s_+`.

S0 need not commute with H or its energy bin. The exact identity

`(I-P_I)S0 chi=-(I-P_I)(I-S0)chi`

proves `||(I-P_I)psi||^2/s<=epsilon^2/s_-`. Because P_I commutes
with H, this upper bound on the conditional probability outside the
original total-energy bin holds at every subsequent source time.
It does not make the post-position state exactly supported in that bin.
The earlier successful energy record is a retained preparation result.

## Arrival at the same local charged atom

For each of the six actual links incident to zero, let A_e require:
the selected link nonblank, every other link blank, true local matter
ground at every cell except its remote endpoint y, and the complete
local h spectral projector at energy epsilon_r at y. The reference
factor is identity. This is a full-carrier spectral effect. The other
local charged atoms and all unselected configurations belong to the
complete complement `I-sum_e A_e`.

The six effects are orthogonal because their nonblank link patterns
differ. They preserve the joint source/reference Gauss constraints.
They are subeffects of the charged arrivals in
[interacting_charged_transport](../coherent_defect_wave/chiral_joint/interacting_charged_transport/PROOF.md).
The source charge is at y, the origin matter returns to its local
ground, and the fixed origin reference is compensated by the retained
flux. Both parallel routes to each of cells 1, 2 and 4 are counted.

S0 commutes with H0. Every arrival requires a nonblank link, hence

`A_e exp[-it(H0+4N_ref)] S0=0`

for every input and every real t. This proves exact zero arrival at
t=0 and under the zero-hopping control. Only edge e can contribute to
`A_e T S0`. Its full norm is at most 15. This remains true between
arbitrary atomic evolution factors, since those factors preserve S0
and A_e. The bound fifteen is the inherited norm of the entire native
edge term on the unrestricted CAR/scalar/link carrier. It includes
every family and both spin components, and is unchanged upon tensoring
with the complete reference and preparation environment. Contracting
with S0 and A_e cannot increase it. The estimate therefore applies to
arbitrary origin matter and arbitrary entanglement with the reference
and environment in S0, including the error vector `psi-phi`. No source
complement is truncated.

On the invariant one-reference-particle sector populated by pump
success, the complete atomic Hamiltonian acts as `a+1` on A_e: seven
local ground factors, the same charged atom at y, one unit link
excitation and the explicit reference energy four. This operator
identity holds on the entire local spectral eigenspace at epsilon_r,
including every degeneracy and any other local representation at
that exact energy. It needs no rank-one target assumption. The full
reference Fock space is retained; its other occupations have zero
weight for this specified successful pump. The exact input-to-arrival
detuning is therefore one. It is the declared link energy, with no
fitted effective hopping or shifted source coefficient.

## Exact atom-resolved first-hop coefficient

Let Q_r be the complete local spectral projector at epsilon_r. The
native gauge and spin symmetries give, for coordinates in family r,

`P_g c_(alpha,a)^dagger Q_r c_(s,b) P_g`
`=delta_ab [w delta_alpha,s I+c sum_i sigma_i,s alpha S_i]`.

The scalar w is exactly the inherited pump atom weight. The real
coefficient c is the residue of the signed complete charged resolvent,
computed from the difference of the two native spin transfer matrices.
The positive resolvent uses their sum. Both retain all target states,
with the exact physical Gram adjoints. The selected root is identified
by its isolating interval in the complete target annihilator.

`motion.atom_first_hop` proves a polynomial remainder identity on that
root's irreducible factor. It gives

`c_ec=-2 w_ec/3`, `c_L=2 w_L`.

The signed residue is also independently enclosed by outward interval
arithmetic before applying these identities. The pump normalization is
`1/sqrt(2d_r)`, the actual hopping coefficient has squared magnitude
`1/16`, and the retained blank-link Gram matrix is `I/d_r`.
No reference coordinate or flux orientation is discarded. Native
hypercharge and gauge covariance exclude other-family refilling and
origin pair annihilation. CAR grading contributes a common phase to
each selected route and therefore does not remove an interference term.

Contracting these factors gives the exact two-ground-spin operator

`M_axis=(w^3/16+3w c^2/64)I`
`       +c^2(2w-c) S_0 dot R_axis S_y/16`.

The code verifies this identity in all three 4-by-4 native Pauli
contractions. In the native-coordinate cube ground,
`<S_0 dot R_axis S_y>=x/48`, where `x=4 epsilon_cube` is the
selected root of the inherited cube characteristic polynomial. Thus

`m_r=||A_e T phi||^2=w^3/16+3w c^2/64+c^2(2w-c)x/768`,

or, using the exact residue identities,

`m_ec=w_ec^3(1/12+x/648)`, `m_L=w_L^3/4`.

The L coefficient is independent of the cube correlation because its
two-spin effect is scalar. Both coefficients are strictly positive.
All six routes have this same comparison moment; equal exact
interacting route probabilities are not assumed.

## Finite-time bound on the complete clock

Write `U0(t)=exp[-it(H0+4N_ref)]`. The first Dyson term on phi has
the exact selected amplitude

`A_e D1(t)phi=lambda e^(-iat)(e^(-it)-1) A_e T phi`.

Its norm is `2lambda |sin(t/2)|sqrt(m_r)`. The detuning and the
phase are explicit. For the actual input psi, edge selection and
unitarity of the atomic factors give

`||A_e D1(t)(psi-phi)||<=15lambda |t|epsilon`.

The second-order Duhamel remainder, using the full H and H0 unitaries,
has norm at most `(lambda tau t)^2/2`, with `tau=360` and
`||psi||<=1`. This estimate has no factor `||H0|| t^2` and no
truncated charged generator. Combining it with the exact zero-order
kernel bounds each actual route probability, conditional only on the
physical preparation, by

`lambda^2 [max(0,2|sin(t/2)|sqrt(m_-)-E(t))]^2`
` <= r_e(t) <= lambda^2 [2|sin(t/2)|sqrt(m_+)+E(t)]^2`,

`E(t)=15|t|epsilon+lambda tau^2 t^2/2`.

The receipt evaluates rational times 0, 1, 3 and 6. Alternating Taylor
bounds enclose `2sin(t/2)` on this interval with rational endpoints;
the first omitted tail is positive and its terms decrease. Every
nonzero-time lower bound is positive for both native channels. For
each family the original-input lower bound at time three exceeds
the upper bound at time six, including the physical preparation's
success probability. These are separate terminal time settings of
the same preparation and clock. No intervening arrival measurement
is inserted into the time-six experiment.

This is an off-resonant flux excursion with a finite return contrast.
The comparison first-order response has its unit-detuning oscillation;
the theorem bounds complete same-H finite-time outcomes. It does not
infer a charged propagation band or macroscopic velocity.

## Complete original-input records

Let q be the physical preparation success probability. The eleven
coarse records at each terminal clock setting have probabilities

`1-q; q p_failure; q p_outside; q(p_line-s);`
`q r_e(t) (six routes); q[s-sum_e r_e(t)]`.

They sum exactly to one because
`p_failure+p_outside+p_line=1`. The complete nonarrival has a positive
lower bound. The preparation failure, pump failure, energy complement
and position failure retain their actual source states. Every source
copy, reference occupation, scalar/link memory and preparation
environment is retained. All branches preserve their stated joint
Gauss constraints.

The result concerns a line-conditioned position preparation. It does
not assert positive motion of the original line without that
intervention, or an arrival probability conditioned on a final
full-H energy measurement. The latter projector need not commute
with A_e. The energy-disturbance bound is larger than the arrival
probabilities and cannot locate those arrivals inside the original
energy bin. The fixed compensator and its flux cost are part of the
construction. No free charged particle, empirical species identity,
relativistic dispersion or continuum limit is identified here.

## Finite repetition counts

The simulated arrival probabilities are small. Define l3 and u3 as the
original-input lower and upper bounds for the union of six routes at
time three, and u6 as its upper bound at time six. The exact contrast
`delta=l3-u6` is positive. Suppose each repetition independently
prepares the entire source and apparatus from its declared physical
input, retaining all failure records. A preparation failure contributes
zero to that repetition's arrival indicator; it is not removed from
the sample count.

For `N_detect=ceil(3/l3)` independent time-three repetitions, the
probability of no arrival is at most

`(1-l3)^N_detect <= exp(-N_detect l3) <= exp(-3) < 1/20`.

The last strict inequality follows from the exact rational partial
sum `sum_(k=0)^8 3^k/k!>20`.

For an independent comparison using N repetitions at each of times
three and six, let D be the difference of the two empirical arrival
frequencies. Its mean is at least delta and its variance is at most
`(u3+u6)/N`. Chebyshev's inequality therefore gives

`Pr(D<=0) <= (u3+u6)/(N delta^2) <= 1/20`

when `N=ceil[20(u3+u6)/delta^2]`. The receipt records both integer
budgets and the total `2N` original inputs in the comparison. These
counts assume independent complete preparations and retain failed
preparations in the sample budget. They give no physical preparation
rate or feasible laboratory runtime.
