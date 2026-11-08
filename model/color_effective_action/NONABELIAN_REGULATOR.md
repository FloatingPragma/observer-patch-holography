# Nonabelian coordinate identities and representation growth

The native current fixes the compact group and the positive electric
metric used here. The existing [coordinate
attachment](../particle_gauge_band/QUANTUM_LINK_PROOF.md) realizes its
matrix coefficients in L²(G,Haar), for G=S(U(3) x U(2)). The following
fixed-graph regulator theorem applies to that stated attachment. The
sixty-state raising-link implementation does not select it through a
capacity parameter; the exact tests below distinguish those sources.

## The low-energy identity a coordinate family must supply

Let C be the positive endpoint Casimir with the source normalization

    C=2C_3+C_2/2+q²/12.                                    (1)

For irreducible representations r,s and a constituent t of r tensor s,

    sqrt(C_t)<=sqrt(C_r)+sqrt(C_s).                         (2)

To prove this, write C=sum_a X_a² in its positive orthonormal current
frame. On the tensor product the current is X_a tensor I+I tensor X_a.
The triangle inequality in the direct sum over a bounds its norm by
(sqrt(C_r)+sqrt(C_s))||v||. Restricting to t proves (2). Dual
representations have the same Casimir. Frobenius reciprocity then also
gives |sqrt(C_t)-sqrt(C_s)|<=sqrt(C_r).

Let P_Lambda contain every matrix-coefficient block with C_s<=Lambda,
including its conjugate. For a unitary matrix-valued coordinate D^r(g),
put U_Lambda=P_Lambda M(D^r)P_Lambda. If psi is supported on C<=Lambda0
and

    Lambda >= (sqrt(Lambda0)+sqrt(C_r))²,                   (3)

then multiplication by D^r and by its adjoint stays inside the retained
space. Therefore both U_Lambda*U_Lambda psi=psi and
U_Lambda U_Lambda*psi=psi. This is a domain-specific identity; the
operator norm of the boundary unitarity defect need not tend to zero.

For any normalized vector of mean electric Casimir at most Ebar,
spectral Markov bounds and (3) instead give

    <I-U_Lambda*U_Lambda> <= Ebar/(sqrt(Lambda)-sqrt(C_r))²,
    <I-U_Lambda U_Lambda*> <= Ebar/(sqrt(Lambda)-sqrt(C_r))², (4)

with either right side capped at one. The norms of the two unitarity
defect operators applied to psi are at most the square root of this bound,
since each defect is a positive contraction. For a retained input,
||U_Lambda psi-M(D^r)psi|| has the same square-root bound. For an arbitrary
full-space mean-energy input the latter comparison has the safe bound
twice that square root. This is the controlled low-energy coordinate
premise that a source-derived family could test.

## Exact failure and required growth of the retained native menu

The five native charge/dimension pairs are (-4,3),(-3,2),(1,6),(2,3),(6,1).
None of their conjugates occurs. Each normalized native raising matrix
V=U/sqrt(d) therefore has V*|external,0>=0. Its opposite unitarity defect
I-VV* is exactly the identity on the zero-electric-energy blank. The
finite coordinate-menu completion has the same defect there because the
conjugate matrix coefficients are also absent. This does not satisfy (4)
at any mean-energy bound Ebar=0. Tensoring spectators or changing the
dimensionless capacity label of the unchanged local operator cannot
remove this matrix identity.

Adding all conjugates would require at least

    1+2(3²+2²+6²+3²+1²)=119                               (5)

link states. That is only the first closure condition. Arbitrary products
of the native Q=(3,2)_1 include the highest-weight block

    (Sym^n 3,Sym^n 2)_n,
    d_n=(n+1)²(n+2)/2,
    C_n=n(7n+18)/8.                                        (6)

The source-normalized fundamental generators satisfy the exact Fierz
identity sum_a T_a tensor T_a=swap/2-I/(2d), for d=3 and d=2. The code
checks every tensor entry on the actual source blocks. On a symmetric
n-fold tensor this gives C=n(d²-1)/(2d)+n(n-1)(d-1)/(2d).
Thus the color highest-weight Casimir is n(n+3)/3, the weak Casimir is
n(n+2)/4 and the charge is n. Substitution in (1) proves (6).
It agrees with the literal source Q block at n=1. At n=2 its dimension
is 18 and its Casimir is 8; this block is outside the native menu.
Thus any coordinate family containing every finite polynomial in the
faithful native representation needs unbounded representation growth.
Finite algebra covariance supplies no particular dependence of that
growth on observer capacity.

## Fixed-graph action convergence under the stated completion

On the same two-vertex graph, conditionally replace each raising entry by
sqrt(d_r) times full coordinate multiplication, exactly as prescribed by
the existing attachment. Keep all given family controls, scalar cutoff,
reference and Hamiltonian coefficients. Colored fermion occupation is zero,
as in the source's declared charged-lepton sector. The electric part is e sum C,
e>0. All other terms form a bounded operator V. At the source defaults,
an explicit bound on the complete fermion/scalar space is

    ||V|| <= 108 + 168 + 24 + 1 = 301.                      (7)

For the magnetic term, a two-edge r-loop has multiplier d_r chi_r, with
norm at most d_r². Summing the three colored blocks and the adjoint gives
2(9+36+9)=108. For hopping, each edge has three diagonal families and
two spin components, four weak-doublet entries of norm sqrt(2) and one
singlet entry of norm one. Both edges and adjoints give
24(4sqrt(2)+1)<=168. The two Higgs components, two spin signs and three
families, with adjoints, give 24. The scalar number has norm one.
These bounds retain every complete state; restricting to the Gauss
subspace cannot enlarge them.

Every bounded insertion changes sqrt(C) on an affected link by at most
two, since the largest native C_r is four. For an input with per-link
Casimir radius at most rho0 and a product cutoff rho=sqrt(Lambda),
r=floor((rho-rho0)/2)+1 insertions are necessary to leave the cutoff.
A path returning to the input spectral subspace needs at least 2r.
For mu>301, the same diagonal-resolvent series as in the
[implemented rotor limit](REGULATOR_PROOF.md) gives

    ||(R-R_Lambda)P_input|| <= 2(301/mu)^r/(mu-301),
    ||P_input(R-R_Lambda)P_input||
       <= 2(301/mu)^(2r)/(mu-301).                          (8)

Here R=(mu+e sum C+V)^(-1). The nonabelian tensor-Casimir bound replaces
the integer flux-step bound. Fixed-input real-time propagators likewise
differ by at most 2 sum_(n>=r)(301|t|)^n/n!, by the interaction-picture
expansion. Complete flux instruments preserve the Casimir support, so
the bound extends to their full records at fixed total unitary time,
including references and all boundary outcomes.

This proves strong-resolvent and complete finite-history convergence on
the dense finite-representation domain for the stated fixed-graph
coordinate completion. It controls the entire effective resolvent action,
including induced operators, without fitting their coefficients separately.
The construction does not suppress those operators or prove that they
reduce to a single Yang–Mills coupling under spatial coarse-graining.

## The coordinate completion changes the same blank experiment

The full coordinate Hamiltonian cannot be attached to the original finite
source by declaring its additional representation states unobserved.
For the same empty input chi_0, the native magnetic term creates one loop
in each of its three colored blocks, with total squared norm
3²+6²+3²=54. Its adjoint annihilates chi_0. The full coordinate adjoint
instead creates the conjugate loops. They are orthogonal to the original
three loops and have the same total norm. Thus the full coordinate model
has loop coefficient 108, while the unchanged Higgs creation coefficient
is 12 in both configurations. Complete measurements of nonvacuum link
flux and scalar occupation have leading probability ratio 9/2 and 9,
respectively. A common positive rescaling of the Hamiltonian and clock
cannot change those leading ratios.

This distinction also has a finite-time certificate. The native source
satisfies ||H||<24 and ||H chi_0||=sqrt(66)<9, hence
||H² chi_0||<216. In the full coordinate model,
||H chi_0||=sqrt(120)<11. Every vector in H chi_0 has total electric
Casimir at most eight, so (7) gives ||H² chi_0||<309 times 11<3400.
Both chi_0 lie in the corresponding H² domains. The unitary integral
remainder therefore bounds the state after its linear term by
108t² and 1700t² respectively. For a projector annihilating chi_0,
with leading amplitude norm at most a, the probability error from its
leading coefficient is at most

    2 a R |t|³ + R² t⁴,                                   (9)

where R is 108 or 1700. The exact receipt uses t=2^-20, a=8 or 11 for
flux, and a=4 for Higgs occupation. The resulting complete flux/Higgs
probability-ratio intervals are disjoint. The failure of direct source
identity is an analytic statement about all outcomes. It does not decide
whether different collective low-energy limits could agree.

The flux and scalar-occupation projectors commute. Their products with
both complements form a complete four-outcome instrument in either source.
The probability ratio uses two sums of that complete record and its
strictly positive denominator bound; it does not condition away an outcome.

## Attachment to a spatial field family

A spatial refinement could reuse this theorem if its actual specification
supplied the graph, coefficient functions, representation cutoffs and
state embeddings. Its finite-volume bound B(a) must replace 301, and its
Casimir radius must grow quickly enough that the explicit tails at total
time T tend to zero. The number of links and the prescribed coefficients
enter B(a); a fixed graph bound is not a uniform continuum estimate.
Gauge covariance of the embeddings, agreement of the electric metric and
coordinate multiplication on each finite-representation core, and a
common state/clock are required before transferring this regulator limit.

The inspected source implementations distinguish the available families.
The gauge–Higgs constructor fixes two vertices and two edges and exposes
no spatial or representation-capacity refinement. The native tree source
has a uniform occupation-projector electric law and no magnetic cycles.
The cubic staggered source takes supplied classical unitary backgrounds.
The U(1) rotor implements its own flux cutoff, as evaluated in the adjacent
proof. None of these specifications identifies a capacity-selected
nonabelian spatial coordinate family satisfying (3)–(8). The finite
raising model has the explicit vacuum defect above; the required growing
coordinate family is unspecified rather than disproved.

Consequently this result establishes a conditional nonabelian regulator
correspondence with source-normalized coefficients and a quantitative
test for its missing source attachment. It gives no physical reference
scale, renormalized QCD coupling, interacting spatial continuum or particle
mass. Those quantities are not assigned by a dimensional count or by the
bare coefficient ratio.
