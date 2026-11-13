# Complete native rest bands and degree counterterms

## Declared source

The carrier, colour currents, seven rational coefficients and electron unit
are those of `openwave_common_hyperfine_spectrum`. The physical path carrier
uses the zero-boundary-flux Gauss construction in
`openwave_strong_step_scaling`. All occupations and edge irreducible
representations are retained. Write the local number as N_x, spin as S_x,
light-quark isospin as I_x, colour Casimir as C_x² and edge electric Casimir
as E_e². The colour term is

    (U−κ/2) sum_x C_x² + κ sum_e E_e².

The path action subtracts

    sum_x (deg(x)−1)[λ S_x² + η I_x² + t N_x]

on paths of at least two cells. The selected action also adds

    sum_x (deg(x)−1)[ΔU C_x² + Δm N_x].

Both sums vanish identically on two cells, in every Fock sector. The second
pair of coefficients is specified in `inputs.json`. It uses the inherited
neutron and axial targets, so agreement is calibration. Its source is a
degree-dependent action, without an assumption that architectural choices
are unique. The path extension and its infinite bulk share one polynomial.
The auxiliary evolution parameter is the same source clock convention as
the parent. The GeV conversion retains the calibrated electron anchor.

## Stability without a local occupation cutoff

For nq fundamental and na antifundamental colour modes, every highest weight
(p,q) obeys p+q≤N=nq+na. Hence C²≤N(N+3)/3 for arbitrary family capacity.
Gauss gives sum C_x²≤4 sum E_e². The positive quark-antiquark spin interaction
K=3 nq na/2+2 S_q·S_a satisfies

    K ≥ nq na−min(nq,na).

This follows from angular-momentum addition and S_q≤nq/2, S_a≤na/2.
Put d=nq−na, c=max(0,κ/4−U), h=pair_triplet and m=mbar−|δ|/2. After
discarding the other positive current squares, the local lower bound is

    αN²+βd²+γN+(h/2)|d|,
    α=h/4−c/3, β=A/9−h/4, γ=m−h/2−c.

The degree subtraction has magnitude at most
(λ+η)N(N+2)/4+tN because the maximum degree is two. Its effect on α,γ is
included exactly in `positivity.py`. Negative ΔU costs at most
|ΔU|N(N+3)/3; positive Δm may be discarded. All remaining coefficient
inequalities are rational. For N=1 use the exact fundamental Casimir4/3;
for N≥2 use α>0, β>0 and the increasing αN+γ. This proves a positive
linear number floor on every open chain and at arbitrary family capacity.

The same file checks the permitted strong charge-conjugation completion
cost −3(Jw+η)min(Nq,Na)/4. Its global minimum function is nonlocal; no
range-one locality is claimed for that term. A range-three nuclear term
conditioned on three successive local occupations N_x=3 pays its negative
bound only at those cells. Each cell belongs to at most three windows, so
the sum of window projectors is bounded by the number of N_x=3 cells.
The rational N=1,2,3 and N≥4 cases prove joint stability without subtracting
an unnecessarily loose nuclear penalty from every particle.

## Exact compact reduction

Fix three quarks, total S=I=1/2 and magnetic projections Sz=1/2, Iz=−1/2.
The native colour epsilon gives the complete symmetric orbital carrier.
A co-located triple contributes one state per cell; a two-plus-one cluster
on an edge contributes two spin/isospin couplings in each orientation.
Let P retain these five states per translation cell. Q contains every
remaining relative configuration, with occupied span at least two.

`exact_native.py` constructs rational vectors in the unnormalized native
colour basis, verifies both Casimir intertwiners and its exact Gram matrix,
and contracts the complete four-cell operator. Its Gram ranks are16 in the
104-dimensional fixed-projection carrier. This is an exact coordinate
reduction, without replacement of the high source by a selected surrogate.
The support of each hopping term extends by one cell. Thus the complete
four-cell identity fixes the translation-covariant infinite-chain operator.

In the order (co-located, two left-oriented states, two right-oriented
states), the compact Bloch matrix has a scalar a, an edge block B, and row
couplings v_l+exp(−ik)v_r. With t=1/8,

    B = [[c,0,−t/2,−3t/2], [0,b,−3t/2,−t/2],
         [−t/2,−3t/2,c,0], [−3t/2,−t/2,0,b]],
    v_l=−sqrt(3/2)t(1,1,0,0),
    v_r=−sqrt(3/2)t(0,0,1,1).

The exact rational expressions for a,b,c are in `analytic.constants`.
At k=0 the even three-state block and odd two-state block retain all five
states. After scaling the co-located coordinate by sqrt(3), the even
matrix and its Gram matrix are rational. Deterministic rational inertia
bisection certifies every eigenvalue used here. Frozen rational vectors are
admitted by exact residual/variance identities, rather than by numerical
eigensolver accuracy.

## Complete complement and current transfer

For a separated two-plus-one configuration of span r≥2, the native colour
energy is 8U/3+4κ(r−1)/3. For three occupied cells of span r≥2, it is
4U−2κ+4κr/3. Retain the local density and spin costs, discard positive edge
squares and bound the degree subtraction. The two resulting rational
floors are the expressions `split` and `three` in `analytic.full_bounds`.
Their minimum C is a lower operator bound on QHQ, including its hopping:
the original positive graph Laplacian is discarded before the number
counterterm is bounded. In the selected bulk action U and mbar include
ΔU and Δm. The exact κ is the declared source coefficient, rather than an
endpoint chosen with an inappropriate inequality direction.

The complete native contraction proves

    (QHP)*QHP = T0+exp(ik)T1+exp(−ik)T1*,
    T0=3t² I4,
    T1=(t²/2) [[0,0,0,0],[0,0,0,0],[1,3,0,0],[3,1,0,0]]

on the edge coordinates and zero on the co-located coordinate. Its
eigenvalues are t²,2t²,4t²,5t² for every real k. Thus ||QHP||²=5/64.
No output configuration is discarded in this multiplication.

If a_top bounds PHP above and C>a_top, the Schur complement gives five
isolated low energies, each in [e_j−ε,e_j], where
ε=(5/64)/(C−a_top). The next spectral value is at least C. This theorem
applies to the countable relative-position fibre: its diagonal confining
operator is self-adjoint on its multiplication domain, and hopping is
bounded. The finite-rank Schur argument does not require finite Q.

Let e_0 be the compact ground energy, g its certified compact spectral gap,
and b²=5/64. The normalized full ground has Q probability at most
b²/(C−e_0)². Its normalized P component differs from the compact ground by
sine at most ε/(g−ε). The actual axial observable is the compressed
I_minus J_plus on this fixed S,I carrier; its norm is at most three.
It preserves occupation span, so P J Q=0. If σ is its certified standard
deviation in the compact trial, the axial expectation error is at most

    2σs+6s²+6q,

where s is the P-vector sine bound and q the Q probability bound. The
trial-vector error is included separately by its exact rational residual.
Native I_plus gives the vector matrix element one and transports the full
neutron eigenstate to the proton eigenstate with energy shift −δ.

## Records and physical boundary

Finite physical input is I/d on the complete Gauss singlet carrier. The
compact filter and its orthogonal complement are both recorded; their
probabilities are rank(P)/d and1−rank(P)/d. No raw occupation poststate
outside Gauss is used. Opposite-charge reference registers balance every
input charge sector. Complete spectral readout has all five low outcomes
and the entire high complement, commutes with the source and has zero
mean readout work. The operator bounds hold with arbitrary quantum
references. Infinite fibres are normalized state or wave-packet domains;
no uniform state on an infinite carrier is asserted.

The displayed k=0 energies are one-dimensional source rest bands with an
explicit unit dictionary. Calibration of these bands and preservation of
the two-cell action do not establish all other bulk hadron masses, a
three-dimensional Lorentz limit, scattering poles or a physical QCD
identification. The fixed-coefficient control is disjoint from the selected
neutron/axial targets and is retained as a source test.
