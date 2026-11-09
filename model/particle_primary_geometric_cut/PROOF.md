# Covariant regional restriction and boundary entropy

## An actual cut of the selected source

Use the [primary contract](../particle_common_contract/PROOF.md), including
all sixty fermion modes per cell, the open Wilson Hamiltonian, product
link cutoff, Gauss projection, positive polynomial state, idle canonical
reference and source time. On physical support the outer ghost links are
constant. The remaining graph consists of the interior sites and their
nearest-neighbor links. No link is divided into added edge registers.

For a nonempty site set R choose the regional factors to be its matter
modes and all links whose two endpoints are in R. All crossing links
belong to the complement. This is a declared regional observation
convention on the unchanged source. It replaces neither its state nor
its dynamics. The region has no asserted proper radius or volume.

Order the regional fermion modes before the complement, with the exact
exterior-power permutation signs. Permute the whole link factors in the
same way. The resulting unitary J_R defines

    sigma_R = Tr_out(J_R rho J_R*),
    j_R(A) = J_R*(A tensor I)J_R.                           (1)

These are the trace-dual restriction and injection of
[the finite regional theorem](../particle_spatial_modular_response/PROOF.md).
Its completely positive restriction, multiplicative injection and tracial
expectation instantiate the same finite fields of
[`CPRestrictionNet`](../../Lean/QFT/CPRestrictionNet.lean). The expectation
need not preserve the interacting state. A physical Hilbert-space tensor
factorization is not assumed: (1) acts on the ambient factors. Its physical
observable injection uses the gauge-invariant regional subalgebra, which
preserves the Gauss support. Compressing an arbitrary ambient operator
does not give a multiplicative map on all regional matrices. Every
operation retains the exterior and the canonical reference in the full
experiment.

The [source cube action](../particle_primary_channel_labels/PROOF.md)
maps an internal link to an internal link of the rotated region. Link
reversal is the source coordinate-inversion unitary, which preserves the
word cutoff. Second quantization of the site permutation and spin action
factors after J_R and J_(gR), since it maps the two one-body summands to
the corresponding summands without mixing them. Thus

    J_(gR) U_g J_R* = U_(g,R) tensor U_(g,out),
    sigma_(gR) = U_(g,R) sigma_R U_(g,R)*.                 (2)

The selected rho is cube invariant. Equation (2) proves covariance of
the partial trace, its support, spectrum, entropy and supported modular
flow. The central binary element is parity; even observables carry the
24 proper spatial rotations. The owner-cell cut is a different regional
algebra and is not silently identified with (1). The literal source
charge conjugation also factors: its sixty-mode site transformations
are even and its link transformations act factor by factor. It therefore
preserves sigma_R and conjugates the regional charge sectors exactly.

## Boundary representations and a source-derived bound

Let b_R be the number of crossing interior links. At an inside boundary
vertex x let k_x be its number of crossing links. Only those factors
carry the complementary action of the gauge group at x. After imposing
Gauss constraints away from this boundary, the actual finite regional
representation has the Schur decomposition from the regional theorem:

    sigma_R = direct_sum_lambda sigma_lambda tensor I_dlambda/d_lambda,
    S(sigma_R) = -sum_lambda Tr(sigma_lambda log sigma_lambda)
                    + sum_lambda p_lambda log d_lambda,       (3)
    p_lambda = Tr sigma_lambda.

All admitted weights and multiplicities are those of the actual partial
trace; no maximally mixed boundary state replaces them. The global
faithfulness bound rho >= c P_G implies
sigma_R >= c Tr_out(J_R P_G J_R*) on its exact support. It supplies finite
supported logarithms, without a volume-uniform faithfulness assertion.

Each native irreducible gauge block has dimension at most six. A word
of degree at most ell belongs to a tensor product of at most ell native
or conjugate representations. Every irreducible summand has dimension
at most 6^ell. Direct sums and multiplicities do not increase this
maximum. Applying this to all crossing links and then decomposing under
the product of inside boundary groups gives

    d_lambda <= 6^(ell b_R),
    0 <= sum_lambda p_lambda log d_lambda <= ell b_R log 6.     (4)

This is a boundary-size upper bound on the representation-dimension
term alone. It is neither a total entropy bound nor an area equality.
With the declared lattice face area A_cut=b_R a^2, the coefficient in
(4) is ell log(6)/a^2, which diverges on the selected refinement.
No finite Newton coefficient follows from this upper bound.

The integer-normalized native hypercharges have absolute value at most
six. Gauss matching consequently gives |Q_(6,x)|<=6 ell k_x on the
physical regional support. Here Q_(6,x) is the regional boundary gauge
generator, including matter and the internal-link endpoints at x;
it is not the matter charge alone. Integer spectral projectors over this finite
range give complete boundary hypercharge records there. Nonabelian
central projectors resolve further sectors. These projectors commute
with the full gauge action; their apparatus is physical. Unsupported
outcomes have zero weight and are retained. An arbitrary regional
matrix outside the physical support is not assigned this charge bound.

For the concrete cut R={0,e1}, one link of an elementary plaquette is
inside and the other three are outside. Each normalized native Wilson
character L_r=chi_r(A B) has the exact Schmidt decomposition

    L_r = (1/d_r) sum_(i,j) [sqrt(d_r) R_r(A)_ij]
                                  [sqrt(d_r) R_r(B)_ji].       (5)

Product Haar measure makes each displayed family orthonormal. There
are d_r^2 equal Schmidt terms, with regional eigenvalue 1/d_r^2 and
opposite boundary charges at the two endpoints. The physical character
belongs to every selected cutoff. Thus its boundary sector has selected
weight at least c and sigma_R >= (c/d_r^2) Pi_r on its Schmidt subspace.
The blank and the conjugate characters give the corresponding zero and
opposite-charge witnesses. This constructs actual populated sectors,
without computing the full regional density or its large multiplicities.
The ten nonzero hypercharge tuples are distinct. Summing their separate
positive contributions in (3) gives the finite lower bound
S_edge >= 4 c log(3*2*6*3*1)=4 c log108. The code authenticates the exact
factored c=2^(-N) from the selected preparation. This strictly positive
bound can deteriorate with capacity; it is not a uniform area coefficient.

## Complete center readout and entropy domains

Let Pi_lambda be the orthogonal central projectors. The complete neutral
apparatus is V psi=sum_lambda Pi_lambda psi tensor |lambda>. A controlled
cyclic shift extends it to a unitary; no outcome is discarded. Since the
regional density is block diagonal, its classical-record marginal is

    chi_(record,R) = sum_lambda |lambda><lambda|
                                     tensor Pi_lambda sigma_R Pi_lambda.

Its nonzero eigenvalues are precisely those of sigma_R, hence

    S(chi_(record,R)) = S(sigma_R).                         (6)

There is no independent coin entropy in (6): the label is an existing
central block. The complete source-reference input is pure; after a
classical dephasing of its recorded outcomes, its full entropy is H(p).
The source can have coherences between matching global sectors, so the
full source measurement is not asserted to be nondisturbing. A coherent
record isometry before classical dephasing instead preserves purity.
These are distinct entropy functionals of the retained experiment.

On the positive support the extended modular generator is
direct_sum_lambda(-log sigma_lambda+log d_lambda). Its central log d term
affects entropy changes of the weights. For supported trace-zero
tangents the inherited finite first law and controlled remainder apply
to (1). Cube rotations carry the entire experiment and its reference to
the corresponding experiment in gR. Charge-conjugation-odd tangents
have zero first-order regional entropy; this does not imply zero
quadratic response or a charged-particle energy.

## Energy and geometry have separate domains

The [covariant energy density](../particle_primary_covariant_stress/PROOF.md)
gives H'_R=sum_(x in R) h'_x and the exact finite Ward current on the
same H. In general H'_R includes crossing hopping and plaquette terms,
so it is not contained in the induced regional algebra (1). One cannot
set K_R=beta H'_R merely by assigning the same region label. The
finite entropy formula, energy conservation and rotation covariance
are compatible source results with this explicit support distinction.

There is a source-defined regional energy in that algebra. Let H_int(R)
sum the whole onsite terms in R, whole interior-link electric and hopping
terms and whole plaquettes whose four vertices are in R. Every term is
physical, and H_int(R) is a covariant gauge-invariant regional operator.
Its exact difference from H'_R is half of each crossing electric and
hopping term plus |V_p intersect R|/4 of each partly contained plaquette.
Every partial plaquette has at least two crossing edges; a crossing edge
belongs to at most four plaquettes. The actual full-carrier source norms
therefore give

    ||H'_R-H_int(R)|| <= (b_R/a) [609 ell^2/4+15+3/14].       (7)

This is an explicit boundary comparison, not a vanishing continuum error.
H_int(R) supplies physical regional phase interventions, with their source
work governed by the existing polynomial-state passivity theorem. The
interacting marginal modular generator is not identified with H_int(R).

The source cut and (3) supply a covariant regional entropy and its full
boundary convention. Proper geometric balls, the Lorentz diamond
modular kernel, fixed-proper-volume area response, generalized-entropy
balance and errors o(ell_radius^4) are additional hypotheses of the
[matched Einstein theorem](../matched_einstein_limit/PROOF.md). No such
physical identification is claimed here.

The entropy decomposition is consistent with the finite gauge-region
analysis of [Donnelly](https://arxiv.org/abs/1109.0036) and the algebra
dependence discussed by [Casini, Huerta and Rosabal](https://arxiv.org/abs/1312.1183).
Those results do not supply the missing proper geometry or source area
coefficient. The derivation above uses the actual native representation
and cutoff rather than importing a continuum entropy law.
