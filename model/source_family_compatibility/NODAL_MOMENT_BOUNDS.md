# Local CAR constraints on the complete momentum response

The centered full-CAR source on the proper-cubic ground average obeys

    0 <= G(k) <= Gamma c^2 I,
    0 <= R*R(k) <= Gamma c^2 I,
    Gamma=(31+sqrt(577))/4 <14.                         (1)

This is a constraint on the actual force and residual Gram, at every
magnetic momentum. It yields an original fixed-energy weight ceiling
`Gamma/(16+Gamma)` on the ballistic family `c=v/a,U=8c`. At the kinetic
node it also bounds the original zero-energy atom by the same number.
Neither assertion supplies a positive atom or a particle assignment.

## Source orbitals and the complete occupation space

Use the action, magnetic momentum, ground average and retarded clock of
`MOMENTUM_RESPONSE.md`. All expectations in this proof are on that one
state. In the signed source chart define eight canonical orbitals:

    a_0s,
    r_s=6^(-1/2) sum_(|y|_1=1) a_ys,
    a_zs for z in {+/-2e1,+/-2e2,+/-2e3}.             (2)

These are orthonormal one-body combinations on thirteen physical sites.
Their CAR follow from disjoint supports and the normalization of r. The
five orthogonal unused combinations of nearest-neighbour orbitals remain
in the source. A real orthogonal change among those six odd-sublattice
orbitals preserves the eta pair form. Thus r has eta sign -1, whereas
the origin and all six axial orbitals have eta sign +1. This is the
coherent neighbour orbital in the source force identity, with no added
mode or spin-only projection.

For each of these orbitals, the full four-dimensional Fock space splits
as a spin doublet and an eta doublet. With

    P_spin=n_up+n_down-2n_up n_down, P_eta=I-P_spin,
    S+=a_up* a_down, S3=(n_up-n_down)/2,
    eta+=epsilon a_up* a_down*, eta3=(n_up+n_down-I)/2,

the exact identities are `S^2=(3/4)P_spin`,
`eta^2=(3/4)P_eta`, and every product of a spin generator and an eta
generator at the same orbital is zero. Each sign of epsilon gives the
same spin-half algebra in its eta doublet. Local changes of basis remove
these signs within the doublets without changing any spectrum.

Let `Q_xy=S_x.S_y+eta_x.eta_y` and set

    H_sigma=12Q_0r+sigma sum_z Q_0z, sigma=+1 or -1.   (3)

Every local spin/eta type projector commutes with (3). In a fixed type
block the center couples only to leaves of its own type. The other leaves
are inactive doublets. This accounts for every one of the 65536 states of
the eight-orbital CAR algebra, including empty and doubly occupied states.

The minimum of an inactive-leaf block cannot be less than the minimum of
the all-active spin star with the same signed couplings. To see this,
extend any density matrix of its active spins by independent maximally
mixed spins at the inactive leaves. Their one-point spin expectations
vanish, so the expectation of the all-active Hamiltonian equals that of
the smaller active star. The variational principle proves the assertion.
The argument is independent of the signs of the omitted couplings. It
also applies when the strong leaf is inactive. The all-active type block
occurs in the full CAR space, so its minimum is the full-CAR minimum.

## Exact spectra of the two signed stars

The all-active spin Hamiltonian is

    H_sigma=12 S0.S1+sigma S0.J,
    J=sum_(l=2)^7 Sl.                                 (4)

The six weak spins have total spin j=0,1,2,3 with multiplicities5,9,5,1.
For j>=1, couple the strong pair S0,S1 to its singlet or triplet. In the
total-spin F=j block, its two-by-two matrix is

    [[-9, sigma sqrt(j(j+1))/2],
     [sigma sqrt(j(j+1))/2, 3-sigma/2]].               (5)

This follows by writing the triplet restriction of S0 as half the strong
pair spin and using angular-momentum addition. The singlet expectation
of S0.J is zero; its singlet-triplet matrix element has squared modulus
`j(j+1)/4`. The other triplet branches F=j+1 and F=j-1 have energies
`3+sigma j/2` and `3-sigma(j+1)/2`. For j=0, the actual energies are -9
and 3; no F=0 triplet is introduced.

The smaller eigenvalues of (5) are respectively

    (-13-sqrt(529+4j(j+1)))/4,
    (-11-sqrt(625+4j(j+1)))/4.                        (6)

They decrease with j. The other branches and the j=0 states have higher
energies. Hence, on the full CAR space,

    H_plus  >=(-13-sqrt(577))I/4,
    H_minus >=(-11-sqrt(673))I/4.                     (7)

The producer constructs the literal spin matrices in all nine
magnetization sectors of dimensions1,8,28,56,70,56,28,8,1 and verifies
their exact integer characteristic polynomials. It separately constructs
all64 Fock states of a three-orbital signed star, checking all eight local
type blocks against the active/inactive spin matrices for both coupling
signs. The local doublet identities supply the reduction of the complete
eight-orbital Fock space, without a dense65536-state diagonalization.

## The same source force at every momentum

The exact source identities in `MOMENTUM_RESPONSE.md` give

    G(k)=[g0+2g2 sum_j cos(2k_j)] I,
    g0=c^2(9/2-12qstar), g2=c^2 q2,
    qstar=omega(Q_0r), q2=omega(Q_0,2e1),
    S(k)=sum_j cos(k_j)^2 in [0,3].                   (8)

The specified24 proper-cubic lifts fix the origin, preserve r, and
permute the six axial orbitals in (2) with phase+1. Their source average
therefore gives `sum_z omega(Q_0z)=6q2`. Consequently the two endpoints
of the affine function in S are

    G(S=0)/c^2 =9/2-omega(H_plus),
    G(S=3)/c^2 =9/2-omega(H_minus).                   (9)

Equation (7) bounds these endpoints by

    Gamma_plus=(31+sqrt(577))/4,
    Gamma_minus=(29+sqrt(673))/4 <Gamma_plus.          (10)

For the strict inequality, `sqrt(577)>23` implies
`(sqrt(577)+2)^2>673`. Also `sqrt(577)<25`, so Gamma_plus<14.
Convex interpolation in S proves the bound on G in (1). The exact
orthogonal response identity `R*R=G-C^2`, with `C(k)=kappa h(k)`, gives
the residual bound because C is selfadjoint and `R*R>=0`. This uses the
same full ground, full occupation algebra and clock as (8). It assumes
neither antiferromagnetic order nor a numerical ground approximation.
The local spectral minimum in (7) need not be attained by that ground.

## Complete original probabilities and the nodal moment bound

For any normalized source coefficient packet f, the actual retarded
embeddings obey

    v=V0 f, u=V1 f, ||v||=||u||=1,
    L u=alpha v+Y f, <v,Yf>=0, alpha=U/2,
    ||Yf||^2<=Gamma c^2.                              (11)

The original signed measure has total mass one, retaining addition and
removal. Let `w=||1_[-E,E](L)v||^2`. The projection argument of
`LOCAL_CHARGE_ESCAPE.md` gives

    alpha sqrt(w)<=E+c sqrt(Gamma) sqrt(1-w).         (12)

On `c=v_physical/a,U=8c`, divide by c and let a tend to zero at fixed
physical E. Uniformly in normalized packets, including a-dependent
packets, this yields

    limsup w <=chi=(31+sqrt(577))/(95+sqrt(577))<7/15,
    liminf (1-w)>=1-chi>8/15.                        (13)

These are combined original weights. The separate addition and removal
weights of an arbitrary packet need not be one half. Equation (13) does
not derive a growing-box state approximation or any positive retained
low-energy weight.

At the kinetic node the source symmetry theorem `NODAL_CHANNELS.md`
makes the original measure one even scalar probability measure for all16
response components. Its moments are

    M0=1, M2=alpha^2, M4=alpha^4+alpha^2 r2,
    r2=G(node)<=Gamma c^2.                            (14)

Here r2 denotes a squared residual norm, distinct from the orbital r in
(2). For any real b, the nonnegative polynomial `(1+bE^2)^2` equals one
at E=0. Its integral therefore bounds the zero atom z0 from above.
Minimizing the integral gives

    z0 <= min_b[1+2bM2+b^2M4]
        =r2/(alpha^2+r2)
        <=Gamma c^2/(alpha^2+Gamma c^2).              (15)

This is the degree-two Christoffel bound, valid also at r2=0. At U/c=8
it gives z0<=chi. The stronger degree-four optimization needs the actual
sixth and eighth moments; (15) does not supply them or a positive atom.
The nodal necessity and conditional locality results use this upper bound
with the same fixed normalized reference state. An arbitrary change of
reference state is outside those nodal comparisons.

## Exact local readout and executable checks

The actual13-site operators H_plus and H_minus have Hermitian Majorana
expansions with252 terms each and coefficient absolute sum117/2. The
complete Pauli score instrument of `MOMENTUM_RESPONSE.md` chooses term l
with probability `|a_l|/(117/2)`, retains both signs of its measurement,
and reports `(117/2)sign(a_l)s`. Its mean is the corresponding expectation
in (9), on every full occupation input and arbitrary reference. A
source/reference half-trace error eta changes the mean by at most117eta.
Both the positive and negative weak-coupling observables are constructed
from the same source operators. Coherent measurement compilation uses the
declared services and retains their failures; no large native tape is
executed by this finite algebra check.

`nodal_moment_bounds.certificate()` reconstructs the source orbitals,
doublet identities, both signed complete spectra, Pauli scores and outward
rational square-root bounds. It keeps physical masses and species
unassigned. No spectral atom, sharp dispersion, scattering state or
relativistic particle is inferred from a response moment ceiling.
