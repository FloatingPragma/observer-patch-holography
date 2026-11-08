# Source-generated composite fermions and their spectral remainder

The centered signed Hubbard source supplies an exact two-block retarded
response at every volume. Its second block is determined by source bond
correlations, and its omitted coupling has a positive computable Gram
matrix. The two-site response closes; a full three-site ground gives a
strictly nonzero remainder. These are fields and correlations of the same
occupation source, not additional independently postulated fermions.

## Canonical operators on the original occupation space

Use precisely the matched source chart, action and full invariant ground
from `CHARGED_SPECTRAL_BANDS.md`. Let a_(x,s) be its annihilators,
N_xs=a_xs* a_xs, alpha=U/2, and

    b_(x,s)=2(N_(x,opposite)-1/2) a_(x,s),
    D=product_x exp(i pi N_(x,up)N_(x,down)).

The exact identity b_i=-D a_i D* proves that the b-family satisfies the
canonical anticommutation relations on the whole Fock space. The a- and
b-families are **not jointly independent canonical families**:
`{a_i,b_i*}=2N_(x,opposite)-1` is a nonzero operator, and their mixed
same-site opposite-spin anticommutator is a spin-flip operator. No modes,
flavours or species are appended to the source.

D preserves total number and hidden spin. Hence b has the same internal
charge and spin transformation as a. The source contact compiler supplies
the local phase `diag(1,1,1,-1)` on all four onsite inputs. The executable
contracts its actual ten-pulse word, retaining all 36 first-export failure
slots. Conjugating a declared complete odd addition/removal instrument
by this phase gives the b instrument. Since b* b=a* a, its original
addition/removal weights on the invariant ground are both 1/2. The
fermionic transfer/reference service required for that odd instrument
is supplied; the phase calculation does not derive it from unread repair.

## Orthogonal response spaces on the same ground and clock

For a finite full-ground density rho, use its Hilbert--Schmidt
purification and the signed retarded generator

    Lret=(K-E0) direct-sum (-conjugate(K-E0)),
    V0_i=(a_i* sqrt(rho)) direct-sum conjugate(a_i sqrt(rho)),
    V1_i=(b_i* sqrt(rho)) direct-sum conjugate(b_i sqrt(rho)).

The scalar product is the original expectation of the CAR
anticommutator. Consequently V0*V0=V1*V1=I. Particle-hole and hidden-spin
invariance give the actual onsite one-body density I/2, so V0*V1=0.
This is orthogonality of cyclic response vectors in a specified state,
not an operator assertion that the mixed CAR vanish. The two retarded
components retain the addition and removal weights separately; the
negative component represents the reversed removal time, not a negative
physical excitation energy of K-E0.

Direct source differentiation, without a state approximation, yields

    Lret V0 = V0 h + alpha V1.                             (1)

Here h is the exact signed one-body hopping matrix. Define

    C=V1* Lret V1,
    R=(I-V0V0*-V1V1*) Lret V1,
    W=[V0,V1],    A=[[h,alpha I],[alpha I,C]].

Then

    Lret W = W A + [0,R],
    R*R=V1* Lret^2 V1-alpha^2 I-C^2 >=0.                 (2)

All quantities in (1)-(2) come from the same actual action, ground and
clock. The source moments also give the independent formulas

    M2=h^2+alpha^2 I,
    C=(M3-h^3-2alpha^2 h)/alpha^2,
    V1*Lret^2V1=(M4-h M3-M3 h+h M2 h)/alpha^2,
    Mk=V0*Lret^kV0.                                      (3)

For U>0 these agree exactly with the direct b-operator construction.
The U=0 free source is its separate elementary response and does not
require division by alpha.

Equivalently the original fourth moment resolves the exact memory Gram:
`M4-[I,0]A^4[I;0]=alpha^2 R*R`. Thus a positive remainder is visible in
the actual original response moments, not only in an auxiliary operator
description.

## The second block is an observable bond correlation

The full finite centered action commutes with hidden spin SU(2) and the
staggered eta SU(2) from `field_symmetry`. Their generators are

    S_x^+=a_(x,up)* a_(x,down),
    S_x^3=(N_(x,up)-N_(x,down))/2,
    eta_x^+=epsilon_x a_(x,up)* a_(x,down)*,
    eta_x^3=(N_x-1)/2.

The normalized full ground projector is invariant under both groups.
Tensoring the tracial exterior, taking the stated translation averages
and the same weak-star subsequence preserves these local symmetry
identities. This invokes the source's full ground selection; an arbitrary
chosen spin-polarized state would not have these identities.

In this state the exact compression is spin diagonal, has zero onsite
and non-neighbour entries, and on each neighbouring bond obeys

    C_(x,s),(y,t) = delta_st 4 h_xy
                      omega(S_x.S_y+eta_x.eta_y).          (4)

For clarity, before using spin/eta invariance the same-flavour coefficient
divided by 4h_xy is

    omega[(N_(x,opposite)-1/2)(N_(y,opposite)-1/2)
      - a_(y,s)* a_(x,s)
          (a_(x,opposite)* a_(y,opposite)
            + a_(y,opposite)* a_(x,opposite))].

The first term is the sum of the spin-3 and eta-3 correlations after
the mixed-vector expectations vanish. The two other terms supply twice
each corresponding correlation. The source signs epsilon_x epsilon_y=-1
on a bond are essential for the eta term. Separate rotational invariance
turns the three equal components into (4). The onsite U contribution
vanishes because it is the mixed V0/V1 overlap. Gauge, spin and
particle-hole symmetry eliminate the onsite kinetic expectation.

The local operators S_x.S_y and eta_x.eta_y act on mutually exclusive
two-site occupation blocks; their sum has norm at most 3/4. On maximum
degree d this proves

    ||C|| <= 3 c d/2.                                    (5)

Each coefficient in (4) is an even two-site source observable, accessible
through the declared complete local coherent readout. No scalar
replacement for its measured expectation is inserted. Homogeneous or
isotropic coefficients require the corresponding symmetries of the
chosen state; (4) itself holds bond by bond.

## A volume-independent residual bound and controlled spectral error

Remove the exact onsite term from the second derivative:
`Y=Lret V1-alpha V0`. Its column is the retarded vector of the odd
operator `[b_i,T]`, supported at the site and its neighbours, and

    ||[b_i,T]|| <= 5 c d/2.

This follows by differentiating `2(N_opposite-1/2)a`: the hopping
annihilator term costs cd/2, and the opposite-flavour number-current
term costs at most 2cd. Two such odd operators anticommute if their
supports are disjoint. Each site has at most 1+d^2 sites within distance
two, with two flavours per site. The row-sum bound on their positive
anticommutator Gram therefore gives

    ||Y|| <= 5 c d sqrt(1+d^2),   ||R|| <= ||Y||.           (6)

The bound is independent of volume and U. Together with bounded h and
C it makes the off-diagonal compression bounded even in the ground
representation of the infinite lattice. Subtracting this bounded
selfadjoint off-diagonal operator shows that the orthogonal compressed
generator Lperp is selfadjoint on its natural domain. This justifies the
same Schur formula on that representation, without a global clock-norm
bound.

For nonreal z, put eta=|Im z|. Exact elimination gives

    W*(z-Lret)^-1 W
      = (z-A-diag(0,R*(z-Lperp)^-1 R))^-1,
    ||W*(z-Lret)^-1 W-(z-A)^-1|| <= ||R||^2/eta^3.         (7)

The memory term has norm at most ||R||^2/eta, and both exterior
resolvents have norm at most 1/eta. Equation (7) is a controlled
broadening-dependent approximation to the actual source response.
Its imaginary part bounds the Lorentzian-smoothed spectral density by
||R||^2/(pi eta^3). It is not a claim that unsmeared spectral atoms or
mass poles are close when eta tends to zero.

There is also a better short-time bound for the original a response:

    ||exp(-itLret)V0-W exp(-itA)[I;0]||
        <= ||R|| min(|t|,alpha t^2/2).                    (8)

Duhamel gives the integral of `R [0,I]exp(-isA)[I;0]`.
That transition block is bounded both by 1 and by alpha|s|, by comparison
with diag(h,C). Integrating either bound proves (8). Its norm bound
controls the original retarded vectors and their bounded matrix elements;
comparison of bounded quadratic forms costs at most twice that bound.
The conjugated, reversed removal component does not itself define a
physical preparation channel for arbitrary unknown coefficient inputs
and references. No complete physical-channel trace bound or native
implementation of the compressed A dynamics is inferred. The separate
actual addition/removal instruments retain their original probabilities
and complementary outcomes, without division by success probability.

Equations (2), (7) and (8) feed the existing bounded source-response
hierarchy with actual operators and measured residuals. Further columns
come from R itself. Setting R to zero without such a calculation would
alter the source dynamics.

## Exact complete-ground controls

`composite.py` reconstructs every number block of the signed full source.
For each control it certifies K-E0 positive by rational elimination,
checks the full nullity, and constructs the normalized entire ground
projector by exact kernel arithmetic. No floating eigenvector is trusted.
All operator CAR, first derivatives, spin/eta symmetries, response Grams
and bond readouts are checked on the full occupation matrices.

For the two-site control c=3/2,U=4, the full ground is unique with
E0=-9/2. Its density has integer denominator 20. The compression is
C=-3h and R*R=0, agreeing with the general dimer formula in
`CHARGED_SPECTRAL_BANDS.md`. The actual dimer closure has no -h lower block.

For the three-site open chain at the explicitly declared coefficients
c=2,U=3/7, the full 64-state ground has energy E0=-45/14 and degeneracy
two. Its normalized projector has denominator 1380. Exact arithmetic gives

    C=-(179/115)h,
    R*R on each spin flavour = (1/13225)
        [[48804,0,-30546],[0,173508,0],[-30546,0,48804]],
    tr(R*R)=542232/13225,
    ||R||^2 <= 173508/13225.                              (9)

This R*R is positive definite; its displayed eigenvalues per flavour are
18258/13225,79350/13225,173508/13225. The two-block source response is
therefore not dynamically reducing. The witness uses specified rational
couplings for an exact independent full-ground control, rather than
relabeling it as the readable-default experiment or extrapolating its
finite poles to a continuum. The operator identities and bounds above
hold at the readable defaults as well.

The construction supplies a canonical composite operator, its actual
state-dependent kinetic block and controlled memory. A relativistic
continuum requires the source momentum-resolved response and the
appropriate original-weight spectral limit of this same hierarchy.
Neither two internally charged CAR descriptions nor their off-diagonal
U/2 coefficient identifies additional species or a relativistic mass.
