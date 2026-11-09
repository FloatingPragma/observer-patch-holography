# Local source readouts of the sixth and eighth nodal moments

The centered cubic source at c=1,U=8 supplies exact finite local observables
for its original sixth and eighth nodal spectral moments. The executable
reduces them to 493 and 27,965 Hermitian Majorana words using the specified
ground's actual magnetic translations and spatial symmetries. Their
selected-ground expectations are not evaluated. The coefficient bounds
give quantitative bounds on the same source's omitted-memory moments;
they decide neither the zero-energy atom nor a physical particle spectrum.

## Source and finite local polynomial

Use the full occupation source, invariant ground prescription and physical
cell Bloch chart of `MOMENTUM_RESPONSE.md`, `NODAL_CHANNELS.md` and
`SPECTRAL_RESIDUE.md`. All quantities use its normalized clock c=1,U=8.
The literal hopping magnitude is 1/2. The inherited inversion and proper
cubic average give one even scalar original measure mu at
k_*=(pi/2,pi/2,pi/2), with total addition-plus-removal weight one.
No spin-source ground or isolated finite cube replaces this reference.

Let delta O=[O,K], a_x=a_(x,up), b_x=(2n_(x,down)-1)a_x, alpha=4.
The actual nodal cyclic representation has L a=alpha b. Self-adjointness
of L and the source exponential moment bounds therefore give

    M_(2r)=alpha^2 <b,L^(2r-2)b>,  r>=1.               (1)

The brackets denote original retarded scalar products. They retain both
physical charge branches, with the removal energy reversed in this signed
representation. Equation (1) is an identity in the same nodal cyclic
space. It does not prepare an infinite plane wave.

For n=2r-2 the right side is the expectation of the even local observable

    A_(2r)=16 Re sum_(z in (2Z)^3, |z|_1<=n)
                  (-1)^(sum(z)/2) {delta^n b_0,b_z*}.  (2)

Here Re means the Hermitian part. The anticommutator vanishes outside the
displayed ball by odd CAR locality. The half-integer momentum phases on
even sites are exactly the signs in (2). The finite-range moment extension
in `SPECTRAL_RESIDUE.md` identifies this local expectation with the
canonical nodal measure, including its specified chart gluing. The order
n commutator only uses bonds in the radius-n source ball. Boundary bonds
outside it cannot occur in these n commutators.

For a site with Majoranas gamma_0,...,gamma_3,

    a=(gamma_0+i gamma_1)/2,
    b=(i gamma_0 gamma_2 gamma_3-gamma_1 gamma_2 gamma_3)/2.

The source numerator4K, up to its scalar shift, consists of the onsite
word -8 gamma_0 gamma_1 gamma_2 gamma_3 and, on each signed bond, the two
words i epsilon_xy(gamma_(x,2s)gamma_(y,2s+1)
-gamma_(x,2s+1)gamma_(y,2s)). These are the actual real signed hopping
coefficients. The algorithm uses Gaussian integers and a common power-of-two
denominator. A Hamiltonian word contributes to a commutator precisely when
its Majorana overlap with the odd input has odd cardinality. No large Fock
matrix or rounded coefficient enters this calculation.

The six successive derivatives of b contain exactly
38,302,2246,19130,164288,1442246 words. Contracting the second, fourth and
sixth derivatives in (2) gives the following exact readouts.

| Moment | Hermitian words before quotient | After magnetic translations | After spatial quotient | Sum of absolute coefficients |
|---|---:|---:|---:|---:|
| M4 | 211 | 181 | 19 | 1264 |
| M6 | 11857 | 9781 | 493 | 277408 |
| M8 | 997010 | 843537 | 27965 | 133578340 |

The M6 and M8 words occupy at most four and five sites respectively. Their
supports can be spatially separated within the local ball. The constant
coefficients are 328,8326,278062. These constants equal the corresponding
tracial control values; **they are not the selected-ground moments**.

## Symmetry quotient and complete measurement

For a magnetic translation t, each Majorana at x acquires the sign

    chi_t(x)=(-1)^(sum_j t_j sum_(i<j) x_i).

The harmless projective constant cancels on an even word. For a signed
axis permutation (pi,s), the source lift carries the sign

    xi_pi(x)=(-1)^(sum_(i<j, pi(i)>pi(j)) x_i x_j).

Both signs are checked against the actual signed links. The permutation
of Majoranas also contributes its fermionic sorting sign. Pure inversion
and the proper cubic lifts generate the 48 spatial transformations used
in the reduction. Unit shifts, including odd shifts, are justified by
`MOMENTUM_RESPONSE.md`, "Source, state and translation convention": the
original-Weyl ordinary translation transports to U_j times the constant
hidden-spin rotation i sigma_j. Full-ground hidden-spin invariance turns
the same ordinary-translation average into magnetic-unit-translation
invariance, which the proper cubic average preserves. Inversion is proved
for that same prescription in `NODAL_CHANNELS.md`. No extra state average
is performed.

The executable applies those signs, sorts the transformed word, and
translates its first site to the origin. Orbit representatives combine
equal expectations. A stabilizer with negative sign would force its
expectation to zero. **The reduced polynomial need not equal the original
operator. Its expectation is equal on this invariant reference.** Applying
it to an arbitrary state or to a cut ground requires its own justification.
Independent occupation-matrix controls check the gauge and fermionic
signs before using any invariant-state equality.

For a reduced polynomial sum_j c_j Gamma_j, each Gamma_j is the Hermitian
norm-one word i^(d(d-1)/2) gamma_(j1)...gamma_(jd). Put C=sum_j |c_j|.
Choose setting j with probability |c_j|/C and retain both effects

    E_(j,+/-)=|c_j|/C (I+/-Gamma_j)/2.

Their sum is I. Assign score +/-C sign(c_j). Its original mean is the
desired moment. The identity word has a zero-weight negative record;
that record is retained. An admitted complete-output half-trace error e
changes this mean by at most2Ce. Ground-marginal and native-service errors
must therefore be supplied explicitly. No large native measurement tape
is executed by this polynomial calculation, and no success probability
is divided out.

## Actual memory bounds and the atom test

In the common scalar nodal response write A=alpha^2=16 and use the exact
source decomposition

    L=[[0,alpha,0],[alpha,0,R*],[0,R,D]],
    y=R, s=||y||^2, t=||Dy||^2, u=||D^2y||^2, H=A+s.

The memory measure is even. Multiplication on these actual cyclic vectors
gives

    M2=A, M4=AH, M6=A(H^2+t),
    M8=A(H^3+2Ht+u).                                  (3)

For example L^3 a=alpha H b+alpha Dy and
L^4 a=alpha^2 H a+alpha(Hy+D^2y). Thus the readouts in (2), together with
the source M4 readout, determine s,t,u. These are measured ground
correlations of the source, not freely chosen Jacobi parameters.

The norm bounds in the table and s,t>=0 imply the unconditional nodal
source bounds

    t<=277408/16-16^2=17082,
    u<=133578340/16-16^3=33378201/4.                   (4)

The same ground has s>36/125 by `NODAL_MEMORY.md`. Jensen applied to the
actual positive memory measure consequently gives, for eta>0,

    M_eta=<y,eta^2(D^2+eta^2)^(-1)y>
       >=eta^2(36/125)^2/[eta^2(36/125)+17082].        (5)

This quantitative source bound is O(eta^2). It supplies no positive limit
as eta decreases to zero.

If u>0, source Cauchy gives su>=t^2 and Hu-t^2>=Au>0.
Minimizing the norm of 1+bE^2+dE^4 using (3) gives

    b=-(Ht+u)/(Hu-t^2), d=t/(Hu-t^2),
    Z=mu({0})<=lambda2=(su-t^2)/(Hu-t^2).              (6)

Its improvement over s/H is A t^2/[H(Hu-t^2)]. For an actual memory
measure, u>0 implies t>0. If s=0 then y=0 and Z=0. If s>0,u=0,
self-adjointness gives Dy=0,t=0 and the cyclic three-dimensional block
has Z=s/H exactly; no zero-over-zero quotient is used. If u>0 and
su=t^2, (6) is zero. None of these selected-source equalities is supplied
by the polynomial generator. A positive finite upper bound does not
prove a positive atom. The coefficients in (6) require certified source
moment values or joint enclosures, including their conditioning.

## Onsite identity and independent full-source controls

For Y=[b_0,T], every bond term is odd on one endpoint and an even spin/eta
operator on the other. The onsite action W therefore obeys

    [ [Y,W],W ]=16Y.                                  (7)

The code verifies (7) on the exact seven-site CAR polynomial. This is an
operator identity. The interacting ground need not commute with W alone,
so its W derivation need not be symmetric in the retarded scalar product.
In particular, (7) does not establish ||[Y,W]||_ground=4||Y||_ground or
a many-body charge gap. Cross terms with the actual kinetic action must
be controlled separately. At c=1,U=8 the difference between the two
anticommutator-norm operators is an explicitly nonzero300-word polynomial.

The full64-state source control uses the specified three-site action
c=2,U=3/7 and its entire exact ground projector, of energy -45/14 and
rank two. The coefficient vector (1,0,-1)/sqrt2 lies in the kernel of its
literal one-body hopping. Integer occupation matrices independently
contract the original moments through degree eight and the b moments
through degree six, verifying (1), (3) and (6). This finite zero-hopping
mode is not a continuum momentum label. The control neither replaces
the c=1,U=8 selected ground nor extrapolates its finite spectral atoms.
Its center-site force also satisfies (7) with multiplier alpha^2=(3/14)^2,
while its actual-ground squared-norm difference is -216/805. Thus even
the source's full-ground norm can fail the inferred onsite norm equality.
