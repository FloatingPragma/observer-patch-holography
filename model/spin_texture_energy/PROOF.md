# Exact texture energy in the ordered quantum source

## Source and state

Use the same cubic spin-1/2 Hamiltonian and selected ground as
[`spin_canonical_energy`](../spin_canonical_energy/README.md):

```
H=J sum_<xy>(S_x.S_y-1/4),  J>0,
omega(S_x)=(-1)^x m e1,     m>1/1000.
```

The ground is the zero-field limit selected by an infinitesimal positive
staggered field along internal axis 1. Its existence and order are the
matched infrared/ground-selection results and exact integral certificate
in [`spin_goldstone_sector`](../spin_goldstone_sector/README.md). The
spin system is the fixed-spacing source limit; the microscopic centered
contact source has no additional prescribed spin potential in this proof.

Residual rotations about axis 1, reality in the S3 product basis and
cubic bond symmetry give the exact nearest-neighbor covariance

```
C_ab=omega(S_x^a S_y^b)=diag(a,b,b),
kappa=-(a+b),          m^2/2<=kappa<=1/2.                 (1)
```

In particular, the possible antisymmetric transverse covariance contains
one imaginary S2 factor and vanishes in the real state. Equation (1),
including its positive bound, is the inherited short-distance identity.
It does not give numerical values of a or b or a measured stiffness.

For a finite collection of onsite rotations choose spin unitaries U_x
whose action on the spin vector is R_x in SO(3), and set R_x=I elsewhere.
The rotated state has mean direction `n_x=R_x e1` and order magnitude m.
Write its energy excess, summed over every affected bond, as Delta E.
Ground positivity makes Delta E nonnegative. A rotation by pi about
axis 1 at one site changes each of its six bonds by `-4Jb`; hence

```
b<=0.                                                    (2)
```

This uses the exact ground inequality for a local unitary, without a
mean-field ansatz or a bound on a source charge gap.

## The full rotation identity and its lower bound

Directly transforming both spins in each bond gives

```
Delta E=J sum_<xy> Tr[C(R_x^T R_y-I)].                    (3)
```

Let `R=R_x^T R_y`. Since `n_x.n_y=R_11`, the summand equals

```
J kappa(1-n_x.n_y)+J(-b)(1+R_11-R_22-R_33).              (4)
```

For a Rodrigues rotation of angle theta about unit axis u,
`1+R_11-R_22-R_33=2(1-cos(theta))u_1^2>=0`.
Equations (1)--(4) imply the exact bounds

```
Delta E >=(J kappa/2) sum_<xy>|n_x-n_y|^2
        >=(J m^2/4) sum_<xy>|n_x-n_y|^2.                 (5)
```

For rotations about either transverse axis the second term in (4) is
zero. This recovers the inherited exact coplanar result
`Delta E=J kappa sum[1-cos(theta_x-theta_y)]` at arbitrary amplitude.
For general textures, rotations of a local frame about its order
direction can change the state and energy through quantum correlations.
Only a uniform rotation in the ground's stabilizer acts trivially on
the state. Dropping the second term in (4) is an inequality, not a
local gauge identification.

All these statements concern the explicit local-unitary orbit of omega.
They are not lower bounds on all quantum states with the same one-point
magnetization, and do not equal an energy minimized over those states.

## A static long-wavelength limit with error control

Take a nonconstant `C^3` map `R:R^3->SO(3)` equal to I outside a cube
`[-M,M]^3`, and set `R_x=R(x/L)` on the integer lattice, with L>=1.
Let `A_i=R^T partial_i R` and write `A_i v=omega_i cross v`.
Each A_i and its derivative are antisymmetric. Taylor's theorem gives

```
R(x)^T R(x+h e_i)=I+h A_i+(h^2/2)(partial_i A_i+A_i^2)+O(h^3),
Tr(C A_i)=Tr(C partial_i A_i)=0,
Tr(C A_i^2)=kappa[(omega_i^2)^2+(omega_i^3)^2]-2b(omega_i^1)^2.
```

The last expression is nonnegative by (1)--(2). There are O(L^3)
affected bonds. Dividing their summed Taylor error by L gives O(1/L).
The remaining Riemann sum has the same error order. Thus

```
Delta E/L -> E_cont(R)
 = (J/2) integral sum_i {kappa[(omega_i^2)^2+(omega_i^3)^2]
                         -2b(omega_i^1)^2} dx
 >=(J kappa/2) integral |grad(R e1)|^2 dx.                (6)
```

For an explicit conservative error, let M>=1 and K>=1 bound the operator
norms of every first, second and third partial derivative of R. Since
`||C||_1<=3/4`, each bond's third-order trace remainder is at most
`J K/(8L^3)`. At most `3(2M+4)^3 L^3` edges need counting. Each axis
contributes `f_i=(J/2)Tr(C R^T partial_i^2 R)` to the quadratic density.
For any coordinate j, `|partial_j f_i|<=3J(K^2+K)/8<=3JK^2/4`.
A derivative of the total density `sum_i f_i` is therefore at most
`9J K^2/4` in each coordinate,
so a cubic Riemann-sum error is bounded by
`7J K^2(2M+4)^3/L`. Therefore the total error obeys

```
|Delta E/L-E_cont(R)| <=10 J K^2(2M+4)^3/L.              (7)
```

This coefficient is evaluated rationally in `energy.py`. The limit is
static energy on the chosen orbit. It is not a derived time-dependent
nonlinear sigma model or equality of a limiting quantum field energy.

## The ordered orbit and topological textures

Uniform spin rotations generate an S2 orbit of ground states. Rotations
fixing e1 preserve the selected state, whereas any other rotation changes
its nonzero order vector. Its stabilizer is therefore exactly SO(2).
This describes the constructed orbit, without classifying every ground
phase. A translation taking n to -n is a global source symmetry, not a
gauge equivalence that replaces S2 by the projective plane.

Point-defect boundary maps `S2->S2` have integer degree. Compact spatial
textures with fixed direction at infinity are maps `S3->S2`; the Hopf
fibration gives their integer Hopf invariant. These are the usual
homotopy invariants of this explicit ordered orbit, not operators proved
to commute with the microscopic Hamiltonian. On a finite unconstrained
lattice a rotation assignment is continuously deformable to the identity;
no exact microscopic topological superselection rule is inferred.
The sphere degree, Hopf fibration and Hopf invariant conventions are those
of [Hatcher, Algebraic Topology, Sections 4.2 and 4.B](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf).

An explicit compact Hopf family is obtained from

```
q(x)=cos f(r) I+i sin f(r) xhat.sigma,
f(r)=pi for r<=1/4,  f(r)=0 for r>=1,
```

with a smooth monotone interpolation. Its degree has absolute value one:
`-(2/pi) integral_0^infinity f'(r) sin^2 f(r) dr=1`.
Let R be its adjoint SO(3) rotation and `n=R e1`. Composition with the
Hopf projection gives Hopf invariant of absolute value one. R is smooth
and equals I inside the inner ball and outside the outer ball, so (6)--(7)
apply. Since n is nonconstant, `E_cont(R)>0`, and

```
Delta E(R(x/L))=L E_cont(R)+O(1).                        (8)
```

This explicit macroscopic Hopf family's energy expectation diverges in
the fixed-J source clock. Under the existing soft spin clock `K/(Jp)`,
p=1/L, its scaled expectation diverges quadratically. A divergent first
moment alone does not exclude a finite-energy weak spectral component;
no such concentration statement is proved here. The construction does
not exclude a distinct bound state or a microscopic stabilized defect.
Shrinking a continuum Dirichlet texture multiplies its static energy by
its length scale; extrapolating that formula below the lattice spacing
would exceed the domain of (6).

For an isolated radial hedgehog, require the direction
`n_x=x/|x|` on `1<=|x|<=L`, with arbitrary onsite lifts R_x and any
completion outside. Equation (4) makes every bond contribution nonnegative.
On this radial region, `sum_i |partial_i(x/|x|)|^2=2/|x|^2`.
Taylor comparison outside a fixed core gives a summable per-radius error
of order 1/r, hence O(log L) total; the leading radial integral yields

```
sum_(bonds within 1<=|x|<=L) [1-n_x.n_y]=4pi L+O(log L),
Delta E>=4pi J kappa L+o(L).                            (9)
```

Boundary shells change only O(1) at leading order: they contain O(L^2)
edges of size O(L^-2). Equivalently, omit a core radius epsilon L,
take the ordinary smooth Riemann limit on that annulus, and then send
epsilon to zero to obtain the lower limit in (9). This latter argument
alone suffices for divergence and does not need the sharper error rate.
The value at the core is unrestricted. An isolated hedgehog of this
rotated-ground family therefore has infinite energy on the infinite
lattice. For a smooth continuum annulus with degree q on each sphere,
the area-Jacobian inequality gives
`integral_S2 |grad_S2 n|^2>=8pi|q|`; its Dirichlet energy has the same
linear lower bound with factor |q|. This is a statement about the stated
texture class, not every quantum configuration carrying an assigned label.

## Physical boundary

No texture here supplies a defect eigenstate, a sharp mass shell,
scattering or a finite Higgs-source energy. A Dirac operator in a
prescribed hedgehog mass background has a separate spectral problem;
its eigenvalues do not establish that this source creates that background
with finite energy. The fixed-spacing quantum ground, the order-parameter
orbit, and the energy of these explicit textures are specified without
identifying their internal spin with physical Lorentz spin.
