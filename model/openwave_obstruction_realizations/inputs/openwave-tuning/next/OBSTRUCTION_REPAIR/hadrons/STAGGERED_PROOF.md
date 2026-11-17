# A finite mass from an interacting colour vacuum

This alternative replaces the old positive number/Laplacian catalogue by a
signed staggered hopping law and its actual filled ground. It demonstrates
a finite mass-generation mechanism. Its three executed coefficients all
miss the proton comparison, and no spin/flavour or physical nucleon
identification is claimed.

There are two sites and three fermionic colours: the full CAR carrier has
64 states. On a single SU(3) tree link with trivial external flux, Gauss
dressing identifies the physical carrier with the complete global colour
singlet. Its particle-number sectors N=0,3,6 have dimensions 1,4,1. The
code constructs their six integer epsilon columns R, checks all eight
total colour generators annihilate R, and proves completeness: the
Cartan-zero subspace has dimension ten and the positive sum of constraint
squares has rank four there. The Gram is diag(1,1,3,3,1,1). Thus this is an
exact physical reduction, not a selected low-energy truncation.

The source is

    H = m(N_left-N_right)
        - t sum_colour(c_left† c_right + c_right† c_left)
        + kappa C_left²,
    t=1/8, kappa in {1/2,1,2}.

The staggered mass m is the actual finite 432-register light coefficient.
The source keeps the two sites and this clock as declared architecture.
Its link Casimir is computed from the literal local colour currents,
C_left²=(2/3) N_left(3-N_left). On the complete physical three-particle
space the normalized matrix is tridiagonal, indexed by n=N_left:

    H_nn = m(2n-3) + (2 kappa/3)n(3-n),
    |H_(n+1,n)| = t sqrt((n+1)(3-n)).

The source code constructs the full 64-dimensional CAR hopping and verifies
its exact intertwiner with this rational Gram representation. Every
physical eigenvalue is enclosed by rational congruence inertia. The lowest
state lies uniquely in N=3. Define B=(N-3)/3. The N=0 and N=6 sectors have
B=-1,+1 and exactly zero source energy, so their actual excitation energy
is -E_ground. This is a genuine difference between eigenvalues of one
Hamiltonian. No energy origin or per-species subtraction changes it.

At m=0 the even-parity vacuum block has diagonal 0 and
D=4 kappa/3-2t with off-diagonal magnitude sqrt(3)t. Hence

    mass_gap = (sqrt(D²+12t²)-D)/2.

For t=1/8 and kappa=1, the square root is exactly 7/6 and the gap is
exactly 1/24. At the actual small m the certified gap is approximately
0.041669169943. The local epsilon creator c_left,0† c_left,1† c_left,2†
and its adjoint give nonzero source-ground records in the B=+1 and B=-1
sectors, respectively. The code evaluates their actual CAR maps and retains
their complementary failure effects. An abstract level is therefore
connected to an explicit local source probe.

A same-clock independent electron rest factor e f†f may use the actual
432 analyzer e. This is an explicit finite tensor construction; it does
not prove a shared spatial Dirac field or a Lorentz limit. With that
conditional electron anchor, the three unchanged kappa members give
0.7807242178, 0.3531114635 and 0.1631551735 GeV. None reproduces 0.93827208943
GeV. The source has no physical proton/neutron pair, no spin-half theorem,
no complete M6 preparation/observation execution and no spatial refinement
claim. Choosing a new kappa using this discrepancy would be a calibration.

The old positive-kinetic obstruction cannot be carried onto this source:
its kinetic term is signed and its vacuum is a filled interacting state.
Conversely, passing its algebra does not validate its physical dictionary.

## The spin-free positive source remains excluded

`kinetic_floor.py` tests the natural intermediate alternative J=0 in the
old source. Let P contain all colocated three-quark configurations and Q
its complete complement. On any path of at least two cells,

    PHP >= 3/8 P,
    QHQ >= 8/3 Q,
    ||QHP||² <= 3/32.

The first is the positive hopping degree, the second the least local
colour cost of a separated three-quark pattern. Each colocated column's
one-hop images are disjoint from those of every other column and have
squared norm at most 6t²=3/32; the code verifies the exact finite source
columns on L=2,3,4,5. The same argument is uniform in path length.
Completing the block square with eta=1/8 gives H>=1/4. With the actual
432/3600 electron coefficients this already exceeds 2.120/1.696 GeV.
The complete exact two-cell ground gives 3.016/2.413 GeV. Thus removing the
installed spin cost alone fails; changing the vacuum/kinetic construction
is a substantive change, as the staggered example makes explicit.
