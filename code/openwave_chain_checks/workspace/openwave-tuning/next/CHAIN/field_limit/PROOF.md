# Controlled fields from native current records

The construction attaches a local arithmetic field readout to the frozen
two-copy native source. Its energy-norm error has an explicit bound that
includes integer rounding and every source move. A separate static
refinement theorem compares cell-averaged fields with a smoothed Coulomb
potential on the unit three-torus. The native source is classical. No
transfer from the 384-state coherent charge-energy source is used.

## Source, preparation and clocks

The source is `auditor/native_field_bridge/FINAL_MANIFEST.json`, SHA256
`4597dc9c226d203d890bc3d3a929bbe1be80e2dc57ec03898c77da4fe5b954e1`.
It consists of two mobile one-load copies on the declared cubic routing
of side two. Each copy has 96 ports and 288 seams. Their signs are readout
labels. A known pair of initial positions and every actual current record
determine the signed charge

    rho_i = rho_0 - B (F_plus,i - F_minus,i).

The orientation of B is +1 at the smaller port, -1 at the larger port.
The implementation calls the literal `native` endpoint routine while
replaying each source record, including zero currents. No terminal source
state is an input to the field readout.

The source uses its frozen joint Poisson clock and supplied electric
scheduler. The readout performs K local arithmetic rounds after preparation
and after each observed attempt. These rounds have their own work counter.
They neither change the frozen source rates nor consume an unreported
physical time. A claim that this readout is a distributed field evolving
with finite latency would require an explicit clock and scheduling service.
This result is a deterministic observation instrument attached to the
complete native record. Its dependence on preparation and clock is exact.

The primary finite-cap output includes the `cap_exhausted` flag. The field
at that flag is the field of the recorded charge at the stopping event;
it is not relabelled as the position at the nominal horizon. In the
uncoupled two-copy case, attempts have rate 576. Thus the exact probability
of a conclusive output at time t and attempt cap C is

    exp(-576 t) sum_(j=0)^(C-1) (576 t)^j/j!.

The complementary probability is a retained field/cap output. The
development checks consume the frozen small-cap laws and compare this
formula directly. Dropping the cap flag loses positive probability.

## Existing statements and the adapter

`N20_native_green_field/n20lib.py` supplies the original cubic routing and
the exact acoustic coefficient 1/8. Its proof establishes the native
current identity and distinguishes a full current tape from its gradient
part. `model/native_spatial_refinement/PROOF.md` supplies a complete
twelve-branch Bloch reduction, a quartic acoustic error bound and an
optical-weight bound. `Lean/Screen/DiscreteCoulombGreen.lean` establishes
the finite Coulomb/Thomson characterization on its twelve-port domain;
its formal theorem is not relabelled as a theorem on every cubic graph.
The graph-independent linear algebra used below is proved here.

The frozen field producer supplies an exact rational side-two Fourier
certificate for lambda_2 >= 1/2 and the degree bound lambda_max <= 12.
Its real-valued local stencil has a spectral truncation certificate but
no general library-rounding theorem. The integer stencil below supplies
a finite arithmetic implementation and a deterministic rounding bound.

The integer z register represents the refined real quantity z/2^b. Its
field need not be integer-valued. This controller is an explicit added
arithmetic architecture, distinct from minimizing energy over native
integer face loads. Real Thomson minimization does not imply that an
integer current tape or integer face repair attains the fractional Green
field in a single history.

The N20 identity for the *mean accumulated current* requires uniform
native scheduling. It is not transferred to the state-dependent electric
scheduler. The pathwise charge reconstruction and the local present-field
bound require no uniform-scheduling premise, so they apply to both
the uncoupled and the reciprocal branches.

## Local convergence, moving sources and rounding

Let L=BB*, G=L+ and R=I-L/12. For neutral rho the unique mean-zero
solution is phi*=G rho. Define the field norm by

    ||v||_L = ||B* v||_2.

The spectral certificate gives ||R v||_L <= theta ||v||_L on potentials
modulo constants, with theta=23/24. Every stencil round is

    phi'(x) = [6 phi(x) + sum_(y~x) phi(y) + rho(x)]/12.

The implementation stores integers z_x=2^b phi_x and rounds the numerator
divided by 12 to the nearest integer, with ties toward positive infinity.
The scalar error is at most 1/(2*2^b). Therefore its field-norm error is
at most sqrt(12*96)/(2*2^b) <= 17/2^b. Constants cause no difficulty:
the field norm and the residual ignore a potential's constant part.

For one block of K rounds, alpha=theta^K and

    e_after <= alpha e_before + (408/2^b)(1-alpha).

A native move changes rho by one signed edge column j B_e. Its target
field changes by |j| sqrt(B_e*G B_e). The edge resistance is at most one,
because the unit edge itself is an admissible flow and the gradient
solution minimizes flow energy. Consequently each source attempt obeys

    e_i <= alpha (e_(i-1) + |j_i|)
                   + (408/2^b)(1-alpha).                 (1)

Here |j_i| is zero or one and the sign/copy label is retained in rho.
With phi=0 initially, e_initial^2=rho_0*G rho_0 <= 4. The code starts
from the valid upper bound two, iterates (1) using exact rational
arithmetic, and rounds each displayed upper bound outward. K=32 and
b=40 are development settings. The worst-case steady bound for a move
at every attempt is alpha/(1-alpha)+408/2^b; actual holding records
reduce it. The displayed source-dependent certificate uses every mark.

The Gauss error and electric-energy error follow without a Green oracle:

    ||L phi-rho||_2 <= sqrt(12) e,
    |rho*phi/2 - rho*G rho/2| <= e.

The second inequality uses ||rho||_2<=sqrt(2) and lambda_2>=1/2. A
second certificate is computed from the exact integer residual:

    e^2 <= 2 ||rho-L phi||_2^2.                           (2)

This follows by diagonalizing L on the neutral subspace and using
1/lambda <= 2. The residual in (2) has denominator 2^b and is evaluated
entirely in integers. Floating Green inversion is used only to report
diagnostic errors; it supplies neither (1) nor (2).

## Native transient on the uniform-clock branch

At kappa=0, each seam of each copy has attempt rate nu=1, with the fair
native move/hold coin. Hence the signed first moment obeys

    E rho_t = exp(-t L/2) rho_0,
    E F_signed,t = B*G (rho_0-exp(-t L/2)rho_0).

The second equality follows by integrating the conditional mean current
B*E rho_t/2 and using the continuity identity. It is the continuous-time
form of the N20 result. Its approach to the static initial-charge field
has the certified transient error

    ||E F_signed,t-B*G rho_0||_2 <= 2 exp(-t/4).         (3)

A zero-start K-round local solver on rho_0-E rho_t adds at most
4 theta^K+(408/2^b)(1-theta^K). Equation (3) is an uncapped expectation.
At finite cap the stopped event-count law must be used; the cap comparison
above prevents replacement of that law by the uncapped semigroup.
Equation (3) describes the field of displaced initial charge. The
pathwise solver (1) describes present signed charge. Those are distinct
observables with explicit sources.

## A quantitative static three-dimensional limit

Let the periodic cubic routing have N cells per axis and h=1/N. The
continuum domain is the unit torus with the ordinary Euclidean Laplacian.
For each nonzero integer mode m, put xi=2 pi m, k=h xi, r=||k||_infinity,
q=|k|^2/8, and R4=3 r^4/8. The proved domain is r<=1/2 and q>R4.
The full twelve-dimensional graph symbol H(k) has

    |lambda_0-q| <= R4,       lambda_1 >= beta=4-r,
    1-|<u,v_0>|^2 <= d=r^2/(4-r-r^2/2)^2,
    u=(1,...,1)/sqrt(12).

These are the pinned source theorem's bounds. They are spatial Laplacian
statements, so their use for a classical Poisson equation transfers no
coherent time evolution. Spectral resolution gives

    |u*H(k)^(-1)u - 1/q|
       <= E(k):=R4/[q(q-R4)] + d/(q-R4) + d/beta.       (4)

The contribution from each optical mode is included. With Q=I-uu*, the
block inverse and ||Q(H(k)-H(0))u||<=r give

    ||u*H(k)^(-1)Q|| <= r/[beta(q-R4)].                 (5)

Indeed QH(k)Q>=beta Q and the scalar inverse block is bounded by
1/(q-R4). Equations (4)-(5) control the full static Green response;
an acoustic eigenvalue alone would be insufficient for a port-localized
source.

For a signed native pair, define the unnormalized Fourier charge vector
rhohat_m by summing its actual charge records times exp(-i xi x) at
their physical cell coordinates. Its total charge coefficient is
s_m=1* rhohat_m. The physical cell-averaged graph potential coefficient,
rescaled by 1/h, is exactly

    Phi_h,m = (h^2/sqrt(12)) u*H(h xi)^(-1)rhohat_m.

The powers of h follow from the inverse discrete Fourier factor h^3 and
the potential rescaling 1/h. The limiting coefficient is

    Phi_m = (2/3) s_m/|xi|^2.                           (6)

Thus -Delta Phi=(2/3)rho on the torus, with zero mean. In infinite
Euclidean space the corresponding fundamental solution is 1/(6 pi r).
Equation (6) is a periodic smoothed Coulomb statement; no finite-torus
fit is substituted for an infinite-volume point-source theorem.

For any signed pair, |u*rhohat_m|<=2/sqrt(12) and
||Q rhohat_m||<=2. Combining (4)-(5) gives the uniform error

    |Phi_h,m-Phi_m|
      <= h^2 [E(k)/6 + (3/5) r/[beta(q-R4)]].           (7)

The rational 3/5 bounds 1/sqrt(3). For fixed nonzero modes, E(k)=O(1)
and the cross term is O(1/h). The result is O(h) for arbitrary internal
ports. If the source is explicitly averaged over all twelve internal
ports, Q rhohat_m=0 and the sharper bound is h^2 E(k)/6=O(h^2).
That internal averaging is a separate comparison readout and is labelled
as such; it is not assumed in the primary native pair.

The executable uses the positive finite filter

    K(x)=product_(a=1)^3 (1+cos(2 pi x_a)).

Its Fourier weights are w_m=2^(-number of nonzero coordinates of m),
for |m|_infinity<=1. The neutral pair removes the zero mode. Summing (7)
against these 26 nonnegative weights gives a sup-norm bound for every
observation position, not only the 64 diagnostic sample positions.
The smooth electric readout is E_h=-grad Phi_h, with Fourier multiplier
-i xi. Its Euclidean vector error is bounded by the same sum with each
term multiplied by |xi|; the implementation uses the rational upper
bound |xi|<=2 pi sum_a |m_a|. This is an explicit coarse observation of
the scalar potential, separate from the microscopic seam values B*phi.
The calculation uses 12 by 12 symbols at N=16,32,64; it constructs no
large lattice and executes no held-out N20 acceptance configuration.

The development source is an actual terminal signed-charge record of the
frozen native process. Its coarse cell coordinates lie at multiples of
1/2. Static refinement places each charge at the same physical coordinate
and retains its internal port label for every even N. This defines the
charge-lifting map precisely. A coarse intercell jump does not become a
single fine-lattice jump under this map. No dynamic trajectory or clock
intertwiner between the different carriers is asserted.

## Observable failures and physical boundary

The field and cap comparisons retain quantitative failure modes. Removing
z-dependence from the symbol violates (4) on a mode with nonzero z
momentum. Substituting an optical eigenvalue removes the acoustic pole
and violates the same bound. Dropping the finite-cap output violates
the complete probability law. Current and rounding controls test the
actual source and arithmetic paths.

The construction supplies a controlled classical present-field readout,
an uncoupled transient theorem, and a static smoothed three-dimensional
Coulomb limit. The frozen reciprocal source supplies a dissipative
electric response with installed signs, kappa, beta and energy-dependent
rates. The field readout here does not feed its delayed approximation
back into those rates. A physical common dynamics requires a specified
distributed field update clock, a source-executed response using that
field, and a dynamic refinement theorem retaining its complete output
law. Transfer to the C10 interacting quantum source, transverse waves,
Lorentz symmetry, mechanical inertia, and the electromagnetic coupling
are not established by these static or dissipative results.
