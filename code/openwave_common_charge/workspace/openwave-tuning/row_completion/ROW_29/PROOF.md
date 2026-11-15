# Maxwell fields and native terminal radiative records

The candidate is `native-local12-dirac-noncompact-quadratic-v4-proposed`.
The literal row criterion is “Maxwell recovered; transverse waves at c”.
Here Maxwell means the vacuum Maxwell equations on the radiative sector,
with two polarizations in every direction and speed one in the same length
and completed-source-time units as the admitted Dirac matter construction.
This is a declared noncompact quadratic architecture, with a=1/2 and b=3.
It does not assert a compact-cosine transfer or an empirical coupling fit.
The effective canonical electric coefficient a rho=1/3 is retained.
The microscopic inputs are the actual oriented LOCAL12 face/incidence lists,
the noncompact quadratic energy, those coefficients, the spectral-code
preparation and the native terminal scheduler. No continuum Maxwell stencil
is installed. The oriented curl identity and rho are derived from the
actual face terms. The coefficients were selected to match the matter speed;
that matching is a checked architectural consequence, not an independent
prediction of c or of an empirical electromagnetic strength. The reduced
native word uses the derived limiting frequency with a certified finite-
source error. Its records verify that operational approximation; the full
size-three word separately uses all frequencies of the actual face source.

## Actual oriented complex and the field dictionary

`maxwell.py` reconstructs the actual 36 by 26 face symbol A(k), the 12 by
36 incidence symbol B(k), and the integer translated face terms. The
frozen provider proves B(k)A(k)=0. At k=0, let Ap=A0^+, P=I-A0 Ap and
Q=I-Ap A0. Let H be the orthogonal projection onto the common nullspace
of B0 and A0^T, and let C put one on each of the two inter-cell edges in
each direction. Define W=(2/3)HC. The new exact rational calculation gives

    W^T W = rho I, rho=2/3,
    V^T V = I,
    P D_i Q = W cross(e_i) V^T,

where A(k)=A0+i sum_i k_i D_i+O(|k|^2). In the code `V_transpose` stores
V^T; V is 26 by 3. Both W and V are fixed rational matrices, not fitted
polarizations. In particular A0 V=0, A0^T W=0 and B0 W=0. The exact
oriented factorization fixes the relative electric/magnetic orientation;
the older squared-symbol dispersion identity alone would not do so.
The executable checks the complete polynomial in three independent
momentum variables, using expanded coefficients.

The physical coordinate is x=h m, the source macro-time is h and the
periodic physical box has side one, with L=1/h cells in each direction.
Use the cell observables

    A_phys = W^T theta / h,
    E_phys = -a W^T E / h^2,
    B_phys = V^T A(k)^* theta / h^2.

The A_phys notation is a vector potential and is distinct from the face
matrix. The microscopic Hamiltonian is H=(a/2)||E||^2+(b/2)||A^*theta||^2
on ker B, with the same retained harmonic factors as the frozen source.
Its exact equations in physical time are

    theta_dot = a E/h,     E_dot = -b A A^* theta/h.

On the low radiative subspace the factorization gives, with curl represented
by i cross(p),

    B_phys = curl A_phys + O(h),
    E_phys = -d_t A_phys,
    d_t B_phys = -curl E_phys + O(h),
    d_t E_phys = a b rho curl B_phys + O(h).

Since a b rho=1 and cross(p)p=0, the limit has both divergence constraints
and both Maxwell evolution equations. There are exactly two positive
low singular values, so no longitudinal propagating polarization is
hidden in the statement. The gauge/longitudinal constraints and the three
zero-momentum harmonic pairs remain in the source. Harmonic constant
fields obey the zero-momentum Maxwell equations. They have normalized
free Gaussian references, not a fictitious normalizable free-particle
vacuum. The row experiment observes the nonzero-momentum radiative algebra.

## Quantitative all-direction control

The provider establishes ||Ap||<4, ||A0||<10,
||A(k)-A0||<=eta=11r, and the Taylor remainder <=17r^2/2, r=|k|.
Its left/right skew block rotations have norms <=4eta and remove the
first-order high/low off-diagonal blocks. Their second-order remainder is

    D(r)=8eta^2+32eta^2(10+eta)+17r^2/2.

After these rotations, the comparison has a high block with singular
values >=1/4-eta and the low block iW cross(k)V^T. The low block has two
positive singular values s=sqrt(rho)r and all remaining values zero.
Rank B(k)=12 and B(k)A(k)=0 give rank A(k)<=24 for nonzero sufficiently
small k. Perturbation preserves 22 high values and two positive low ones;
therefore there are no further positive low singular values.

For completeness, the polarization control uses the Hermitian dilation
[0,A;A^*,0], not an eigenvector estimate for A A^* divided by r^2. The
comparison's positive low cluster at s is separated from zero and the
negative cluster by at least (4/5)r. Its separation from the high block is
larger on the admitted domain. For each of the two orthonormal exact positive low eigenvectors,
projection of its eigenvalue equation onto the comparison complement gives
an unwanted component <=D/((4/5)r-D). Combining two columns costs sqrt(2);
the polar-overlap frame bound and extraction of left/right halves together
cost at most another 2sqrt(2). The resulting bound is
4D/((4/5)r-D). Restoring the two skew rotations and their
left/right coordinates is bounded conservatively by 32eta. Thus paired
left and right low frames can be chosen with error at most

    zeta = 32eta + 4D/((4/5)r-D).

`bounds` requires eta<1/16, r+D<1/4-eta, D<r/4 and zeta<1/8.
The paired comparison frames for a unit transverse e are
W e/sqrt(rho) and iV(cross(k/r)e). The positive dilation eigenspace
supplies orthonormal left and right frames even when its two frequencies
split; its frequency operator Omega_h need not be diagonal in this frame.
The eigenvalue bound implies

    ||Omega_h-|p|I|| <= dw=(5/4)D(h|p|)/h.

All constants are outward rational bounds. No finite-lattice eigensolver
is called at L=2^34,2^36,2^38. Those sizes enter a theorem about the
actual finite periodic symbols; the executed native terminal instruments
below have only the observed spectral registers and their clock.

Fourier normalization on the L^3-cell box is L^-3/2=h^3/2. Combining this
factor, the microscopic oscillator amplitudes sqrt(a/(2h|p|)) and
sqrt(h|p|/(2a)), the coarse W,V maps and the displayed h powers gives

    A amplitude = sqrt((1/3)/(2|p|)),
    E and B amplitude = sqrt((1/3)|p|/2).

Paired left/right frame error zeta controls B without dividing an arbitrary
vector error by |k|. The source and limit amplitudes and orientations then
differ, relative to these amplitudes, by at most
64(zeta+dw/|p|), provided dw<|p|/4. This loose bound covers the projection
maps, the sqrt/inverse-sqrt frequency factors and all three observables.
For every fixed bounded nonzero momentum set it tends to zero with h.
Finite sums of such Fourier modes therefore recover Maxwell wave packets;
the sum of their individual weighted bounds controls the sum. This is a
controlled radiative limit, not an interacting or ultraviolet limit.

## Exact free cap, preparation and spectators

The observed algebra contains the two traveling polarizations at one
allowed momentum p=2pi(1,0,0). The opposite momentum, optical modes, other
momenta and harmonic coordinates are spectator factors. The selected
oscillator annihilators define its Fock encoding. A source basis state
|0>, |y>, |z> maps to native basis 00,01,10; 11 is retained padding.
The exact source's observed Hamiltonian on total occupation at most one is

    H_obs = Tr(Omega_h)/2 + dGamma(Omega_h).

The low-frame polar overlap fixes this dictionary. Native preparation is
a retained-failure reset to 00 followed by the emitted native Hadamard on
the selected polarization wire. Under the dictionary it prepares
(|vac>+|one_photon_e>)/sqrt(2). This is a spectral-code architecture;
no unimplemented encoder from an old integer edge register is asserted.
The comparison Hamiltonian is diag(p,2p,2p,0), including its visible
vacuum scalar and the specified padding action. Its distance from the
source Hamiltonian is <=2dw on the physical code, so terminal evolution
differs by at most 2T dw. The comparison word is compiled from outward
intervals for p, with all scalar Walsh terms retained.

A total-occupation cap changes canonical commutators. It does not change
the free Heisenberg equations: P_cap commutes with H, hence
[PHP,POP]=P[H,O]P for each linear field O. Equivalently
[H,a_j]=-sum_l Omega_jl a_l holds exactly after compression. The declared
Hamiltonian is the compression of the full Hamiltonian, including its
zero-point energy; it is not obtained by squaring truncated field matrices.
`source_contract.py` checks the matrix-unit identity and shows that the
incorrect squared-truncation prescription fails. No Weyl kick or charged
interaction is admitted on this cap.

Tracing spectators is an exact observable reduction. For every normalized
trace-class rho on observed modes, spectators and an arbitrary reference,

    Tr_s[(U_obs tensor U_s) rho (U_obs^* tensor U_s^*)]
      = U_obs (Tr_s rho) U_obs^*.

It holds also for subnormalized states, correlated spectators and each
observed instrument branch. Spectator zero-point phases cancel in this
identity; no infinite-energy subtraction is being used. All finite-box
coexact vacua and the three Gaussian harmonic references are normalized,
so the product reference is a valid instance. The packet never calls
those factors a flat or delta-normalized state. A full local field has
spectator contributions: the reported means are its selected radiative
Fourier component, or the full mean when spectators have zero first
moments. Spectator fluctuations are not claimed to vanish. The explicit
finite correlation control tests the trace identity and rejects replacing
a correlated joint law by the product of its marginals.

## Native measurements, common clock and costs

`radiative_limit.py` emits and independently executes every native pulse
for T=1/8,1/4; both polarizations; positions x=0,1/4; both field X/Y
settings; and both clock X/Y settings. The position phase is an apparatus
setting for the chosen physical field quadrature. For e=y or z its sign
observable has ideal means X=cos(2pi(x-T)) and Y=sin(2pi(x-T)). Thus
A=e amplitude_A X, E=-e amplitude_E Y and
B=-(e_x cross e) amplitude_B Y. These are the compressed source field
first moments under the dictionary, with the stated source errors.
The source-free Maxwell equations above derive their dispersion and
orientation; the quadrature labels alone are not used as proof.

The same phase-clock attachment advances by h/8 after each source macro.
For T/h macro powers, its terminal phase is T/8. Exact group composition
coalesces the free source and clock rotations into the emitted terminal
word. This is a fresh terminal-only implementation with its own native
error certificate. It does not transfer an old per-macro pulse count or
claim to execute T/h individual failures and event reads. Intermediate
observations or moving-current interventions are outside its contract.
An actual split-clock word and the coalesced word agree in the control.

On the three-dimensional physical code, the measured bit X/Y compresses
to the appropriate vacuum/one-photon field matrix. The other polarization
has zero field mean and is mapped to a fair sign law by this two-outcome
POVM. Its dilation can populate the padding output during analysis; that
output is retained rather than discarded. This implements a field first-
moment measurement, not a claim to reproduce the full continuous electric
field spectral distribution.

All eight basis outcomes are retained for each native C/E pointer law,
including padding. Each of the three one-bit fine native effects is
stored, and their tensor product with the unnormalized success vector
specifies every fine joint output. First export complements and reset
timeouts are additional outcomes under the frozen complete instrument
contract. Three blank and three pointer-ready reset failures are priced.
No probability is divided by success. The exact complete-record TV bound
is the sum of word/pointer error, reset failure and 2T dw. A normalized
field sign has expectation error at most twice this TV bound, to which
the physical-field dictionary error above is added. Native precision is
fixed in these executions and leaves a nonzero floor; the continuum
source limit itself has error tending to zero. The primitive compiler's
precision refinement can be used for a joint zero-error limit, but it is
not claimed to have been executed here.

The experiments recover the same speed one from the field phase and
clock quadratures. They distinguish doubling b, which gives speed sqrt(2).
Other controls reject a reversed face entry (BA fails), a longitudinal
three-scalar substitute and deletion of photon frequency. The size-three
full cap-one implementation separately retains all 646 coexact modes,
377 padding states and 2048 field/clock outcomes. It is a finite native
member with its own exact frame certificate, not evidence that a huge
spatial tensor product was explicitly allocated.

The supplied coefficients, spectral-code preparation, radiative readout,
terminal schedule and native clock are architecture inputs. What is
derived is their actual oriented Maxwell limit, its physical field
normalization, two transverse branches, common speed and bounded native
records. Independent full-criterion reviews remain required for a row
validation decision. No producer self-approval or official-row change is
made by this packet.
