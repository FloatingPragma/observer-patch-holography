# A pixel-dependent mass law realized in complete native memory

## Readable design and source domain

The retained minimal complete public response (MCPR) construction specifies
eight register dimensions, source monomials, response signs, a phase offset,
an amplitude and a determinant normalization. These are simulator design
inputs. Their selection was informed by charged-lepton data. The execution
in this package reads no charged-lepton reference values and performs no
fit to observed masses. This separates runtime derivation from historical
calibration; it does not turn that calibration into a prospective prediction.

The two supplied source-coordinate pairs are read from the forward pixel
certificate and the measured-endpoint projection certificate. Their decimal
values define exact conditional inputs to this construction. The interval
proof here encloses the resulting formulas; it does not prove that either
coordinate pair is selected by a complete physical endpoint calculation.
In particular, the measured-endpoint branch inherits measured-alpha input.

Write the conditional masses as t_i=m_i/v>0 and lambda_i=t_i/sum(t).
The producer calls the retained MCPR response and readout functions directly.
The rational checker independently encloses their normalized shape:

    chi=(alpha_U^2/512)/(1-alpha_U/77),
    zeta=((P-phi)^2/(21*pi))/(1-alpha_U/27-alpha_U^2/135),
    r_j=1+sqrt(2)*exp(-chi)*cos(2/9+zeta+2*pi*j/3),
    lambda_i=sorted(r_j^2)/sum(r_j^2).

The absolute normalization is

    kappa=(alpha_U/50)/(1+alpha_U/31+alpha_U^2/310),
    t_i=exp(kappa)*r_i^2/(2*6^14*product(r_j^2))^(1/3).

Exact rational Taylor bounds prove positivity, ordering and sum(t_i)<1
at both declared coordinate pairs. The last result follows by bounding
the cube of the displayed sum, so it uses no numerical cube-root test.

## Exact inverse of the complete memory coefficient

Let P2,P3,P10 be the actual native rank-one family projectors. Direct exact
matrix multiplication gives tr(Pi Pj)=1/5 for distinct indices in this set,
and their spanning-vector Gram determinant is

    d=2/5-2*sqrt(5)/25 > 0.

For Y=u P2+v P3+w P10, its three elementary eigenvalue invariants are

    tr(Y)=u+v+w,
    e2(Y)=4*(uv+uw+vw)/5,
    det(Y)=d*u*v*w.

Consequently a positive triple (u,v,w) consisting of all three roots of

    z^3-z^2+(5/4)*e2(lambda)*z-e3(lambda)/d=0

makes Y positive with exactly the eigenvalues lambda_i. The ordering is
chosen as w<u<v, with u assigned to P2, v to P3 and w to P10. Three
disjoint rational intervals have alternating polynomial endpoint signs
for both source-coordinate pairs. The degree is three, so these intervals
exhaust the roots, including multiplicity. This is an existence and
exhaustion proof, not a numerical eigenvalue match.

For the complete reset2 controller put a=p(1-q), b=pq and D=1-a^2+ab.
Its unconditional C-duty coefficient is

    Z_C=[(1-a)P2+a(1-a)P3+ab P10]/D-(1-p)P2.

The subtracted term accounts for all disabled records. It does not
condition on a successful event. Its trace is p. Define

    A=v(v+w),  B=v(1-v+w)+(v+w)w,
    -A*p^2+B*p-w=0,
    q=1-v/[1-p(v+w)].

Reduction of all three coefficient differences by this quadratic proves
the exact identity Z_C=pY. Rational endpoint signs give two distinct
solutions satisfying 0<p<1/12 and 0<q<1/6. The denominators are positive
on the certified boxes. The two admissible branches are retained.

The native left and right charged-lepton support Casimirs are 9/8 and 3.
The complete current probes with independent durations

    theta_L=asin(sqrt(12p))/sqrt(9/8),
    theta_R=asin(sqrt(6q))/sqrt(3)

therefore implement these rates. Each probe retains thirteen outcomes.
There are 169 fine current pairs times 128 native repair records per old
memory input, or 512 records after the gauge-covariant binary grouping.
All twelve old memory inputs and all failed/disabled records are retained.
An independent stationary solve explicitly sums this complete binary law.

Choose the feedback gain g=sum(t_i)/p. Then

    g Z_C=sum(t_i)Y,
    spec(g Z_C)={t_1,t_2,t_3}.

The spectral theorem supplies a unitary family frame taking this matrix
to diag(t_i). This frame acts on the gauge-trivial family factor and hence
commutes with all twelve actual native gauge currents. It preserves the
canonical kinetic form. This establishes the same finite mass operator
up to that allowed frame, rather than equality of three labels alone.
It does not establish an additional alternating-group family action.

## Correlated weak-pulse dynamics

On each C-enabled record apply exp(-i h g X tensor P_new); on every other
record apply identity. This pulse service is part of the chosen simulator.
Compilation of that service into a prescribed smaller native gate alphabet
is not asserted. Let T be the complete classical transition and Psi_h the
complete classical/quantum channel. On operator-valued memory blocks use
the sum of trace norms; an arbitrary untouched finite reference is allowed.

Since a failed C event resets every input to port2,
T=(1-p)R+pS with S stochastic. For every zero-sum block vector e,
||Te||_1<=p||e||_1. Let Pi z=pi sum(z), where pi is the stationary
distribution; ||Pi||<=1 and ||I-Pi||<=2. The controlled-unitary bound gives
||Psi_h-T||<=2gh. For an arbitrary normalized initial classical/quantum
state z_0, including correlations with the reference, set
e_k=(I-Pi)z_k. Thus

    ||e_k||_1 <= 2p^k+4gh/(1-p).

The reduced step on a stationary product input has first derivative
-i[g X tensor Z_C, .]. Its distance from the unitary step with that
generator is at most 4g^2 h^2: each Taylor remainder is at most 2g^2 h^2,
using the double commutator and unitary trace-norm contraction. The
additional reduced error on e_k is at most 2gh||e_k||_1 because sum(Te_k)=0.
Telescoping N steps, with h=T_phys/N where T_phys denotes dimensionless
model time, yields

    ||sum(Psi_h^N z_0)-U_T(sum z_0)||_1
    <= 4g*T_phys/[N(1-p)]
       +(8/(1-p)+4)*g^2*T_phys^2/N.

For z_0=pi tensor rho the first term vanishes. The bound tends to zero
uniformly on states and finite references, so the actual correlated
pulse service realizes the stated Hamiltonian in this weak-pulse limit.
At finite N the unread channel need not be unitary; no matrix logarithm
of the full memory channel is used. The time and gain are dimensionless
simulator coordinates, not a derived clock or a physical energy unit.

## Physical boundary

This construction joins the conditional pixel-dependent charged coefficient
to native current probes and complete memory dynamics. It supplies no
identification of physical interacting poles, no renormalization map to
charged-lepton measurements, no complete quark or neutrino spectrum, and
no theorem that an occupied bound state creates topological winding.
It does not derive a relativistic continuum or an Einstein limit. Those
claims require the same source to supply their states, observables and
controlled limiting operations.
