# Source locality and bounded physical-energy windows

For the centered full-CAR source, a positive nodal zero-energy atom is
necessary and sufficient for positive **original probability in a bounded
physical-energy window** for packets in a fixed shrinking nodal neighbourhood.
This equivalence uses one fixed normalized selected reference. The theorem
does not determine the atom, identify a sharp dispersion, or construct a
relativistic field limit.

## 1. Actual source and complete response

Use the centered source, magnetic translations and selected cubic-averaged
reference of [MOMENTUM_RESPONSE.md](MOMENTUM_RESPONSE.md). Its complete
addition/removal response has the doubled GNS representation

    C_t(k) = integral exp(-itE) dnu_k(E)
           = V0(k)* exp(-it L(k)) V0(k),     V0*V0=I16.       (1)

Consequently `||C_t(k)||<=1`. The signed removal energy reverses its physical
time; the combined physical-energy measure is the pushforward by `E -> |E|`.
No restriction to the one-body dynamics, Gaussian state or particle-number
sector enters (1). The canonical every-momentum measure of
[SPECTRAL_RESIDUE.md](SPECTRAL_RESIDUE.md) includes the magnetic node.

The full-CAR cut estimate in
[full_fock_field/README.md](../full_fock_field/README.md), instantiated for
this centered action in [LOCAL_MATTER_READOUT.md](LOCAL_MATTER_READOUT.md), is

    ||tau_t(a_x)-tau_t^(R)(a_x)||
      <= B_R(t)
      <= 6c(2R+1)^2 |t| 2^(-R) (exp(48c|t|)-1).             (2)

The onsite interaction is removed in a product onsite interaction picture.
The even boundary bonds retain norm at most `c`, uniformly in `U`. The cut
field is odd. Its anticommutator with a creation field outside the cut cube
is exactly zero. Thus (2) bounds the original distant anticommutator by
`2 B_R(t)` before taking the state expectation. Ground-state clustering is
not a premise.

## 2. Uniform matrix Fourier tail

Use the periodic two-step cell frame, with `q=2k`, cell index `m in Z^3`,
parity positions `s in {0,1}^3` and two hidden-spin components. For
`n=||m||_infinity>=2`, any two sites `2m+s,s'` are separated by at least
`2n-1`. Cutting at radius `2n-2` in (2) bounds each scalar response entry by

    768 c |t| n^2 4^(-n) exp(48c|t|).

The norm of a 16-by-16 matrix is at most 16 times its largest entry. There
are `24n^2+2<=26n^2` cells in the `n`-shell. Since `n^4<=21*2^n` for every
positive integer `n`, summing shells gives, for integer `N>=1`,

    sum_(||m||infinity>N) ||C_t[m]||
       <= A0 c|t| exp(48c|t|) 2^(-N)
       <= A0 exp(49c|t|) 2^(-N),       A0=6709248.           (3)

For the elementary bound on `n^4`, check `n=1,...,6` and use
`(1+1/n)^4/2<=2401/2592<1` thereafter. The last inequality in (3) uses
`c|t|<=exp(c|t|)`. At `t=0` the tail is zero by CAR.

This absolute Fourier convergence defines a continuous matrix kernel.
The direct-integral contraction holds almost everywhere, hence everywhere
by continuity. Its fixed-time kernel agrees with the canonical measure in
(1), since both are continuous in each Bloch chart and coincide almost
everywhere. This justifies use at the node without assigning an arbitrary
exceptional fibre.

## 3. Bernstein estimate and Bloch gluing

Let `0<epsilon<=1`. Truncation to `||m||infinity<=N`, where

    N = ceil(log2(A0/epsilon)+49c|t|/ln(2)),

gives a matrix trigonometric polynomial `P_N(q)` with uniform error at most
`epsilon` and norm at most `1+epsilon`. The latter estimate uses the
contraction of the original response, rather than a sum of Fourier norms.

For any scalar trigonometric polynomial of degree at most `N`, the exact
interpolation formula is

    T'(x) = sum_(j=1)^(2N) w_j T(x+theta_j),
    theta_j = (2j-1)pi/(2N),
    w_j = (-1)^(j+1)/(4N sin(theta_j/2)^2),
    sum_j |w_j| = N.

The formula follows by evaluating the Fourier modes `|ell|<=N`; the weight
sum follows from the finite cotangent identity. Applying it entrywise gives
`||partial_(q_j) P_N||<=N(1+epsilon)` in matrix norm. With
`d=||k-k_*||_1`, its periodic-frame change is at most `2Nd(1+epsilon)`.

The physical-cell convention includes
`D(k)=diag_s exp(-ik.s)`, tensored with hidden spin. Since
`||D(k)-D(k_*)||<=d`, conjugating the original contraction costs at most
`2d`. Therefore, on a local lift through any Brillouin-zone face,

    ||C_t(k)-C_t(k_*)|| <= 2epsilon+2d[(1+epsilon)N+1].       (4)

The source's reciprocal-cell gluing transports both packet and response.
The physical `cos(k_j) Gamma_j` symbol is not globally periodic in that
frame. The periodic polynomial and the explicit `D(k)` in (4) account for
this distinction.

## 4. Exact outward arithmetic for the joint limit

The original Lorentzian effect is

    Q_eta(k) = integral eta^2/(eta^2+E^2) dnu_k(E)
             = eta integral_0^infinity exp(-eta t) Re C_t(k) dt
             = -eta Im G00(i eta,k),       eta>0.           (5)

Here `Re` is the Hermitian real part. Integrating (4) gives a bound on the
original effect, without an inverse-memory expansion. In real arithmetic,
for `0<d<=1`,

    ||Q_eta(k)-Q_eta(k_*)||
      <= min(1, d[8+4 log2(A0/d)+196c/(eta ln(2))]).         (6)

The minimum uses `0<=Q_eta<=I`. The difference is zero when `d=0`.

The executable bound uses no floating logarithms. Let `L(epsilon)` be the
smallest integer with `A0/epsilon<=2^L`. The exact series
`ln(2)=2(1/3+(1/3)^3/3+...)` gives `ln(2)>2/3`. Thus the integer cutoff

    N(t,epsilon) = L(epsilon)+ceil(147c|t|/2)

has tail at most `epsilon`. Integrating its upper bound gives

    error <= min(1,
      2epsilon+2d[(1+epsilon)(L(epsilon)+1+147c/(2eta))+1]). (7)

In particular, with `epsilon=d<=1`, the weaker rational expression
`d[8+4 ceil(log2(A0/d))+294c/eta]` also bounds the difference.
`nodal_locality.py` returns (7). Its integer logarithm handles exact powers
of two and rational values on either side of them.

## 5. Conditional finite physical-window prescription

Fix one normalized `c=1,U=8` selected reference and scale its Hamiltonian by
`c_ph=v/a`, with `v>0` fixed. Arbitrarily changing normalized ground states
with `a` requires a separate uniform atom premise. Suppose the actual
scalar node measure from [NODAL_CHANNELS.md](NODAL_CHANNELS.md) satisfies

    z0 = mu({0}) = m0/(16+m0) >= z > 0.                    (8)

The paired spin/eta stars of
[NODAL_MOMENT_BOUNDS.md](NODAL_MOMENT_BOUNDS.md) give
`z0<=chi=(31+sqrt(577))/(95+sqrt(577))<7/15`. The code rejects a larger
proposed lower bound using exact rational comparisons: for `0<z<1`, test
`(95z-31)^2<=577(1-z)^2` only when `95z-31>0`. Equation (8) is a
hypothesis, not a computed value.

For every unit packet supported in `||k-k_*||_1<=aK`, set

    tau = max(1,2352vK/z),       eta=a tau/v.

The rational dynamical term in (7) is at most `z/8`. Choose an integer
`j>=0` such that

    2^(-j)[12+4(ceil(log2 A0)+j)] <= z/8,
    a <= a_max = 2^(-j)/K.                                (9)

This finite dyadic search terminates. The added `4` in (9) controls jumps
of the ceiling function for every smaller spacing: for `0<d<=D=2^(-j)`,

    ceil(log2(A0/d)) <= ceil(log2(A0/D))+ceil(log2(D/d)),
    (d/D) ceil(log2(D/d)) <= 1.

Thus the logarithmic term is also at most `z/8`. Since the zero atom
contributes exactly `z0` to every `Q_eta(k_*)`, no convergence rate of that
nodal score is needed. Equations (7)--(9) imply the original packet score
is at least `3z/4` for every `a<=a_max`.

Choose integer `N>=1` with `t_N=1/(1+N^2)<=z/4`. For any complete unit
measure,

    Q_eta <= w(|E|<=N eta)+t_N[1-w(|E|<=N eta)].

Consequently the physical window `[0,N tau]` has original probability at
least `(3z/4-t_N)/(1-t_N)>=z/2`. A vanishing norm-squared packet tail outside
the prescribed cap changes the asymptotic conclusion by a vanishing amount.

The complete readout of
[SPECTRAL_WINDOW_ACCESS.md](SPECTRAL_WINDOW_ACCESS.md) can be joined to
this budget if its total signed-score error is at most `z/8`. Then its
observed score is at least `5z/8`; subtracting that readout error once more
certifies original window probability at least `z/4`. All complementary
records retain their original probabilities. The theorem supplies no
regulator-uniform state preparation or claim that the growing native tape
was executed.

Conversely, the fixed-reference atom-necessity corollary in
[SPECTRAL_RESIDUE.md](SPECTRAL_RESIDUE.md) bounds the limsup original weight
of every bounded physical window by `z0`. This proves the stated equivalence
for bounded-window survival. Neither implication isolates positive energy
away from zero, a spectral graph, a mass shell, or a scattering field.

## 6. Finite source controls

The producer reconstructs the actual 16-row interacting dimer at
`c=1/8,U=1`, its full ground and all energy projectors. Every original
addition/removal component is included. On all four CAR modes it verifies

    M0=I4,       M1=h,       M2-M1^2=I4/4,
    Q_1=q I4,   q in [0.7994456731,0.7994456732].

The positive moment defect rules out treating this original response as a
unitary one-body group. Each local charge has original probability `1/2`;
its omitted complement cannot be restored by conditional normalization.
These are finite source controls of the contraction and probability
conventions. The dimer ground is not substituted for the thermodynamic
reference, and its finite spectral gap does not decide (8). Tests also
check the actual magnetic Gamma matrices and their reciprocal-cell gluing,
the Bernstein interpolation identity, and complete source spectral scores.
