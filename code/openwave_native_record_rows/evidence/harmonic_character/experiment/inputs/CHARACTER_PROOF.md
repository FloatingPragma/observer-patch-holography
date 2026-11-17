# A connected native boundary counter and a positive electric character sector

This is a theorem-only N+ construction proposal, with no new source
execution. It follows the compatible-flux conditional theorem and its
finite Maxwell attachment. It supplies an explicit record connection for
the harmonic boundary reference and a positive limiting charged sector.
It does not equate a real Euclidean affine shift with a Lorentzian state,
claim full finite-word Maxwell reflection positivity, or supply target
inertial mass/G. Section7 proves a stronger, separate finite positivity
theorem for the harmonic counter block transfer itself.
No source experiment may run before a separate task and preregistration.

## 1. The same shell harmonic coordinate

Fix a finite relative cubical shell and the source-derived rational values
u=1+h, M_xy=u_xu_y, v=1/u. Here h is the exact finite-graph harmonic hit
law obtained from the pinned native branch kernel, not a measured value
fitted to a desired stress. Set epsilon=dv and kappa=epsilon^T M epsilon.
Assume the nontrivial two-boundary voltage, so kappa>0. The preceding
attachment proves C epsilon=0 and G^T M epsilon=0; the relative harmonic
space has dimension one. The gauge quotient splits into that coordinate
and the positive transverse oscillator complement.

In the ideal compatible-flux Gaussian limit the harmonic increment on
time link t is

    ell_t=kappa^(-1) epsilon^T M e_t.                     (1)

It is the change in x for a harmonic spatial potential x epsilon. The
harmonic action is (kappa/2) sum_t ell_t^2. The variables ell_t are
independent centered Gaussians with variance1/kappa, independent of the
transverse modes. This follows by the M-orthogonal decomposition of the
same quadratic action; it does not import a second scalar source.
Temporal scalar-potential gradients drop from(1) because G^T M epsilon=0.

On each finite member all coefficients in(1) are rational. The finite
source Gaussian readout has e_t=sqrt(8/n) Z_t with integer cochains Z_t.
Choose an explicitly computed common denominator D and integer vector a
such that kappa^(-1)epsilon^T M=a^T/D on the electric edges. Therefore

    ell_t=delta_n a^T Z_t,       delta_n=sqrt(8/n)/D.       (2)

No rounding of individual increments is needed in this rational case.
If an implementation instead uses approximate coefficients, its complete
phase error must be bounded by |p| times the accumulated increment error;
that is an additional budget, not silently zero.

## 2. A uniform counter tied to every accepted source increment

Choose an integer Q>=2 and a cyclic counter c in Z/QZ. Initialize c0
uniformly using native current bits from a disclosed persistent one-unit
two-port pair. Draw b=ceil(log2 Q) bits, reject integers>=Q, and repeat.
The probability of accepting one draw isQ/2^b>1/2, so expected bit calls
are below2b. All rejected initialization words are retained. If capped at
R attempts, failure is retained and has probability at most2^-R; conditional
on success the counter is exactly uniform. There are no fresh per-bit
unit preparations, and the endpoint-law mutation remains visible in these
movement bits.

An entire compatible-flux source word contains all its time links. After
its integer compatibility test accepts, update

    c_(t+1)=c_t+a^T Z_t mod Q                             (3)

for every one of its time links, retaining the pre/post counters and exact
integer sums. A rejected incompatible source word leaves the counter at
its prior value and is still retained in the archive. Such rejected words
are not miscounted as samples from the conditional field law. Initialization
is independent of the source increments; every subsequent counter position
is connected to the actual accepted source history by(3). It is not a
disconnected Gaussian boundary variable.

For multiple independent field-word samples, either disclose fresh uniform
counter initializations and count them, or retain the counter between
samples: addition of an independent increment to a uniform cyclic value
preserves uniformity. In the latter construction conditional counter values
need not be independent between samples, but all character differences
below depend only on the new increments. No unproved IID counter assumption
is needed.

The physical circumference is B=Q delta_n, and the counter position is
x=delta_n c modulo B. The native words thus sample a wrapped version of
the same harmonic increments. Choosing Q as n changes is part of the
observer's declared boundary apparatus, not a modification of native.py.

## 3. Positive limiting transfer and the character sector

For the Gaussian harmonic limit on the circle R/BZ the wrapped one-step
kernel is

 k_B(y-x)=sum_(j in Z) sqrt(kappa/(2pi))
                    exp[-kappa(y-x+jB)^2/2].              (4)

It integrates to one. The normalized characters

 chi_k(x)=B^(-1/2) exp(i p_k x),       p_k=2pi k/B

are eigenvectors, with exactly positive eigenvalues

    lambda_k=exp[-p_k^2/(2kappa)].                         (5)

This follows by integrating the wrapped Gaussian against the character
and unwrapping the real line. Every eigenvalue is positive, so the full
wrapped convolution is a positive self-adjoint transfer. The rank-one
projection onto chi_k commutes with it. Tensoring with the independently
reconstructed positive transverse transfer gives a positive quantum state
in that harmonic character sector. No normalizable free-particle vacuum
on the real line and no global-time zero-mode deletion is required.

The physical canonical momentum and energy in this sector are exactly

    <p>=p_k, Var(p)=0, H_harm=p_k^2/(2kappa).

Consequently the real Lorentzian harmonic electric field has

    <E_harm>=(p_k/kappa)epsilon,       B_harm=0.            (6)

This is a source-normalized state target. The circle/sector labels are
supplied preparation choices; their existence does not derive a physical
charged source or its quantization law from OPH.

Choose integer N->infinity, Q=nearest integer to2pi N/delta_n, and
k=nearest integer to kappa N, with deterministic tie rules. Then
|B-2pi N|<=delta_n/2 and

 |p_k-kappa| <= (pi+kappa delta_n/2)/B =: e_p.              (7)

For sufficiently large n,N, Q>=2 and B>0. The mean-field error is at most
(e_p/kappa)||epsilon|| and the harmonic energy error from kappa/2 is
at most e_p+e_p^2/(2kappa). These are arithmetic preparation errors, fixed
before observations, not stress or G fits. Taking n->infinity at fixed N
first, then N->infinity, is a valid diagonal; the compatible-word and
native cap error budgets must be recomputed for every member.

## 4. The boundary character is an observable of the connected records

For J time links, the normalized boundary factor is

 w_k=exp[-2pi i k(c_J-c_0)/Q]
     =exp[-i p_k sum_t ell_t].                            (8)

The equality is exact modulo the counter period. Its real and imaginary
parts are bounded public-record readouts. The initial uniform reference
cancels from(8), while the recorded update law proves that it is tied to
the same flux history. The Gaussian-limit normalization is

    alpha_k=E[w_k]=lambda_k^J>0.                          (9)

For a bounded observable of the complete source word, its phase-weighted
ratio E[w_k O]/alpha_k is the Euclidean boundary matrix element in the
character state. The corresponding quantum interpretation applies to
observables with the usual OS ordering and reflection prescription; an
arbitrary complex weighted classical variable is not thereby a positive
measurement probability. Positivity follows from the limiting transfer
and its spectral projector, not from pretending w_k is an acceptance flag.

In particular E[w_k ell_t]/alpha_k=-i p_k/kappa. Multiplication by the
canonical Wick factor i gives the real physical mean in(6). For t!=s,
E[w_k ell_t ell_s]/alpha_k=-p_k^2/kappa^2; its Wick factor gives the positive
physical electric product. At t=s there is the extra Gaussian contact
1/kappa. Separated-time products or the exact transfer energy(5) avoid
misreading that contact as negative physical electric energy. The phase
preparation is therefore distinct from the real affine Euclidean mean
used in the earlier static Dirichlet theorem.

At finite n the counter and phases are actual functions of retained
native records. The conditional full-moment and exponential-moment limit
proves convergence of(8) and polynomial phase-weighted numerators. To see
uniform integrability directly, the phase has modulus one, and the
conditional Gaussian tail from the compatible-word theorem controls all
polynomial tails. On compact sets the phase is continuous; counter
wrapping causes no discontinuity in its character. Finite-word or arbitrary
finite-mode reflection positivity is not claimed. The positive state is
the controlled Gaussian-transfer limit.

## 5. Error and cost accounting for the complex normalization

Let eta_den bound the complete error in the estimated complex denominator
relative to alpha_k: finite Gaussian-word approximation, native caps and
resets after compatibility conditioning, counter initialization failure,
coefficient/circumference/phase arithmetic, and sampling. Let eta_num
bound the corresponding numerator error. If eta_den<alpha_k, then

 |nu_hat/alpha_hat-nu/alpha_k|
   <= eta_num/(alpha_k-eta_den)
       + |nu| eta_den/[alpha_k(alpha_k-eta_den)].          (10)

This is a direct difference of ratios, valid for complex denominators.
For a record observable |O|<=B_O, |nu|<=B_O suffices. Requiring
eta_den<alpha_k/2 is a useful fixed margin. Small phase normalization is
not free: alpha_k=exp[-J p_k^2/(2kappa)] can make this cost large.

For independent accepted field words, empirical real/imaginary phase
means obey Hoeffding bounds since their values lie in[-1,1]. A component
error epsilon has failure probability <=2 exp(-L epsilon^2/2) for L
samples. A bounded numerator of magnitude B_O replaces the denominator2
in this exponent by2 B_O^2. Union bounds and(10) therefore give an explicit
finite sample budget. Unbounded polynomial insertions first use a frozen
cutoff and the proven conditional moment-tail bound; they do not inherit
a bounded-observable cost for free.

If a fully native binary display of a real phase component is desired,
compute a certified rational approximation q to(1+cos theta)/2, generate
l fair native current bits on the persistent auxiliary pair, and compare
the resulting l-bit integer with floor(2^l q). The mean of2*flag-1
approximates cos theta within2 times the rational/truncation error. The
sine component is analogous. All l native calls and both outcomes are
retained. This optional signed-estimator compiler is not a physical
postselection onto a complex probability.

The compatible-word acceptance cost remains of order n^(r/2) attempts
per accepted word on a fixed complex, each attempt paying its complete
native coordinate-instrument cost. Multiply that by the number of words
required for(10), plus counter initialization and optional phase-display
calls. Complete-word native TV error delta_n is amplified by compatibility
acceptance exactly as in the parent theorem, and again by alpha_k in(10).
The archive error for all attempted words is accounted separately.

## 6. What this adds and what still needs execution

This construction gives a positive limiting electric character sector
whose boundary reference is connected to the same native compatible-flux
records. Its real canonical electric mean, finite energy and normalization
are explicit, including Wick/contact terms and circle approximation.
It closes the disconnected-reference ambiguity at the theorem level.

An executed packet still needs a fixed shell and rational harmonic basis,
all native coordinate gates and caps, integer compatibility and counter
logs, finite phase arithmetic, normalization lower bounds and a complete
independent reader. A practical phase-normalization budget may be costly
and must be demonstrated. A reconstructed spectral sector is not silently
called a primitive native unitary or an already performed charged-state
intervention. Local continuum stress, curved vacuum contributions, target
inertial energy and Newton's constant remain separate joins.

## 7. Stronger result: finite compatible-word character positivity

The earlier sections deliberately used only the proved Gaussian limit.
For the counter's own block transfer a stronger finite result follows
directly from the native ternary law. It does not establish full local-time
reflection positivity of the finite Maxwell word.

Let A be any integer k-by-m compatibility matrix, with dependent rows
allowed, and Z_i independent sums of n ternary variables with probabilities
p_i<=1/2. For phi in R^m, exact torus Fourier inversion of the integer
constraint gives

 E[exp(i phi dot Z) 1_(AZ=0)]
   =(2pi)^(-k) integral_[-pi,pi]^k
        product_i f_i(phi_i+(A^T theta)_i)^n dtheta,
 f_i(s)=1-p_i+p_i cos s >=0.                              (11)

The constraint probability is the same integral at phi=0 and is strictly
positive, for example because all ternary variables can be zero. The
ratio in(11) is therefore real and nonnegative for EVERY phi. The proof
does not infer this from a covariance: it uses the exact source
characteristic functions and the full integer compatibility filter.

Take a in Z^m and the conditional counter increment a dot Z mod Q. The
finite Q-by-Q stochastic convolution P has Fourier eigenvalues

 lambda_(n,k)=E[exp(2pi i k a dot Z/Q) | AZ=0] >=0.         (12)

Thus P is exactly positive self-adjoint in the ideal finite source law.
All eigenvalues lie in[0,1], the constant eigenvalue is1, and each character
is a positive invariant sector. If all p_i<1/2, the integrand is strictly
positive everywhere, so every lambda_(n,k)>0. More generally, if the
readout uses electric coordinates only, those have p_E<1/2, and all
coordinates with p_i=1/2 have phi_i=0. The integrand is then strictly
positive near theta=0, again proving lambda_(n,k)>0. The proposed curved
Maxwell probabilities satisfy exactly this electric-only condition.

Independent accepted words concatenate to this Markov counter even if
the records inside each word have nontrivial correlations. Their fresh
source preparations or conditional reset-error budget must be explicit.
The counter's block time is the word index. It is not automatically the
local time of every plaquette inside a spacetime word. For a J-link word,
the limiting harmonic block generator is J times the one-link generator;
this factor must be retained in clock comparisons.

For a strictly positive eigenvalue define the exact finite block energy
H_n(p)=-log E[exp(-i p ell)|AZ=0], at the available counter characters.
Its interpolation as a real function of p is positive and analytic in
any neighborhood on which the characteristic stays positive. Its exact
source phase insertion is

 i E[ell exp(-i p ell)|AZ=0] / E[exp(-i p ell)|AZ=0]
     = dH_n(p)/dp.                                       (13)

At finite n, H_n is generally NOT p^2/(2kappa), and(13) must not be
replaced by p/kappa. The all-moment Gaussian limit gives H_n(p)->p^2/(2kappa)
and its first two derivatives on compact p intervals bounded away from
vanishing normalization. This follows from phase-weighted moment
convergence and the ratio bound(10). For J links use Jp^2/(2kappa).

Equation(12) supplies exact finite positive counter-sector reconstruction
from the conditional native distribution, with energy read directly from
its character response. The Gaussian step supplies its Maxwell harmonic
interpretation and canonical low-frequency normalization. Fully capped
persistent native implementations still need their source-law error bound;
a small TV perturbation does not automatically preserve positivity of all
unresolved tiny eigenvalues. It suffices for a finite chosen character
menu to bound every perturbation below its explicitly positive eigenvalue
floor. A globally positive compiled finite kernel needs a separate proof
or exact independent preparation.
