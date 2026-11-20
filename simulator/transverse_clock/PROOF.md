# When does the scalar cone limit a message?

## Free propagation

Bernhard Mueller's eight-slot update is $U=SC$. The local coin is $C=2P-I$ with $P=|s\rangle\langle s|$ and
$s=(1,\ldots,1)/\sqrt8$. In a Fourier fiber the flip-flop shift $S(k)$ moves each outgoing amplitude by one
signed step and exchanges the two slots on its axis. Thus $S^2=I$. Three coordinates serve as spatial dials
and the fourth is private. The proof uses lattice units on the infinite lift. A periodic measurement uses it
only before a periodic image can reach the detector. [Free update](../law.py)

For $v=Ss$ and $a=\langle s,Ss\rangle=\tfrac14\sum_j\cos k_j$, one has $Us=v$ and $Uv=2av-s$.
The dispersive frequencies obey $\cos E=a$. Their group velocity has components $\sin k_j/(4\sin E)$.
Since $\sum_j\cos^2 k_j\ge4a^2$, its norm is at most $1/2$, approached by long waves. The orthogonal
complement of $s,Ss$ evolves by $-S$ with frequencies zero or $\pi$. Limits cover coincident bands.
A bound on group velocity alone does not control finite-time communication.

## All-slot exponential bound

Let a source occupy $\hat n\cdot x\le0$ and a detector occupy $\hat n\cdot x\ge d$, for a unit vector $\hat n$.
Weight each position amplitude by $\exp(\theta\cdot x)$, where $\theta=\eta\hat n$ and $\eta>0$.
The shifted Fourier operator contains $S(k+i\theta)$, whose norm is $D=\exp(\max_j|\theta_j|)$ and whose square is the identity. Define

$$
 A(\theta)=\operatorname{arcosh}\left(\tfrac14\sum_j\cosh\theta_j\right).
$$

The complex $a$ lies in the ellipse with semiaxes $\cosh A,\sinh A$. In the real planar direction of angle
$\phi$ the support function of the $j$-th cosine ellipse is $\sqrt{\cosh^2\theta_j-\sin^2\phi}$.
The function $x\mapsto\sqrt{x^2-\sin^2\phi}$ is concave on $x\ge1$. Averaging proves the inclusion.
The two roots of $z+z^{-1}=2a$ therefore have modulus between $e^{-A}$ and $e^A$.

The row space spanned by $\langle s|$ and $\langle s|S$ is invariant. With second-kind Chebyshev polynomials
$V_j$, the row $g_j=\langle s|U^j$ equals $V_j(a)\langle s|-V_{j-1}(a)\langle s|S$, where $V_{-1}=0$.
The Laurent sum gives $|V_j(a)|\le(j+1)e^{jA}$, including repeated roots. Expand $U^t$ about $(-S)^t$
using $U+S=2SP$. Its $j$-th additional term is $2(-S)^{t-1-j}S|s\rangle g_j$. The factor left of
the ket has norm at most $2D$. Plancherel and removal of the position weight give

$$
 \|\Pi_{\rm det} U^t\Pi_{\rm src}\|
 \le D[1+(1+D)t(t+1)]\exp(t A(\theta)-\eta d).
$$

This covers every initial slot polarization, both flat bands and an untouched quantum reference.
Since $A(\eta\hat n)=\eta/2+O(\eta^3)$, fixed $d/t>1/2$ gives exponential suppression. For $B$ noninteracting
signal particles initially confined to the source half-space, the expected receiver occupation is at most $B$
times the square of this bound, by applying the estimate to each one-particle marginal. The probability of
detecting any particle is no larger. A polynomial number of independent occupied-rail trials cannot overcome
this suppression. Pre-existing remote carriers, intermediate writes and active receivers that transmit during
the flight require their own analysis.

## Sharp cost of a faster passive signal

For an equal-slot state localized at the origin, the squared norm of the weighted wave is
$M_t(\theta)=\mathbb E[\exp(2\theta\cdot X_t)]$ under its complete position measurement.
The norm bound before removal of the weight gives $\limsup t^{-1}\log M_t\le2A(\theta)$.
At $k=0$ in the shifted Fourier integral the equal-slot return is $T_t(\cosh A)=\cosh(tA)$.
For $\theta\ne0$, in an open neighborhood the larger root has modulus exceeding $\exp(A-\epsilon)>1$.
Its reciprocal cannot cancel it in $T_t=(z^t+z^{-t})/2$ at large $t$. Integrating the squared return over this
fixed positive volume proves the reverse limit. This argument also works for the equal-slot detector alone
and for the axial, transverse-uniform preparation, using its one-dimensional Fourier integral.

On an axis $A(\eta)=2\operatorname{asinh}[\sinh(\eta/2)/2]$. For $1/2<u<1$ take the unique $\eta>0$
with $A'(\eta)=u$. Under the measure tilted by $\exp(2\eta X_t)/M_t(\eta)$, $X_t/t$ concentrates on $u$.
Apply exponential Markov bounds on either side using $M_t(\eta\pm h)/M_t(\eta)$.
Differentiability makes both bounds decay exponentially for sufficiently small $h>0$. Removing the tilt on
$|X_t/t-u|<\delta$ gives the lower bound for that window. For the half-line take a window just above $u$,
then shrink it towards $u$. The upper bound is the original exponential Markov inequality. The same proof
normalizes the tilted equal-slot subprobability by its own moment. Consequently

$$
 -\lim_{t\to\infty}t^{-1}\log\Pr[X_t\ge ut]
 =2(\eta u-A)=4[u\operatorname{asinh} b-\operatorname{asinh}(b/2)],
 \qquad b=\sqrt{(4u^2-1)/(1-u^2)} .
$$

At $u=1$ only the straight path arrives. The equal-slot launch has probability $(1/8)16^{-(t-1)}$ there
over all slots. A directed launch $C|{+x}\rangle$ has probability $16^{-(t-1)}$. Both give rate $\log16$.
Beyond that front the probability is zero. Within the cone the best uniform exponential bound is trivial.
For the transverse-uniform equal-slot experiment put $P_t=2\,4^t T_t((6+z+z^{-1})/8)$.
Its integer Laurent coefficients follow $P_{t+1}=(6+z+z^{-1})P_t-16P_{t-1}$ from $P_0=2$ and
$P_1=6+z+z^{-1}$. Squared coefficients divided by $4\,16^t$ give the equal-slot position law on every site.
Other slot outcomes carry its missing normalization. [Exact finite tail numerators](records/result.json)

## A message with local regeneration

On an axial segment of length $r$, a directed launch reaches the first allowed arrival with probability
$p=16^{-(r-1)}$. Each forward step after the first uses the off-diagonal coin entry $1/4$. On an oblique
free route all shortest paths ending in the same slot add coherently. The forecast counts permutations by
first and last step. Acquisition scatters the actual eight amplitudes at each reachable site. Pruning outside
the monotone box is exact at the shortest arrival. Each receiver reads both rails of each bit.

With $B$ copies per bit the segment delivers a bit with probability $1-(1-p)^B$. For a $b$-bit frame
and $m$ equal segments each way, the ideal authenticated return probability is $[1-(1-p)^B]^{2mb}$.
With one service tick per regeneration or remote reply the round trip takes $2mr+2m-1$ ticks. Initial
preparation and terminal read are boundary services. A sufficient budget for failure at most $\epsilon$
is $B\ge\log(2mb/\epsilon)/[-\log(1-p)]$ for $p<1$. For fixed $r>1$ the required parallel copies grow
logarithmically with distance. The asymptotic axial speed is $r/(r+1)>1/2$. Taking $r$ of order
$\log_{16}d$ and using order $d\log d$ copies per bit per segment approaches speed one with polynomial
total resources. This statement charges the declared parallel service. The circuit depth of reading and
regenerating an arbitrarily large bank in serial hardware is outside it.

The stations have reserved, initially empty lanes for each trial, hop, bit and direction. A source uses
fresh preparations for each occupied rail. A relay prepares its next frame only from a complete local read.
Other branches and failed preparations need no coherent reset and are included in the counts. The passive
comparison uses the same parallel budget per launch and no relay. The addressed comparison carries the
payload explicitly by the supplied local courier. The wave runs use the courier only to calibrate the
chart and deliver terminal reports. The message between wave stations travels through detections.

The sender chooses a fresh nonce after the forecasts are sealed. Distinct request and acknowledgment tags
bind the paired endpoint key, trial identifier and payload. The home station checks the returned nonce.
Keys and sender entropy are disclosed for offline replay. This is a trusted laboratory transaction with
explicit setup. The scalar deadline uses a flat Euclidean interpretation of the three reciprocal chart
ranges and the independently measured wave speed.

## Two pulse clocks and their finite resources

The original clock combines adjacent private Fourier orders. Here each prepared mode has private order zero.
An independently enrolled pulse circulates in a *public transverse direction*. At real momentum let
$w=(Ss-as)/\sin E$. In the basis $s,w$ the update is a planar rotation. The preparation
$(s+iw)/\sqrt2$ has eigenvalue $e^{-iE}$. The two-dimensional packet superposes five Fourier orders on
each axis with the registered binomial weights. It is uniform in the other two coordinates. After one
free update per service, local counters measure all position bins, summing over slots. Combining
those reports gives the first circular position moment. Its phase measures the pulse's displacement.

Separate preparations supply the reference speed $c$, longitudinal speed $v_x$, and transverse rest
lap rate. They predict the moving-to-rest ratio $\sqrt{1-(v_x/c)^2}$. For an infinitesimally narrow,
long-wave pulse this follows from the magnitude $1/2$ of the continuum group velocity. A finite-width
grating need not attain that limit by increasing the grid size alone. The complete finite-gate forecast
uses the actual prepared modes and dispersion. Reported passes use the registered tolerances.

For $M$ equally spaced counter thresholds the cumulative probability error at any cut is at most $1/(2M)$.
Summation by parts around the circle bounds the error in its complex moment by $\pi/M$. The angular
uncertainty at a nonzero measured moment $z$ is at most $\arcsin[\pi/(M|z|)]$. The moment couples adjacent
orders, whose frequencies differ by at most $\pi/L$ by the group-velocity bound. The sum of their absolute
overlap weights is at most one by Cauchy–Schwarz. Between enrolled samples spaced by $L/8$, its path stays
within $\pi/8$ of its starting value. At the registered visibility threshold 0.45 this excludes an
unobserved full turn. A uniform transverse translation erases the clock hand while preserving its longitudinal marginal.

The local challenge adds three reversible tag bits. On crossing the transverse seam at the edge of the
launch-centered chart, the tag changes by the signed crossing. Inverse Fourier transformation of the
prepared wave, injection into tag zero and a Fourier transform on the eightfold lift preserve its norm.
At $y=0$ the eight lifted addresses separate the tag outcomes. An aperture fixes $x,y,z$ and the private
address and accepts every slot. Both uniform spectator coordinates contribute the probability factor
$1/L^2$. The experiment uses $2^{16}L^4$ midpoint preparations per aperture and wait. It includes every
miss. For waits below $6L$ an initially zero-tag amplitude ending at $y=0$ and tag one has precisely one
positive net cut crossing. A path with winding $-7$ or $9$ would exceed the nearest-neighbor support.
The initial preparation has finite width. A tagged click identifies a winding of that amplitude and
a sampled return peak, without defining a point worldline.

The sealed local forecast propagates the two dispersive components by a rotation and the orthogonal
complement by $(-S)^t$, including both flat bands excited by writing the tag. Acquisition independently
iterates the full eight-slot gate. Suppressing the write leaves the different windings coherent in
the original torus. Summing measured tags after the write instead sums their probabilities.
