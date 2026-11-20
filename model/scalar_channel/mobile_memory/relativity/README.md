# Common dispersion and clock records on the mobile-memory carrier

The mobile-memory source supplies charge motion, cage transfer and a retained seam current. Its
kinetic terms give incompatible speeds in the shared cell-Laplacian momentum convention, and its
clock records follow the supplied source time. A separate incidence transport variant shares a
single spectral speed for charge and cage and admits a two-read clock with the corresponding
mass-derivative rate. Both calculations use the stored carrier. The reference dynamics and their
record labels are preserved. [Mobile-memory source](../README.md)

## Dispersion test of the mobile-memory source

Let $L$ be the Laplacian of the cubic cell dual with unit link weights. Use its eigenvalues
$\lambda$ as the squared graph momentum for both objects, in a common energy unit. When all
scalar links are written, the one-cage generator above the filled reference is
$H_{\rm cage}=M I+L$, with $M=2\sqrt 6$. At a sender and receiver two seams apart, $L_{ji}=0$ and
$(L^2)_{ji}=a>0$, where $a$ is the number of two-step paths. It follows that
$(H_{\rm cage}^2-M^2 I)_{ji}=a$. Every proposed $c^2 L$ has zero at that entry. This is an
integer obstruction to the squared shell for every $c$ under this convention.

Lift one cell amplitude uniformly to its twelve ports by the isometry $P$. The actual charge
operator $Q$ satisfies $P^*(Q-\mu I)P=(\eta/12)L$. Min–max bounds the ordered low charged
energies by $\mu\le E_j\le\mu+\eta\lambda_j/12$. For $0<\lambda_j\le6$ the candidate
parameters $\mu=1$ and $\eta=1/8$ give

$$
\frac{E_j^2-\mu^2}{\lambda_j}\le \frac{11}{512},\qquad
\frac{(M+\lambda_j)^2-M^2}{\lambda_j}\ge 4\sqrt{6}.
$$

The two intervals are disjoint on the 20-, 80- and 320-cell carriers. A common positive speed
would order both spectra by the same $\lambda_j$ and give the same ratio. A common rescaling
of time or graph momentum preserves their separation.

The receiver test compares the native cage motion with the spectral benchmark
$\sqrt{M^2+2ML}$, matched at the rest energy and its first slope in $\lambda$. On the native
graph a distance-two arrival starts at order $t^4$ in probability. The benchmark has a
distance-two matrix element of magnitude greater than $a/36$. The square-root
subordination formula weights the off-diagonal heat kernel by
$e^{-M^2 s}s^{-3/2}/(2\sqrt\pi)$. Cubicity and the adjacency series give a kernel lower
bound $e^{-6Ms}(2Ms)^2 a/2$. Integrating gives
$a(2M)^2/[8(M^2+6M)^{3/2}]>a/36$. After removal of the rest phase both operator norms are at most
six. At $t=1/128$ the native arrival is at most $(at^2/2+36t^3)^2$ and the benchmark arrival
is at least $(at/36-36t^3)^2$. Their difference is positive for every integer $a\ge1$.
The replay keeps every cell outcome and a parked control. The benchmark is a nonlocal spectral
comparison for this diagnostic. [Dispersion replay](unchanged_dispersion.py)

## Clock test of the mobile-memory source

The addressed cage clock starts from a sharp localized cage with a written link. Its ideal
arrival is $\sin^2(gt)$. The finite source-on fuel read has the same frequency $2g$ with an
apparatus-dependent offset and contrast. At equally spaced waits take first differences of
the full-fuel probability. They obey the recurrence
$d_{k+1}+d_{k-1}=2\cos(2g\Delta)d_k$. The clock reader uses early records to recover the
frequency and checks later waits. The blank-link and no-sender references are stationary.

With the record fixed throughout the wait, the charge source and cage clock act in disjoint
factors. The joint probability of a charge-cell record and a clock result factors. The rate
therefore stays the same in the active and parked charge arms, including reads conditioned
on detecting a charge at the receiver. Direct sparse propagation of the tensor source
independently checks the ideal joint occupations. On 20 cells the active arm raises the charge
receiver probability by 0.157928 and the two clock readings are both one radian per source tick.

An internal tag with a same-kinetic mass split $\delta$ gives the same exact obstruction for
a carried probe. The two single-charge arms differ by $\delta I$, so their relative phase is
$\delta t$ for every spatial preparation, including a received packet. The rest-normalized
rate is one. A positive square-root shell at fixed nonzero momentum requires the
rate $(2\mu+\delta)/(E_{\mu+\delta}+E_\mu)$, which is strictly below one. The tag is a diagnostic of
the additive mass architecture in the one-charge sector. The native converter changes charged
occupation in its two-excitation sector. [Clock audit](unchanged_clock.py)

## Incidence transport variant

Use a cell mode for each object and add a mode on each actual external seam. The two charged
signs use the same spatial generator. The scalar cell mode denotes the selected intact cage
of that cell, and the auxiliary scalar seam mode carries the entire excitation label.
Orient the seams, with weights $w_s$ read from held link controls, and put $\sqrt{w_s}$ and
$-\sqrt{w_s}$ into the incidence column $B_s$ at its endpoints. The local generator is

$$
D_m=\begin{pmatrix}mI&cB\\cB^*&-mI\end{pmatrix}.
$$

The supplied coefficient is $c=1$ in the common graph units. The mass settings are those
of the two rest energies above. The exact square has blocks $m^2I+c^2BB^*$ and
$m^2I+c^2B^*B$. Thus both positive branches obey $E_m^2=m^2+c^2\lambda$ for the same
incidence and momentum convention. This shell is a property of the supplied operator $D_m$. The positive-band preparation dresses cell modes into
cell and seam amplitudes. The raw operator also has negative energies, including $-m$ on
$\ker B$. The declared single-object preparation and filter admit the positive branch.
An interacting vacuum model and the opposite-charge converter across the added seam modes
are further premises for a full source.

In the addressed experiment, the same charge or cage is prepared for the written and
blocked arms. Every cell and seam outcome is retained. In a joint trial the measured
signed current through the actual charge half-link enables the cage link when its stored
value is positive. The local current is the continuity commutator with seam occupation.
The positive-band filter after that read has an explicit failure record. During the later
wait the charge and cage evolve with the held pointer. The matched shuffled arm preserves
the current and charge records and the total cage marginal. Their joint records differ.
[Incidence source and detector replay](incidence.py) · [Record-controlled cage](feedback.py)

## A clock started by its first read

Prepare a definite mass label in the positive graph mode. Its first $X$ read supplies the
recorded sign of the relative amplitude. The phase originates in that event from a state
with no relative mass-tag phase. A declared positive-mode filter for each mass keeps both
admitted branches and reports every rejection. Identical terminal position detectors read
the two quadratures of the tag. From the parity of the first and last outcomes they
recover the phase advance and its frequency. A separate parked-reference read recovers
the energies used to infer $c^2$.

The two masses are $m\pm\delta/2$, with $\delta=1/8$. At fixed graph mode the observed
clock frequency is $E_+-E_-$. Divide by the clock reading for the rest mode to obtain

$$
r_\delta=\frac{2m}{E_++E_-},\qquad
\lim_{\delta\to0}r_\delta=\frac m{E_m}.
$$

The finite difference deviates from $m/E_m$ by at most
$\delta^2/[8(m-\delta/2)^2]$ by the third-derivative remainder. An independent cell-versus-seam
occupation read gives $m/E_m$ from the difference of its two probabilities. With declared spectral
momentum $p=\sqrt\lambda$, the shell derivative is $v_g=dE_m/dp=c^2p/E_m$ and the normalized
infinitesimal rate obeys $r_0^2+(v_g/c)^2=1$. The equal-mass control carries zero clock phase
advance. The filter, read bases and waits are supplied and the first read is a marked trigger.
[Two-read clock replay](clock_read.py)

On the actual addressed carrier loop with the remaining links blocked, supply equal spacing
$a$. The momentum label $p$ is then the Fourier wave number in those spatial units, with
$E^2=m^2+(4c^2/a^2)\sin^2(ap/2)$. Prepare the circulating Fourier mode for each object.
Read the continuity current at every successive half-link and multiply the sum of its
means by the half-link length $a/2$. This directly gives
$v_g=c^2\sin(ap)/(aE)$ while retaining all current outcomes. In that same mode use the
first-read mass clock and the cell-versus-seam occupation read. The spatial phase of the
circulating preparation is given, and the initial mass tag has no relative phase.
Both resulting group speeds are bounded by $c$. The exact finite-spacing defect
$1-r_0^2-(v_g/c)^2$ is

$$
\frac{4c^2\sin^4(ap/2)}{a^2E^2}\le (ap/2)^2.
$$

For an equally spaced ring family with fixed $p$ and $a\to0$ this gives the relativistic
clock-rate identity of the supplied shell. The energy and group rate have limits
$\sqrt{m^2+c^2p^2}$ and $c^2p/\sqrt{m^2+c^2p^2}$ by the sine limits.
The stored sphere refinements and a physical distance calibration
require their own observer attachment. The bounded replays concern the declared graph or
ring momentum, stationary controls during each free interval and a first-read clock
trigger. They leave autonomous timing and the full interacting proper-time law as further
requirements.
[Circulating clock and current records](ring_clock.py)

## Reproduction

With the dependencies in `model/requirements.txt`, run from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.scalar_channel.mobile_memory.relativity
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m pytest -q model/scalar_channel/mobile_memory/relativity
```
