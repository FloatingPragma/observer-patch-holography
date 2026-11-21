# A shortest preparation of the four-corner clock

A head-on encounter of two ordinary constituents can write a matched pair at
one address. Conditioned on that wire record, 26 native services suffice to
place the pair in the signed circulation around a plaquette. Twenty-five
services cannot do so with deterministic on-site controls. The count starts
at the write and includes the first compulsory newborn hop. The off-contact
launch that brings the two constituents to the vertex costs one further
service. The preparation preserves an unknown superposition of the two rails
and its correlations with an untouched reference. [Native update](../law.py) ·
[Ordinary launch and sequential 28-service preparation](clock_loader.py) ·
[26-service circuit](clock_loader_optimal.py)

The carrier is a four-axis torus with every period at least three. The target
has two occupied matched slots at each of four corners, with amplitudes
$+1/\sqrt8$ and $-1/\sqrt8$. Let $F$ exchange the ends of every routed
edge and let $\mu_r=\varphi_r/2$ on rail $r$. The target has $F=+1$ and is
orthogonal to the local equal-slot vector at every site. Each subsequent
service on the protected branch therefore multiplies it by
$-e^{i\mu_r}$, at any opening age. A local rail analyzer can read the
resulting relative phase. The signed circulation and the conserved bright
span of an unmodified birth are source properties. The operations that
transfer weight between them are supplied parts of this apparatus.
[Plaquette state and invariant phase](../pair_persistence.py) ·
[Bright birth subspace](clock_modes.py)

## Why 25 services cannot suffice

The controller may apply arbitrary position-preserving unitaries between
complete native services, with local registers and no injection of other
material. It must load the entire selected birth deterministically.
Group the material positions by graph distance from the root. Write $x,y,z$
for the state norms on the root, its neighbors and the entire remaining
carrier. The initial triple is $(1,0,0)$. The target has norms
$(1/2,1/\sqrt2,1/2)$. Local unitaries, including the contact coin, preserve
each norm. Dropping a global phase, each partial hop follows $e^{isF}$ for
an angular interval $\theta_j$. Set $T=\sum_j\theta_j$.

Almost everywhere along the interpolated hops,

$$
  \dot x^2+\dot z^2\le y^2,\qquad
  |\dot y|\le\sqrt{1-y^2}.
$$

For the first inequality, the derivative of the root norm is bounded by
the norm of the neighbor slots routed into the root. The derivative of
the last shell's norm is bounded by the neighbor slots routed into that
shell. These slot sets are disjoint because $F$ is an edge matching.
Terms internal to a shell contribute zero to its probability derivative.
For the second inequality, apply Cauchy–Schwarz to the coupling between
the neighbor shell and its complement using $\|F\|=1$. The inequalities
hold for the norm derivatives almost everywhere, including intervals with
zero population. Excursions away from the target plaquette are included.

The angle $\arcsin y$ changes at speed at most one. Its endpoints are
$0$ and $\pi/4$, so $T\ge\pi/4$. If $T<5\pi/12$, its upper envelope
$\min(s,T-s+\pi/4)$ lies below $\pi/2$. With $h=(T+\pi/4)/2$ the chord
between the endpoints of $(x,z)$ obeys

$$
  \frac1{\sqrt2}\le\int_0^T y(s)\,ds
  \le \int_0^h\sin s\,ds+
       \int_h^T\sin(T-s+\pi/4)\,ds
  =1+\frac1{\sqrt2}-2\cos h .
$$

This forces $h\ge\pi/3$, contradicting $T<5\pi/12$. For the first 25
cosine-opening services, $1-\cos v\le v^2/2$ and $\pi<22/7$ give

$$
 \frac1\pi\sum_{j=0}^{24}\theta_j
 \le\frac{\pi^2\sum_{k=1}^{25}k^2}{8\cdot128^2}
 <\frac{668525}{1605632}<\frac5{12}.
$$

## A circuit that attains 26

No gate acts before the first service. Write $c_j=\cos\theta_j$ and
$s_j=\sin\theta_j$. Apart from the known first rail phase and common hop
phase, each neighbor receives $i d$ with $d=s_0/\sqrt8$. At age one
orient the root amplitude along the six unwanted edges with real
amplitude $-m$ each, where $m=d\,c_1/s_1$. Put
$g=\sqrt{(c_0^2-6m^2)/2}$ on each of the positive $x$ and $y$ edges.
The six unwanted neighbor outputs vanish, since $i c_1d-i s_1m=0$.
The two desired neighbors each hold $i(s_1g+c_1d)$.

From age two the circuit uses only the four plaquette edges. Let $a,c$
be the root and diagonal norms and $b$ the combined norm of the two
middle sites. Orient the root equally along both edges and the diagonal
equally with the opposite sign. At each middle site put
$i b\cos\beta_j/\sqrt2$ toward the root and $i b\sin\beta_j/\sqrt2$
toward the diagonal. The exact hop gives

$$
\begin{aligned}
a'&=c_j a-s_jb\cos\beta_j, &
c'&=c_jc+s_jb\sin\beta_j,\\
(b')^2&=(s_ja+c_jb\cos\beta_j)^2+
         (c_jb\sin\beta_j-s_jc)^2 .
\end{aligned}
$$

Use $\beta_j=(1-u_j)B+u_jC$, $u_j=((j-2)/23)^2$, for
$2\le j\le25$. Solving $a_{26}=c_{26}=1/2$ then gives each corner
norm $1/2$. The terminal on-site rotations map those four rays to the
specified signed slots. Before each intermediate native coin a local
unitary maps its input ray to the chosen orientation followed by the
inverse rail coin. The subsequent native coin cancels that inverse.
The table uses 108 on-site gates, including the terminal rotations. Each preserves
the rail label. The only residual factors at age 26 are
$e^{-i\sum_j\theta_j}e^{-i\mu_r}$.

The split angles use only the native opening and the target corner
norms. For $G(B,C)=(a_{26}-1/2,c_{26}-1/2)$ the executable encloses a
zero in a box of radius $10^{-30}$ about
$(0.1878769362445742693320310536844898,\,
1.1136838729257317506159348841672850)$. Outward interval arithmetic
propagates the shell amplitudes and their two derivatives at each service.
It checks all three post-hop amplitudes stay in the positive orientation used by the
circuit. For the supplied invertible rational matrix $M$, it verifies
$\|I-MDG\|_\infty<10^{-5}$ throughout the box and that
$(B,C)-MG(B,C)$ maps the box into its interior. Contraction gives an
exact real solution. The gate table uses floating-point evaluations of
the enclosed angles. A separate two-constituent calculation checks
the gate table against the full target, including both rail phases, on
a centered and a wrapped plaquette. Its largest amplitude error at
age 26 is $1.11\times10^{-15}$ or less on both carriers. The two local
corner records are checked at ages 26 and 133. [Interval certificate
and compiled gates](clock_loader_optimal.py) ·
[Independent literal checks](test_clock_loader_optimal_literal.py)

## A uniform circuit for births at unknown locations

The 26-service minimum concerns one prescribed root. A separate 35-service
circuit uses the same eight-slot matrices at every site. It sends a bright
birth at *any* root into a coherent sum of two protected plaquette states.
The matrices are selected by the pair's local age and conserved rail.
One installed table therefore also serves births at different updates.
[Uniform circuit](clock_factory.py) · [Apparatus and observation window](clock_factory_protocol.json)

Let $\ell_x^+$ be the normalized signed circulation with root $x$. Multiply
its two middle-corner amplitudes by $-1$ to obtain $\ell_x^-$. Both have
zero equal-slot amplitude at every site. The reciprocal route gives
$F\ell_x^\pm=\pm\ell_x^\pm$. Define
$w_x=(\ell_x^++\ell_x^-)/\sqrt2$. At the root $w_x$ has amplitudes
$(+1,-1)/2$ in slots $(0,2)$, and at the diagonal it has the same two
amplitudes in slots $(1,3)$. The two slot pairs are orthogonal. Hence
$W|x\rangle=w_x$ is an isometry even when the source is coherent across
roots and entangled with another register.

**Single-parity success has a sharp bound.** On the material space write
$L_\pm|x\rangle=\ell_x^\pm$. The four axial neighbors of a given translate
overlap it by $\pm1/4$ each. All other distinct translates have zero
overlap. In spatial Fourier coordinates the Gram operators are

$$
 G_\pm(k)=(L_\pm^\dagger L_\pm)(k)
 =1\mathbin{\pm}\tfrac12(\cos k_x+\cos k_y),\qquad
 G_++G_-=2I .
$$

A coherent success record that gives the *same* ancillary output for all
origins and maps each localized input to $\sqrt p\,\ell_x^+$ has Kraus
operator $\sqrt p\,L_+$. Contractivity at $k=0$ forces $2p\le1$. Projecting
the uniform circuit's output onto $F=+1$ attains $p=1/2$ at every individual
root. The complementary parity carries the other half. For a general
coherent input $\psi$ the two probabilities are
$\langle\psi|G_\pm|\psi\rangle/2$. The uniform plane goes entirely to the
even outcome. The uniform classical mixture of orthogonally marked roots
gives one half to each. The bound concerns a common coherent material
output. Reading the origin and preparing a site-dependent state changes the
channel and can discard the original spatial coherence.

**A compact single-parity channel on four unwrapped axes.** Suppose every
translated birth maps isometrically to a translate of one finite material
packet, with the same final apparatus state. In routing parity
$\sigma=\pm1$, write $a_j(k)$ for the Fourier amplitude on positive slot
$j$ and choose $z_j=e^{-ik_j}$. Its partner amplitude is $\sigma z_j a_j$.
Local darkness gives $\sum_{j=1}^4(1+\sigma z_j)a_j=0$. At
$z_1=\cdots=z_4=-\sigma$, vary only $z_j$. The other three terms vanish,
so $a_j=0$ along that line except at its center. Continuity gives
$a_j=0$ at the center too, for every $j$. This contradicts the isometry
condition $2\sum_j|a_j(k)|^2=1$. The proof also applies on any torus
with all four periods greater than $4R+1$ for a packet supported in
$[-R,R]^4$. The parity, darkness and norm identities are Laurent
polynomials of degree span less than each period and extend from the
sampled momenta. A fixed short private period can invalidate the
extension. For private period one, the on-site antisymmetric private-slot
pair gives orthogonal translates in the dark $\sigma=-1$ band. On either
carrier the dark coin contributes $-e^{i\mu_r}$. The hop contributes
$1$ in positive parity and $e^{-2i\theta_j}$ in negative parity.
At maturity the two service eigenvalues have opposite signs and the
same inter-rail phase increment. [Finite-carrier checks](test_clock_factory_parity.py) ·
[Rhim and Yang, *Physical Review B* **99**, 045107 (2019), on singular flat bands](https://doi.org/10.1103/PhysRevB.99.045107)

**Compiling the 35 services.** A uniform even/odd slot phase plate is a
two-port phase gate on every routed edge, since $F$ exchanges its ends.
After two native hop groups of angles $a,b$, a phase between them can set
the transfer probability anywhere between $\sin^2(a-b)$ and
$\sin^2(a+b)$. This follows by adding the amplitudes of hopping in just the
first or just the second group. The next plate aligns their phases for the
third group of angle $c$. Choosing transfer $\sin^2 c$ permits complete
refocusing. Choosing $\cos^2 c$ permits complete transfer. Unitarity fixes
the orthogonal input's output ray, and a final diagonal plate corrects its
phase as well. The three compiled blocks give identity, $F$, and $F$,
respectively, after removal of the common hop phases.

| Ages at group boundaries | Lower accessible power | Chosen power | Upper accessible power |
| --- | ---: | ---: | ---: |
| 0, 3, 5, 6 | 0.00004068 | 0.00007226 | 0.00016897 |
| 6, 14, 22, 28 | 0.172628 | 0.389758 | 0.576962 |
| 28, 30, 32, 35 | 0.002751 | 0.511209 | 0.555247 |

After refocusing, a local plate divides the input between slots $+x,+z$.
The first reciprocal swap carries them to their respective neighbors.
A local slot interchange directs the first branch toward $+y$. The
second branch points back along $-z$. After the second swap, a plate
maps the two terminal slots into the root and diagonal rays of $w_x$.
The complete 35-service output on rail $r$ is
$e^{-i(\sum_{j=0}^{34}\theta_j+\mu_r)}w_x$.
There are 35 pairs of rail-indexed eight-slot matrices in the table.
It is installed at all $V$ addresses, with a local control on the
13-bit age. Each cohort traverses 35 native services and 35 on-site
layers, including its terminal plate. A wall-clock table fired only
once misses a birth arriving after those layers. Its passive bright
state has zero overlap with $w_x$. [Full two-constituent and reference checks](test_clock_factory_literal.py)

**One stationary receiver.** After loading, both parities acquire the
same contact factor $-e^{i\mu_r}$ per service. Define
$\Phi=\sum_{j=35}^{a-1}\theta_j$ at age $a\ge35$. The spatial ray is
$e^{-i\Phi}(\cos\Phi\,w_x+i\sin\Phi\,Fw_x)$. If the native wires retain
orthogonal records of the origins and $h_x$ is their birth probability
at one common age, tracing those records gives the local incidence at $z$,

$$
 \frac{\cos^2\Phi}{2}(h_z+h_{z-\hat x-\hat y})
 +\frac{\sin^2\Phi}{2}(h_{z-\hat x}+h_{z-\hat y}).
$$

Uniform unknown locations give incidence $h/V$ per original launch,
where $h=\sum_x h_x$. One detected pair then costs $V/h$ launches in
expectation. A balanced coherent rail input at a common age gives the
same hand wherever this incidence is positive. For different ages the
unnormalized local hands must be added with their actual birth weights.
Two equally likely cohorts at mature ages 150 and 4246 have identical
local incidences and opposite rail hands because
$4096(\mu_1-\mu_0)=\pi$. Discarding their age labels erases the hand.
This comparison requires a nonreturning wire through the later read.
The observations use ages below 8192, before the age-controlled sequence
can restart on a surviving pair. The measured time is the stationary
address counter's service count. [Local records and delayed-birth control](test_clock_factory_literal.py)

## How many launches distinguish two elapsed services?

On the enrolled mature cohort the signed rail hand advances by
$-\delta=-\pi/4096$ per service at a fixed station when comparing
equal parities. A separate local optical return checks the parity.
Two local four-port count words, taken on independent launches, identify
an even lapse within the supplied range $0,\ldots,8191$. The calibrated
rate and range are inputs to elapsed-time inference. The complementary
rate inference uses the optical spacings and retains the per-service
$\pi$ alias. Both readings include the original launches that missed the
station. [Optical counter and primary counts](clock_comparison.py) ·
[Count-only readers](clock_inference.py)

### A necessary terminal budget

Consider two equally likely candidate lapses separated by two services.
Let $s_j$ be the local incidence per launch at station $j$. Its incidence
and the classical probabilities of the recorded miss outcomes are the
same under both candidates. The receiver has no access to material at
another site on a miss. Even if the starting phase is known, a successful
local read sees two rail states differing by a relative phase $2\delta$.
After removing a scalar phase their unitary difference has eigenvalues
$e^{\pm i\delta}$. The root fidelity on a hit is at least
$|\operatorname{tr}(\rho U)|\ge\cos\delta$, by applying the unitary to
one side of a purification. This holds with unresolved positions,
mixed states and untouched references. A balanced pure rail state
attains equality.

Suppose $N_j$ independent terminal launches occur at station $j$.
The root fidelity including an orthogonal miss flag is at least
$1-s_j+s_j\cos\delta$ per launch. The flags give a sharper comparison.
Condition on the total number $K$ of detected pairs. The Helstrom
binary error on that pattern is at least $f(\cos^{2K}\delta)$, where
$f(x)=(1-\sqrt{1-x})/2$. The expectation of this bound is attainable
for pure balanced rails with known phases and an optimal joint receiver.
Convexity of $f$ gives the following necessary condition for total
binary error at most $\epsilon$,

$$
 \prod_j(1-s_j\sin^2\delta)^{N_j}\le 4\epsilon(1-\epsilon).
$$

At $\epsilon=10^{-6}$ even unit incidence requires 21,128,268 terminal
launches. Equal allocation across the eight registered station and
carrier settings requires at least 50,432,354 per setting, or
403,458,832 in total. Concentrating on the setting with the greatest
incidence requires at least 176,693,041. These lower bounds allow the
initial phase to be known without acquisition cost. A dephased rail can
give identical records for every lapse. A finite sufficient budget needs
a positive visibility bound. [Registered incidences and reproducible bounds](clock_budget.py)

### A sufficient budget at one mouth

In each fixture choose the second ordinary launch mouth $(1,0,0,0)$
before drawing a separate cohort. Let $f_A,f_B$ denote the incoming
fractions at the first mouth. Just the two first-service writes give
the following lower bound on incidence at the selected detector,

$$
 s \ \ge\ \tfrac14\big[f_A f_B\sin^2\Phi+
                  (1-f_A)(1-f_B)\cos^2\Phi\big].
$$

Here $\Phi$ is the accumulated partial-hop angle after loading at the
first appointment. Its squared cosine is 0.25671903. The same weights
recur at every booked even lapse. These two contributions alone give
$s>0.065$ on both fixtures. Within the 16-service enrollment all birth
phases span at most $15\delta$. Their normalized sum has magnitude at
least $\cos(15\delta/2)>0.99998$. These source bounds use neither the
observed contrasts nor a chosen subset of held count words.

Each of 54 settings receives $2^{34}$ launches. With failure probability
below $10^{-110}$ every active positive/negative axis subtotal is at
least 557,787,402. Conditional on its subtotal, the positive count on
that axis is binomial. Give the two axes simultaneous Hoeffding
intervals, allowing $2\times10^{-11}$ error per original-start input
probability. With family coverage at least $1-5\times10^{-7}$ the full
angular width for a difference of two such local reads is below
0.001527 radians. Adjacent even candidates are 0.00153398 radians
apart. On these events the elapsed decoder retains exactly the true
candidate. The sample-size certificate uses outward interval arithmetic.
[Budget certificate and decoder](clock_budget.py) · [Interval tests](test_clock_budget.py)

| Independent original launches | Eight-station comparison | One mouth on each fixture |
| --- | ---: | ---: |
| First active count word across booked settings | 549,755,813,888 | 34,359,738,368 |
| Each of eight later active words | 549,755,813,888 | 34,359,738,368 |
| All nine active words and both control menus | 14,843,406,974,976 | 927,712,935,936 |

The two individual receivers each decoded all eight hidden lapses at
the reduced budget, including 2 and 7,926 services. The passive control
gave a unique correct interval on zero of its 16 comparisons. All 16
dephased comparisons retained the full 4,096 even candidates. The
second acquisition uses the authenticated source probabilities and an
independent multinomial seed. Every seven-outcome histogram includes
the miss counts. [Reduced count archive](clock_budget_counts.npz) ·
[Archive checksum](clock_budget_counts.sha256)

## Transporting the loaded clock

The preparation uses two independent ordinary feeds focused on the selected
mouth with probabilities \(3/4\) and \(1/2\). The first marker accepts
\(3/16\) of ordinary starts there. The neighbouring mouth takes \(1/16\)
and the unmarked outlet takes \(3/4\). The 28-service loader prepares the
signed loop on the selected branch. It rests until age 128. The transport
plates on both rails are identical. They act before or after a complete
native service at the occupied corners. Laboratory counters enroll each
appointment, with the same phase convention at each station. [Source and
gates](clock_transport_protocol.json) · [Routing and readers](clock_transport.py)

### How much can move in one service?

Let \(s\) be the required outgoing slot, \(D_r\) the contact coin on rail
\(r\), \(w_r=D_r^\dagger|s\rangle\) and
\(\delta=(\varphi_1-\varphi_0)/2\). For a translation by one positive
axial link, each leading corner of the target has exactly one predecessor
in the original loop. That corner must supply all its weight through slot
\(s\) for a unit-probability translation. A rail-blind local preparation,
including one entangled with a local auxiliary system, supplies a reduced
slot density \(\rho\). Its two routed probabilities obey

\[
 \min_r \langle w_r|\rho|w_r\rangle
 \le \left\|\frac{|w_0\rangle\langle w_0|+|w_1\rangle\langle w_1|}{2}\right\|
 =p,\qquad
 p=\frac{1+|d|}{2},\qquad
 d=\frac{e^{i\delta}+7e^{-i\delta}}8 .
\]

Here \(d=\langle w_0,w_1\rangle\) follows from the equal-slot weight \(1/8\)
and its complement \(7/8\). Thus \(p<1\). The normalized vector
\((w_0+e^{-i\arg d}w_1)/\sqrt{2+2|d|}\) attains \(p\) on both rails.
At each corner a common on-site plate maps its signed loop ray to this
vector. The native update delivers the selected slot into its partner
at the translated address. A common output plate maps the partner to
the signed loop ray there. The successful material operation has rail
amplitudes \(\sqrt p\) and \(\sqrt p\,e^{-i\alpha}\), where
\(\alpha=\arg d=-\arctan[(3/4)\tan\delta]\). This action also gives the
success probability and relative phase for a rail qubit entangled with
an untouched reference. [One-hop basis and flag tests](test_clock_transport.py)

Local absorbers record each rejected address and slot at its first failed
hop. The four transmitted corners leave identical null flags, which
preserves their relative phases. The plates run on a fixed schedule.
Comparing the stored flags after the journey gives total survival
\(p^J\) after \(J\) successful translations. On original starts, the
first-failure probability at hop \(j\) is \((3/16)p^{j-1}(1-p)\).
The transport probability per step is \(0.9999999356575108\). A
transport step moves the hand by \(-0.000575242845\) radians. A parked
service moves it by \(-0.000766990394\) radians. With one translation
every \(s\) address services the signed rate per service is
\(-\delta+(\alpha+\delta)/s\).

### Local reception and reunion

The booked route crosses a periodic seam, moves two links along \(x\)
and one along \(y\), then returns on the other sides of the rectangle.
A second route traverses the rectangle in reverse. Fresh preparations
supply the four rail counters and every complementary outcome at each
read. At spacing four, after three hops the displaced station receives
\(0.046874990952\) per original launch. Its corner clicks disappear when
the loop stays home. The home station receives zero at this appointment
on the transported branch. The displaced station reads its own count
141. Three spatial links separate it from home, so the first causally
possible home comparison is count 144. The experiment compares station
logs after acquisition. It measures the transported material at actual
booked stations and supplies their clocks and phase standards.

| Services per link | Forward reunion difference (rad) | Reverse reunion difference (rad) |
| ---: | ---: | ---: |
| 1 | 0.00114983 | 0.00115219 |
| 4 | 0.00114686 | 0.00114758 |
| 8 | 0.00115238 | 0.00115075 |

Each entry subtracts the separately parked reading at the same elapsed
service count. All six differences agree with the finite forecast
\(6(\alpha+\delta)=0.001150485295\) radians within their simultaneous
half-widths of at most \(0.000188984\) radians. Each terminal word uses
\(2^{44}\) virtual multinomial exposures and keeps the four rail ports,
the missed corner, both natal refusals and every first-failure bin.
The familywise error bound is \(0.001\) across 1,167,156 bin comparisons.
Suppressing the absorbers allows rejected amplitudes to interfere again.
On the forward spacing-four return, \(5.88\times10^{-8}\) per original
start lies outside the four target corner rays, and the hand differs
from the screened reading by \(1.17\times10^{-5}\) radians. The service
rate depends on the driven transport and its local filters. A free
material proper-time rate requires an independent construction.
[Complete likelihoods](clock_transport_likelihoods.npz) · [Original
counts](clock_transport_counts.npz) · [All arrivals and comparison
intervals](clock_transport_records.json.gz)

A separate moving ordinary clock selects two positive contact lines of
an exchange-odd ordinary pair at matched modal current. Its reported
moving and stationary local rates are 0.173920838552 and 0.176797604897
radians per local tick. The driven loop here occupies the marked,
matched sector and uses a directed slot filter at each move.

## A freely dispersing hand on a moving receiver path

A sharp ordinary head-on approach reaches the origin in one service. Its
first capture writes the local marker with probability \(1/2\), leaving the
equal-slot ray \(S\) with equal coherent rail amplitudes. After capture the
pair takes every opening and mature service freely. On the axial strip
\((1,1,16385,1)\) the enrolled stations are centered at \(z=2n/7\), where
\(n=3136,3920,\ldots,7840\) counts mature services after the 128 opening
services. Each station supplies a coherent 25-address aperture and the same
four rail comparisons. These are destructive readings on fresh launches.
Reconciling them requires the booked locations, their relative service
labels and a common analyzer phase convention. No spatial or age winding
enters any record.

For \(\mu_r=\varphi_r/2\) the bright dispersion is
\(\cos E_r(k)=(3+\cos k)\cos\mu_r/4\). At \(v=2/7\), \(E_r'(k)=v\) has
two positive roots \(k_{r,L},k_{r,H}\). Each gives a right-moving component
at \((k,+E)\) and a second at \((-k,-E)\). The two rails have slightly
different roots. Applying rail zero's low momentum to rail one would give
velocity 0.284210, a shear of 0.00150445 sites per service.

Write \(u(k)_j=e^{ijk}\) for \(j=-12,\ldots,12\). The low receiver projects
\(u((k_{0,L}+k_{1,L})/2)\) orthogonally off the six unwanted steering
vectors from the two rails, then normalizes the resulting ket \(q\).
The high receiver interchanges \(L,H\). The common response \(q^\dagger
u(k)\) vanishes at each rejected stationary point. In the unwrapped mature
propagator, stationary phase gives the selected rail amplitude

\[
h_r(n)=\frac{q^\dagger u(k_r)\,\langle S|P_{+,r}(k_r)\psi_{128,r}(k_r)\rangle}
 {\sqrt{2\pi n |E_r''(k_r)|}}\,
 e^{-in(E_r(k_r)-vk_r)-i\pi\,\mathrm{sgn}(E_r''(k_r))/4}
 +O(n^{-3/2}).
\]

All stationary points are nondegenerate. The exact finite opening
\(\psi_{128,r}\) is a trigonometric polynomial before projection and the
gapped spectral projectors are smooth. The rejected terms therefore start
at order \(n^{-3/2}\). At finite time the aperture's trigonometric
polynomial also admits momenta near each of its isolated zeros.
The admitted intensity decays as \(1/n\).
The signed hand \(\arg(h_0\overline{h_1})\) has asymptotic rate
\((E_1-vk_1)-(E_0-vk_0)\) per booked service.

A single-address slot filter cannot make the same selection. With fixed
independent slot vectors \(S,T,O\), the routing satisfies
\(F(k)S=S+(\cos k-1)T+i\sin k\,O\). At a stationary point of sign \(s\),
the eigenray is proportional to \(e^{-i\mu}F(k)S-e^{isE}S\), and
\(s\sin E=\cos\mu\sin k/(4v)\). Its coefficients are an invertible
linear image of \((1,\cos k-1,\sin k)\), with diagonal entries
\(-i\sin\mu,e^{-i\mu},ie^{-i\mu}\). Any three distinct points on the
unit circle are noncollinear. The three unwanted, populated eigenrays
thus span the entire accessible one-address slot space. A positive local
effect which annihilates all three also annihilates the wanted ray.

| Common receiving port | Asymptotic advance in 4,704 services | Exact finite advance | Simultaneous count interval | Admission at last visit per start |
| --- | ---: | ---: | ---: | ---: |
| Low positive branch | 2.951209 | 2.948391 | [2.940842, 2.958570] | 0.00012949 |
| High positive branch | 0.291419 | 0.291437 | [0.287350, 0.295535] | 0.00061687 |
| Single address, equal slots | multiple components | 0.736800 | [0.713360, 0.757805] | 0.00001652 |
| Low aperture, rail phase washed | undefined | undefined | full circle | 0.00012949 |

The advances are radians modulo \(2\pi\). Equal cadence admits per-service
frequency aliases spaced by \(2\pi/784\). The two positive rates above
use the native branch assignment. Every signed low and high advance passed
the 0.05-radian gate, including count uncertainty at family coverage
\(1-10^{-6}\). On the low branch the leading complex amplitude had at most
3.39% relative error. **The high branch failed the 10% complex-amplitude
gate at all seven appointments on both rails.** Its largest error was
32.93%, and the last exceeded 13.59%. The high branch's relative rail
phase is considerably more accurate than either individual amplitude.

The 28 simulated words contain \(481{,}036{,}337{,}152\) original starts,
including the four rail outcomes, marked misses and capture refusals.
The receiving apparatus moves along a specified path through a spreading
field. A bound, localized specimen carrying an observer or a reciprocal
timing standard requires additional construction. [Full receiver and
source](clock_ballistic.py) · [Original counts](clock_ballistic_counts.npz) ·
[Independent literal opening and readout](test_clock_ballistic_literal.py)

A separate matched-velocity compact experiment supplies narrow incident
bright bands centered at the two momenta. Its traveling ordinary-slot
experiment reads opposite energy signs with a real aperture. Here one
ordinarily formed pair supplies the occupied bright continuum and the
complex receiver selects the same energy sign on both rails.

### How long can a commissioned packet survive?

The native age has \(L=8192\) states. A pair with birth service \(b\) has
age \((n-b)\bmod L\) at service \(n\). Write \(F\) for reflected routing and
\(P_\pm=(I\pm F)/2\). Both the loaded wave and its routed image are
orthogonal to the equal-slot inlet at every site. The contact coin thus
acts on rail \(r\) by \(-e^{i\mu_r}\). The partial hop acts by 1 on
\(P_+\) and by \(e^{-2i\theta_j}\) on \(P_-\). These commuting actions
preserve darkness even when the opening starts over. The first 128 angles
sum to \(129\pi/4\), since \(\sum_{k=1}^{128}\cos(k\pi/128)=-1\).
The remaining 8,064 full hops have even parity. The rail phases complete
128 and 129 turns, respectively. One cycle therefore acts on either rail
by \(P_+-iP_-\); four cycles act by the identity, including on a rail
entangled with an untouched reference. An opening held permanently at
its mature value gives fidelity one half on the equal-parity loaded state
after the first cycle.

A sharp ordinary head-on launch supplies a first write with probability
one half. After loading, a corner receiver at age 128 has conditional
admission 0.3716405. At the same age on successive cycles its admissions
are 0.0315885, 0.1283595, 0.4684115 and 0.3716405. All four rail ports
and the complement are evaluated at each independent terminal appointment.
Per original launch each conditional weight is multiplied by one half
and the refused first write carries the remaining half. These are exact-law
probabilities evaluated numerically. An independent real-space calculation
follows the original source through all 32,896 services and checks both rail
columns against the parity formula.

The simulator specification enrolls a natal wire longer than its run. To
exercise reverse capture we also run the specified vertex on shorter natal
wires with the periods given below. Reverse capture requires the occupied
wire head and age zero. For a
marker on its natal periodic wire of width \(W\), the first joint return
occurs at \(\operatorname{lcm}(L,W)\). Widths 64 and 4,096 admit every age
wrap; width 65 first admits service 532,480. At every eligible site
the loaded state has zero equal-slot amplitude. The same vertex releases
one half of a bright equal-slot input. This argument covers a routed
marker at any site where its head is actually occupied. A general
bright pair requires knowing that site. The material revival makes no
claim about revival of the whole marker or controller.
[Age, routing and capture](../SPEC.md)

An age-only plate table starts again on each wrap and produces 0.1266511
capture-active weight after its first restart. A finite controller can
schedule the first loading and four passive cycles. Give each incoming
enrolled pair a blank five-state program digit, carry it with the matched
slots, increment it modulo five at the native age wrap, and enable the
local plate only at digit zero. This reversible address update commutes
with translations and preserves spatial and rail coherence. The digit
is extra apparatus. Constituent A carries it in the ordinary sector.
The eligible vertex rotates the ordinary amplitude against the matched
equal-slot inlet independently in every digit, rail and reference state.
Its inverse recovers an arbitrary coherent input on all these sectors,
including after bright discharge. The other wire configurations retain
the native spectator action. On the dark
branch all four passive wraps retain zero bright weight. On the fifth,
the table reenables and the bright weight returns to 0.1266511.

Five states suffice and are necessary for an isolated reversible program
clock that supplies those exact enable words. At a common native age let
\(V\) advance the controller one lap and \(S\) be its prepared support.
Perfectly readable enabled and disabled words give \(S\perp V^dS\) for
\(d=1,\ldots,4\). Unitarity makes all five translates pairwise orthogonal.
The same argument for unbounded \(d\) rules out permanent one-shot
disarming in such a finite closed program clock. An external synchronization
record, a changing material input, or a reset carries additional resources.
The terminal records here resolve a finite recurrence; they supply no
absolute time standard. [Cycle and program](clock_cycle.py) ·
[Literal field propagation and inverse program](test_clock_cycle_literal.py)

### Reading the lap and the age at one corner

Consider the known natal root of the balanced, sharply formed two-rail
packet. Let the native age be even, between 128 and 8,190 inclusive. At
these ages each mature two-service hop returns the spatial ray. At the
root, \(w_x\) has squared norm one half and \(Fw_x\) vanishes.
Writing \(\Phi=\sum_{j=35}^{127}\theta_j\), the admission per original
launch on lap \(h\) is

\[
 s_h=\frac{1+\Re[(-i)^h e^{-2i\Phi}]}8 .
\]

The four values in order are 0.1858202, 0.0157942, 0.0641798 and
0.2342058. The quadrature phase advances by \(-\pi/4096\) per native
service relative to the age-128 reference, independent of the lap.
Consequently the six-outcome law distinguishes all 16,128 enrolled
epochs modulo 32,768. The denominator includes every launch. Conditioning
on a click erases the admission, as does an unobserved uniform mixture
of natal sites for this lap comparison. This reader uses one known
birthplace, sharp first-write incidence and a coherent terminal rail
reference.

The instrument enrolls four age-128 count words, one per lap, and accepts
only counts as terminal input. Across three menus, 36 words of \(2^{38}\)
simulated launches cost 9,895,604,649,984 original starts. All eight
blinded acquisitions on the native opening gave a single correct address
class, including an independent reading four full laps after the reference.
The saturated-opening control gave four candidates per reading and the
dephased control 4,032. Each control has its own four calibration words.

For any enrollment of this size drawn from the declared domain, simultaneous
Hoeffding intervals cover the incidence and rail phase at failure probability
at most \(10^{-6}\). The four incidences differ by more than 0.048 and the
smallest exceeds 0.0157. Except on an axis-count outage of probability
below \(10^{-100}\), each native phase axis has at least two billion counts.
At that floor, each possible phase differs from its true value by at most
\(\arcsin(2\sqrt2 q)\), where \(q\) bounds the single-coordinate error of
the estimated quadrature including the dyadic source allowance. The sum of
the two phase errors is below the even-age spacing \(\pi/2048\).
The interval calculation checks both separations with outward rounding.
The uniform budget covers any booked set of these visits. The acquisition
samples eight terminal words per menu.
[Executable source and budget](clock_epoch.py) ·
[Archived counts](clock_epoch_counts.npz) · [Independent native formula](test_clock_epoch.py)
