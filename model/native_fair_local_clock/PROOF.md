# One counted clock for native local experiments

The fair federation programme counts every fixed repair used for preparation
or readout and every uniform background repair. The construction below joins
these operations on that counter, including their complete quantum
instruments. Its limiting clock is an observation convention. It supplies no
physical energy, proper time or autonomous unitary evolution.

## 1. Literal programme and observation

Let the declared fair source have a finite graph with $M$ seams and maximum
degree $\Delta$. Each repair uses the original C or E record partition, or
their fixed mixture, fair placement coins and repair probability one. Choose
a rational partition $0=t_0<\cdots<t_b=t$ with $Mt_i$ integral. Bin $i$
contains a fixed native word $W_i$ of length $L_i$ followed by
$U_i=M(t_i-t_{i-1})-L_i\geq0$ independently uniform repairs. A fixed word
consists only of literal fixed-seam repairs. Its endpoints lie in a specified
finite root set $S$.

The programme is the concatenation of these exact commands. It uses
`export_oldest`, capacity one, and the actual source resource calculation.
Its total attempt count is

$$T=\sum_i(U_i+L_i)=Mt,\qquad U=\sum_iU_i,\quad L=\sum_iL_i. \tag{1}$$

The code rejects a word longer than its bin and authenticates every derived
field against recompilation. The native admission limits are 4096 command
nodes, $2^{32}$ total attempts, and the source record-storage and description
limits. Repeating a fixed *sequence* requires its literal command nodes;
the grammar's repetition field repeats one command. No nested word-loop
instruction is introduced.

Observe all fixed-word records, uniform records only on designated seams
whose endpoints lie in $S$, and final quantum ports in $S$. Retain bin labels,
the fixed-word ordinal, native selector tag, law, coin, raw value, decoded
current, and order within each bin. Trace the exact global tick addresses
and every unselected record value. The source retains its original full
FIFO and exported archive. This observation map does not erase them from
the programme. It also retains a genuinely untouched external reference $R$.
Remote graph ports that the programme repairs are not untouched references.

Every underlying fixed operation has all 25 original slots, including zero
hold and zero-weight law/record slots. Uniform operations retain all
$1+24M$ slots. No successful outcome is normalized or postselected. The
word query in `history.py` invokes the owning kernel and authenticates the
input operator, raw slots, records, source digest and FIFO partition. It is
a finite CP-word query, distinct from a claim that an exponentially expanded
history tree completed on the host.

## 2. Complete-instrument clock comparison

Let $B$ contain the graph radius-$r$ ball around $S$, with $K$ internal seams.
The target process runs the *same complete fixed instrument* $W_i$ at the
left boundary of each bin, followed by independent rate-one Poisson repair
clocks on its $K$ seams for the whole bin. It exports the observation just
defined. The target is a mathematical comparison process, not a replacement
native scheduler.

For any joint input on source and $R$, the half-diamond distance between the
native selected-local instrument and this window target is at most

$$a+c+d,\quad
a=\min\!\left(1,|S|{U\choose r+1}(\Delta/M)^{r+1}\right),\quad
c=\min(1,UK^2/M^2),\quad d=\min(1,KL/M). \tag{2}$$

The final sum can be capped at one. If the initial joint window/reference
states differ by half-trace distance $\epsilon$, add $\epsilon$.

**Proof.** Couple the native source to the process that discards every
outside-window repair and its record. Seed the backward cone by the entire
set $S$, including all observed and fixed-operation endpoints. A fixed
operation cannot expand this cone outside $S$. Any exit from $B$ therefore
requires at least $r+1$ chronological uniform seams. There are at most
$|S|\Delta^{r+1}$ candidate oriented edge paths and ${U\choose r+1}$ time
choices. Each independent selector supplies a prescribed seam with
probability $1/M$. The union bound gives $a$. On the complementary event,
the complete observed CP map agrees, including all off-diagonal source and
reference matrix units. Unread disjoint trace-preserving instruments cancel
under partial trace. This is the owning local-history cone argument with
the fixed root-supported instruments retained.

For every one of the $U$ background ticks, compare its window selector with
Poisson clocks over an interval of length $1/M$. With $\lambda=K/M$, the
marked distributions have exact total variation
$\lambda(1-e^{-\lambda})\leq\lambda^2$. Coupled marks use the same complete
native instrument. Multiple Poisson events and every failure/zero outcome
are retained. Tensor-reference extension and composition contract the
half-trace norm, so their full instrument distance is at most $c$.

This intermediate process has Poisson duration $U_i/M$ after $W_i$ in bin
$i$. Insert an independent Poisson interval of length $L_i/M$ immediately
after that word and before its background interval. This gives exactly the
target duration $t_i-t_{i-1}$. If all inserted intervals contain no event,
the joint maps and observed records agree. Their combined no-event
probability is $e^{-KL/M}$; hence the additional distance is at most
$1-e^{-KL/M}\leq d$. The observation discards exact times, so moving the
fixed word's native tick positions to the boundary changes no exported
label. Finally apply the triangle inequality and initial-state
contractivity. This proves (2). The sharper one-tick coupling and reservation
probabilities are not needed for the rational certificate. $\square$

Omitting $d$ would assign the background time $T/M$ to only $T-L$ background
repairs. For $U=0$, the native window has just its fixed word while the
target adds a positive-duration repair process. The reservation term is a
real finite-duration correction.

At the population-response level, the fair update is
$P=I-G/(2M)$. Consequently each uniform epoch's exchange time is
$U_i/(2M)$, half its per-seam attempt time. The existing endpoint-response
comparison applies separately on those intervals. It does not replace the
complete-history estimate (2), nor does it convert a fixed preparation word
into free evolution.

## 3. Finite rooted transport and initial states

Use an authenticated whole-cell corner ball of radius $R$ from
`native_rooted_refinement`, and let $S$ be its concentric cell ball of radius
$s\leq R$, retaining every port of those cells. The window contains the
port-graph radius $R-s$ ball around $S$: a single port edge crosses at most
one cell edge. Thus `attach_window` supplies the premise of (2) with
$r=R-s$, $|S|$ equal to the actual inner port count and $\Delta=6$.
It verifies that every fixed endpoint, observed seam endpoint and final
port lies in $S$. Observing records incident to a smaller detector region
requires including their outside endpoints in $S$.

The address map is the authenticated native map
$((f,d),p)\mapsto12(f4^\ell+d)+p$. It preserves ordered seam endpoints,
and therefore the coin and current conventions. Two matching routed
windows, fixed words, uniform laws, partitions and observation maps have one
common target process. Their output distance is bounded by the sum of
their two bounds (2), including their joint initial errors. The code rejects
a changed law, clock, observation map or fixed-word transport. Both initial
errors must explicitly name the same joint reference. This identifier
records the theorem's premise; it does not prove an unspecified input's
closeness to that reference.

For the declared fair canonical global state with $V=240\,4^\ell$ ports
and $q$ occupations, its $v$-port reduction differs from the Bernoulli
product of declared density $p$ by at most

$$\epsilon\leq\min\left(1,\frac{v(v-1)}{2V}
                           +v|q/V-p|\right). \tag{3}$$

This specialization uses a common independent external reference. An
arbitrary initially correlated reference is covered by (2) only when the
actual *joint* window/reference error is supplied. The numerical value of
$p$ is input, not a selected cosmological density.

An exact preparation instrument can be included as $W_i$, with its entire
archive on both sides of (2). A preparation theorem comparing its unread
marginal to a reference state applies only after the specified preparation
archive has also been traced. It cannot replace that archive by an
independent copy. A subsequent complete repair or readout channel contracts
such a valid marginal error. The same distinction applies to a readout
decoder's failure/erasure bound: it may be added only for its proved input
state and its declared output map, with the failure flag retained.

## 4. Admission and service boundaries

Literal root-cell fixed words are admitted at growing levels subject to
the resource limits. At levels above two, the native `_seam` validator
rejects fixed intercell addresses. A uniform graph selector remains
admitted. A multicell detector can therefore use a fixed uniform E epoch
and coverage decoding when that separate theorem's complete incident-record
hypotheses hold; it cannot silently substitute a fixed intercell tree.

The generic host executes source geometry only through level two, at most
64 expanded steps and 65,536 expanded branches. At a larger level it
returns `RESOURCE_STOP_UNEVALUATED_CARRIER` before validating the supplied
operator or constructing the graph. These terminal outcomes are preserved.
The kernel and finite rooted evaluator supply their explicitly stated
domains; equation (2) does not make a stopped host history complete.
If a finite experiment gives all unsuccessful host outcomes a separate
flag, any comparison to a success-only target must include their actual
probability, rather than dropping or renormalizing them.

For a positive fixed per-seam duration, $T=630\,4^\ell t$ eventually
exceeds the $2^{32}$ programme ceiling. The finite rooted evaluator
certifies selected windows through level fourteen; exact all-radius
oriented subdivision does not itself certify floating port routing at all
levels. A unique infinite routed process requires a specified, validated
continuation. No all-capacity executable or arbitrary-accuracy conclusion
is obtained by taking a symbolic limit of this bounded programme.

The separate proper-code source-word service accepts a supplied work
register and has its own pulse counter. The native repair command grammar
does not itself contain that operation. The complete twelve-bit phase
experiment supplies a concrete separately declared preparation/control
recipe and retained failures. These facts do not prove that coherent
services are impossible on a growing carrier. They also do not produce a
growing-carrier load/intertwining implementation or identify one pulse with
one repair attempt. Any such service join must bind the actual encoding,
initialization, full export and timing map. This theorem uses only the
admitted repair commands.

With such a separately proved attachment, native hold records can bracket
a service and supply its order relative to repair events. The operational
labels are the native attempt count and the independent service pulse
ordinal. Bracketing alone supplies no duration or conversion rate between
the counters, and does not determine concurrent uniform repairs.

The certificate compiles the owning 1931-attempt local preparation followed
by uniform epochs on the shared partition $(0,4,12)$. The final eight units
are a uniform E readout interval. The programme has 1933 nodes; source
levels zero, two and nine pass the actual admission checks, while level ten
exceeds the repeat and total-attempt ceilings. A separate short-time
comparison retains that same complete 1931-record fixed preparation on
certified 348-port windows at levels six, ten and fourteen. At level
fourteen, $t=1/630$ and declared canonical density $1/3$ with the same
independent external reference, (2) plus (3) is
less than $1/80000$. This is a bound for the stated mathematical local
instrument, not an expanded level-fourteen execution. The long readout
interval's coverage bound and the short-time spatial bound describe their
respective programmes; they are not silently combined into one error.

Exact global time records are outside the comparison: the native event
times are atomic and Poisson times are continuous. Full global record
lengths also distinguish the two schedulers. Equation (2) concerns the
specified binned selected-local records, their disturbance and $R$.
