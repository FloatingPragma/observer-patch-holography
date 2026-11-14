# Common completed-word clock and locality

The source is the frozen proposed local twelve-slot walk W=SR. The geometry,
native event dilation, retained record recoder, preparation and all padding
are those of `auditor/native_walk_bridge`. The clock attachment has its own
identifier, `local12-completed-word-clock-v1-proposed`. The attachment adds
an operational clock to this source; it does not change the source step or
identify it with the C10 Hamiltonian.

## Exact classical comparison

Let the connected six-regular graph have M seams. In the one-load classical
domain, an independent uniformly selected seam touches the load with
probability p=6/M, independent of the occupied vertex. Conditional on contact,
the six seams are uniform and the placement coin is fair. A contact can hold
the load; a moving-current counter is therefore a different clock.

Nonincident choices leave the load fixed. They retain the seam, coin and native
C or E label. For a fixed law, every specified n-letter global path has
probability (2M)^(-n). The first-contact attempt has probability
p(1-p)^(n-1). A cap K has the separate timeout probability (1-p)^K. These
probabilities sum to one before any conditioning. Given the current load,
sampling the discarded-length variable and every nonincident letter rebuilds
the same complete classical native path. The executable does not discard them.
For the operational contact-count replay, the initial unit position is known.
The observer reconstructs each subsequent position from the signed-current
tape; simulator endpoint occupations are audit-only. A zero C current by
itself cannot classify contact when the initial position is unknown.

The transition matrices satisfy

    T_global = (1-6/M)I + (6/M)T_local,
    T_local  = I-L/12.

Assigning every seam the same Poisson attempt rate r gives total rate Mr and
local contact rate 6r. Both generators equal -r L/2. Thus a per-seam rate is
volume consistent; assigning a fixed total global rate is a different
scheduler. The checks use existing development sizes two, three and four,
including 55,296 actual native source columns at size two.

This is a classical diagonal-state and marked-path identity. It is not an
identity between the full recorded quantum channels, and it does not turn
successive coherent W steps into repeated fresh native coins. Starting from
one vertex with a uniform event state, the two-step return probabilities are
25/72 for W and 7/24 for T_local squared. Their difference is exactly 1/18.

## One operational phase coordinate

The clock uses a fresh proper-code qubit prepared in the plus state. After
each completed source word, one emitted native Z pulse applies exp(-i h Z/16)
to that clock. The base member has h=1. The ideal joint step is
W tensor exp(-i h Z/16). The clock coherence is exp(i n h/8) after n completed
steps. Both the preparation and the
phase increment are explicit architecture choices. A different phase
increment defines a different clock attachment.
The source convention is exp(i theta Z). Here theta=-h/16, while the measured
Bloch coherence 2 rho_10 has relative phase -2 theta=h/8. The pulse angle and
readout phase differ by this factor of two. Under t=n tau0 h, the measured
angular clock frequency is 1/(8 tau0).

Configured X and Y analysis words, followed by the actual native C or E
pointer, produce contrasts cos(n h/8) and sin(n h/8). The executable retains every
native fine record, both binary outcomes, source-export complements and reset
timeouts. It executes base steps zero through eight and refined members with
h=1/16,1/32,1/64, reaching 512 clock ticks. The domain 0<=n h<=8 has phase in
[0,1] radians and a unique monotone lift. Scaling the clock increment with h
allows a fixed effective-time window while the source step count grows.
The reported sign observable
can assign zero to a failure flag. It does not divide its expectation by the
success probability.

The source control certificate supplies a complete W Stinespring/operator
error at most 6.559285619e-11 and a source-preparation error at most
1.111382776e-14. The clock's own words and native pointer have outward rational
error bounds from the existing source-control certificate. The joint bound
adds all twenty-three common blank resets, one clock blank reset and one
ready-pointer reset. Each reset is bounded on all seventy-eight native input
columns. With at most 512 W steps, the complete joint error is below
4.410372e-7 in the stated total-variation/half-diamond comparison. No source
or apparatus failure is removed.
This joint numerical bound uses the frozen size-two W compiler. Changing h
in the one-qubit clock does not compile a larger spatial graph. The clock
family can attach to a larger ideal W_L, but its implemented joint bound
requires that member's own source word and accumulated error certificate.

Each unconditional sign expectation changes by at most twice this error.
A disk of radius three times the error encloses the complex two-quadrature
contrast. Its phase error is bounded by radius/(1-radius), below 1.323114e-6
radians. Adjacent clock phases remain separated. The certification assumes
the declared ideal coherent load/export services and exact retained source
response; it supplies no uncontrolled laboratory-noise bound.

The default clock wire is 23, outside the common work wires zero through 22.
A joint experiment with another apparatus using wire 23 must remap one of
them to a distinct address. The full retained-current construction uses 898
wires, so its clock and probe require fresh addresses such as 898 and 899,
with further detector auxiliaries separate. The emitted clock-tick API checks
the complete declared source register size. Its bare-W error and duration
cannot be reused for W_F without adding the full current words, preparation
and reset costs. A physical controlled-W_F probe must also control those
current additions; controlling bare W alone can leave the Gauss code.
Independent clock preparation is stated; a correlated input requires its full
joint state rather than this product calibration.

## Refinement time and implementation duration

For cell spacing h and one positive fixed scale tau0, define

    t(n,h) = n tau0 h.

`effective_time` implements this map with exact rational arguments. Every
particle, field and apparatus observable in one declared refinement uses the
same h and tau0. The source dispersion cos(theta)=1-lambda/12 and
lambda(k)=|k|^2/8+O(|k|^4) then give acoustic speed squared 1/(48 tau0^2).
For tau0=1, the speed is 1/sqrt(48). The time coordinate can be read from the
native clock phase as 8 tau0 times its unwrapped phase, within the declared
finite calibration window.

This coordinate is not derived physical seconds. A source pulse is specified
by twelve integrated port coefficients c_p. Under a declared rectangular
control profile with per-port amplitude cap nu, its shortest duration is
max_p |c_p|/nu. The executable sums this quantity over all 432,425 base-W
pulses, then adds 864,850 explicitly supplied load/export latencies and any
explicit native-attempt latencies. Literal binary64 coefficients are retained
as exact dyadic rationals. A separate outward rational Q(sqrt(5),sqrt(2))
enclosure records the ideal-coefficient duration and its difference from the
literal source.

The profile is a conditional implementation model. Doubling nu while halving
the service and attempt latencies halves its duration without changing the
source word. Extra idle time can lengthen it. The bare native rule therefore
does not select that duration. The calculation covers the base W word;
clock-tick, detector, initialization and retained-current extension words add
their own pulses and latencies. Their costs cannot be replaced by the base
word count.

A refinement requiring duration tau0 h must supply corresponding bandwidth,
service-latency or distributed-execution scaling. With fixed serial bandwidth
and positive service latency, a word has positive duration independent of h;
its speed in length units h tends to zero as h decreases. Choosing the
effective clock does not prove the implementation can meet it. The supplied
size-two compiler is not a compiler family with a proved uniform resource
bound over all graph sizes. For strict implementation convergence, the
accumulated bound n times the per-W error must vanish as h decreases. Fixed
binary64 source precision does not establish this requirement. The finite
512-step bound and the ideal-source refinement theorem have distinct scopes.

## Stroboscopic locality and encoded processor states

R changes only event slots at one position. S either holds that position or
crosses one actual graph edge. Consequently a matrix element of W^n between
positions at graph distance greater than n vanishes exactly, including for
arbitrary retained references. A current shift attached to the same local
transport preserves this position-support statement. The exact numerator
replay checks complete position input columns for one through four steps.

The serial address-encoded compiler has a different intermediate-state
interpretation. Its fifth reversible gate, a bit-five flip, takes the
encoded position zero with event zero to encoded position three. Those ports
have graph distance three. Other prefixes occupy source-padding addresses.
All final S columns return to clean work and move at most one graph edge.
Thus the compiled word proves the graph light cone at completed-step
boundaries; its individual processor pulses do not establish a distributed
graph-local interaction law.

A concrete local implementation uses twelve event modes at each graph site,
on-site R operations, and matching-mode transport only between adjacent sites.
Its native code services and many-body statistics require a separate
implementation proof. A global shared event coin in a many-body extension
also requires a cluster-locality check; the one-load light cone does not
establish it.

The receipts provide an executable common clock interface, native phase
records, exact classical thinning and explicit implementation conditions.
They make no claim of a measured light speed, Lorentz proper time, autonomous
mass clock, interacting continuum field or official OpenWave cell award.
