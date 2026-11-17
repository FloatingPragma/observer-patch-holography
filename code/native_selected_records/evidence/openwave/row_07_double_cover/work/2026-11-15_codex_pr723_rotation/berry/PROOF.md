# Native optical holonomy with a complete retained instrument

## Source transfer follows from native endpoint records

There are three one-unit source states on a triangle. Ticket e addresses one
of its seams with probability p_e, and the native fair coin holds or moves the
unit when incident. Exact branch enumeration gives

`P = I - (1/2) sum_e p_e b_e b_e^T`, `sum p_e=1`.

Here b_e is the oriented incidence column; no target wave or spin matrix is
evolved. Every matrix entry is reconstructed from the native current table.
P preserves the constant ray and its orthogonal, zero-sum plane. The optical
input f=(1,1,-2)/sqrt6 belongs exactly to that plane.

In the orthonormal chart u=(1,-1,0)/sqrt2, v=(1,1,-2)/sqrt6, the derived matrix is

`P_opt = I/2 + [[(1-3p0)/4, sqrt3(p1-p2)/4],
                [sqrt3(p1-p2)/4, -(1-3p0)/4]]`.

The chart is used for stable computation of an invariant subspace. Applying
it in code is not a physical projection between events. Products in an
unnormalized rational chart give the same answer exactly, independently.

Uniform weights give the degenerate optical eigenvalue1/2. A positive weight
loop encircling that point transports the slow optical eigenline. The source
experiment applies the actual P products rather than replacing them by
instantaneous eigenspace projectors. The finite sign is established by exact
rational multiplication; no adiabatic-limit inference is needed for the
recorded numbers. Resolution affects the deviation from a returned ray.

## Complete quantum memory and eraser

Six native vacuum coin events create the declared coherent64-ticket memory.
Each body branch has amplitude1/sqrt128. Its retained flag is(ticket,coin,j),
where j=-1,0,1, giving384 labels. For a fixed flag, current continuity recovers
the unique input from the output, so the isometry's columns are orthonormal.
This is verified across all source columns, not only the chosen input ray.

Projecting the flag onto its uniform vector gives

`A = sum_flags K_flag/sqrt384 = P/sqrt3`.

The other383 outcomes are the rows of a stated Householder orthogonal basis.
Together their effects sum to identity. This eraser is a declared coherent
service. It is not a claim that merely declining to read a classical record
implements erasure. Committing flags first gives a different CP map: on the
initial ray, coherent success probability25/192 changes to1/384, and the
resulting state is mixed. The coherence is an essential premise.

For the two interferometer arms let a=P_loop f and b=P_ref f after N source
events. Original joint-success probability and fringe visibility are

`p_success = 3^-N (||a||²+||b||²)/2`,
`V = 2 Re<a,b> / (||a||²+||b||²)`, `p_minus|success=(1-V)/2`.

All removed vector norms in numerical calculations remain in the log of
p_success. ALL_OUTCOMES.json and ROTOR_OUTCOMES.json retain every first-
rejection flag at every stated stage, its reach probability, and the final
success. Rejected mass plus success sums to1. Failed trials are not omitted.

## Conditional polar-analyzer interpretation

An equilateral triangle relative to an external polar analyzer n(theta) gives
the declared ticket coupling p_e=(1+cos(theta-2pi e/3)/2)/3 before64-ticket
rounding. Substitution into the derived P_opt yields, before rounding,

`P_opt(theta)=I/2 + [[-cos(theta), sin(theta)],
                    [sin(theta), cos(theta)]]/8`.

Its slow eigenline has coordinates(sin(theta/2),cos(theta/2)). This half-angle
is derived from the native transfer and declared polar coupling; it is not
inserted as the source update. The actual experiment uses the integer ticket
tables and finite native transfer products. One96-step turn gives overlap
-0.998664740954 and visibility-0.997122765729. Repeating the exact same recorded
program twice gives overlap+0.998664740954 and visibility+0.991866004768.
The apolar cos(2theta) coupling gives positive overlap after one external
turn. Thus the analyzer-to-scheduler coupling is essential and still supplied.

The native process has not realized a spatially rotating analyzer, established
that ticket law from physical motion, or derived an SO(3) spin representation
for a persistent charged object. The construction is a conditional source-
record mechanism for relative holonomy, not a physical spin identification.

## Costs and robustness

At the96-event polar turn, original success is10^-85.0181555 despite a strong
conditional fringe. A constant component of the prepared source is invariant
under every P and eventually dominates the decaying optical part. This
requires an exactly specified optical preparation; the threshold calculations
for small leakage are reported, and no physical preparation precision is
claimed. There is no practical experiment asserted at these success costs.
Native source laws outside this selected ray, apparatus and schedule remain
unscored. No conclusion about the physical energy dictionary follows.
