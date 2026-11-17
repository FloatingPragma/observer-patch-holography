# Row 1: repaired observer interfaces

This is an additive replacement for two interfaces of the frozen `../row1/`
candidate, manifest `85f9eefcf4c2292c1a2f47f522149cf1e28852705376d3abc2a65024243c7f3a`.
That version and the negative [Fable-B review](../../../../openwave-push/reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B/REVIEW.md)
remain intact. The review is the source of the defects and of the closing
difference-equation construction. These are component repairs, not an award
of row 1 or an identification of physical electric charge.

## A supplied field, a resolution record and a conditional reading

`phase_reader.read_phase` requires a separate field argument. The integer
record F sets n_e=8|F_e|+8 before any sample is taken. A separate producer
`declared_field` can encode an independent integer F'. The reader never
reads F' or constructs its field from F. It checks the carrier's oriented
closure, unit-circle samples, field gluing at every corner and closing
point, every principal increment against pi/2, and each cell's integral
degree to tolerance 1e-9. Failed checks return a reason and no degree.

For the declared family, the phase is

    phi_e(t)=g_a+(g_b-g_a)t+2 pi F'_e t+A_e sin(2 pi t),
    |g_v|<=0.31, |A_e|<=0.23.

If |F'_e|<=2|F_e|+1, then

    sup|phi'_e|/n_e <= [2 pi (2|F_e|+1)+0.62+0.46 pi]/[8|F_e|+8] < pi/2.

The last strict margin has numerator 1.54 pi-0.62>0. Each principal
increment is therefore the unique actual increment; summing telescopes.
Vertex g terms cancel on a closed contour, giving its actual field winding.
For F'=F the stronger pi/4 bound from the prior review still holds. The
independent field can disagree with the record, and the comparison gate
then rejects. Reversing a nonzero winding is an executed rejecting control.

The default and maximum registered resource budget is 100000 sample values.
A caller may reduce it, but cannot expand the registered experiment.
The required count is computed in exact integer arithmetic before any
callback or large-integer float conversion. Budget exhaustion is retained
as a resource stop, not a zero charge or failed physical prediction.

The source must justify the field's regularity. Samples alone cannot do
that: adding n_e turns on an edge gives exactly the same sampled U(1)
values. The verification deliberately retains this alias counterexample.
Every returned successful reading says it is conditional on regularity;
`regularity_verified_from_samples` is always false. This is a mathematical
limit of finite sampling, not a remaining programming defect to hide by
another numerical threshold. The fields used in the positive controls
meet the bound analytically; the alias fields do not.

## The closing records determine charge without a remembered start

`closing_reader.observe_closing` receives exactly thirty E records and a
region. There is no scalar filling, preparation label, earlier current
history, initial per-port array, or terminal snapshot. Each seam occurs
once, in the order actually recorded. A record contains its seam (a,b),
coin c, current j and pre-event contrast d. The unchanged native endpoint
law gives

    j=[d-(2c-1)(|d| mod 2)]/2.

Let delta_v(k) be cumulative changes from the first k closing currents,
with delta_a reduced by j and delta_b increased by j. If x is the unknown
state at sweep start, each recorded contrast supplies

    x_a-x_b=d_k-delta_a(k)+delta_b(k).

The icosahedral seam graph is connected, so these equations determine x
up to one common additive constant, if and only if their cycle constraints
are consistent. The observer solves them with temporary gauge x_0=0 and
subtracts their mean. This gives the unique initial deviation vector
y_v=x_v-(sum x)/12. Uniform shifts of the physical filling have no effect.
The observer checks that y is integral and is either zero or exactly one
+1 and one -1. It replays all thirty native events from y and checks the
contrast equations and sector at every step. This sector is invariant
under the native law: differences 0 do nothing, differences of magnitude
1 move or hold an excess/deficit, and an adjacent excess/deficit annihilates.

The observed regional charge is sum_R[y_v+delta_v(30)], in {-1,0,+1}.
The pointer command is its negative, following the existing orientation
convention. The absolute filling is neither reconstructed nor required.
A nonnegative physical representative exists at every uniform filling
q>=1. The known neutral-pair preparation remains a way to prepare this
sector, but is not an input to the new observer or its measured charge.

This uses the same emitted E contrasts as the prior experiment. It does
not show that a continuum physical field emits those integer contrasts,
nor that the preparation occurs spontaneously. It reads the state at the
last recorded event. Unrecorded movement afterward is outside that
measurement time and is explicitly retained as an indistinguishability
control, not claimed to be detected.

## The complete conditional instrument

The instrument begins at the start of the closing sweep. Its source
space has 133 possible relative states; its pointer has twelve possible
unknown initial ports. The closing records fix exactly one source column.
Every thirty-coin closing history has probability 2^-30, including all
coins which produce zero current. The source-prefix preparation and its
probability are not part of this newly defined conditional experiment.

`admit_closing` checks the existing native pointer policy for all twelve
initial columns, reconstructing its terminal loads from its actual
coin/current record. A supplied `final_loads` is ignored. If the pointer
record contains m native attempts, the fine effect is diagonal: it is
2^-(30+m) on each compatible pointer column of the uniquely measured
source column, and zero elsewhere. Routing, discovery and measurement
timeouts remain distinct outcomes with no readout. No success renormalization
is performed.

The corresponding branch operation is explicit. For the measured source
start y_i, its reconstructed finish y_f, compatible pointer starts x and
their reconstructed finishes p_f(x), take

    K_t=2^(-(30+m)/2) |y_f><y_i| tensor sum_x |p_f(x)><x|.

A fixed current tape adds the same known load-change vector to every
compatible pointer start, so different starts have different finishes.
Thus K_t^dagger K_t is precisely the reported diagonal effect. The map
rho -> K_t rho K_t^dagger retains the branch's actual terminal states;
no ideal unitary pointer or unrecorded reset replaces the native events.

Completeness follows by summing the two fair-coin branches at each bounded
native step. `EXACT_LAWS.json` additionally executes every source-sector
branch and sums all 2^30 coin histories for each of 133 starts by exact
integer dynamic programming. The pointer aggregate probabilities are
calculated independently from its discovery, routing and cycle stages.
For retries R=2, put h=2^-R, r=1-h. Each port has five incident seams.
Discovery timeout is h^5. First discovery on incident seam k has weight
h^k r; if the reached port has route length L to the cycle start, its
routing timeout, measurement timeout and success weights are respectively

    h^k r (1-r^L),
    h^k r r^L (1-r^(6|degree|)),
    h^k r r^(L+6|degree|).

Summing over k and adding discovery timeout gives exactly one for every
pointer start and each of the three commands. All fine histories are
retained in this law; the finite sample tapes only test the executing API.

## Source and status boundaries

The phase reader and closing observer are distinct component interfaces.
No new theorem identifying their encodings with electromagnetism is claimed.
The unchanged source law, geometry, pointer and strict validator are pinned
in the manifest. The observer policy is changed, so this is a new candidate.
The six old filing passages identified by Fable-B still describe the old
candidate; this directory does not rewrite a filing or move its public pin.
Independent review of these new bytes is requested through the sync file.
