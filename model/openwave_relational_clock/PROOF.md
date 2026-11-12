# Operational native clocks and wave arrival

## One fixed construction and a complete clock instrument

The source is exactly the cap-four construction in
`openwave_geometry_refinement/source.py`, with all80000 source coordinates.
The left native cell is prepared in level i=0 or1. The right native cell is
prepared in (|0>+|1>)/sqrt(2). The geometry is in its vacuum. A retained
preparation flag supplies a half-weight field-vacuum experiment and a
half-weight exact packet experiment. These are different preparations of
one source; no state-transfer assumption identifies them.

The right-cell clock read has a fair setting bit. The X setting reads the
eigenstates(|0>±|1>)/sqrt(2); the Y setting reads
(|0>±i|1>)/sqrt(2). Each of the other six native levels is an explicit outcome
in either setting. The sixteen setting/native effects sum exactly to I₈,
as checked with symbolic algebra. The setting choice is retained, and a
simultaneous sharp X and Y outcome is not assigned. The clock phase origin
is the declared equal-phase preparation at event zero; its four active
probabilities are1/2,1/4,0,1/4.

Let v₀,v₁ be the environment amplitudes of the two native clock branches,
with the external1/sqrt(2) already accounted for. In active order
X+,Y+,X−,Y−, the conditional-field probabilities at every environment
coordinate are

q_phi = |v₀+exp(−i phi)v₁|²/8,
phi=0, pi/2, pi, 3pi/2.

They sum to (||v₀||²+||v₁||²)/2. Their sign differences reconstruct
Z=<v₀,v₁>. The receipt supplies the full geometry/site arrays; it multiplies
these arrays by1/2 for the original field-mixture preparation. Native
complements have explicit zero weights on this input. The phase-setting
apparatus is a specified finite measurement service, with its source work
retained. No implementation by a narrower passive repair scheduler is claimed.

## The event word and its partial order

The experiment uses sixteen copies of the fixed source unitary
U=exp(−iH/2). The label one-half specifies this gate, rather than a fundamental
global-time variable. The event register is finite. Following each U, a
receiver effect reads x≥5. The no-hit branch continues coherently. A hit
terminates that branch and reads the native setting/sign, geometry number
and field position. All timeout records are kept, including the field vacuum.

Within one terminal layer these three effects act on distinct tensor
factors and commute. Symbolic Kronecker algebra checks this on all eight
clock levels, including the complement. The executed arrays also check
clock-first and receiver-first probabilities. Exchanging adjacent commuting
readouts connects every topological ordering of that layer; the complete
Kraus products, conditional quantum outputs and probabilities agree.
Evolution gates and readout layers are not freely reordered.

Completeness of the clock effects and receiver/no-hit effects gives V*V=I
for the whole history word on the entire source carrier, also after tensoring
an untouched reference. An orthonormal extension supplies a finite unitary
dilation with setting, history and observer registers. Copying an orthogonal
record |r> to |r>|r> provides two observer records. This operation does not
copy an arbitrary quantum state. The finite completed-instrument identity
is the same algebraic shortcut used in
`finite_observer_readout_closure/PROOF.md`; its native scheduler and source
probabilities are not transferred to this construction.

The phase read occurs only on the terminal hit or timeout branch. The
vacuum controls for individual calibration events use independent complete
preparations. Repeated sharp clock observations are not inserted invisibly
into the evolving no-hit branch.

## Calibration from actual native records

For every calibration event n, the vacuum experiment produces the complex
clock coherence Zvac(n). Its phase is enclosed from the certified source
state errors. The intervals have positive adjacent increments below pi,
starting from the event-zero phase. This proves a unique increasing phase
lift over the finite event domain. Divide that lift by2pi to define a
reference clock coordinate in cycles.

The same background-zero event-to-cycle dictionary is used for both
source-on and source-off wave-arrival laws. Let p_i(n) be their actual
clock-superposition arrival weights, including clock backreaction on the
wave, and c(n) the calibrated cycle coordinate. The reported statistic is

C_i = sum_n c(n)p_i(n) / sum_n p_i(n).

It uses observed clock-phase and arrival probabilities. It does not insert
a static ratio or convert a different wave preparation by a mean lapse.
The two hit probabilities are approximately0.167957103 and0.167730113.
Complete finite propagation, calibration and denominator bounds give

C_1−C_0 in [0.0013506,0.0013547] native clock cycles.

An additional directly joint observable is Zarr_i=sum_n Zhit_i(n), reconstructed
from the four aggregate arrival/sign frequencies. Its visibility is about
0.524858 and0.526111 for the two backgrounds. The relative phase

−arg(Zarr_1 conjugate(Zarr_0))/(2pi)

lies in [0.0004679,0.0004687] cycles. This circular coherence statistic differs
from the mean arrival coordinate. A common constant phase-reference change
cancels from the relative product. The calibration phase is measured relative
to its event-zero preparation, so its arbitrary coordinate origin also cancels.

The displayed clock cycles are simulation-generated observables of this
specified observer experiment. The construction supplies no physical-second
calibration, Einstein proper-time theorem or relativistic continuum dictionary.

## Errors, work and finite sampling

The source, positive-root and Chebyshev error certificates are those of the
same fixed geometry package. A sixteen-step history with vector errors
δ₀,δ₁ has clock-overlap error at mostδ₀+δ₁+δ₀δ₁, plus the bounded complex-dot
rounding. For a complex disk centered at z of radius e<|z|, its argument differs
from arg(z) by at most e/(|z|−e). Interval atan2 and exact rational norm bounds
certify the phase enclosures, including branch choices. Ratio bounds retain
the original hit probability and timeout; they do not discard failure weight.

The receipt also gives the mean source work of each terminal native-clock
and position read, together with its finite numerical error. Since the source
spectrum lies in[0,32], every complete read has absolute mean work at most32.
No physical apparatus cost is inferred from this source bound.

For a finite sampling version, use N=2^40 independent prepared trials per
background and retain every original timeout, setting and field flag. Let
eta=1/200000. Simultaneously controlling the eight aggregate unconditional
arrival/sign frequencies within eta has failure probability at most

16 exp(−2N eta²) <10⁻²⁰.

Because the original field mixture contributes a factor1/2,
Z=4[(q_X+−q_X−)+i(q_Y+−q_Y−)]. Thus frequency errors bounded by eta imply
an overlap-disk radius at most12eta. The enlarged phase enclosure still has
strictly positive lower endpoint. This proves a finite statistical
separation for the declared preparation and setting service. The large
repetition experiment is not executed; it is an analytic resource bound.
