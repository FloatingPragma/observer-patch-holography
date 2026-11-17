# A field emitted by the same finite reciprocal-work event

This packet observes the unchanged R-work law. Its new claim is an exact
finite characteristic-field equation, local budget continuity, and a source
contact → field → reciprocal detector join. Physical energy and electromagnetic
identity remain unproved. The separate native checkpoint experiment supplies
storage/readout, not a hidden change to R-work's declared transport services.

## 1. Public source and local field

The public reference specifies twelve signed loads and two slots per ring
observer. Native records change the two addressed loads by (-j,+j). An `emit`
record additionally places `(e,c,j)` in both blank slots. An `absorb` record
removes matching slots and applies the inverse load translation. A swap record
exchanges its two addressed slots. These records reconstruct all later loads
and slots without reading a hidden source state. Complete detector failures,
silent native events and explicit interventions remain in that reconstruction.

Let u_i and v_i be zero for empty slots and the signed current j otherwise.
Define E_i=u_i+v_i and B_i=u_i-v_i. This is an invertible real observation of
the two signed amplitudes: u=(E+B)/2, v=(E-B)/2. It does not reconstruct the
seam/coin metadata, which the actual reciprocal detector also requires.

The identity

    ε_i := (E_i²+B_i²)/2 = u_i²+v_i²

is exact. Since each occupied slot has j=±1, ε_i equals its occupied-slot count.
Thus this field norm equals the slot term of the old, independently declared
budget W=Σ_i,a x_(i,a)²+Σ_i occupied_slots_i. No new scale or fitted energy
functional was chosen after observing a propagation curve.

## 2. Field dynamics derived from the executed swaps

R-work first swaps the two slots at every site, then swaps lane1 at i with
lane0 at i+1. Composing these two operations gives, for arbitrary slot contents,

    u_i(t+1)=u_(i-1)(t),     v_i(t+1)=v_(i+1)(t).                 (1)

This identity follows from the actual swap law. The reader does not evolve
its own wave equation or apply a fitted propagator. Taking sum and difference
in (1) yields

    2E'_i=E_(i-1)+E_(i+1)+B_(i-1)-B_(i+1),
    2B'_i=B_(i-1)+B_(i+1)+E_(i-1)-E_(i+1).                    (2)

Each characteristic propagates one site per complete cycle, so in a region
without emission, absorption or intervention both E and B satisfy

    F_i(t+1)+F_i(t-1)=F_(i-1)(t)+F_(i+1)(t).                  (3)

On the linear span of observed impulses, the actual one-cycle matrix T obeys
T^T T=I and T²-ST+I=0, where S is the spatial nearest-neighbor sum. Its Fourier
phases are exp(-ik) on u and exp(+ik) on v for k=2πm/L. Equivalently, in the
(E,B) basis the mode matrix is cos(k)I-i sin(k)σ_x. The finite phase is determined
by the two local swap layers, not an oscillator frequency supplied separately.
Zero and Nyquist modes are retained, including their degeneracies.

The formal smooth, long-wavelength identification with one transverse field
component in one spatial dimension is ∂_t E=-∂_x B, ∂_t B=-∂_x E, with one edge
per cycle as the declared unit. This identifies a mathematical wave structure.
The experiment does not supply a smooth physical continuum, a laboratory clock,
an electromagnetic detector/current normalization, or a quantum energy-frequency
relation. A classical Fourier phase is not by itself a photon energy.

## 3. Local continuity and the actual source event

Define the bond flux during one transport cycle by

    J_(i+1/2)(t)=u_i(t)²-v_(i+1)(t)².

Using (1),

    ε_i(t+1)-ε_i(t)=J_(i-1/2)(t)-J_(i+1/2)(t).                (4)

Summing (4) on the ring gives exact conservation of total field norm. For a
right-only occupied slot E=B=j and the outgoing norm moves right; for a
left-only slot E=-B=j it moves left. The continuum expression EB=u²-v² is the
same-site difference of characteristic intensities; (4) uses the precise
finite bond flux rather than replacing it by a continuum expression.

For an initial contact dipole x_a=j, x_b=-j, both native coins give current j
and zero output loads. R-work then writes the two actual committed records.
The field at that site is E=2j, B=0 and ε=2. Therefore the same event changes
load-square by -2 and field norm by +2. The released field is tied to the
actual recorded annihilation current, including its sign; it is not injected
as an unrelated planted pulse in the contact experiments.

The two records travel to the antipodal observer and first coincide after
L/2 cycles. Its complete detector checks equality of both records and zero
recipient endpoints. It consumes the records and recreates x_a=j,x_b=-j,
changing field norm by -2 and load-square by +2. The same reference defines
all three terms of this source → field → detector budget. Metadata mismatch,
single arrival, busy endpoints and disabled detectors remain explicit outcomes.

This is stronger than passive bookkeeping: a transported output changes a
remote observer, and a tag intervention can preserve field norm while blocking
the inverse response. It remains a dimensionless conserved resource of an
explicitly enlarged source. The finite equality does not calibrate rest energy,
frequency energy or laboratory radiation energy.

## 4. Native checkpoint implementation and its limits

Every checkpoint slot is assigned the preregistered symbol 0 or 1+4e+2c+(j+1)/2.
The native service encodes that symbol in two binary total-two twelve-port
registers, with a supplied base66 codebook. Full78 current-only reset, native
mixed-edge writes and current-only census are actual integer `native()` calls.
The controller receives current records, not the private initial register
configuration. Complete native fine words make the arithmetic independently
replayable. Fresh checkpoint-specific coin tapes and initial states are used.

Both native carriers have total load2 at every stage, so every encoded symbol
has exactly the same total load4 as the blank symbol. Relative to that stated
blank, the encoded information has sharp native charge difference0. This is
an implemented relative-charge statement; it does not identify the carrier's
absolute native load with electric charge or establish a neutrino species.

For a conclusive pair of symbol services, the measured packet currents give
the same E/B as the original symbolic record reconstruction. Inconclusive
services carry no guessed field value. The existing proven service bounds
compose without independence: per symbol, failure ≤35/536870912 and at most
10652 native attempts. The 960 normal services use a union bound; the 240
disabled-write services have a separate denominator. A forced write timeout
and a forced census seed timeout expose both failure flags and retained work.

These checkpoint instruments prepare fresh native representations of public
records. They certify the declared encoding/readout map. They do not show the
same native registers physically executing the entire symbolic world history.
Their reset/write/read work is additional apparatus work. Their bare native
load-square is constant for all encoded symbols and is not the field norm ε.

## 5. Evidence and interpretation

The frozen DESIGN/TASK/code precede the evidence. Corrected execution passes
all nine frozen gates. Its complete receipt and compressed records retain:

| Quantity | Result |
| --- | --- |
| Source/control cases | 493: 360 contact flights, 12 controls, 112 impulse columns, 9 planted local readers |
| Symbolic records, including public references and observation marks | 355,205 |
| Retained literal-native attempts in the R-work source | 86,751 |
| Contact launch and antipodal absorption | 360/360 correct, with budget 2 throughout |
| Exact three-time wave equalities | 147,840/147,840 |
| Complete finite mode responses | 112/112; maximum phase residual 2.6686×10^-14 |
| Native normal symbol services | 960 conclusive and correct; 0 inconclusive |
| Native disabled-write services | 240 conclusive blank; 0 inconclusive |
| Native E/B checkpoint pairs | 480 conclusive and correct |
| Retained native memory attempts / stages | 1,341,601 / 6,722 |
| Forced write / census failures | 32-attempt write timeout; 120-attempt seed timeout; neither guesses a symbol |

All field, current, guard, budget and native replay defects are zero. Every
control discriminates as preregistered. The union failure bound for the 960
normal services is 525/8388608; for all 1,200 normal/disabled services it is
2625/33554432. These are analytic worst-case bounds, not rates inferred from
this deterministic fair-bit realization. The native-call counts above count
retained source attempts; cached controller-table computation is not presented
as additional measured histories.

`headline.py` reconstructs fields and rechecks native fine records; it does
not trust stored pass flags as evidence. The invalid recorder-alias attempt
and its original source execution are preserved separately and are not used
in these numbers. No gate or source parameter was changed after that failure.

This packet does not merge the row29 occupation-Maxwell source or the RER
installed radiative-pair Hamiltonian with the R-work source. Those remain
different constructions. The neutral-even PMNS comparison remains failed;
this field has no flavor law and provides no row10 pass. The existing local
outgoing-shell proposal and physical radiation/energy identification are distinct
remaining row8 tests. No new universal prerequisite is imposed on that platform
criterion, and positronium branching remains its explicit follow-on.
