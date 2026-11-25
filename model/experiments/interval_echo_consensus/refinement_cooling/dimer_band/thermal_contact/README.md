# Measuring equilibrium with a third body

A thermal contact made from Interval Echo exchanges has an exact equilibrium criterion. Under this contact, the joint state of two separately prepared bodies stays stationary exactly when both density matrices are Gibbs states at the same inverse temperature. This holds for positive definite product preparations of any two nonempty bodies made from the guarded dimer registers. Each body can initially contain arbitrary internal correlations and coherences. A second contact compares the same thermometer to a transition with twice the energy gap. The instruments require a supplied state-independent calendar, sharp energy and position readers, and a controlled contact window. [Protocol](PROTOCOL.md) · [Native contact and complete simulated records](measure.py) · [Two-gap protocol](DOUBLE_PROTOCOL.md)

## Physical contact and fixed points

Each register holds one double on a native pair. Its two stationary modes have energies 32 and 33 in exchange units. At each appointment a random tag tells the controller to idle, to open one dimer for a free half-period, or to exchange one entire dimer of each body. Every option has positive probability. Independent fair material half-contact readings index the options, with spare tags assigned to idle. A half-period implements the sign operator Z in the energy basis. To exchange registers the two occupied rail addresses undergo two transpositions made from native π exchanges, each followed by the reverse route needed to restore intermediate ports. On blank memory the full π pulse equals the load swap on all nine amplitude columns including their phases. The path check covers every pair of the six deployed dimers on the 20-owner board. A separate four-port matrix check verifies the logical SWAP on every energy-basis column. Consequently every completed contact conserves the sum of the measured dimer energies. Timing, routing and the lamps that record the energy signs are specified apparatus. [Operator and routing tests](test_contact.py)

For the mixture T of these unitaries, if Tρ=ρ then every positive-weight unitary fixes ρ under conjugation: the weighted sum of the squared Hilbert–Schmidt distances from their mean is zero, since each conjugate has the same squared norm as ρ and the mean is ρ. The converse follows by linearity. Invariance under every local Z makes the state diagonal in the joint energy basis. Let a faithful independent input then have joint probability A(x)B(y), with no factorization inside either body assumed. Exchanging any two unlike signs gives the identity

\[
\frac{A(x_i=1,u)}{A(x_i=0,u)}
  =\frac{B(y_j=1,v)}{B(y_j=0,v)} .
\]

The contexts u and v on the remaining registers vary independently for every i and j. All single-site upper-to-lower odds must equal the same positive number r. Starting at the all-lower word and flipping sites therefore gives A(x) proportional to r raised to the number of upper signs in x, and likewise for B. Conversely these two distributions have joint weight depending only on the total number of excitations, so every offered operation fixes them. Even for initially correlated bodies, invariance under all the phase operations removes every off-diagonal entry, and the cross-body swaps generate every permutation of sites. The remaining joint weights are exactly arbitrary functions of the total excitation number. Requiring the two bodies to be independent restricts these functions to the common Gibbs odds derived above. The contact therefore singles out the Gibbs preparation on the whole faithful independent-input domain. The unique inverse temperature in the measured unit gap is β=−log r. It can be zero or negative on this bounded spectrum. If A and B are each stationary with the same faithful thermometer C, their odds agree with C’s and the A:B joint state is stationary. The proof applies to arbitrary finite body sizes when the contact can visit all participating registers. On a coherent input the contact is applied before any sign read. The two-time energy-word assay uses already diagonal preparations. An absolute kelvin scale is unnecessary for this temperature comparison.

## Independent thermometers and correlation controls

The three bodies contained one, two and three dimers. A and C came from direct filters. B came from a separately constructed bounded-retry repair. Independent local calibrations measured upper populations **0.33219, 0.33444 and 0.33317**. Of 524,288 original starts in each contact bank, 390,464, 219,332 and 165,630 A:B, B:C and A:C pairs were prepared. All first refusals are stored by owner. Each complete energy word before and after contact, with its tag, entered the census. The maximum observed change in any word probability was 0.00053, 0.00081 and 0.00082 respectively. The A:C test used a fresh cohort. Heating and cooling C changed the mean energy entering A to +0.06268 and −0.08550 per accepted appointment, against +1/16 and −1/12. When the kinetic exchange was disabled in the contact, the flow was exactly zero despite the measured population difference.

A more stringent control prepared B with the same local temperature readings as A but with a hidden positive correlation. On a fresh fair flag, the source read both sharp energy labels and copied the first to the second through an addressed private phase when they differed. This is a copy of a classical read. The observed mean flow into A was −0.00033. The probability that both B signs were upper fell by **0.02730**; the exact predicted drop is 1/36. Merely measuring zero mean flow would have misclassified this contact. In a separate phase test all four input dimers occupied the first sharp physical rail, which is a coherent superposition of the energies. Equal-register swaps alone preserve it. With the free half-periods in the menu 262,375 of 524,288 launches changed their recorded position; parking those phases gave zero. The tests use simultaneous 99% Hoeffding allowances propagated from the independent thermometer calibrations. [All calibration, control and original-launch bins](records.json)

## A material thermometer with twice the gap

Two dimers whose signs agree have energies 64 and 66. Position correlations after free flight, with two local quadrature settings, isolate the double-frequency term cos(2t). The two underlying native gaps and both time settings are measured in the same units as the single-gap clock. In the contact the two independent B signs are read. If they agree, a fair offer swaps both dimers with C's equal pair. Other offers leave the complete word untouched. This sign-reading transaction applies to the dephased energy records; the preceding theorem supplies the coherent-input analysis for the single-gap contact. Both directions of the material transfer move two units of energy. Balance requires

\[
\frac{p_C}{1-p_C}=\left(\frac{p_B}{1-p_B}\right)^2 .
\]

Thus an upper population of 1/3 in the single-gap source predicts **1/5** in the double-gap thermometer at thermal equilibrium. In 162,449 accepted pairs out of 524,288 original launches, the latter gave **0.20073**, with net energy transfer +0.00261 per accepted appointment. A deliberately prepared double-gap source with upper fraction 0.33319 instead transferred **+0.07287** into B, against +2/27. Parking exchange gave zero. The calibration, all joint words, no-ops, invalid double-gap preparations and refusals are retained. [Clock and contact counts](double.json) · [Two-gap code](double_gap.py) · [Independent comparisons](test_double.py)

This construction establishes equality and transitivity of an operational temperature in the stated material contact family. It concerns sealed registers and a finite contact measured at its endpoints. Identifying physical heat exchange in freely interacting unshuttered matter, or attaching an energy unit or a kelvin to the readouts, requires additional evidence.

From the RER root, with numerical-library thread counts set to one:

```sh
pkg=model.experiments.interval_echo_consensus.refinement_cooling.dimer_band.thermal_contact
.venv/bin/python -B -W error -m "$pkg.measure"
.venv/bin/python -B -W error -m "$pkg.double_gap"
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/refinement_cooling/dimer_band/thermal_contact
```
