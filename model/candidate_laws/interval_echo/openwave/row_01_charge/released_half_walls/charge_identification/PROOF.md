# When a local calibrated meter must report half a unit

Consider a ring of sharp two-stock readings $s_j\in\{0,1\}$. The occupied stock has calibrated load $A=1$ or $2$. A detector can use an arbitrary real function $f$ of $w$ consecutive port readings. The same function is used at each translated placement, and an extended reading sums the outputs of its local placements. The input is the complete observer word. A detector sensitive to an independent field or an unread material variable has a different input. For a noisy detector one can use its calibrated conditional mean on each word; conclusions then concern expected signals.

Three assumptions on this apparatus can be tested separately. Its local excess vanishes on both alternating $w$-site words. Its integrated output on known sharp stocks agrees with the added load $A(\sum s_j-N/2)$. Two defects whose entire detector neighbourhoods fit disjointly on the measured word can be assigned distinct patches. The stock calibration below checks the second assumption without fitting any separated wall measurements.

## Compact stocks determine the integral on every sharp word

Set $m=w-1$. A binary $w$-site window is a directed edge $e:u\to v$ between its initial and final $m$-bit words. A periodic ternary-port record on the binary stock gives a closed walk in this graph. Use the symmetric density $A(s_j+s_{j+1}-1)/2$ as a mathematical reference. Its integral counts $A(\sum s_j-N/2)$ on every word. Write $d(e)$ for the local excess of $f$ over this reference.

Prepare each of the $2^w$ binary strings of length $w$, with zeros on all other ports. The ring has $N\ge 2m+1$ sites. The stock with no occupied port fixes the residual $d_0$ on the zero self-loop. From vertex zero, append the bits of $u$ to reach $u$ in $m$ steps. Let $h(u)$ be the residual accumulated on this prefix. Appending $m$ zeros at any vertex $v$ gives a return path with residual $z(v)$. The compact stock spelling $e$ has residual
$$
T_e=h(u)+d(e)+z(v)+(N-2m-1)d_0 .
$$
The compact stock spelling $0v$ has residual
$$
T_{0v}=h(v)+z(v)+(N-2m)d_0 .
$$
Consequently $d(e)-[h(v)-h(u)]=T_e-T_{0v}+d_0$. When the compact stocks pass exactly, $d$ is a difference of endpoint values; it sums to zero on every closed walk. The construction only uses sharp stock integrals with known total insertions and removals. It accommodates nonlinear local responses.

For a finite-resolution instrument, suppose the absolute integrated error is at most $\delta$ on every compact calibration. Then $|d_0|\le\delta/N$, and the residual after subtracting the endpoint difference is at most $(2+1/N)\delta$ per window. The integrated error on any held ring word is at most $(2N+1)\delta$. The bound allows calibration errors of either sign; it makes no independence assumption. The [Lean file](LocalCharge.lean) proves the padded-edge identity, its error inequality and cancellation on an arbitrary finite cycle. The [integer certificate](census.py) checks the padding identity for every window at the two registered sizes, by comparing counts of all possible local words independently of their assigned outputs.

## Separate readings on the emitted defects

The accepted whole-ring word contains two occupied-equal bonds following addition, or two vacant-equal bonds following removal. Every other bond alternates. The patch of a particular wall consists of the $w-1$ window starts that contain that bond. When the two patches are disjoint, each window in the first patch has a translated copy in the second: its local word is set by the sign and its offset from the only equal bond inside it. The two histograms are equal. Every remaining window is one of the two locally calibrated vacua. Thus the two patch outputs are equal, their sum is the added load, and **each is $+A/2$ or $-A/2$** on that individual sharp record. The proof applies to the integrated result of any calibrated detector in this class.

If individual vacuum responses are bounded by $\eta$, there are $N-2(w-1)$ remaining windows. If the two integrated patch calibrations can differ by at most $\xi$, the error per inferred half is bounded by
$$
\tfrac12\big[(2N+1)\delta+(N-2(w-1))\eta+\xi\big].
$$
For example, on the fourteen-port ring, an exact translation calibration, exact vacuum zeros and $\delta=0.01$ ordinary units give a conservative radius $0.145$ ordinary units around each signed half. This is a required accuracy example; no uncertainty is assigned to an unperformed external-probe calibration.

The census checks the two histograms and both vacuum patterns on the sharp words from every retained cohort. Across 32 settings, seven release times and 4,194,304 original simulated trials per setting and time, all 77,669,791 gated outcomes pass for each $w=2,3,4,5,6$. These are repeated cohorts. For a positive single-load injection on the native 14-port level-zero ring at time three, 1,789,126 of the original 4,194,304 trials qualify. Turning off the flight contrast gives 66,814. These numbers reanalyze the [complete archived port records](../words.npz); the [compact reduction](counts.json) retains the widths, totals, all times and controls. No fresh acquisition is used for the census.

A local raw two-ended contrast $A(s_j+s_{j+1}-1)$ gives a unit response on a positive single-stock wall. One point-like increment changes two adjacent bond responses by a unit each. This meter doubles the calibrated integral. For six-site lookups it fails all 64 compact-stock normalizations on both rings. Dividing by two gives the globally normalized response and a half at the wall. A separate foil adds the gradient of the last bit to the local density. It passes every integrated calibration, while separate local readings of the two alternating backgrounds expose it. A gradient supported on two consecutive occupied bits instead passes both tests and redistributes the density between adjacent windows. The integrated two-patch rule survives. A sharp flat port window uses an integer host subtraction and its own boundary convention. Electromagnetic identification requires a separate calibration of a physical electric field and charge on this source.

Run the deterministic concentration-of-records and calibration checks from the repository root:

~~~sh
p=model/candidate_laws/interval_echo/openwave/row_01_charge/released_half_walls/charge_identification
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B "$p/census.py"
(cd Lean && lake env lean -j1 -DwarningAsError=true "../$p/LocalCharge.lean")
~~~
