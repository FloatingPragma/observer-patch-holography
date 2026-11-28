# What returns after one material circuit at 1045

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

In this finite written-field model, the receiver can count a complete circuit of a signed load and then read a quantum record at the returning charge. At a local wait of π, **0.40645% of all original launches** return to station 1045 after one net positive circuit and **0.28915%** after one net negative circuit. In both cohorts, comparing nine pairs of source and receiver readouts excludes every separable joint record. Two reference cohorts with equal initial word populations and equal real coherence subsequently prefer opposite circuit directions. Closing one ring flip removes this distinction even though the material still completes circuits. [Design and selection history](DESIGN.md) · [All circuit probabilities and refusals](records.json) · [Simulated detector counts](counts.json)

## Counting a real material circuit

A five-face wheel in the same receiving cell supports a stationary current at the original rates. Five addressed native pulses move the arriving charge and its two-word record into the wheel. After selecting a degenerate energy, all five side counters agree on the current, and a receiver two rim steps away can still verify its quantum correlation with the source. [Stationary five-face measurements](five_faces/README.md)

Use the same fresh reference write at 98, clock-bypassed four-fan flight, and sharp charge and unresolved face-0 pass at 1045. Under the adopted extended written-field flight rate the late gate at 30.8125 admits **1.58437% of original launches**. After this pass, shutter every other contact as a whole Laplacian term. Leave active the three material sides of the actual face, the flip of its three heads and the original degree potentials. A supplied register coupled to the 1049→1053 material hop counts +1 on that crossing and −1 on the reverse. The load address, both end words and this meter remain readable after free local evolution. The meter and shutters are part of the apparatus. [Original arrivals and nine-setting witness](../WITNESS.md) · [Microscopic checks](triangle.py)

| Position along the positive circuit | 1045 | 1049 | 1053 | 1045 |
| --- | ---: | ---: | ---: | ---: |
| Excess head word | 365 | 364 | 332 | 334 |
| Vacancy head word | 146 | 147 | 179 | 177 |
| Physical degree | 6 | 5 | 5 | 6 |

On either sign chart the three material exchanges have rate $t=4$. The two configurations at 1045 also communicate at ring rate $g=1/2$. This last transition transfers no material. The four configurations close exactly under the admitted contacts for all 16 upstream word combinations. Complementing every head maps every term of the full incoming flight and both detector bases between the charge signs. The source gauge here is a common minus sign. The counter follows material displacement in either chart. The signed transported charge has the additional factor $q=+1$ for an excess and $q=-1$ for a vacancy. Distinct registered windings become distinguishable and can change the interference observed by a reader who discards the counter.

For a hub return, conservation at each rim vertex forces the algebraic crossings of all three material sides to agree on every contributing walk. Thus a single registered winding requires traversing all three sides net once. Coherent paths can include reversals and unobserved extra loops. A resolved intermediate trajectory would require a different instrument. The 65-position cyclic meter has a possible integer alias. A different winding congruent to ±1 requires at least 192 material hops between hub endpoints. The resulting amplitude allowance is smaller than $2.76\times10^{-54}$ at every scheduled wait. The receipt includes every register and final word for each reference report and for the unread reference, all source refusals, both signed charts and the 129-position check.

## What the retained source reference distinguishes

The source analyzer direction $(-0.062372,0.938206,-0.340411)$ was fixed using the existing late arrival measurements before computing circuit probabilities. For the two source reports the conditional receiver records have equal $X=0.475443$ and $Z=-0.031678$ means. Their $Y$ means are $0.428521$ and $-0.100010$. To reconstruct the axis, subtract the product of the marginal means from each joint Pauli row, then take the unit vector perpendicular to the corrected X and Z rows. Its sign makes the Y difference positive. The material detector is independent of this calibration. The earlier appointment uses the same axis. [Arrival tomography](../word_reads.json)

| Native local wait | Positive circuit per original | Negative circuit per original | Directed score, source + | Directed score, source − |
| --- | ---: | ---: | ---: | ---: |
| π/8 | 0.04914% | 0.04419% | +0.004646 | +0.001263 |
| 3π/8 | 0.14210% | 0.11086% | +0.038341 | −0.003060 |
| π | 0.40645% | 0.28915% | +0.165657 | −0.037986 |

The first two columns include every original launch and ignore the source analyzer. The score is the positive-minus-negative single-circuit return probability conditioned on the initial face and the specified reference reading. The late numerical guard is $6.25\times10^{-7}$ per score or initial Bloch component. At π the difference is **0.203643**. A closed ring gives the same score **0.021633** in both cohorts. Halving the ring rate gives **0.172333** and **−0.024033**. These comparisons change no incoming preparation or material rate. Selecting a reference report partitions shared records and leaves all unread terminal probabilities unchanged.

For each return, the receipt gives the two-by-two amplitude map $K_m$ from the incoming to the returning pair and all receiver–source Pauli moments after that map. Applying $K_m$ to the local side sends the partially transposed initial joint record $R$ to $(K_m\otimes I)R(K_m^\dagger\otimes I)$. For an invertible $K_m$, evaluating the Rayleigh quotient on the inversely transformed negative eigenvector gives an upper bound $\lambda_{\min}(R)s_{\min}(K_m)^2<0$. At the late appointment and a local wait of π both directions have $s_{\min}>0.15093$ and an upper bound below $-9.60\times10^{-6}$ per original. The single-cycle maps have full rank at all three tested positive waits. This gives a finite implication from the measured incoming correlations and calibrated single-circuit transfer to the two observers’ returning measurement statistics. [Separability and partial transposition](https://arxiv.org/abs/quant-ph/9605038)

The direct detector check reads the two returning words in one of three bases and compares them with three independent source settings. Its witnesses were chosen after the model probabilities were available. For the positive and negative circuits the conditional witness intervals from simulated counts are respectively **[−0.004137, −0.003904]** and **[−0.008739, −0.008470]** at simultaneous 99% coverage across all six tested cases. A separable record has nonnegative expected score. Both early-arrival controls and both late ring-cut controls give positive intervals with these same instruments. The count receipt retains original-start denominators. [Witness definition and selection](COUNT_TEST.md) · [Count banks and intervals](counts.json) · [Independent count replay](witness.py)

The small-time mechanism can also be checked exactly. With chain ends ordered as the table and $u$ the local wait, the four positive-circuit entries first arise from three, four or five exchanges. The unique minimal paths give weights $t^3/3!$, $g t^3/4!$ for each diagonal, and $g^2t^3/5!$, with the common corresponding powers of $iu$. Therefore

$$\det K_{+1}=\det K_{-1}=\frac{g^2t^6}{2880}u^8+O(u^9).$$

The replay independently gets the native coefficient $16/45$ from integer matrix powers. This gives invertibility for all sufficiently small positive waits when the native transfer and ring flip are open. With $g=0$, a positive circuit can only connect word 365 to 334 in the excess chart, and a negative circuit can only connect 334 to 365. Conditioning on either direction therefore leaves a single returning word; its joint record with the reference factors exactly. A material circuit alone is insufficient to preserve the record.

## Can the current persist?

At the **native rates**, every stationary ensemble on this shuttered triangle has **zero material current on all three sides**, even if an experimenter assigns arbitrary static scalar potentials to all three physical ports. This extends the finite wait tests to every energy and arbitrarily long stationary observations on this aperture. For the native degree potentials the four energies are $8-4\sqrt2,4,8+4\sqrt2,21$. [Exact stationary threshold and characteristic polynomial](StaticThreshold.lean)

To see the threshold, subtract an energy $E$ from the Hamiltonian. Write the resulting diagonal in chain order as $(d,e,f,d)$. If the two hub amplitudes are $x,y$, the rim amplitudes must be $(dx-gy)/t$ and $(dy-gx)/t$. Two independent modes at the same energy therefore require independent hub pairs. Substitution into both rim equations gives

$$ed=fd=t(t-g),\qquad eg=fg=-td,\qquad d^2=g(g-t).$$

For $0<g<t$ this is impossible. At $g=0$ the irreducible real chain also has simple spectrum, since the amplitude at its first site recursively determines all the others. Conversely, for $g\ge t$ one may set $d=\sqrt{g(g-t)}$, $e=f=-td/g$, and the two free hub amplitudes give the degeneracy. The native ratio is $g/t=1/8$. A real symmetric Hamiltonian with simple levels admits real stationary modes. Their imaginary bond currents vanish. The upstream spectator copies do not change the measured current because each material instrument acts identically and leaves the spectator word unchanged.

As a sensitivity comparison, keep $t=4$ and the actual ports, and **tune** $g=(20+16\sqrt2)/7=6.089631$ to double the line at $E_*=8+4\sqrt2$. For each of the other two energies $E_j$, a binary controlled-wait meter has pass amplitude $(I+e^{-i\pi(H-E_*)/|E_j-E_*|})/2$. Passing both meters selects the complete degenerate line. A separate copy at each side equalizes the analyzer's site voltages, takes a quarter swap and counts the two addresses. Four times destination-minus-source probability reads the corresponding material current. The same counts can be made after any later free wait. On the native energy lines the two addresses have equal probability on each of these detector copies.

At the tuned rate **0.468486%** of original late launches pass the two energy filters with the source reference unread. The sustained spatial current, in native transfer units conditioned on the pass, is **[0.411560, 0.411566]** on each side. The two source reports separately give **[0.924696, 0.924706]** and **[−0.215816, −0.215805]**. Multiplying by the measured charge sign gives signed current. At an externally supplied rim phase of magnitude 1/64, the two adjacent energy lines split by $2|\phi|\times1.971187$. The stationary current extrema have magnitude 1.971197. The unordered spectrum is the same for both signs of the applied phase. This check demonstrates the sensitivity of the apparatus and the consequence of supplying the degeneracy. An emitted magnetic field, angular momentum in physical units, $g\simeq2$, and a wound spin frame require further premises. A many-face excitation can have a different spectrum.

## Verification and scope

Each Fourier mode is a four-site Hermitian matrix. Summing the counter amplitudes recovers evolution with no meter shift. A separate local recurrence propagates both input words with 32 Taylor substeps of order 24 on the 65-bin register. At the three counted rates $\|\,|H-12I|\,\|_2\le12$ by its maximum absolute row sum. Each substep has $z=12u/32<1.18$ and Taylor tail $z^{25}/25!$. Allowing 16 binary64 roundoffs per term step and 25 additions gives a per-substep vector error below $(25\times16+25)2^{-53}e^{1.18}+1.18^{25}/25!$. The 32 substeps, the final global phase, two input columns and the measured Fourier difference below $2\times10^{-13}$ fit inside the $3\times10^{-10}$ transfer-norm allowance. The Hermitian integral remainder bounded by $(12u)^{192}/192!$ applies to both lifted and finite evolutions with matching short paths. The earlier incoming model propagation has its own certified moment errors. For an initial Pauli error $\epsilon$, orthogonality gives joint trace-norm error at most $4\epsilon$. The output filter and ratio guards carry this through to the reported signs. [Source precision bounds](../README.md#precision-and-reproduction) · [Calculation](measure.py)

Run from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.electron_field.fan_path.local_launch.material_loop
.venv/bin/python -B -W error -m "$pkg.measure"
.venv/bin/python -B -W error -m "$pkg.witness"
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/electron_field/fan_path/local_launch/material_loop/StaticThreshold.lean)
~~~

With the phase convention $H_{1053,1049}=-t e^{i\phi}$, positive spatial flow from 1049 to 1053 equals $-\partial H/\partial\phi$ at zero phase. The energy-splitting comparison uses its magnitude. The local circuit and correlated return bear on the material field and angular response in OpenWave rows 5–7. A native stationary magnetic source, a spin-½ response and fermionic exchange of the identified carrier remain open under the published tests. [OpenWave criteria at the checked revision](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#simplest-test-per-criterion)
