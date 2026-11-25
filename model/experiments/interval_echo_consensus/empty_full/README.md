# Empty–full interference after one seam repair

**On registered Interval Echo, the empty–full multiplier is exactly one for a fresh seam.** Prepare two loads in a coherent superposition of both being inside or both being outside the target seam. After one repair, an observer can recover the original bright port with certainty or extinguish it by changing a local phase by π. A reversible change to the stored record preserves every sharp load probability and changes the bright probability to one half at both phases. In the historical port arrangement, the native bright counts were 1,048,576 and 0 out of 1,048,576 simulated starts. The changed record gave 524,738 and 524,671. [Native instrument and control](measure.py) · [Literal twelve-port acquisition](historical.py)

The status board asks for a second preparation that reads the multiplier k and distinguishes the parity record E from the seam-total record T. The original total-two experiment with an empty and a single seam measures a different multiplier, g. E and T both predict a phase-independent 3/8 in that experiment. The archived paper program also gives the ideal answer for empty versus full. The construction here connects that conditional formula to a prepared source, held count banks and a finite detector under the active law. The externally measured choice of record law remains open. The measured outputs here are material occupation patterns following controlled recombination. A free spatial intensity pattern behind a double slit is a separate experiment. [Status board](https://status.floatingpragma.ai/) · [Original record-law derivation and fringes](https://github.com/FloatingPragma/observer-patch-holography/blob/f2f264de2743d4c5d0194ff93bd792ec6e7708e2/model/record_laws/README.md#theorem-3-what-observations-read) · [Protocol](PROTOCOL.md)

## Preparing the historical source

Use target seam (0,1) and keep total load two throughout. Call the sharp state with one load at each of ports 2 and 3 X. Call the sharp state with one at each of 0 and 1 W. The input and the measured bright ray are both (X + W)/√2. Other seam contacts are shuttered. The apparatus supplies pulse timing, local number phases, per-contact addressing, fresh memories and terminal load counts. Every preparation enters the denominator.

The short source uses the two native exchange rails (0,5) and (1,7), with separated sharp unit loads on 5 and 7. On a single-occupancy rail, a π/2 hold makes a balanced coherent exchange. A local number phase removes the known relative factor i. On the target seam, insert a π number phase on port 0 between two native half-holds and undo it after the second half. In the odd sector this reverses the off-diagonal exchange and retains the diagonal energy κ+1/2. In the equal-occupation configurations 00 and 11 both holds have zero energy. Setting the total duration to π/[2(κ+1/2)] supplies phase −i on the odd terms. One more balanced exchange on (1,7) cancels those terms. The retained state has amplitude 1/√2 on each of the two desired words, up to one common phase. No measurement participates in this preparation. At κ=8 every inverse contact in the reciprocal bright reader can use a positive native hold complementary to 2π. [Local circuit](measure.py)

For the literal X at ports 2 and 3, conjugate the same circuit by native π swaps along 5–4–2 and 7–6–3. At integer contrast a π hold swaps the endpoint loads without an address-dependent phase. The routing leaves the target untouched. It is executed on all 78 words of the twelve-port total-two cell; no terminal word is dropped. After the fresh write apply the indicated number phase on port 0, undo preparation and count returns to the initial word. The bright event is the jointly recorded terminal word. [Routing and full-word census](historical.py) · [All original-port counts](historical.json)

## Which record preserves the phase?

For source endpoint loads (a,b) the registered write attaches the coin c and the contrast d=a−b. If p=(a+b) mod 2 and j is the loss from the first endpoint, the output loads give d=2j+(2c−1)p. Conversely p=d mod 2 and j is fixed by c and d. Thus recording (c,d) has exactly the same record coincidences as recording (c,j,p), the published E instrument. In particular 00 and 11 attach the same normalized coin memory and preserve their off-diagonal load amplitude. The original T instrument records the full seam total in place of p and makes these memories orthogonal. Reading after recombination gives

$$P_{\rm bright}(\theta)=\tfrac12+\tfrac12\operatorname{Re}(k e^{-i\theta}).$$

The changed-write control uses the same 99-dimensional local space. Keep the native columns except for input 11. Its output stays 11, with normalized memory (|0,2⟩−|1,2⟩−|0,−2⟩+|1,−2⟩)/2, written in the (coin,contrast) basis. This vector is invariant under simultaneous coin and orientation reversal. It is orthogonal to the empty-seam memory and to both symmetric coin records for incoming collisions. Completing the isometry with the native echo formula gives a Hermitian unitary and preserves every sharp endpoint load distribution. On the original total-two sector its **unread load channel** equals T. The private memory alphabet and basis measurement need not reproduce T's named fine record.

Using integer coefficient comparisons on every one of the 6,084 ordered pairs of cell words gives 6,084 agreements of the registered map with E and 6,084 of the control with T. The two reduced maps differ on 110 units: the unique full word paired in either order with one of 55 empty-seam words. They agree on every sharp input. The label identity for E extends to all nine local qutrit keys; the T control was constructed for the total-two sector. A bright fringe by itself determines k and admits other record laws sharing that value. The matrix-unit comparison establishes the stronger channel identity. It requires a fresh input memory. For a sharp single load, two writes with a fresh memory at each step return its initial position with probability 1/2. Reusing the same unread memory gives one by the involution of the echo. [Channel and memory checks](measure.py) · [Independent circuit and trace tests](test_fringe.py)

## Counts and held predictions

Each entry gives the bright returns among 1,048,576 original starts at the historical addresses. The single-arm cases physically retain only the selected branch; the archive includes all rejected starts.

| Seam event or intervention | Phase 0 | Phase π/2 | Phase π | Phase π/3 |
| --- | ---: | ---: | ---: | ---: |
| Bypass | 1,048,576 | 525,152 | 0 | 786,108 |
| Native echo | 1,048,576 | 523,628 | 0 | 785,310 |
| Total-sensitive echo | 524,738 | 524,987 | 524,671 | 523,129 |
| Resolve the branch | 524,169 | 524,149 | 524,999 | 524,391 |
| Pass empty only | 262,381 | 261,992 | 261,657 | 261,969 |
| Pass full only | 262,376 | 262,047 | 261,734 | 261,645 |

The commissioning bank at ports 5 and 7 uses phases 0 and π/2 to estimate two quadratures. Their sine and cosine forecast the four held phases without refitting, including phases of opposite signs. The same estimates forecast all six phases at the literal historical addresses. For example the native π/3 forecast was 0.750314. Its held twelve-port rate was 0.748930. All 24 original held returns and all 36 cross-port returns agree with their joint allowances. The 99% family radius across every outcome bin and rejection of both campaigns is 0.002527 per proportion, rounded outward; each forecast also propagates both training errors. Synthetic shot noise and probability rounding (13 decimal places before normalization) are distinct from empirical measurement of a physical seam. [Commissioning counts and forecasts](results.json) · [Independent held bank](historical.json)

From the repository root, use its NumPy and SciPy environment.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -B -W error -m model.experiments.interval_echo_consensus.empty_full.measure
python -B -W error -m model.experiments.interval_echo_consensus.empty_full.historical
python -B -W error -m pytest -q model/experiments/interval_echo_consensus/empty_full
```
