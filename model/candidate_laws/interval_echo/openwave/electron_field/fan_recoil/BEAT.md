# Is the filtered beat a common rest clock?

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

Use the vacancy on the three- and four-fan prefixes of the same 320-cell path. Keep the two lines passed by the coarse parked-reference filter in the [path experiment](../fan_path/PROTOCOL.md). Prepare the six Pauli eigenstates in the archived orbital code with a fair recorded source draw and a vacuum in each spectator. At the terminal hub, require the charge, two-dimensional code and every spectator vacuum. On that joint acceptance, independently read X, Y and Z, retaining both signs and all refusals. The rare-event budget is set by a numerical sensitivity check on the archived code. The agreement test pools every prepared input and both receiver signs. No setting is selected by its observed contrast.

Use generator dwell settings 0 and 3 to acquire the beat ratio, with held reads at 1 and 2, and retain the two archived appointments 27.0625 and 30.8125 as additional observations. Apply the same schedules at both lengths. For each route compare a coherent coarse filter and a fine-energy dephased control, which retains both lines and erases their relative phase. Use the same detector without any rephasing across lengths. The neutral clock is a separately read outer-ring coherence on the same dwell channel. Ring clocks are acquired on independent trials at each route. Use 2^28 trials per route, phase treatment, dwell and receiver axis, and the same shot budget per clock quadrature. Report simultaneous 99% intervals for the compared records, with a conservative union bound over 64 marginal confidence intervals. Treat a pooled gap ratio disagreement as resolved only if the two intervals are disjoint.

## Frequency from accepted clicks

The detector projects onto actual superpositions of the two local head-word parities at the receiving hub. With local code vectors $d_1,d_2$, the X and Y outcomes project onto $(d_1\pm d_2)/\sqrt2$ and $(d_1\pm i d_2)/\sqrt2$. The Z outcomes are the two parity labels. The same code defines the source polarities, so these preparations and projections have an orbital-word implementation in the model. They are supplied coherent operations in the simulator. The ideal spectral filter and the fine-energy dephasing are likewise supplied controls.

Let $T_j$ be the accepted transfer amplitude associated with line $j$, with the source and receiver containing their corresponding spectator vacua. The probabilities on the fair source mixture satisfy
$$
p_{\rm coherent}(t)-p_{\rm erased}(t)
=\operatorname{tr}(T_0^\mathsf{T}T_1)\cos[(E_1-E_0)t] .
$$
The transfer matrices are real before their dwell phases. Erasing the fine label removes the cross term and leaves the filter acceptance unchanged. Summing the two receiver signs for each axis gives the same full accepted projector. A fair recorded source draw therefore gives an ordinary Bernoulli joint click at each setting when the source label and port sign are pooled.

For nonzero measured contrast at dwell zero, the contrast ratio at dwell three is $\cos[3(E_1-E_0)]$. There is no alias in this inverse. The parked threshold is one neutral-ring unit above the parked charged ground, and the released generator adds a positive bridge Laplacian. The gap between any two accepted energies is consequently below one in the generator's units. The angle for three dwell units lies in $[0,3]$. Exact inertia at both ends of the inherited threshold bracket certifies that two lines pass. Eigensolver output is used to generate detector probabilities and validate the forward model. The readout uses the count ratio and an independent ring clock. It receives no energy values.

The confidence intervals use the binomial law of the marginal counts and propagate through the signed contrast ratio and the independent clock phase. A guard from the dyadic eigenvector and eigenvalue enclosures pads the modeled probabilities. The phase readings at dwells one and two test the prediction from the zero and three dwell counts. The remaining scheduled observations and all terminal refusal categories are kept in the receipt.

## Clock comparison

| Joint read | Three fans | Four fans |
| --- | ---: | ---: |
| Beat in independent neutral-ring clock units | 0.58551 | 0.36144 |
| Simultaneous 99% interval | [0.57611, 0.59532] | [0.34370, 0.37785] |
| Coherent accepted counts, dwell 0 → 3 | 51,945 → 189,429 | 19,469 → 49,702 |
| Phase-erased accepted counts, dwell 0 → 3 | 168,522 → 167,458 | 76,091 → 76,366 |

Every count in the table uses 805,306,368 simulated launches pooling the three receiver axes. The difference of clock ratios has simultaneous 99% interval **[0.19827, 0.25162]**. The prediction intervals at held dwells one and two overlap their independently acquired intervals on both lengths. All phase-erased dwell readings agree within their joint uncertainty with the constant reference acquired at zero and three. The coherent terminal Y sign difference at held dwell one is unresolved with the fixed budget for the last odd input. Its joint probability difference has intervals $[-1.21,0.95]\times10^{-5}$ on three fans and $[-0.73,0.55]\times10^{-5}$ on four. The count-only gap uses the total same joint click across both port signs. [Full records](beat.json)

The two filtered lines produce a readout beat on both paths. They give different frequencies relative to the same independent ring standard. Their relative phase survives the specified coherent filter and joint read. A common scalar shift of the two admitted energies only multiplies the accepted amplitude by a common phase and leaves all these records invariant, provided the same two lines pass the filter. This detector probes the splitting. An absolute rest-energy inference would require a neutral or compensator reference on the same extended preparation, plus a physical calibration. The finite path's two-line splitting fails as a length-independent rest clock on the tested fixtures.

Run from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_recoil.beat
~~~
