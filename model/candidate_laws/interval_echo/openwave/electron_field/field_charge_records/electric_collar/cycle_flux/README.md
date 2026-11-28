# Testing the finite Coulomb field on the 1045 arrival

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

Move the signed arrival from 1045 to 1051 or 1052 and compare it with a copy held at 1045. The measured head changes have the correct discrete divergence for the measured dipole. They also have nonzero circulation after free face evolution and after randomizing the native face pulses. The finite Green theorem requires zero circulation. A separately programmed set of reachable head words gives exactly the Green field in the ensemble mean. [Protocol](PROTOCOL.md) · [Original-launch count banks](results.json) · [Whole-cell preparation](../../../dressed_exchange/intact_path/face_turn/FRAME.md)

### The premise observers can test

Use the 30.8125 arrival in both charge charts, with the original face pair and its supplied 21 additional internal heads. The full-cell arrival passes on about 1.58437% of all starts. An enabled material contact acts on its two words as $4(I-X)$, so exposure $\pi/8$ exchanges them. The two common hops preserve the ordered face amplitudes by flipping just the traversed head. Each comparison uses one arm with the source at 1045 and one with the source relocated. Their twelve material-port records give the neutral load difference $\rho=q(\delta_b-\delta_{1045})$. The head observer records whether the head names the larger endpoint of each seam. Its unit is one head changed by an isolated complete native material hop. With incidence matrix $B$ that credits the larger endpoint, the independently measured difference of head means $F$ obeys $BF=\rho$. Every reached sharp word has the required port counts.

The actual icosahedral seam graph maps isomorphically to the committed Lean twelve-port table. The executable checks the integer Green identity after relabeling. The theorem `thomson_unique_orthogonal` says that the measured neutral $ρ$ and a divergence-matched, cycle-orthogonal $F$ imply the unique finite field $F_G=B^T Gρ$. Here the twenty oriented triangular reads form a matrix $C$ with cycle rank 19. They span all source-free cycles, so measuring $CF=0$ supplies the exact remaining finite premise for this head-field interpretation. Exact positive pivots give a lower bound of $3/4$ on $CC^T$ on the space orthogonal to constants. Thus the same comparison with uncertainty gives

$$
   \frac{\|CF\|_{\infty}}{\sqrt{3}}\;\le\;\|F-F_G\|_2\;\le\;\frac{2}{\sqrt{3}}\|CF\|_2 .
$$

The lower bound uses three seams per triangle. The upper bound follows because $F-F_G$ is a cycle and $CF_G=0$. Both apply to real population means as well as rational histograms. [Finite Green and uniqueness theorems](../../../../../../../../Lean/Screen/DiscreteCoulombGreen.lean)

### Equal divergence and different circulation

Every exact sharp-word comparison has an integer on the moved seam, while $F_G$ needs $q\,11/30$ there in the outward orientation. The sharp comparison therefore fails for every pair of reachable words, at every wait. Means can be fractional. The component contains 6,610 fixed-source words and 5,258 moved-source words at either destination. Each random face command permutes the head records, and the idle weight makes the connected random walk converge to uniform. For both routes, the limiting distributions miss the Green prediction by about 0.01238747 on one seam and retain a maximum triangle circulation of exactly $155/5258$.

For the 1051 route and the positive load, these simulated banks give the following circulation intervals. Each bank has $2^{34}$ original launches including refusal. Both signs and both routes reject zero circulation at each enrolled native wait and after 4,096 random face commands. The uncertainty is simultaneous 99% over every seam and face comparison and includes the inherited source and propagation bounds. A triangle compares the single-seam reading with the sum along the other two seams, so its circulation is an observable difference between two path estimates of one potential drop.

| Exposure and field schedule | Triangle | Circulation interval |
| --- | --- | ---: |
| Free wait 1 | 1045–1051–1052 | $[0.72229,0.72417]$ |
| Free wait 4 | 1045–1051–1052 | $[0.33156,0.33345]$ |
| All faces shut, wait 4 | 1044–1045–1051 | $[-1.00094,-0.99906]$ |
| 4,096 random native face pulses | 1044–1045–1051 | $[0.02851,0.03039]$ |

The independently programmed passing control uses 21 reached words per route across its two arms. All the routes to those words use native face swaps with the measured material position parked. The rational probabilities are chosen to satisfy the Green target. The control passes every triangle count interval, and its exact ensemble means satisfy the theorem premise. The free-wait result covers the four enrolled waits. The all-time sharp obstruction and the exact randomization limit have separate proofs. These comparisons distinguish conserved material transfer from the canonical finite electric-field candidate. A source law that produces the requisite circulation-free field, an independent electromagnetic calibration and a physical range refinement are required for rows 1 and 20.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.electric_collar.cycle_flux.measure
~~~
