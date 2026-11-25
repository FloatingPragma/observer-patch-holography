# Fluctuations seen by the four-port observer

The observer records a four-bit word before, between and after two conditional repairs. The observed entropy production can have either sign. In the forward experiment 16,271 of 65,536 simulated histories had negative production (24.83%). Separate calibrations of the reference and endpoint populations resolve the sign of all 16,271. The predicted fraction is exactly one quarter. Over all starts the mean production was positive, 0.20931 nats per history. The independently derived value is 0.20678 nats. [Frozen protocol](FLUCTUATION_PROTOCOL.md) · [All eight cohorts](fluctuation.json)

## The quantity the observer can measure

Use the existing eight-port path with four material units. Its calibrated distribution on the four visible bits is R. Prepare another ensemble P with odd parity on the first two bits at frequency one quarter, using the existing fixed-time Gibbs source. Refresh port 0 and then port 2. Each refresh holds the other three visible values and mixes the selected port with the same four hidden ports. A new ensemble prepared this way supplies an independent calibration of the endpoint law Q. In the forward cohort every triple (x,m,y) is kept. In the reverse cohort independently generate another Q and refresh port 2 before port 0. The comparison is to the reverse triple (y,m,x). The address clock and sharp occupation instruments are supplied. All measurements concern fixed public appointments.

For a given trajectory define

\[
\Sigma(x,m,y)=\log\frac{P(x)}{R(x)}-\log\frac{Q(y)}{R(y)}.
\]

The reference is material. Its values count inaccessible arrangements of the four hidden ports. R, P, and Q can each be read out before comparing any forward and reverse trajectory counts. A decrease in stochastic entropy production here means a negative value of this dimensionless observable. The positive ensemble mean measures the loss of relative entropy to that bath. A different experiment would be needed to express it as heat in energy units. [Reference and Gibbs preparation](FIXED_TIME.md) · [Repository definition and finite proof](../../../../Lean/Thermodynamics/FluctuationTheorems.lean)

## Why reversal predicts an independent measurement

For refresh address j, write u for the other three visible bits. The conserved material and the uniform five-site shuffle give the full visible kernel

\[
 K_j(x,y)=\mathbf{1}_{u(x)=u(y)}R(y)/R(u(x)).
\]

Consequently R(x)K_j(x,y)=R(y)K_j(y,x) for all words. Each of the four addresses has an explicit neighbouring-swap route into and out of the same bath. Exact arithmetic and native gate columns check the induced visible kernel on every one of the 70 microscopic inputs. In particular, hidden preparation within an allowed four-port read cannot change its subsequent visible transition probabilities. Publishing the visible word after each appointment therefore preserves this description for the next appointment.

Multiply the identities along any finite fixed address sequence. Every internal R factor cancels. For any positive visible initial law P, let Q be its actual final distribution. The resulting forward path weight is its reverse weight under Q and the reversed address schedule multiplied by exp(Σ). Self-hold has positive probability, so Q stays positive, and the paired path supports coincide. This proves the pointwise and level-set odds relations on the complete visible history space. Summing gives the integral identity E[exp(−Σ)]=1. Taking the mean of the endpoint formula gives D(P‖R)−D(Q‖R), which is nonnegative by stationarity and data processing. Each refresh is reversible on the declared records. Its hidden sequence of measured microscopic contacts and discarded controller actions is a distinct record space.

## Counts and controls

There are 64 feasible three-word trajectories in the two-step experiment. The exact test also evaluates every one of the 4,032 impossible triples. Each cohort uses 65,536 fresh simulated launches. All three marginal calibration histograms and all five trajectory cohorts are included in the receipt. The histogram indices use base 16 in chronological order with zeros implicit.

| Quantity | Derived | Measured |
| --- | ---: | ---: |
| Fraction of forward paths with negative production | 0.25 | 0.24828 |
| Mean exp(−Σ) | 1 | 0.99622 |
| Mean Σ (nats) | 0.20678 | 0.20931 |
| Forward to reverse odds, least multiplier | 25/64 | 3,781 / 10,193 |
| Forward to reverse odds, greatest multiplier | 100/49 | 10,523 / 5,149 |

The 64 direct count-ratio intervals overlap the 64 independent endpoint/reference odds intervals. All eight distinct level-set multipliers agree with their independent forward and reverse frequencies under the registered allowances. The equilibrium source gives zero production on each supported path. Reversing the order incorrectly puts 46,164 of 65,536 histories outside the correct reverse support. Disabling material exchange on the two public repairs leaves every word unchanged. Replacing the correlated reference with a flat visible distribution gives a signed aggregate discrepancy of 0.12055, outside the simultaneous allowance 0.02344. The elementary counting model and the material simulator share the carrier and declared Born readout; the visible forecasts use inventory counts, while the launch simulation samples each scheduled contact separately. [Producer](fluctuation.py) · [Exact kernels and control tests](test_fluctuation.py)

The simultaneous failure allowance is 1% over the two specified families. One half goes to 176 Clopper–Pearson binomial intervals on marginal and path probabilities (Bonferroni). The other half goes to at most 1,024 bounded-mean comparisons (Hoeffding). For two independent frequencies both single-frequency errors are propagated. The tables include ordinary sample fluctuations.

Replay from the RER root.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_consensus.conditional_bath.fluctuation
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/conditional_bath/test_fluctuation.py
~~~
