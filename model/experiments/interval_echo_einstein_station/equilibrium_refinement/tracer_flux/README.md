# Reading a joint rim correlation from an exterior receiver

A separately launched signal detects a change that every individual rim occupancy counter misses. Two exposure settings, together with a reversible randomized preparation step, allow a single exterior receiver to recover the full distribution of rim occupancy on this finite source. The distinction between transmission and proper area can also be tested. A single transmission calibrated on the ordinary transfer predicts almost twice the independently counted change on the held transfer. [Full simulated counts and responses](results.json) · [Native dynamics and acquisition](measure.py)

## What the observers do

Use the [five-station rim and enrolled sources](../CLOCK_BALL.md) of the scheduled round-trip clock. Its two-hop ball has six inner and five rim addresses. A twelfth enrolled station is outside the timed ball and inside the parked one-cell apparatus of the 20-cell host. The two residents each carry load two. The independent tag carries load one. All three move during transmission. The schedule supplies the shutter. Two separately prepared material-clock cohorts read each exposure using the native dimer hand with two quadrature settings. Their reconstructed phases are 0.25011, 0.50000 and 1.00001 radians. This calibration refers to the supplied isolated hand.

* **Inward.** Start the tag at the outside receiver. At the shutter classify the receiver as empty, holding a resident, or holding the tag. The first event is the inward hit, with probability $F$. The terminal can be read against load-zero and load-one standards, using the local echo's agreement flag and undo. Routed endpoint swaps bring the terminal to the reader. These two commuting flags distinguish all three alternatives. Probabilities sum over every history giving the terminal flag.
* **Outward.** On a fresh copy with the same initial resident preparation, launch the tag uniformly from a vacant slot inside the ball. Sampling all eleven sites and heralding vacancy admits nine in every word. Score every admitted launch, recording a hit if the tag reaches the outside receiver. Call the hit probability $R$. The outward acquisition has its own simulated counts.

The generator has 660 sharp words for one tag and two residents. Each contact uses the nine native blank-memory columns, including exchanges of a resident with the tag, and the same host diagonal as the regional experiment. All contacts to the host outside the enrolled region are parked as in that experiment. The basis labels enumerate the probability of a receiver report. An observer only needs to read the terminal and the clock. For the uniform resident reference, the inward and reverse expected hits satisfy $F=9R$. Indeed the real symmetric generator gives a symmetric transition-probability matrix. Summing over the 55 possible outside-tag starts and the 495 possible inside-tag starts transposes the same block. This reciprocity is a check independent of the entropy law.

## Prediction after ordinary calibration

The registered ordinary transfer moves probability from word $(0,1)$ to $(0,2)$, sending one resident from the interior to the rim. The held transfer changes an equal mixture of $(0,2),(1,9)$ to an equal mixture of $(0,1),(2,9)$. The latter two stocks agree on all single-site marginals, the incoming native energy distribution, and the marker volume. At the uniform monitored reference the interior entropy differential is minus the differential of the mean exterior completion score $C=\langle\log\binom{5}{k}\rangle$, where $k$ counts rim residents. The ordinary change in $C$ per unit transfer is $\log5$. The held change is $\tfrac12\log(2/5)=-0.458145$ nats. [Independent interior and rim acquisitions](../clock_ball.json)

| Half-radian shutter | Ordinary change | Held change | Held prediction from the ordinary gain (nats) | Simultaneous 99% interval from synthetic counts |
| --- | ---: | ---: | ---: | ---: |
| Inward hit probability $F$ | −0.020730 | +0.011462 | −0.88993 | [−0.9092, −0.8771] |
| Scaled outward probability $9R$ | −0.019835 | +0.010206 | −0.82810 | [−0.8739, −0.7764] |

The ratio is $\log 5$ times held over ordinary response. The same ratio follows using opacity $1-F$. The intervals combine eight binomial Clopper–Pearson intervals with a Bonferroni correction, from $2^{26}$ shots per bank. With the contrast term switched off and the alphabet, incidence, shutter and stocks fixed, the inward prediction becomes +0.0780 nats. Thus the held signal depends on the dynamics. Among two-rim words with the identical exterior score $\log 10$, the active inward probability ranges from 0.000597 to 0.015840.

## Two exposures recover the exterior completion score

This is a different, active acquisition. Before the tag arrives, independently permute the six inner sites and five rim sites uniformly, preserving the unknown $k$ on each run. Each permutation is a product of routed endpoint exchanges. All 25 within-group transpositions restore their intervening addresses, and the native sharp contact exchanges their loads. The randomized operation can change energy and erases location information inside each class. The held unshuffled comparison above uses its original stocks. A permutation on each fresh copy of the unknown source is needed for the two-exposure result.

Let $p_k$ be the unknown probabilities of $k=0,1,2$. Class $k$ contains respectively 15, 30 or 10 words. For each class the independent generator or an enrolled pure-class panel supplies the inward transmission coefficient $f_k(t)$ after the shuffle. The deposited matrix uses the finite generator. It includes a two-resident calibration of the instrument's transmission. No entropy or held contrasts enter this response matrix. With exposures $t$ and $2t$ the three recorded equations are normalization, $F(t)=\sum_k p_k f_k(t)$ and $F(2t)=\sum_k p_k f_k(2t)$. At $t=1/2$ their condition number is 596.1. A single signed linear decoder recovers $p$ and hence $C$. The two stocks' observed difference is **−0.45704 nats**, with simultaneous 99% interval **[−0.47226, −0.44182]** containing the independent rim measurement. A separate outward detector at the unused quarter-radian setting checks the inferred distributions:

| Held stock | Outward probability predicted from two inward readings | Separately observed outward frequency |
| --- | ---: | ---: |
| Split | 0.00257592 | 0.00257635 |
| Spread | 0.00161142 | 0.00161160 |

These are seeded, aggregated binomial acquisitions with $2^{34}$ trials per setting, calculated from the transition probabilities. No individual trajectories of that many trials were evolved. A familywise 99% Hoeffding bound covers the six frequencies with the response matrix fixed. The resulting limit on each outward discrepancy is 0.0000161. The unconstrained linear inversion retains small out-of-range statistical fluctuations. The charge estimate and bounds include them.

An independently counted calibration propagates uncertainty in the instrument itself. Three registered pure-class standards at each of the two inward settings and the withheld outward setting provide nine fresh binomial counts. With these counts as the calibration, the same held acquisitions yield **−0.45679 nats**, with simultaneous 99% interval **[−0.48912, −0.42445]**. The predicted split and spread outward rates are 0.00257613 and 0.00161174; both agree with the independently observed rates above. Here the common frequency error $\epsilon=0.00001527$ covers all fifteen banks. For any physical class distribution, its weights are nonnegative and sum to one, so the error introduced by one row of the calibration matrix is at most $\epsilon$. Including the held bank gives $2\epsilon$ per inward setting. For the fitted charge decoder $v$, the two-stock charge error is therefore at most $4\epsilon(|v_1|+|v_2|)$. The observed least singular value, 0.002907, exceeds the simultaneous calibration error bound 0.000038. These bounds describe the simulated protocol with a fixed shutter and source ensemble.

### Exact separation of the three classes

Write $D=2H$, an integer symmetric matrix. On disjoint launch and arrival events the coefficients of $t^2$ and $t^4$ in a transition probability are $D_{ij}^2/4$ and $[3(D^2)_{ij}^2-4D_{ij}(D^3)_{ij}]/192$. Summing integer entries at the exterior and averaging over each registered class gives

$$
a=(5/4,1,3/4),\qquad b=(-17915/288,-2629/48,-23359/192).
$$

Thus the determinant of the three response rows at times $t,2t$ has first nonzero term $12\det[\mathbf 1;a;b]t^6=(42811/192)t^6$. In fact the response is invertible at **every nonzero real algebraic exposure** in generator units. To see this, write each transition probability by the finite spectral decomposition of the rational symmetric matrix $H$. Its eigenvalues are algebraic and its distinct-eigenvalue projectors are polynomials in $H$ with algebraic coefficients. The determinant is consequently a finite sum $\sum_j c_j e^{i\omega_jt}$ with algebraic coefficients and distinct algebraic frequencies after grouping. The displayed nonzero Taylor coefficient ensures some $c_j\ne0$. For algebraic $t\ne0$, Lindemann–Weierstrass gives linear independence over the algebraic numbers of the exponentials of the distinct algebraic numbers $i\omega_jt$. This establishes injectivity of the ideal instrument. Its statistical conditioning is exposure-dependent. [Lindemann–Weierstrass theorem and formal proof](https://isa-afp.org/entries/Hermite_Lindemann.html)

The same reasoning isolates what one exposure can do. For the original held transfer the $t^2$ coefficient cancels and the $t^4$ coefficient is $581/24$, so its inward signal is nonzero at every nonzero algebraic dwell. But the two equal-charge words $(4,9)$ and $(2,9)$ have distinct $t^4$ coefficients $-175$ and $-6559/96$. No nonzero affine gain on this one inward probability represents the completion score for all preparations at any such dwell. The class reconstruction establishes that the correlation can be read through transmission. Identifying a physical area still requires observers to determine transverse proper area or neighboring proper volumes with independent rulers at a fixed proper volume, and to predict its held response at the ordinary coupling.

## Reproduce

From the repository root, use a single numerical thread:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.tracer_flux
python -m "$pkg.measure"
python -m pytest -q --pyargs "$pkg"
~~~
