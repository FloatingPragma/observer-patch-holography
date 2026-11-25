# A transverse coincidence area on the monitored two-load source

The four rim stations can measure how far apart two resident arrivals land on the same triangular detector face. On a deposited stock that the original return clocks scarcely distinguish from the reference, this area increases by **2.38714 squared clock units**. The same apparatus fed by two independent copies gives a decrease of **5.95469**. The joint image therefore resolves a source change that the twelve individual clocks miss. At fixed survey volume, its ordinary calibration fails on a held change of the same source. There is also an exact positive stock on which all the timing and image data used to calculate this area stay fixed while the first-order interior entropy changes. [Native-wave bounds](certificate.json) · [Acquired model counts](counts.json)

## What each launch measures

Use the existing 80-cell board, the same adjacent cells 0 and 3, four cut ports 3, 10, 36 and 44, and both resident double loads. The two loads occupy one of 276 unordered pairs of 24 ports. As in the entropy experiment, park the exterior contacts and monitor the complete sharp word every half-unit. On the resulting record the native transition probability is $P_{xy}=|\langle y|\exp(-iH/2)|x\rangle|^2$. The reference stock is uniform on these words. The instrument uses the same initial word and one ensuing two-load trajectory for a coincidence. The three pixels on the receiving face are existing ports. The registered material clock requests a joint exterior read after a positive geometric number of half-units with mean 4096, obtainable with twelve fresh fair bits per interval. This appointment and each initial preparation are supplied experimental operations. [Original clock protocol](../README.md) · [Joint protocol](measure.py)

Mount at one apex $o$ and start if it is occupied. At the appointment read the three opposite ports together. At each of four apexes the archive counts all eight outcomes: three possible distinct coincident pixel pairs, three exclusive single pixels, no opposite pixel, and failed start. For an opposite pair $a,b$, its probability per preparation is

$$q_{ab}=\sum_{x:o\in x}p_x\left[\frac{P}{4096}\left(I-\frac{4095}{4096}P\right)^{-1}\right]_{x,\{a,b\}}.$$

The control uses **two fresh independent copies**. Both must pass their respective apex checks. At each independently scheduled end, select one of the two occupied ports by a fair coin. Count distinct opposite pixels across copies. If $h_o$ is the apex occupation and $v_a$ is the selected-port probability conditional on launch, the all-attempt pair probability is $2h_o^2v_av_b$. The control stores three pairs, an eligible miss and a failed double launch. At zero flight the same-word coincidence is impossible because the occupied apex leaves only one other resident. The independent control can already produce pairs.

## Why this quantity has area units

Determine the detector triangle's three sides from the six directed clocks for its three edges. Their symmetrized sides fix the triangle area $A_f$ by Heron's formula. Let $y_a$ be its two-dimensional vertex coordinates in any orthonormal detector frame. A joint pair contributes the symmetrically oriented separation vector $\pm(y_a-y_b)$. A miss contributes zero. Write $e$ for an unordered pair and $d_e$ for its separation. The matrix of second moments $M=\sum_e q_e d_e d_e^T$ gives the registered face area

$$\frac{\sqrt3}{2}\sqrt{\det M}=A_f\sqrt{3\sum_{e<e'}q_e q_{e'}}.$$

Indeed each squared two-edge determinant equals $4A_f^2$. The determinant formula follows by Cauchy–Binet. A uniform conditional distribution on the three pairs returns $A_f$. Any affine change of detector coordinates multiplies the image area by the absolute determinant of the change. An image confined to one edge has zero area. Dividing by $\sum q_e$ yields the area conditional on a double hit. Keeping the unconditional $q_e$ includes throughput. We sum the four material face areas. The primary **live** estimate remeasures the three sides of each face on each stock. A separate control freezes the four reference face areas. These are areas of measured coincidence distributions in individually ranged frames. Identifying their sum with the proper surface area of a geodesic ball would require common physical radial distances and a measured angular coverage. The [five-station checks](../README.md) already test the common flat embedding on this setup.

## Ordinary calibration and held transfers

Use the parent experiment's ordinary transfer from zero to one occupied rim port. Independently measured clock ranges set the tetrahedral survey volume. Its compensating transfers move weight between the equal-energy, equal-rim-load sharp words $(3,36)$ and $(3,10)$. Finite arms use their deposited compensation settings. Set $\eta=-\Delta H(I)/\Delta A$ only on the ordinary transfer. The live conditional and throughput population gains are respectively 0.00053475 and 0.00071318 nats per squared clock unit. Independent model cohorts give conditional gain $[0.0005136,0.0005579]$. The table gives held $\Delta H(I)+\eta\Delta A$ per unit source transfer with simultaneous 99% coverage. [Source derivatives](results.json) · [All-attempt acquisition](acquire.py)

| Held stock change at matched survey volume | Conditional image | Throughput image |
| --- | ---: | ---: |
| Zero to two occupied rim ports | [−21.141, −19.495] | [1.6736, 1.6783] |
| Fixed mean native energy and mean rim occupation | [−9.3693, −8.4816] | [0.7351, 0.7368] |

On the parent's deposited 50-bit clock-insensitive stock, the exact live conditional area change is in $[2.38713687,2.38713797]$. The independent-copy image changes by $[-5.95468620,-5.95468607]$. The corresponding count intervals are $[2.1793,2.5945]$ and $[-6.2381,-5.6748]$. This is the actual dyadic stock whose twelve old success probabilities differ from reference by less than $7.48\times10^{-10}$. The seeded simulation samples $2^{60}$ independent starts per setting in multinomial batches. Coverage uses the joint union bound over at most 4,000 observed binomial marginals and event sums, including the bulk word readout. The bands include the certified native probability errors. These are synthetic records on this finite apparatus. [Raw counts and stock tickets](counts.json)

## Which records determine the entropy

Form 36 rows on the 276 word probabilities: normalization, ordinary native energy divided by 256, nine independent initial clock heralds, twelve complete return successes, twelve same-trajectory pair pixels multiplied by 256, and the indicator of two occupied rim ports. The exact native matrix $L$ has a verified dyadic proposed right inverse $B$ with $\|LB-I\|_\infty<0.00002541$. Vary $p=p_*+tB(LB)^{-1}e_{36}$. Both signs stay positive for $|t|<0.00002450$. This preserves **all twelve three-bin clock laws and all twelve pair pixels**, their associated survey volume and every registered same-trajectory image area, including both normalizations. The sum of old initial herald rates is $3P_1+4P_2$, where $P_j$ is the probability of $j$ occupied rim ports. Hence $dP_2/dt=1$, $dP_1/dt=-4/3$, and the uniform-stock interior entropy changes at $(4\log 4-3\log 6)/3>0.05663$ nats per unit $t$. This excludes first-order balance for every differentiable area determined solely by these fixed rates. New single-apex and exclusive-single counts can vary on this curve. [Integer proof](certificate.py)

At $t=-1/65536$ the deposited $2^{50}$-ticket source has old clock-success differences below $1.84\times10^{-14}$ and pair-pixel differences below $2.07\times10^{-16}$. A correction below $3.452\times10^{-9}$ in every word probability puts this ticket stock on the exact constant-record level set. Every corrected probability exceeds 0.001. Its exact entropy loss exceeds **0.00468837 nats**, using interval logarithms and the entropy gradient bound on the joining segment. The acquisition archives the actual ticket weights and measures their entropy loss in $[0.00469417,0.00469591]$ nats. [Certified bounds](certificate.json)

An independent exterior read does capture the completion information. On any stock the four single-apex rates give $\mu=E[k]$. With the original herald sum $h$, $P_2=(3\mu-h)/2$ and $P_1=h-2\mu$. The boundary capacity $K=P_1\log4+P_2\log6$ therefore has its own all-trial measurement. The multiplicities $1,4,6$ of each interior word imply

$$H(p_I)+K(p)=\log276-D(p_I\Vert(p_*)_I).$$

This gives stationarity of the sum at uniform stock. It also explains a simple finite control that preserves **every bin** of all the registered clock, joint and independent-copy experiments. The 190 strictly interior pairs all fail every initial herald. Sort these words by native diagonal energy and pair them within each class. Their five class sizes are 30, 10, 1, 123 and 26, giving 94 disjoint equal-energy pairs. Toss one fair coin to redirect half a donor's uniform weight to its partner. The smallest final word probability is $1/552$. Mean energy, $\mu$, $K$, every clock and every screen law stay fixed. The interior entropy decreases by exactly

$$\frac{94}{276}\left[\tfrac32\log(3/2)+\tfrac12\log(1/2)\right]=0.0891038505686\ldots\ {\rm nats}.$$

This control starts at a stationary point of $H+K$ and its entropy change is second order there. The separate 36-row construction gives nonzero **first-order** entropy response while preserving the clock and coincidence records used in the area estimator. The full-bin control's independently counted finite loss is $[0.0891029,0.0891048]$ nats.

## Verification

The propagator certificate removes a scalar phase, computes the 29th-order matrix exponential Taylor polynomial at $1/128$ scale, and squares seven times. Exact integer dot products on 21-bit limbs enforce the $2^{80}$ fixed-point accuracy at each round. The resulting native transition matrix has maximum row $\ell^1$ error below $1.167\times10^{-19}$. Integer residuals for the discounted return and terminal equations and the maximum resolvent norm 4096 give probability errors below $9.56\times10^{-16}$ and $4.78\times10^{-16}$ respectively. Fresh count simulation uses a faster floating kernel separately checked against these bounds. [Wave calculation](wave.py) · [Regression checks](test_joint.py)

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.area_survey.joint_screen
python -m "$pkg.measure"
python -m "$pkg.certificate"
python -m "$pkg.acquire"
python -m pytest -q model/experiments/interval_echo_einstein_station/equilibrium_refinement/area_survey/joint_screen
```
