# Can a charged and a neutral source give a positive bound quadrupole?

On this Interval Echo fixture, observers can select a below-threshold preparation, read a connected three-site region with charge +1 and a disjoint four-site region with charge zero, and measure a **positive** axial material quadrupole in the same cohort. The two-charge signature occurs in **0.055990%** of bright sharp launches. Most launches instead separate the signs into oppositely charged regions. Two independent cells containing all four excesses and all three vacancies separately have an energy below every seven-defect level of the joined cell. Consequently the measurement establishes a conditional two-region precursor and gives an explicit competing dissociation channel. The OW-16 deuteron remains unestablished. [Prospective choices](PROTOCOL.md) · [Held energy aperture](APERTURE.md) · [Held regional read](FRAGMENTS.md)

The two nominal constituent stocks also admit cheaper breakup on identical isolated cells. The lowest (2,1) energy is 119.18591 against 117.32096 for (2,0) plus (0,1). The lowest (2,2) energy is 142.24260 against 135.05547 for (2,0) plus (0,2). These are eigenvalues of the specified cell generator with the same total inventory and reference energy in each comparison. A common chemical potential cancels. The comparison constrains binding energetics; kinetic trapping can still preserve a prepared stock during a finite flight. [Signed-sector operator and regression](test_signed_pair.py)

## The separation threshold depends on what the two observers prepare

Use row 14's unit-filled reference and its independent 2-excess/1-vacancy (positive) and 2-excess/2-vacancy (neutral) source inventories. An observer opens all 30 internal seams and nine surviving private references per cell. The external gateway contacts stay shuttered. The comparison is between two otherwise identical cells, either with the two source stocks separated, or with all seven defects on one cell and the second left filled. Material, stock size, energy origin and gateway placement match. Host subtraction removes a common additive energy. In exchange units the operator on each signed stock is

$$
H=\tfrac12\sum_{ij}(I-S_{ij})+
8\sum_{ij}(q_i-q_j)^2+
8\sum_{i\text{ private}}[(1+q_i)^2-1],\qquad q_i\in\{-1,0,1\}.
$$

Complete integer word enumeration gives the minimum diagonal potential. Each exchange contribution is positive semidefinite. A sharp word with minimum ordinary diagonal energy supplies the variational upper bound. Independent-cell energies add. The last column comes from separate cohorts of 16,777,216 starts selecting one of the 39 physical energy terms uniformly, reading its outcome and multiplying the mean by 39.

| Stock (excess, vacancy) | Native ground interval | Registered minimum-diagonal source mean | Contact estimate |
| --- | ---: | ---: | ---: |
| Positive (2,1) | [112, 119.5] | 119.5 | 119.573 |
| Neutral (2,2) | [136, 144] | 144 | 144.001 |
| Joined (4,3) | [200, 209] | 209 | See signed banks |
| Excesses alone (4,0) | upper 133 | 133 | 132.980 |
| Vacancies alone (0,3) | upper 52.5 | 52.5 | 52.552 |

The positive and neutral standards straight through the extended material pointer give 16,254,669 and 8,386,714 positive exits respectively, each out of 16,777,216 trials. Their expected calibrated signed responses are 0.938 and zero. The nominated separated sector has floor 248. The alternative upper bound 185.5 is strictly below the joined floor 200. This competitor also exists in the same positive threefold rotation sector. Rotating both independent cells together gives three orthogonal product words. Direct exchange counts give the four-excess compression $133 I_3$, the vacancy compression $52.5 I_3$, and the product compression $185.5 I_3$. The threefold Fourier column of each rotation residue therefore has trial energy 185.5. This excludes absolute binding against all stock-preserving two-cell partitions in this fixture, including in the polarized sector. The time to dissociate in a larger, fully opened federation needs a different experiment. [Generator and census](material.py) · [Raw energy and shape banks](records.json.gz) · [Lean channel comparison](SignedPair.lean)

## Read both local charges and the axial moment on one launch

The source evolves under the native generator with a coherent energy comparison centred at 209. A threefold source coherently compares physical rotations about the gateway axis. Its plus-one residue accepts one in three source preparations. The sharp source bypasses this operation. The terminal unsigned support read commutes with the charge meters. The reader chooses uniformly among connected 3+4 partitions of the occupied sites using only occupancy and adjacency. Each retained target belongs to one of these disjoint regions; one records the three-site charge and the other the four-site charge. They may touch along open contacts. A regional read senses the signed material on each of its ports by the 02-versus-20 refocused ternary contact. Forward sensing on 02 and inverse sensing on 20 cancel the word-dependent common phase. The total phase is $2\pi\,130 Q/31$ on that region. A supplied 31-arm clock with repeated phase kicks reads bin $130Q\bmod31$, injective over the possible charges. Two such independent meters reading +1 and zero identify the retained class. All other regional outcomes, impossible partitions and clock refusals are recorded. Routing, coherent labels, regional assignment, timing and the material references are supplied apparatus. [Full ternary phase check](standards.py) · [Partition and independent readout](fragments.py)

For the shape meter, a uniformly chosen port supplies its material pointer outcome. At matched coordinates an observer records the Bernoulli acceptance proportional to $|3(z\cdot d)^2-|d|^2|$, then applies the sign of the geometric weight. Here $d$ is distance from the common centroid of the seven defect sites, each assigned unit mass. The unconditioned probe banks also use the cell centre and two transverse axes. A unit host has zero differential pointer response. Native phase analysis at 1/64 and 3/64 of a reference pulse period gives opposite signs. Known axial and transverse antipodal pairs give +3.997 and −0.798, against ideals +4 and −0.8. The absent source gives +0.002. The inverted cohort yields −4.009, +0.799 and +0.006 respectively. These measurements use the calibrated signed material interaction with ordinary 0, 1 and 2 loads. [Calibration counts](standards.json)

| Eight-clock cohort, one exchange time | Regional (+1,0) reads / bright clocks | Axial moment on (+1,0) | Inverted analysis on (+1,0) |
| --- | ---: | ---: | ---: |
| Native sharp | 153,537,065 / 274,222,274,788 | +0.06969 ±0.01372 | −0.06889 ±0.01372 |
| Native plus-one rotation | 51,164,290 / 91,407,406,903 | +0.06423 ±0.02368 | −0.06302 ±0.02368 |
| Exchange alone, sharp | 55,176,562,397 / 120,539,169,110 | +0.16472 ±0.00095 | −0.16484 ±0.00095 |
| Exchange alone, plus-one rotation | 18,326,983,737 / 40,179,366,028 | +0.16706 ±0.00148 | −0.16703 ±0.00148 |

Each register includes $2^{38}$ original simulated starts, drawn directly as multinomial counts on a 40-bit probability grid. The clock-bright denominator includes launches with no eligible partition. The first column of fractions uses a separate sign-blind regional cohort. The companion word archive keeps all 27,720 signed words at arrest in the two analysis arms of each of the four runs, with 13 admitted category counts per word. Its axes follow the run order, then forward/inverted analysis, the 12 charge-pair/pointer combinations plus an unassigned support, and the OW-14 sector index. The conditional word draws use a separate seed (161404) and retain exactly the registered category totals. [Full terminal words](joined_words.npy.gz) Each moment uses co-recorded local charge and pointer outcomes; 153,509,458 and 51,173,466 trials survive the two native forward shape selections, including zero-weight responses. The exchange-only generator has a zero-energy uniform ground in every stock and loses the binding inequality. Its larger positive conditioned moment shows the contribution of the geometric selection. [Registered counts, controls and refusals](fragments.json)

Without the regional charge selection, the four-clock sharp cohort gives **−0.65560 ±0.08338**, with **23,643 of 16,738,153** bright words even *admitting* a connected split with the original charges. The plus-one cohort gives −0.66963 ±0.14377. Testing existence after a signed survey admits more words than choosing an unsigned partition ahead of the regional charge reads. The common material charge, pose, source and law stay fixed during each registered comparison. [Full four-clock banks](records.json.gz)

### Uncertainty from the timing aperture and finite samples

The enrolled native sharp and normalized cyclic sources have mean 209 and variance 4.5. The spectral enclosure is [200,513], by nonnegative exchange and the sum of diagonal and off-diagonal absolute row weights. Put $\delta=E-209$. The absolute amplitude of one uniform 32-delay clock is $|\sin(\delta/16)|/[32|\sin(\delta/512)|]$. On $39\le\delta\le16\pi$ the numerator decreases and the denominator increases. Above $16\pi$ and up to 304 the numerator is at most one and the denominator is at least $32\sin(\pi/32)$. The full high-side amplitude is therefore at most 0.318822. Cantelli's bound gives high-energy source probability at most $4.5/(4.5+39^2)$. After eight simultaneous bright readouts its unnormalized weight is below $3.362\times10^{-11}$. A free dwell preserves it.

Let $K$ be the contraction for a successful regional pair, including a register for its observed assignment. For the normalized bright source and its truncation below 248, the sine of the angle between the normalized outputs is bounded by the discarded norm divided by the unnormalized accepted norm of the full source. This follows by projecting orthogonally to the truncated output. The same argument with an explicit environment handles a mixture of terminal assignments. The measured joint yields per successfully prepared source exceed 0.00053 even at simultaneous 99% confidence. On all seven-defect words, $\|Q\|\le7\times 2.524$ by the largest registered geometric weight. Before any regional read, removing the whole spectrum at or above 248 can thus shift either conditional moment by less than **0.00890**. Both native forward signs survive this correction. The corresponding unconditioned four-clock moment shift is below 0.021.

For $n$ admitted trials the signed read lies in [−1,1]. Each empirical signed mean and the enrolled gain have error at most $\sqrt{2\log(2\cdot4096/0.01)/n}$ with familywise 99% sampling coverage over these settings. The ± entries propagate these bounds through the calibrated signed estimate, using the known zero unit-host signal in the regional apparatus. Probability-grid error contributes less than $10^{-6}$ in the reported regional moments. The coverage describes the finite simulated counts. [Bounds and tests](test_signed_pair.py)

## Positive vector polarization on the earlier six-load band

There is an independent sign result on the older six-identical-load cap. In its twelve-dimensional band below 144 the physical vector representation has common energy 132.92397416. Use the positive fivefold rotation about tip zero and the three cap-orbit amplitudes proportional to $z$, $x+i y$ and $x-i y$. A uniform twelve-label coherent source with an appropriate contraction prepares these columns with original success probability 1/12; the additional 8,192-arm native energy port leaves approximately 0.08314 successes per original launch. For the circular polarization $m=+1$ the native material probe gives **+0.69918 ±0.53753** about the equal-load centroid. The $m=0$ preparation gives **−1.43806 ±0.53869** and their isotropic three-component mixture −0.00706 ±0.53653. The old sharp cap gives −3.47195 ±0.15945. The cell-centred reads are consistent with zero. This vector is an orbital representation of six identical charges on the empty-host convention. Its total material charge is six. [Measured moments](records.json.gz) · [Earlier band certificate](../PROOF.md)

The sign also holds on the exact low-energy vector multiplet. The cap at tip $v$ has centroid $(1+\sqrt5)v/6$, cell-centred moment zero, and centroid moment $-(1+\sqrt5)^2[3(z\cdot v)^2-1]/6$. The twelve vertices have second and fourth axial moments $1/3$ and $1/5$ along the tip axis. Weighting by circular intensity $1-(z\cdot v)^2$ gives the exact cap moment $(1+\sqrt5)^2/15$. The axial intensity $(z\cdot v)^2$ gives minus twice that number. Projection from the low band onto the twelve caps is an equivariant isomorphism, since both have dimension twelve and the squared noncap norm is at most $10/(4\cdot11^2)$. The vector irrep appears once in cap functions: a rotation-equivariant vector-to-cap map is fixed by evaluation at one tip, whose fivefold stabilizer fixes just its axial line. Thus the relative cap weights of the dressed vector modes are unchanged. The charge read is diagonal, so the complement enters as a convex weight. Its moment lies in [−6,12] because the six squared distances about their centroid sum to at most six. Even this worst-case dressing preserves the two signs. [Machine-checked final inequalities](SignedPair.lean)

OW-16 asks for a deuteron. A common calibrated length and electric scale, a nuclear spin assignment and an independently identified stable proton–neutron bound source are still needed. The six-charge multiplet does not supply those identities to the seven-defect signed measurement. The checked cheaper separation limits the physical interpretation of the latter.

From the RER root, replay the simulated counts and certificates with one numerical thread.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_16_binding.signed_pair
.venv/bin/python -B -m "$p.measure"
.venv/bin/python -B -m "$p.fragments"
.venv/bin/python -B -m "$p.standards"
.venv/bin/python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_16_binding
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_16_binding/signed_pair/SignedPair.lean)
```
