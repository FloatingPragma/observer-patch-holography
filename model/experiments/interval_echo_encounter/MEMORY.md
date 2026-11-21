# A scheduled native repair with its working memory retained

At duration one, cell 0 applies Interval Echo’s native 99-state echo U to seam (0, 1). The seam begins blank. Its action on that input is the contrast write

$$
W|a,b\rangle = \frac{1}{\sqrt 2} \sum_{c=0,1} |c,a-b\rangle\,|y_c(a+b)\rangle .
$$

The first arm keeps that coherent working memory and lets the federation evolve for four more units. The matched erase arm applies U† immediately after U and then evolves for the same duration. The bypass arm uses the same separated preparation and free evolution. All 36 sharp, initially nonadjacent preparations receive the service regardless of the loads present on (0, 1). The only observer reading in the retained arm is the final seam/quiet measurement on cell 0. Its 21 internal seams disjoint from ports 0 and 1 provide the primary outcome. The free flight, the local service and the reading use the fixed sphere and the registered operators. The service appointment, waiting durations and reader are supplied. [Protocol](MEMORY_PROTOCOL.md) · [Native operators](../../candidate_laws/contrast_memory_spherical_v1/law.py)

## Observer results

| Measurement, all 36 separated preparations | Retained write | Write followed by inverse |
| --- | ---: | ---: |
| Final pair on an internal seam disjoint from the repair endpoints | 0.12672518 | 0.00348463 |
| Weight in the reference bound band immediately after the service, global diagnostic | 0.05241507 | 0.00291210 |

Retaining the repair increases the disjoint-seam yield by 12.324 percentage points and the immediate bound-band weight by 4.950 percentage points. Erasing the write restores the bypass probabilities on every local seam and on quiet, within 10⁻¹⁰. The disjoint terminal detector and unconditional enrollment exclude a coincidence with the written seam and selection on a first detection. The comparison changes both the coherent load-memory state prepared by the service and the subsequent response of the tagged seam. It measures their combined effect. The schedule gives the load-memory system an external intervention. A first record of autonomous encounter and a later repeated reading of that object are separate requirements. [All terminal probabilities](memory.json)

## Evolution with a written seam

The other seams have blank working memory and preserve their blank subspaces. On the addressed seam the matrix is UGU† throughout the flight. The code removes the blank term for that seam from the two-load Hamiltonian, acts with the remaining terms on every memory component, and inserts the native nonblank term on its two actual endpoints. This accounts for the physical effect of leaving that seam written during the flight. [Implementation](memory.py)

Binary input loads have contrast 0 or ±1. The write occupies the symmetric coin combination for zero contrast and two coin states for each parity of the ±1 contrast labels. The calculation checks directly in all 99 local basis states that the write maps the incoming states into this five-dimensional memory space times the binary endpoint occupations, and that the full tagged term maps each of the three component spaces into itself. Blank terms on the other seams preserve binary loads by exchange. A separate one-cell calculation evolves the full eleven-state memory and all 78 possible two-load configurations, including all twelve doubles, and compares the final probabilities and the load-band diagnostic. [Independent full-basis check](test_retained_memory.py)

The load-band diagnostic immediately after the service uses the blank-sector bound projector on the loads and identity on the working memory. This determines which load components the write creates or removes in the reference band at the instant of the service. A subsequent joint Hamiltonian acts on the retained memory. The terminal seam read is available to the patch. The projected spectral quantity uses the state of the full federation.

## A clock formed in the same write as the pair

One sharp double on port 0 and a vacancy on port 1 produce a binary pair and a coherent sign on their shared seam. During free flight the sign stops evolving on the word with both ports occupied. A patch can measure the resulting slow beat. At the registered contrast $\kappa=8$ its frequency is 0.00260641 per radian of a remote material hand. The observed simultaneous 99% interval is [0.002214, 0.002981]. This source tests the coupling between the written record and the matter that created it. The prepared double and the written address fix the object’s birthplace. [Frozen source](dwell_protocol.json) · [All sign, quiet and clock counts](dwell_counts.npz)

The native echo sends $|20\rangle$ and $|02\rangle$ to $|11\rangle$ with signed contrast $+2$ and $-2$, respectively, and the same symmetric coin. The second quadrature uses an inverse half exchange of the double before the echo. At reception the reverse echo restores a double on the first or second port precisely for the tagged pair. Every other binary load word gives quiet. The two double-location effects can be read as commuting sharp occupation flags. The supplied terminal apparatus may route a gateway occupation to a degree-five neighbour by a blank-seam swap and undo the route. A vacancy rail and a native echo expose the load as a contrast flag. The ordinary public reader can commit each orthogonal flag. These formulas use ideal isolated preparation and terminal contacts. The free Hamiltonian acts on all seams between them. [Physical columns and terminal effects](test_dwell.py)

Two disjoint native dimers in the most distant cell supply a stationary timebase. Each contains a single unit. Their two positions give the cosine and sine quadratures in the same trial. Eighteen native writes on empty seams isolate the cell and the two hands. A vacuum write prepares contrast zero, whose free contact vanishes for all nine endpoint words. Both dimers use degree-five ports. Their equal diagonal energies cancel in position measurements, leaving a hand rate of one for every contrast in the protocol. One dimer receives a $3\pi/2$ preparation pulse. The controller supplies four load units in total, the initial addresses, ordered appointments, isolated source and terminal gates, and pooling of the patch records. The 20-cell flight retains all $\binom{228}{2}=25{,}878$ exterior binary words, two sign parities and all four reference outcomes. The remaining twelve ports belong to the timebase cell. The simulated panels draw independent multinomial ensembles from the complete joint probabilities. No quiet trial is discarded. [Station and resource counts](dwell_result.json)

The source write performs work. At the native tag the initial sharp double has energy $24\kappa+d_0/2$ above the empty exterior, with $d_i$ the original port degrees. The written pair has expectation $10\kappa+(d_0+d_1-2)/2$. Their difference, $-14\kappa+(d_1-2)/2$, is the change in the expectation of the free Hamiltonian during preparation, with the same reference dimers in both comparisons. Free binding and phase measurements begin after this service.

### Why the rate is small

Let $P=|p\rangle\langle p|$ project onto the tagged pair and $Q=I-P$. Write $K$ for exterior exchange, $A_{ab}$ for active adjacency and $\ell_a$ for the original port degree minus its active degree. The energy relative to $12\kappa$ and the timebase in sign parity $s=0,1$ is

$$H_s=K-\kappa(\ell_a+\ell_b+2A_{ab})+(4\kappa+s)Q.$$

The diagonal contrast formula retains the original degree references. The last term follows from the 99-state echo: the two signed, symmetric-coin states with loads (1,1) lie in the write image and map to blank memory with zero generator. At every other endpoint word these two signed states are orthogonal to that image. They see the ordinary contrast-square energy and sign reversal. The comparison that always ticks restores the whole $4\kappa+s$ term on $P$. The single-term controls restore just $s$ or just $4\kappa$. Their full Hermitian extensions agree on the other 97 local dimensions. [Native restriction](dwell_source.py) · [Independent contact assembly](dwell_checks.py)

For the positive law the complement lies at least $D_s=19.5+s$ above the sharp tag diagonal on all three released boards. This follows from positivity of $K$, a maximum off-tag saving $\ell_a+\ell_b+2A_{ab}\le3$, and the tag values $\ell_0+\ell_1=2$ and $\langle p|K|p\rangle=4.5$. For reversed contrast a lower bound on the distance above the complement uses $\|K\|\le2\max d_v$. Interlacing and the variational principle then give a unique isolated line in each parity (lowest for positive contrast, highest for reversed). The coupling $v=QK|p\rangle$ has norm 1.5. The eigenvector’s tag weight satisfies $w_s\ge[1+\|v\|^2/D_s^2]^{-1}$ by the complementary resolvent. All other lines together carry at most $1-w_s$ of a sharp launch. Thus the probability of finding the pair on its birth seam is at least $(2w_s-1)^2$ for every exposure. The smaller bound at positive contrast is 0.976609. This proof uses the registered graph incidence at each size and extends to graphs satisfying the same inequalities.

For $0\le s\le1$ the simple eigenvalue has derivative $\|Q\psi_s\|^2$ by Hellmann–Feynman. An excursion off the seam is necessary for its clock to advance. Put $D_{ab}=6-\ell_a-\ell_b-2A_{ab}$ off the tag. The eigenvector’s unnormalized tail is $-(\kappa D+B_s)^{-1}v$ when its tagged amplitude is one. With $d=4.5$ and $u=d-2\kappa$, the low eigenvalue lies in $[u-\|v\|,u]$. The remainder therefore has norm at most $M=\max(d,2\max d_v+1-d+\|v\|)$ using active degrees. For $\rho=M/(\kappa\min D)<1$, the triangle and reverse triangle inequalities applied to $I+D^{-1}B_s/\kappa$ give

$$
\frac{C}{\kappa^2(1+\rho)^2+C}\le E_1-E_0\le
\frac{C}{\kappa^2(1-\rho)^2+C},\qquad C=\|D^{-1}v\|^2.
$$

Four of the nine outgoing moves leave an adjacent pair and five do not. Their numbers $D_{ab}$ are three and five and each exchange amplitude is $-1/2$. This gives $C=29/180$ exactly, without fitting the five spectral measurements. A second solver uses the scalar equation $E=u-v^*(QH_sQ-E)^{-1}v$ and computes the rate from $\langle\psi_0|Q|\psi_1\rangle/\langle\psi_0|\psi_1\rangle$. Values of $\kappa^2(E_1-E_0)$ at 8, 16, 32, 64 and 128 are 0.166810, 0.163912, 0.162500, 0.161802 and 0.161456. Only eight is the candidate parameter. [Exact coefficient and finite enclosures](dwell_scaling_result.json)

### What the observers can distinguish

Write $a_s=\langle p|\exp(-itH_s)|p\rangle$. The probability of a tag hit in either quadrature is $(|a_0|^2+|a_1|^2)/2$. The differences of the two terminal sign probabilities read the real and imaginary parts of $a_0a_1^*$. Independent reference positions give the elapsed phase with the same enrollments. The full signed difference uses the denominator of all preparations. The complete spectra contribute to both amplitudes. The bound $\arcsin((1-w_s)/w_s)$ limits the phase error from every complementary line together. The inference combines this with simultaneous Hoeffding bounds on all signed marginals and the reference duration. Sampling at successive appointments resolves winding because the line separation is at most one and the combined adjacent phase error is below a half turn. Missing visibility or an uncertified native line gives an unresolved rate.

| Tag interaction | Contrast | Lowest sampled hit fraction | Clock rate interval |
| --- | ---: | ---: | --- |
| Native | 8 | 0.99176 | [0.002214, 0.002981] |
| Continuous ticking and contrast energy | 8 | 0.89091 | [0.999675, 1.000300] |
| Native binding, continuous ticking | 8 | 0.99140 | [0.999742, 1.000308] |
| Continuous contrast energy, native ticking | 8 | 0.89342 | Unresolved |
| Native | −8 | 0.99298 | [0.001906, 0.002938] |
| Continuous ticking and contrast energy | −8 | 0.84141 | [0.999664, 1.000308] |
| Native / continuous | 0 | 0.00025 / 0.00014 | Both unresolved |

All eight settings comprise 1,627,389,952 simulated trials over 97 appointments. The six initial settings share one 99% simultaneous budget. The two single-term controls share a second 99% budget and an independent seed. In particular, retaining the binding surface while restoring the unit tick preserves a large hit fraction and recovers rate one. The clock is sensitive to how the law responds to joint occupation. Its rate in this stationary preparation supplies no measurement of relativistic time dilation. [Complete primary acquisition](dwell_result.json) · [Single-term controls](dwell_controls_result.json)

### Can this clock travel with a particle?

For each seam the native free law conserves its blank projector and its written contrast magnitude and sign-reversal parity. The image of the write is invariant under contrast reversal (swap the input loads) and under the contrast square (the stored magnitude is fixed for each input). The echo-conjugated free generators therefore have zero mixed blank/written blocks. On the written block they commute with each of these projectors. Other seams can act on the shared endpoints but commute with a projector on this private memory. The exact integer check exhausts all 99 local states and both generator parts. The sign observable can oscillate because it reads coherence between conserved parity sectors. The location of the written contrast is fixed. Transport of this particular clock mark requires a service or a change of the law; a clock encoded in moving loads has a different source and readout. These results give a stationary dressed excitation and an emergent frequency. Collision-born mobility, a rest-mass assignment and an empirical particle dictionary require separate measurements. [Memory charge decomposition](dwell_scaling.py)

## Reproduce

From the RER repository root with the model dependencies installed, run

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m model.experiments.interval_echo_encounter.memory
python -B -O -W error -m model.experiments.interval_echo_encounter.dwell
python -B -O -W error -m model.experiments.interval_echo_encounter.dwell_scaling
python -B -O -W error -m model.experiments.interval_echo_encounter.dwell_checks
python -m pytest -q -W error model/experiments/interval_echo_encounter
```

For the extended audit, set both BLAS and OpenMP thread counts to one and run the last module with `--extended`. It reconstructs the generator from physical contact columns and propagates a full state by short Taylor steps. The bilinear sum of squares at half exposure returns the sharp-state amplitude at full exposure for a real symmetric Hamiltonian. The audit restores the common phase, encloses Taylor and floating-point errors, and checks every joint bin for every setting.
