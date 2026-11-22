# Three spatial response components read by neighbouring patches

On the intact 20-, 80- and 320-cell screens, three neighbouring one-cell observers can jointly distinguish three independent variations in the sender’s local load distribution. Their recorded counts identify all 144 held patterns across nine sender sites and two exposure settings. At three additional sender sites, a ten-cell patch predicts the neighbours’ measurements at three untrained exposures. Its largest error in an individual detection probability is 6.32 × 10⁻⁷. [Response measurements](result.json) · [Held station forecasts](forecast.json)

## Preparation and readout

The free generator is the registered Interval Echo one-load generator on all ports and seams. A sender prepares one sharp load per trial. Twelve anonymous source ports define the possible preparations. A recipient reads one of its twelve ports or its quiet outcome. A distinct ensemble is prepared for each recipient. The joint reader pools the three local inboxes by reported preparation and exposure. The three recipient readings on a row come from independently prepared ensembles. Each source, neighbourhood and encounter has its own finite aperture. Cohort pooling and estimation are supplied analysis. [Candidate implementation](../../candidate_laws/interval_echo/candidate.py) · [Acquisition and inverse](instrument.py)

At a short calibration dwell the source measures the probability table for its own ports. The gap between two sets of off-diagonal probabilities gives the internal incidence. Its first nonconstant Laplacian eigenspace has dimension three. The source prepares small positive and negative variations of the uniform port probabilities along an orthonormal frame of that measured subspace. The three signs in a held word select a direction between the calibration axes. None of the preparations requires a coherent multiport source. The simulation uses complete Born probabilities for the short calibration and separately draws finite counts for the indicated held tests. No embedding coordinates enter the inference.

The lab supplies patch membership, addressed ports, the mixtures, the native-time exposure and local detection. Between preparation and detection the excitation evolves freely. The elapsed-word clock is idle. In the one-load sector the generator is the federation Laplacian divided by two, plus 6κ times the identity. The constant contributes a common phase. The same measurements therefore apply to the contrast-memory reference and are insensitive to κ in this sector. [One-load reduction and band](../../evaluation/repair_band.py) · [Specification](../../candidate_laws/contrast_memory_spherical_v1/SPEC.md)

## Direct test of three independent responses

At each of the nine sender sites, each recipient's calibration includes the uniform preparation and the six signed single-axis perturbations. Eight different source mixtures are held out for each of two exposure settings. The rank-two control keeps the two singular components of the calibrated receiver response that minimize its squared error. A separate control supplies an unmodulated, uniform source for all intended words. All counts, including quiet, remain in the receiver distributions. [Frozen protocol](PROTOCOL.md) · [Exact probabilities and finite counts](records.npz)

| Native exposure | Sender sites | Smallest / largest singular value of joint response | Relative held-signal RMS error after retaining two components | Decoded patterns |
| --- | ---: | ---: | ---: | ---: |
| 2 | 9 | 0.248–0.499 | 17.7–40.3% | 72/72 |
| 4 | 9 | 0.205–0.286 | 16.3–24.7% | 72/72 |

Each local setting uses 4,194,304 simulated preparations. The sham decodes 16 of 144 words, with 18 expected for uniform guessing. The exact three-component interpolation agrees within 3.4 × 10⁻¹⁶, as expected from the linearity of mixtures. The physical tests are the observed third singular component, its effect on held predictions, and decoding with independent counts. Pooling matters. Only one of the 54 individual recipient panels has a smallest-to-largest singular-value ratio of at least 0.05. All 18 three-recipient panels meet that threshold. A receiver's probability of detecting a load ranges from 0.0372 to 0.0937 over the held settings. The operator fit and the held predictions use the absolute, unconditional counts. [Per-site results](result.json)

## Forecast at three additional sender sites

A ten-cell observer patch consists of the sender and two layers of neighbours. One very short survey at external dwell 1/32 gives the conditional port readings in that patch, including the quiet complement. The reader infers contacts from the off-diagonal readings and total port degrees from the return probabilities. It builds the finite Hamiltonian from these local observations and the single-cell hop coefficient 0.5 inferred at exposure 1/8. The simulator keeps all external links during the actual survey and the subsequent flight. The finite forecast uses the projected Hamiltonian on the surveyed patch. Propagation out of the patch and back can produce a discrepancy. [Held forecast protocol](FORECAST_PROTOCOL.md) · [Forecast calculation](forecast.py) · [Independent calibration](../../evaluation/patch_spatial_readout.json)

An incidence-isomorphism check preserves total port degrees and the sender, inner and outer patch roles. The three sender cells have the same marked geometry. Adding immediate neighbours distinguishes the 20-cell site, while the 80- and 320-cell stars agree. The complete ten-cell surveys distinguish all three sites. Their similar responses reflect shared local structure. [Independent geometry and forecast checks](test_forecast.py)

The three enrolled senders are cell 10 on the 20-cell screen, cell 40 on the 80-cell screen and cell 160 on the 320-cell screen. Each prediction uses only the short survey, preparation, declared form of the walk and exposure. All three receiver cells at each site, every held exposure and all eight words enter the comparison. The exact conditional probabilities used by the survey are a supplied input to this test. Its finite-count check draws 1,048,576 preparations per held setting per recipient.

| Held dwell | Largest absolute error in any port or quiet probability across boards | Largest relative RMS error in the three-component signal | Relative RMS error with two components | Patterns decoded |
| --- | ---: | ---: | ---: | ---: |
| 0.75 | <3 × 10⁻¹² | <2.0 × 10⁻¹⁰ | 18.76% | 24/24 |
| 1.5 | <4.3 × 10⁻⁹ | <2.4 × 10⁻⁷ | 18.50% | 24/24 |
| 3 | <6.32 × 10⁻⁷ | <2.22 × 10⁻⁵ | 18.04% | 24/24 |

[Held station measurements](forecast.json) · [Surveys, forecasts and counts](forecast_records.npz)

The measured directions describe three components of a prepared occupation distribution. Their transmission is detected across one intercell step. Establishing an effective spatial metric requires independent ranging and elapsed times on material patch clocks. A continuum comparison also requires matching physical regions across refinements. Neither experiment measures transport of an arbitrary coherent quantum state in a three-component wave band.

## A clock that runs during free flight

Six disjoint seam memories in every cell contain signed contrast of magnitude two. A single port carries amount one or two and every other port is vacant. A native write of this contrast has output loads (1,1), orthogonal to all the permitted inputs at a clock seam. The echo fixes those signed states for either coin. One hand has initial sign plus. A second has amplitude $(|{\rm even}\rangle+i|{\rm odd}\rangle)/\sqrt2$ in the sign-flip eigenbasis. The other four use the even state. The local hand term is $4\kappa I+b(I-X)/2$. Two sharp sign reads give the cosine and sine of the angle $\theta=bt$. Each detecting cell can read its own hands while the same load moves through their endpoints. This is a stationary, synchronized apparatus. Ordered independent cohorts resolve the hand's winding. The candidate's guided interval counter has a different preparation and timing rule. [Native contact matrices](../../candidate_laws/contrast_memory_spherical_v1/law.py) · [Source and carrier restriction](shell_source.py)

On the blank seams the mobile load sees exchange $b/2$ times the graph Laplacian and contrast $\kappa q^2$ per incident edge. Each port has one signed seam. Its active degree plus its original number of degree references is five. The free Hamiltonian restricted to this sector is the sum of the hands, $5\kappa q^2 I$, and $bK$ with $K=L_{\rm active}/2$. All federation contacts are on during the native flight. Half exchange and reversed contrast give the same load and hand distributions at equal hand angles. Two additional controls disable external exchange or park the load while the hands continue to run. [Flight protocol](shell_protocol.json) · [Joint outcomes](shell_counts.npz)

The apparatus supplies synchronized signed memories, coin zero, a sharp load, addressable preparation and read windows. The native write of (2,0) or (0,2) produces the requested sign on balanced loads with a common coherent coin. Initializing coin zero and clearing the fuel require preparation. A spanning tree of the blank contacts distributes the mobile load by native partial exchanges. Each daughter receives the squared norm of its requested subtree. A local degree-reference pulse corrects the phase left at each port. At degree six a blank-bond swap to a degree-five neighbour and its inverse implement that phase. The native $\pi$ swap has unit transfer amplitude for both amounts at $\kappa=8$. Preparation is held fixed when the flight law changes. Terminal echo, sharp signed-difference read and reverse echo against an empty reference distinguish loads 0, 1 and 2. The same swap gives this read at a gateway. The region detector groups commuting port occupation measurements and reads the two signs at the occupied cell. All the other cells are quiet in the one-position sector. The controller supplies terminal collection and synchronization. [Operator checks and source compilation](shell_source.py) · [Measurement tests](test_shells.py)

## Angular motion and a candidate depth coordinate

Four adjacent children make each coarse face. The acquisition verifies connectivity and the two fine interfaces per coarse adjacency from native seam incidence. The twenty such regions serve as corresponding detector addresses across the three boards. Their identification is an apparatus choice. No measured distance is assigned to a carrier edge. A cell reports which of the two ends of its signed matching received the load. The complete joint datum consists of region, depth bin and two local hand signs. The angular source on a cell is $\max(0,1-hd/2)$ with $h=2^{-\ell}$ and $d$ its contact distance to a marked macroface. The smooth source distributes this amplitude equally over the twelve ports. The polarized control occupies just one end of each pair. Two anchors and all four flight arms use the same decoder. [Three-screen measurements](shell_result.json)

| Cells | Kinetic energy, smooth source | Peak depth TV, anchor 0 / 1 | Region TV at $4/h$, anchor 0 | Region TV at $4/h^2$, anchor 0 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.053571 | 0.025535 / 0.053699 | 0.018751 | 0.018751 |
| 80 | 0.011503 | 0.014552 / 0.016882 | 0.009248 | 0.019305 |
| 320 | 0.002677 | 0.006156 / 0.007432 | 0.005523 | 0.022054 |

TV denotes total variation. The depth comparison uses the same angular marginal and equal population in the two depth bins. The region comparison uses the emitted distribution. These are Born-law distances. The acquisition retains 547,356,672 simulated trials. At each size and launch, the peak over the scheduled exposures separates native transport from each disabled-transport control beyond the simultaneous 99% sampling allowance of 0.00985 in angular TV. The half-rate, reversed-contrast arm stays inside that allowance. The polarized source gives depth TV 0.5 at emission on every board. The largest sampled hand-angle error across all arms is 0.00415 radians or less. [Exact distances, counts and individual controls](shell_result.json)

For any extension with the same cell and external matching, let $J$ lift one amplitude per cell with normalization $1/\sqrt{12}$, and $Q=I-JJ^*$. The first nonconstant icosahedral eigenvalue is $5-\sqrt5$. Removing a matching changes the internal Laplacian by at most two. Thus $K\ge gQ$ with $g=(3-\sqrt5)/2$ (the registered matching gives an internal half-Laplacian gap of one). In a unit state of energy $E$, $\eta=\|Q\psi\|^2\le E/g$ at every exposure. Any balanced six-versus-six depth sign $D$ has $J^*DJ=0$. Choosing its sign separately in each detector region bounds the depth TV by $\sqrt{\eta(1-\eta)}+\eta/2$. Hence *every* such read loses depth contrast at vanishing energy, including alternative relabellings of the six ports.

Could a coherent detector recover depth from port phases? At any exposure with $\eta<1$, the incident state has trace distance $\sqrt\eta$ from its normalized cell-mean projection. Attaching the same clocks and applying a common quantum instrument contracts that distance, even when the instrument keeps its quantum output. Two incident states with identical normalized complex cell means and energies at most $e<g$ therefore differ by at most $2\sqrt{e/g}$ in any terminal record. The premise includes the angular phases. A detector may distinguish waves with different angular phases or information encoded in additional occupied cells or memories.

The slow modes can also be counted. The compressed form $J^*KJ=L_{\rm dual}/24$ follows by summing the three external contacts per cell. Write $E=E_{\rm int}+E_{\rm ext}$. The incidence matrix of an external matching has norm at most $\sqrt2$. The triangle inequality bounds the square root of the compressed form on $J^*\psi$ by $\sqrt{E_{\rm ext}}+\sqrt{E_{\rm int}/g}$. For a spectral subspace below $e<g$, projection to cell means is injective and has squared norm at least $1-e/g$ on unit vectors. The min-max principle sandwiches its dimension between the counting functions of $L_{\rm dual}/24$ at $e$ and $(1+1/g)e/(1-e/g)$. The cell means account for all soft degrees of freedom up to this energy comparison.

For either smooth launch the squared unnormalized cell amplitudes sum to at least one twentieth of the cell count. Adjacent amplitudes differ by at most $h/2$. There are three external ends per cell, giving $E\le5h^2/16$. The spectral bound $\|(e^{-i\theta K}-I)\psi\|^2\le2|\theta|E$ makes every terminal measurement of these waves converge to its parked-signal value at $\theta=O(h^{-1})$, including the synchronized clock read. The finite slower-time data support detectable transport. They determine neither a limiting equation on that scale nor a physical ruler. A freely propagating radial field would require an independently refinable family of distinguishable depth records on the same clock. A higher band, another collective preparation or a depth stored across cells would require its own propagation and detector analysis. [Smooth-description hypotheses](../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf#page=8)

## Do the two moving loads have distinct relativistic masses?

On the 20-cell board the experiment compares the uniform mode with the entire first positive kinetic line. Its eigenvalue in hand units is 0.0199624555 (multiplicity one). A supplied isolated fuel qutrit couples to the common home port through a native half exchange. The data arm prepares the requested mode by the tree above. During free flight the fuel is parked with zero or one active degree reference. The inverse preparation, a supplied fuel phase of 0 or $\pi/2$ and inverse half exchange give three outcomes: fuel, home or return elsewhere. If $A$ is the full spectral return amplitude and $F$ the fuel amplitude, the first two probabilities are $|Fe^{-i\phi}\pm A|^2/4$. The complement is $(1-|A|^2)/2$. Both clock signs are read jointly with each analyzer. The fuel is an auxiliary apparatus port with a switched contact and a common reference convention for both amounts. The supplied preparation, contact windows and readout are outside the free flight. [Return protocol](shell_spectroscopy_protocol.json) · [All 240 eigenvalues and endpoint intervals](shell_spectroscopy_result.json)

The rest frequency relative to $r=0,1$ fuel references is $R_q=(5-r)\kappa q^2/b$. The measured kinetic increment $\epsilon$ is the same for both mobile amounts on the same mode. A relativistic law with a common momentum and speed would give the same squared-frequency shift for both. The difference is exactly $2(R_2-R_1)\epsilon$. It depends on the fuel convention. Equal kinetic increments at fixed nonzero momentum also preclude *distinct* relativistic rest masses for this assignment, since $\sqrt{m^2c^4+p^2c^2}-mc^2$ strictly decreases with $m\ge0$. This test does not supply a spatial momentum calibration or identify a particle. [Two-load return operator](shell_spectroscopy.py)

| Flight law | Active terms on the parked fuel | Exact difference in squared shifts | Simulated beat counts, simultaneous 99% interval |
| --- | ---: | ---: | ---: |
| Native | 0 | 4.7910 | [4.126, 5.554] |
| Native | 1 | 3.8328 | [3.262, 4.403] |
| Half exchange, reversed contrast | 0 | −9.5820 | [−11.049, −8.199] |
| Half exchange, reversed contrast | 1 | −7.6656 | [−8.794, −6.516] |

The endpoint intervals include zero in all four settings. The independent repeat uses the same shot budget and a separate seed. Each time has independent rest and line quadratures. Their complex product after conjugating the rest gives the beat. Equal windows of 128 times at either end have means separated in phase by 385 steps of the beat, for single-frequency returns. Each measured real or imaginary component is a product sum of four independent sample means. Changing one trial changes the window average by at most $2/(128N)$ for $N$ trials per analyzer. Bounded differences controls both window arguments, and separate Hoeffding bounds control point arguments, winding and hand duration. The combined error probability is at most 0.01 over both sectors, laws and reference settings. Each acquisition retains 4,303,355,904 simulated trials, including the missed-return category. Cumulative probabilities are rounded to a $2^{-36}$ grid before sampling. The phase budget covers the resulting mean error and the verified source residual. [Frozen beat test](shell_beat_protocol.json) · [Beat intervals](shell_beat_result.json) · [Independent counts](shell_beat_counts.npz)

Both estimators check the nonnegative integer count in each joint bin and the shot total at each exposure before pooling over clock outcomes. Valid cohorts with insufficient phase visibility or an unreadable hand give an unresolved interval. An independent test constructs the two position spectra from matrix elements of the 99-state contacts. Their analytic detector laws agree with the archived joint bin counts of both seeds within a simultaneous sampling bound. The contact law, carrier and declared preparations are shared inputs to this check. [Independent detector checks](test_shells.py)

## What a different realization must change

For depth stored inside a cell, a detectable pair of states with matched cell means at vanishing kinetic energy requires the gap between a uniform cell state and some internal excitations to close under refinement. Those soft internal states must also exchange signals across independently ranged depth regions. An alternative is to use collective states spanning several cells or evolving memories, where the single-position bound does not apply. In either case the tests require corresponding prepared regions, a radial and angular propagation law on the same hands, and a mapping of the resulting records to the scalar instruments. For two distinct relativistic masses on one spatial dictionary, the moving objects must acquire different kinetic increments at the same nonzero momentum. Changing only the contrast or the parked fuel energy leaves the kinetic operator here identical for both amounts.

## Reproduce

From the repository root, run the following to read the three published acquisitions on any supported platform.

```sh
python -B -O -W error -m model.experiments.interval_echo_space.shell_archive
```

This replay verifies the published count files by SHA-256, calculates each joint detector probability from the frozen protocols, checks every measured bin against a simultaneous 99% sampling bound, and feeds the original counts into the experiment’s readers. It compares every resulting summary with the stored result at publication precision.

The following commands regenerate the acquisitions from their seeds and require exact count equality. On macOS arm64 six cohorts differed from the archived draws by a single trial moving between two bins. The three archival reads reproduce the published statistics on that platform.

```sh
python -m model.experiments.interval_echo_space
python -m model.experiments.interval_echo_space.forecast
python -B -O -W error -m model.experiments.interval_echo_space.shells
python -B -O -W error -m model.experiments.interval_echo_space.shell_spectroscopy
python -B -O -W error -m model.experiments.interval_echo_space.shell_beats
python -m pytest -q model/experiments/interval_echo_space
```

Each command accepts `--output DIRECTORY` to write an independent copy. In `records.npz` keys identify the board size, sender, dwell and record type. The last receiver outcome is quiet. `forecast_records.npz` adds the local survey, model forecast and uniform-source records. The JSON files include the relabellings, detector cells, data selection and protocol hashes.
