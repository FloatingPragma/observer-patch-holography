# One local repair rule for overlapping regions

A cell records its own volume bit and the local endpoint of each shared face. On a four-cell cycle the same stationary preparation gives an exact entropy–aperture first law on every region, including overlapping and disconnected choices. One face calibrates a gain of one bit per open aperture. A separate change involving two faces and their correlation predicts a gain of $8/5$ open apertures and a loss of $8/5$ bits per unit of perturbation on the region containing cells 0 and 2. The observers read aperture settings and word probabilities separately. [Local transition and observer map](regional.py) · [Exact enumeration of 10,000 joint words](regional_results.json)

## Microscopic transition

Each interface stores one pair of sharp endpoint letters. The five legal words, in cycle order, are $(0,0),(1,1),(1,2),(2,2),(2,1)$. On an interface update two fair bits choose the previous word, the same word, the same word or the next word. The update acts on both ends of that interface. It reads only the current word. The row and column sums of the transition matrix are one. Its two-step kernel has every entry at least $1/16$, so its two-update contraction factor in total variation is at most $11/16$, with $\|P^2\mu-P^2\nu\|_{\rm TV}\le(11/16)\|\mu-\nu\|_{\rm TV}$. The unique stationary law assigns probability $1/5$ to each legal pair.

The separate cell volume bit changes from 0 to 1 with probability $1/4$ and from 1 to 0 with probability $1/2$. Detailed balance gives marginal probabilities $(2/3,1/3)$. The interface and cell updates use independent commands, with enough available updates on each register. Their product has a unique global reference even on a graph with cycles, from any initial distribution on legal words. For a point initial condition, the independent outputs converge to the product reference. The same follows for a correlated start by linearity on a finite state set. The receipt also propagates the exact probabilities of a preparation with every word and bit initially zero for 4, 16 and 64 sweeps. This is an added repair rule on the declared records. Appending independent records to other dynamical matter leaves the entropy variations of that matter to be accounted for separately.

## How a region obtains its area and entropy

For a region $R$ consisting of a set of cells, the observer reads both endpoints on each interior face, one endpoint on each face crossing the boundary, and the volume bit at every included cell. Each open cut shutter has a fixed independently specified footprint of one aperture unit. The local endpoint symbol 0 reports a closed shutter. Either symbol 1 or 2 reports open. The two ends agree on the shutter status in every legal word. Write $A_R$ for the number of open cut shutters, $V_R$ for the sum of the visible volume bits, $n=|R|$, $i$ for its interior faces and $c$ for its cut faces.

At equilibrium the fully seen face word has probability $1/5$. A closed cut endpoint has one possible partner. Each specific open cut endpoint has two possible partners. Its marginal probability is accordingly $1/5$ or $2/5$. For every accessible local word $w$,
\[
 P_R(w)=5^{-i-c}3^{-n}\,2^{A_R(w)+n-V_R(w)},\qquad
 -\log P_R(w)=C_R-(\log2)A_R(w)+(\log2)V_R(w).
\]
The multiplicity follows from the transition and compatibility of the endpoint words. The area counter reads the shutter word and includes every closed outcome. All regional marginals restrict to the same records on intersections because they come from one global law. The code independently checks all 65 strict restriction maps on the four-cell cycle for the reference and four signed transfers, including a change in correlation between faces. No tree or conditional independence between overlapping *regions* is required.

Let $Q_R=P_R+\epsilon D_R$ be any differentiable change supported on these words, with $\sum_wD_R(w)=0$. It can arise from a change that correlates arbitrarily many interfaces. Direct differentiation of Shannon entropy gives
\[
 dS_R+(\log2)\,d\langle A_R\rangle=(\log2)\,d\langle V_R\rangle .
\]
At fixed volume this forces first-order stationarity in nats. The volume coefficient and area gain reflect equal registered units in this experiment. For a footprint $a$ and a volume increment $v$ they become $(\log2)/a$ and $(\log2)/v$. For finite changes the exact discrepancy has a sign and magnitude,
\[
 \Delta S_R+(\log2)\Delta\langle A_R\rangle-(\log2)\Delta\langle V_R\rangle
 =-D_{\rm KL}(Q_R\Vert P_R).
\]
This proves strict maximization on the admitted regional distributions at fixed mean volume. It does not assert that a perturbation leaves the reference invariant under the dynamics. The test prepares changed distributions to measure the response.

The Lean proof checks the four repair commands and both endpoint fibers, obtains every cut's probabilities by summing over complete tickets, and differentiates the Shannon entropy of an arbitrary normalized path through this equilibrium. Its fixed-volume theorem gives zero derivative for entropy plus the expected aperture times $\log 2$. The same ticket bank permits different regional selections. [Kernel-checked regional first variation](../../../Lean/ObserverPatchHolography/EinsteinBranch/TileCapacity.lean)

## Independent calibration and two falsifying changes

The cycle comprises $0$–$1$–$2$–$3$–$0$, with its faces numbered from 0 in that order. Set $f_e=5\,1_{\rm open,e}-4$, centered under the uniform face law. The calibration perturbs $P$ to $P(1+\epsilon f_0)$ and measures cell 0 at fixed mean volume. The independent word probabilities give $dS/\log2=-4/5$, while its two apertures give $d\langle A\rangle=4/5$. Freeze this gain. The held preparation uses $f_1+f_2+f_0 f_2$ in place of $f_0$ on the same global ensemble. For the overlapping regions $\{0,1\}$, $\{1,2\}$ and $\{0,2\}$ the respective area slopes are $4/5,4/5,8/5$. All entropy slopes have the opposite signed values in bits. The pure correlation change $f_0 f_2$ holds the means of every shutter and volume bit fixed. On $\{0,2\}$ at $\epsilon=1/16$ the entropy nevertheless falls by 0.0276328 nats. All occupations, misses and signed masses in this exact population calculation come from the finite local rules. [Four transfers and blank-start evolution](regional_check.py)

Two changes distinguish the mechanism. Double the physical footprint of one face while preserving its two hidden partners. It then calibrates one half-bit per aperture unit instead of one. Alternatively bias one fully observed legal word, $(1,1)$, to have twice the weight of each other word by a reversible local Metropolis acceptance. The single-end law becomes $(1/6,1/2,1/3)$. Transfer probability from $(2,2)$ to $(1,1)$. The aperture and volume readings are fixed, while $dS/d\epsilon=\log(2/3)$. A global fit to the mean number of open shutters cannot remove this signal. [Control kernels and unit tests](test_regional.py)

## Refinement and the distinction between first and finite variations

A sequence with one bit per physical face and footprints $a\to0$ would give $(\log 2)/a\to\infty$. A distinct refinement uses $d$ closed letters and $d+1$ open letters at either end, allowing every pair within the same class. Uniform reversible sampling gives each letter $d$ or $d+1$ partners. Its area gain is $\log(1+1/d)/a$. For independently surveyed $a$ and integer $d=1/(\eta a)+O(1)$ this tends to finite $\eta$ with per-face log discrepancy $O(a^2)$. With $O(\ell^2/a)$ faces, the total log discrepancy is $O(\ell^2a)=o(\ell^4)$ when $a=o(\ell^2)$. This is a separate growing local alphabet or multi-cell encoding. It controls the first-law coefficient once the proper footprints and surface cover are supplied.

Finite matter deficits face a stronger test. At fixed membership let $A=\sum_e a_e O_e$. The reference shutters $O_e$ are independent zero/one variables. For any changed law $Q$, the entropy variational inequality and the moment bound of an independent bounded sum give
\[
 D(Q\Vert P)\ \ge\ \sup_t\left\{t\,\Delta\langle A\rangle-t^2\sum_e a_e^2/8\right\}
 =2(\Delta\langle A\rangle)^2/\sum_e a_e^2 .
\]
The new law $Q$ can have arbitrary correlations. With total footprint $O(\ell^2)$ and $\max a_e=o(\ell^2)$ a nonzero finite area signal $\Delta\langle A\rangle\sim c\ell^4$ therefore requires $D/\ell^4\to\infty$. The receipt takes $\ell=1/k$, $a=\ell^4$, $d=k^4$, $k^2$ faces and a uniform increase $\ell^2/8$ in open probability. The area changes by $\ell^4/8$ as the independently selected gain tends to one. The exact lower bounds on $D/\ell^4$ are $1/2,2,8,32$ for $k=4,8,16,32$. Even without independence, an area readout with full range $W=O(\ell^2)$ gives $D\ge2(\Delta\langle A\rangle/W)^2$ by total variation and Pinsker's inequality. Such a fixed affine readout cannot supply a little-$o(\ell^4)$ finite deficit at a nonzero leading signal.

The demonstrated stationarity is therefore a first variation of a specific visible regional state. A nonlinear Einstein inference would need separate control of a changing state family and its regions. Identifying the aperture and volume counters with proper ball measurements in moving frames, and the regional modular variation with stress from a matter action, would also require independent physical measurements. The repair law alone has no such stress or motion calibration.

~~~sh
python -m model.candidate_laws.interval_echo_partition_repair.regional_check
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/test_regional.py
~~~
