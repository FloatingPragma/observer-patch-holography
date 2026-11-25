# Regional entropy at finite source strength

Two executable Interval Echo sources permit an exact comparison between a change of entropy and the first-order modular-energy prediction. On the protected two-word source, the nonlinear error has the same fourth power of the routed radius as the linear term. Selecting the entropy maximum fixes the reference and removes the linear collar signal in this preparation. On the stationary angular source, relative entropy grows with the measured region while the signed entropy derivative changes sign under the same enlargement. The two observations use different native preparations. [Clock and word apparatus](../../../../../model/experiments/interval_echo_einstein_station/collar_scaling/README.md) · [Stationary angular preparation](../angular_energy/README.md)

## An all-radius result on the protected clock and word records

The protected source writes two tagged energy words in distinct, uniformly sampled cells. Each word has nine outcomes. Write \(P_-,P_0,P_+\) for the distributions with opposite, independent and common tags. These are exact integer tables divided by 2592. \(P_0\) is the product of the equal one-word marginals and \(\Delta=P_+-P_0=P_0-P_-\). Mixing the registered preparations realizes the entire family \(P_z=P_0+z\Delta\), \(-1\le z\le1\). The two individual marginals, the full joint distribution of the two ordinary word-energy sums and the protected clock tape are constant. The complete record of which source slots appear in the region includes all misses. [Native source and energy distribution](../../../../../model/experiments/interval_echo_einstein_station/shared_variance.py) · [Timing and admission law](../../../../../model/experiments/interval_echo_einstein_station/collar_scaling/clock_balls.py)

For a region \(B\), let \(q_B\) be the probability that both words appear. The placement law is independent of \(z\). Condition on the occupied slots. The only variable conditional probability table is \(P_z\) when both words appear. Consequently the regional record entropy \(S_B(z)\) and relative entropy \(D_B(z\Vert r)\) satisfy the exact identities
\[
\begin{aligned}
 D_B(z\Vert r)&=q_B D(P_z\Vert P_r),\\
 S_B(z)-S_B(r)&=q_B\big[(z-r)k(r)-D(P_z\Vert P_r)\big],\\
 k(r)&=-\sum_{ij}\Delta_{ij}\log P_{r,ij},\qquad
 S_B(z)-S_B(0)=-q_B D(P_z\Vert P_0).
\end{aligned}
\]
The last identity follows because \(\log P_0\) is a sum of one-word functions and \(\Delta\) annihilates both marginals. Independent ancestor tags give the same equations when \(q_B\) counts just the joint admissions sharing an ancestor and the record resolves the slots. All other conditional word pairs then have law \(P_0\) for every \(z\).

The nonlinear error is uniformly controlled without a lower bound on the probabilities of rare *regional* word outcomes. The conditional Fisher information
\[
 J(s)=\sum_{ij}\frac{\Delta_{ij}^2}{P_{s,ij}},\qquad
 J(0)=\frac1{324},\qquad
 \frac{229}{86580}\le J(s)\le\frac{44}{12285}
\]
comes from the source tables. The lower estimate uses the larger of the two endpoint probabilities term by term. Convexity of \(J\) gives its upper bound from the two endpoints. Since \(k'(s)=-J(s)\), the exact second-order remainder obeys
\[
 \frac{229}{173160}\,q_B(z-r)^2
 \le D_B(z\Vert r)\le
 \frac{22}{12285}\,q_B(z-r)^2 .
\]
The two constants and all three reference slopes are independently computed from the native writer. [Exact fractions and enclosures](written.py)

For the full routed balls, incidence bounds give \(q_B=\Theta(\ell^4)\) when the radius in microscopic swap units tends to infinity and its ratio \(\ell\) to the maximal station range tends to zero. This holds both for the global tag and the twenty-ancestor preparation. Thus the divided first-order error \(D_B/(|z-r|\ell^4)\) tends to zero uniformly with \(|z-r|\). At any fixed nonzero finite transfer \(D_B/\ell^4\) is bounded away from zero. Only the small-source limit eliminates that error at the retained order. [Uniform incidence and dimension bounds](../../../../../model/experiments/interval_echo_einstein_station/collar_scaling/README.md)

**The stationary reference on this family is unique.** If \(q_B>0\) then \(S_B\) has its unique maximum at \(z=0\) and both derivatives of the incident area and volume read from the protected clocks vanish identically. This simultaneously gives first-order fixed-volume generalized-entropy stationarity for every admitted region on this source variation. The measured collar derivative is then exactly zero. The clocks and ordinary energy records agree at every other reference too, but \(S_B'(r)=q_B k(r)\ne0\) there. Even a level-dependent choice \(r_\ell\) cannot preserve a nonzero normalized collar derivative while removing this discrepancy. The Fisher bound implies \(J_{\min}|r_\ell|\le |k(r_\ell)|\le J_{\max}|r_\ell|\). On the resolved full balls the stationarity residual is \(o(\ell^4)\) exactly when \(r_\ell\to0\), which also makes the normalized linear collar response vanish. The clock radius in this calculation belongs to the routed station distance. That distance has the two-dimensional packing bound in the source proof.

The unique maximum also follows for any accessible classical two-word distribution with these fixed marginals and the same placement law. The independent reference is accessible. Its entropy exceeds the joint distribution's entropy by the mutual information, which vanishes exactly at the product. One physical preparation option selects that maximum by refreshing the second tag with a fresh fair tag before writing, with probability \(a>0\) per independent attempt. This sends \(z\) to \((1-a)z\) per attempt. The selection is a property of that preparation service. The protected readout alone imposes no tag refreshing.

### Finite interventions and recorded information

Per joint admission the reversed-to-shared change has linear modular prediction \(+0.00661025507\) nats, nonlinear subtraction \(0.00649321353\) nats, and actual entropy change \(+0.00011704153\) nats. At the product reference the linear term is zero. The finite entropy deficits there are \(0.00161169089\) nats for reversed tags and \(0.00149464936\) nats for shared tags, per joint admission.

As an independent check on the finite relative-entropy cost, score the archived complete records on each of the four hosts with \(\log(P_z/P_0)\) when both words appear, scoring every other outcome zero. The selected \(P_z\) is the source's multinomial law. This is an unbiased log-likelihood read of \(q_B D(P_z\Vert P_0)\); it uses no fitted histogram-entropy estimator. On the finest annulus of the 327,680-station host, \(q_B=21243/26843463680\). The reversed and shared scores give respectively \((1.27271\pm0.00294)\times10^{-9}\) and \((1.18223\pm0.00294)\times10^{-9}\) nats per trial, with one standard error under the specified laws. The exact expectations are \(1.27544\times10^{-9}\) and \(1.18281\times10^{-9}\) nats. These are reanalyses of the deposited simulated counts. [Per-host data and relative entropies](results.json)

## A time-stationary family with a sharp information bound

The 79-word angular reference \(\sigma\) is faithful, invariant under cell rotations and stationary under the isolated native cell Hamiltonian. In its one-load block the first triplet is a single copy of a three-dimensional irreducible representation. Each of its rays has probability \(p\), with \(6.044424856454\times10^{-7}<p<6.044424856456\times10^{-7}\). The exact triplet projector is \(Q_{ij}=\mathbf n_i\cdot\mathbf n_j/4\) and \(p=\operatorname{Tr}(Q\sigma)/3\). It commutes with the native Hamiltonian, so the source's waiting-time dephasing preserves this probability. Projecting the integer bounds of the physical filter into \(Q\) encloses \(p\). [Preparation and native energy weights](../angular_energy/README.md)

Relative to a chosen pole, give its axial triplet ray probability \(p(1+x)\) and the two transverse rays probability \(p(1-x/2)\) each. All other weights stay fixed. The sources \(\sigma_x\), \(-1\le x\le2\), are stationary under the same isolated flight and have the same **full distribution** of total energy. These are physical one-load preparation lotteries in the true stationary reference. Since only commuting probabilities change,
\[
 D(\sigma_x\Vert\sigma)=p\,f(x),\quad
 f(x)=(1+x)\log(1+x)+(2-x)\log(1-x/2)
       =\tfrac34 x^2+O(x^3).
\]
Endpoint expressions use continuity. A single terminal measurement can attain this information. A sequence of thirty calibrated rotations along actual adjacent ports maps the three real orthogonal modes to three sharps. Final occupation distinguishes axial, either transverse, and all other words, including all other load sectors. Their probabilities are \(p(1+x),p(2-x),1-3p\), giving the same relative entropy as the full density matrices. The implementation verifies the two-port rotations against the native exchange and addressed phases, and samples those three outcome probabilities. [Analyzer construction and regional spectra](angular.py)

For any region or complete observation channel \(\Phi\), including a sequence of controlled measurements with every outcome retained, data processing bounds the finite entropy remainder
\[
 0\le x\,\operatorname{Tr}\big[\Phi(\dot\sigma_0)(-\log\Phi(\sigma))\big]
        -[S(\Phi(\sigma_x))-S(\Phi(\sigma))]
     =D(\Phi(\sigma_x)\Vert\Phi(\sigma))
     \le p f(x).
\]
Here \(\Phi\) is a fixed trace-preserving linear channel on states, extended linearly to traceless variations. Common subsequent native evolution and appending independent spectators are allowed. The bound applies on the output support. It is uniform also for growing observer menus. For an inferred sensitivity divided by \(\ell^4\), the generic bound on its finite-difference error tends to zero if \(p|x|/\ell^4\to0\). This source supplies no identification of \(\ell\) with a proper radius. The protected word experiment above admits the stronger bound from its measured admission probability. [Data processing for quantum channels](https://doi.org/10.1007/BF01609396)

At \(x=1/2\), the home-centered native graph balls contain 1, 6, 11 and 12 ports. Their respective relative entropies are \(6.1392\times10^{-9},2.7341\times10^{-8},8.3434\times10^{-8},1.06790\times10^{-7}\) nats, approximately 5.75%, 25.60%, 78.13% and 100% of the full value. A positive spectral sum evaluates each finite remainder directly. The nested entropy derivatives have different signs. At the single pole the occupancy \(w\le2/12\) by rotation and total load at most two. Its derivative is \(p/4\), so \(S_1'(0)=(p/4)\log[(1-w)/w]>0\). On the enclosing six-port star the certified quadrupole coefficient gives \(S_6'(0)/p\in[-0.231148,-0.230975]\). On the full cell \(S_{12}'(0)=0\) by triplet degeneracy. Monotonicity of relative entropy across nested regions therefore supplies no common sign of the entropy variation. Connecting this signed matter response to optical curvature still requires the independent spacetime measurements on the same regions. [Nesting data and counts](results.json)

Run `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.modular_orbit.relative_entropy.measure --write` from the repository root.
