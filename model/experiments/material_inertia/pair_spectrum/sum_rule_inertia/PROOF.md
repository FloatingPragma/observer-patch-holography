# Energy-weighted response and a scalar inertia test

## The tested effective description

Let L be the unweighted Laplacian of the independently specified 20-owner graph. The proposed description is Hₛ = A L + V, A > 0, with arbitrary real diagonal V. Its normalized ground wavefunction ψ can be chosen strictly positive on the connected graph. Its position records pᵤ = ψᵤ² fix ψ without an energy fit. For a real diagonal spatial grating F = diag(f), the ground-energy eigenvalue E and completeness give

~~~text
M₁(f) = Σₙ (Eₙ−E)|<n|F|ψ>|²
      = <Fψ,(Hₛ−E)Fψ>
      = A Σ_{unordered edges uv} ψᵤψᵥ(fᵤ−fᵥ)².
~~~

The last equality follows by substituting Hₛψ = Eψ and pairing each edge's two endpoints. Every on-site potential cancels. The statement permits many lines and a highly nonuniform envelope. The resting state's nodeless ground identification, a common A, the nearest-neighbor kinetic graph and the point-owner density/grating dictionary are essential hypotheses. Removing them removes this identity's application; it does not prove that another massive composite is impossible. With ħ = 1 and one owner step as the declared graph length, the conventional scalar coefficient corresponds to m = 1/(2A).

On the native pair, the coarse population read has effects Cᵤ = (1_{owner(x)=u}+1_{owner(y)=u})/2, and the coherent grating is F = Σᵤ fᵤ Cᵤ. Cᵤ are generally not projection operators on this extended object. The scalar hypothesis identifies their effective records with one point-owner density. No exact transfer of this algebra, spatial metric or scalar mass is presumed. The test challenges that proposed interpretation.

## Complete native observer read

R is the exact lowest-doublet projector and P the exact 24-line moving-band projector, the same source filters used for the recoil measurement. The initial sharp source is projected into R, including its refusal. Let g denote its normalized doublet coordinates, K = PFR, and W = H|P − E₀. A coherent two-arm clock implements the successful amplitudes

~~~text
a = exp(−iWτ) K g,
b = K exp(−i(H|R−E₀)τ) g.
~~~

Applying the signed amplitude grating before or after a common flight, followed by the final band filter, supplies these arms. The supplied Y analyzer gives probabilities ||a−ib||²/4 and ||a+ib||²/4. Their difference is Im<a,b>. The joint final-band refusal is the mean squared norm of the two complementary grating amplitudes; the grating refusal is the remaining norm deficit. All these events and the initial-source refusal form the five terminal categories. Their nonnegative probabilities are calculated from squared norms. Tiny normalization roundoff is recorded before a scalar normalization used for the seeded multinomial draws. No signed event is clipped into an acceptable law.

The occupation read uses an independent fresh launch, reads both sharp addresses and chooses the excess or vacancy label with a fair classical bit. It estimates pᵤ = <Cᵤ>. The count readback computes intervals for the population-corrected edge expression directly from those counts. Native amplitude vectors are derivation and numerical-certification tools, not a substitute observer read.

Let δ be the width of the rest doublet and w the full moving-band width. Replacing its small internal evolution in b by identity changes Im<a,b>/τ by at most δ, since ||K|| ≤ 1 and ||exp(−i(H|R−E₀)τ)−I|| ≤ δτ. The remaining expression is

~~~text
sτ = Σₙ |(Kg)ₙ|² sin((Eₙ−E₀)τ)/τ.
~~~

All gaps in this expression lie in [0,w]. For wτ < π/2, sinc(wτ) M₁ ≤ sτ ≤ M₁, with M₁ = Σₙ |(Kg)ₙ|²(Eₙ−E₀). Thus the two clock settings bound the whole first moment, not a chosen peak. The widths wτ are about 0.0507 and 0.1014 on the native band. This inequality has no phase aliases. δτ is accounted for and no individual-doublet eigenvector resolution is assumed.

For a declared isolated scalar description, all its density-transition lines belong to its finite owner Hilbert space. Native final-band losses remain explicitly measured. An effective model with external internal transitions or a renormalized density operator is a different hypothesis. Likewise, this measurement is not the full native unfiltered high-frequency first moment, whose immediate recoil has the known opposite sign.

## Numerical source and propagation bounds

The source replay recomputes the parent Schur counts, full-space Ritz frame and residual r. The old spectral separation estimates give band projector error εP and doublet projector error εR. Polar alignment changes each numerical frame by at most twice its projector error. If a is the sharp rest-filter admission, the normalized projected source coefficients differ by at most η = min(2,2εR/√a). For owner u, the amplitude of its population effect changes by at most

~~~text
bᵤ = 2εR + ||sqrt(Cᵤ)R_num|| η.
~~~

If the observed population interval is [l,u], enlarge it to [(max(0,√l−bᵤ))², min(1,(√u+bᵤ)²)]. The edge expression is monotone in each nonnegative population, so endpoint substitution bounds it. This avoids assigning a uniform absolute probability error to nearly empty sites.

The bandwidth-restricted response operator avoids amplifying the full formation energy. Polar alignment of K costs at most 2(εP+εR). The divided sine has operator norm at most w and is Lipschitz under self-adjoint generator perturbations by the Duhamel formula. The effective band-Hamiltonian error is bounded by r + 3·80 εP², using the parent shifted-generator norm bound 80. Including both K factors, the anchored energy error and conservative padding, the implementation uses

~~~text
B = 8w(εR+εP) + 4r + 240 εP².
~~~

B is about 1.88 × 10⁻⁹ native energy units. The coefficient error contributes η times the spectral width of the actual two-by-two divided-fringe response operator, since a scalar identity component cancels between normalized states. Add δ for replacing rest-doublet evolution, then the observed signed-count sampling allowance. Divide the upper end by sinc(wτ) and the lower end by one. These are finite-precision residual certificates under the parent's stated double-precision spectral estimates, not interval-arithmetic or kernel-checked spectral proofs.

An independent native 10/01 seam has h = 8.5I − σₓ/2, so its two sharp exits occur with cos²(t/2) and sin²(t/2). Three fresh count banks at commanded appointments 1/4, 1/2 and 3/4 infer its gap from 2arcsin(√p_swap)/t. The first is calibration and the others are held checks. A parked compensating excess preserves total load relative to a host reference. This independent reference response defines the native exchange unit; the protected controller and long-time repetition are supplied. The exact declared law fixes the reference gap at one; the finite records check this clock response and do not add a physical time calibration.

Each count bank enrolls 2⁵⁸ starts. Conditional on the number passing the initial filter, the Y score lies in [−1,+1] and the owner indicators in [0,1]. Hoeffding halfwidths use log(2·200·20/0.01). There are fewer than 200 banks and at most 20 relevant scores per bank, so the union bound gives at least 99% simultaneous sampling coverage for all reported comparisons. This sampling statement is separate from the numerical certificate assumptions. Neither final-band nor grating successes are conditioned out of the fringe score.

## Split inference, controls and falsifier

One native source and shape, launch owner 0 with first-shell projector column 0, calibrate A. All other launch sources, shell shapes and column 10 are held relative to that one interval. Intersect the two clock intervals for each setting. A held interval rejects the common law when it is disjoint from the calibration even after expanding that calibration by an additional 5%. The script records all outcomes, including the other second/third-shell settings. This is an exploratory split with known spectra and a deterministic pilot. It is not a retrospective preregistration.

A uniform scalar reference and a nonuniform reference ψ proportional to 1+0.2f₁ both use A = 2⁻¹⁹. For the latter, Vᵤ = −A(Lψ)ᵤ/ψᵤ fixes a ground energy of zero. Both are explicitly supplied standards; their passing held profiles validate the sum-rule reader. Their scalar coefficient is not derived from Interval Echo.

The native changed law halves the full free-flight generator while preserving alphabet, routing, symmetries and initial source. Its inferred response coefficients halve and reject the original calibration. A recorded which-arm control removes the cross-arm term with native source, grating and timing held fixed; its Y contrast is null. An independently changed held Y bank at fixed total counts fails the count-only replay. Scientific falsification is a common A passing every held native profile; the record alphabet permits it, and the scalar standards produce it. Receipt corruption is only an additional custody check.

At this finite scope the first-shell response transfers across the enrolled sources after independent population correction. The third-shell response does not. The stable neutral band therefore has reproducible kinetic susceptibility, but this particular scalar owner interpretation cannot supply its common mass. No result here prohibits an internal matrix kinetic law, nonlocal transport, a different centre or a long-wavelength mass found on refinements. OW-19's physical mass and electromagnetic dictionary remain unresolved; no criterion grade, official award or filing authority follows from this research packet.
