# A quartic shell and one signed Yukawa line at two stations

Two Interval Echo stations use the finite scalar doublet to read a declared quartic potential and a complex doublet–singlet interaction. At the chosen lower-component vacuum the potential reading is **0.00338**, with simultaneous 99% interval **[−0.07493, 0.08183]**. A second preparation with the same single-station probabilities and perfectly correlated addresses gives **2.98895 [2.92226, 3.05692]**. On a held complex pair of doublets the oriented interaction reads **−0.396 + 0.299i**, against the forecast **−0.4 + 0.3i**. Both measurements use all launches, including the reference outcomes. [Protocol](PROTOCOL.md) · [All terminal counts](results.json)

## The field chart and its material reference

For a finite complex doublet h, prepare one unit load in the three-address state

$$
\frac{(h_u,h_d,1)}{\sqrt{1+|h_u|^2+|h_d|^2}}.
$$

The first two addresses carry the doublet of integer hypercharge 3 and the last holds the coherent reference. This chart admits arbitrary finite amplitudes and a nonzero reference probability. Two native plane rotations supply the three nonnegative magnitudes and two private phases supply the arguments. The source checks the prepared amplitudes for each exposure. The singlet uses the same chart at (χ,0). [Native scalar action and triangle pulses](../README.md)

Cell 0 holds the three ports (0,5,11). Three endpoint exchanges, each made of nearest-neighbour π contacts, map them to (156,161,167) at cell 13. The 127 contacts restore all intermediate loads. All six active ports have their own degree-five phase reference, and the three pairs at each station are incident seams. The three material basis columns have unit amplitude with unchanged references. Their permutation therefore transfers all coherent superpositions. The receiver applies its own analyzer after transfer. Each of the three terminal loads is read by two native ternary echoes against zero and double standards. The six raw clicks distinguish the three load locations and an empty receiver. [Native pulse, route and echo implementation](apparatus.py)

Each joint trial uses separate, independently initialized triangle frames. Scalar-frame readouts are made locally at cell 0 or following delivery to cell 13. Another prepared dimer supplies a gain reading. Upper and lower outcomes of its fair contact give scores 1 and 2 with mean λ = 3/2. A π/3 contact gives signed scores +1 and −1 with respective probabilities 3/4 and 1/4, hence mean y = 1/2. The two gain calibrations contained respectively 131,131 and 196,187 upper-port hits in 262,144 launches. The radius r = 1, the gains, the chart, and the pairings are apparatus choices in dimensionless pulse units.

## A fixed reader for the potential

On a sharp triangle read assign 1 to either scalar address and −r² to its reference. Write N = |hᵤ|² + |h_d|². The single-frame mean is (N−r²)/(1+N), and the reference probability is 1/(1+N). Measure one copy at cell 0 and an independent copy at cell 13. The ratio of the mean gain-weighted product to the two-reference probability is exactly

$$
\frac{\lambda[(N-r^2)/(1+N)]^2}{[1/(1+N)]^2}
=\lambda(N-r^2)^2.
$$

This is a fixed joint observable and a separate positive normalization. It works for every finite doublet on this chart, including h = 0. An empirical interval requires a resolved denominator. For positive λ and r the imported finite theorem identifies the global minimizers with the entire nonzero shell N = r². The material program chooses (0,r) as its vacuum representative. Its separate relative-phase test records all 262,144 bright X outcomes against the reference, and 130,959 bright Y outcomes. [Kernel-checked reader identities](ShellLine.lean) · [Potential and minimizer theorem](../../../../../Lean/Screen/ElectroweakBreakingComposition.lean)

The upper, lower and complex (3/5,4i/5) shell directions are all consistent with zero; the signed estimators can fluctuate on either side. The (0,1/2) reading is 0.84529 [0.80906, 0.88207], and (0,3/2) gives 2.36828 [2.05908, 2.71074]. In the correlation control, the sender measures the first address and prepares the second frame in that same sharp address, then transfers it normally. Both one-frame laws stay fixed. At r = 1 the signed product is always +1 in this arm, while the two-reference probability becomes 1/(1+N). All three shell readings under this control exceed 2.9 at the lower confidence limit.

## The oriented doublet–singlet coefficient

Choose the odd exterior sector's colourless doublet row of q = 3 and singlet row of q = −6 for ψ and χ. The scalar h has q = 3. The two Weyl rows have the same declared chirality, and the sum of the three integer charges vanishes. The weak epsilon contraction gains det U under a common two-dimensional frame change U. It is invariant for the special-unitary block. This is one specified line of the finite field grammar. The reading retains the order of ψ and h. [Exterior and Weyl component dictionary](../../exterior_register/baryon_labels/README.md)

With R the reference and j a chosen component, a native pair analyzer gives the signed population difference for each quadrature of the transition |R⟩⟨j|. The other component has score zero. Half the complex combination of the X and Y means equals zⱼ/(1+|z|²) for an input doublet z with real-positive reference. The singlet gives χ/(1+|χ|²). Measure ψ locally, h after delivery, χ locally and the separate signed gain. For each ordered term the product of its three complex means is recovered by summing all eight X/Y joint score means with factors i to the number of Y settings, divided by eight. Sum the two ordered terms with opposite signs. The three-reference coincidence probability removes the three normalizations, giving

$$
y(\psi_u h_d-\psi_d h_u)\chi.
$$

At (0,r) the line is yr ψᵤχ. The real upper input gives 0.502 − 0.001i. For the imaginary upper input the imaginary reading is 0.497 [0.404, 0.597] and the real interval contains zero. The lower input is consistent with zero. Reading χ = i rotates the line into the imaginary quadrature; removing the singlet amplitude gives zero within resolution. Exchanging the real input doublets gives −0.499 + 0.002i. Adding the two terms symmetrically on those same records gives +0.496 + 0.001i. Parallel complex doublets give zero for the alternating line; the symmetric reading has imaginary part 0.497 [0.404, 0.598]. The held complex case has intervals [−0.492, −0.307] and [0.214, 0.392] on the real and imaginary parts. Replacing the gain exposure by a fair dimer makes the declared coefficient zero and gives both quadratures consistent with zero on all eight preparations. [Acquisition and signed reconstruction](measure.py) · [Independent raw-bin and frame checks](test_shell_line.py)

Every cohort includes 262,144 joint draws from the propagated finite law. The archive stores the counts by ordered terminal address, with the corresponding six-bit material words in the transport record. The blocked comparison sets exchange to zero on every inter-owner link. The receiver then reports only its six-click empty word, including after each X or Y analyzer. All 278 blocked cohorts retain the local and gain outcomes where applicable. The 99% intervals use a joint Hoeffding allowance over at most 1,024 bounded means, with the finite sampling-grid tolerance and denominator uncertainty. The reference quotient concerns absolute joint frequencies measured on fresh cohorts. It does not condition the matter reading on successful delivery.

These experiments implement the declared quartic packet and single real Yukawa coefficient of PR-68 and PR-69 on the registered source. The quartic is a functional of classical field coordinates represented by the prepared frames. The measurements do not supply autonomous vacuum selection, a propagating Higgs, a physical fermion mass, family matrices, or a coupling scale in laboratory units.

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -m model.experiments.interval_echo_gauge.scalar_grammar.shell_line.measure
nice -n 15 .venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/scalar_grammar/shell_line
(cd Lean && nice -n 15 lake env lean -j1 ../model/experiments/interval_echo_gauge/scalar_grammar/shell_line/ShellLine.lean)
```
