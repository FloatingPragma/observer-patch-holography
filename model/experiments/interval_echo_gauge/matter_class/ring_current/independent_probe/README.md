# An independently prepared probe distinguishes two orders of meeting orbital currents

Two owners independently put one unit of load into a winding around five material ports. A third observer prepares a separate unit split between a probe port and a reference. On a single twenty-cell Interval Echo board, matched blank contacts let each source exchange load with the probe. Each owner checks whether its winding returned. The receiver compares the final probe against its reference, keeping every failed encounter in the original denominator. Both source signs, both visit orders and two pairs of marked directions were enrolled before taking counts. [Protocol](PROTOCOL.md) · [Native preparation, contacts and detector](apparatus.py)

**The receiver distinguishes at least thirteen signed response directions.** Reading 66 input/output port pairs gives twelve single-source columns and two order-difference columns. The thirteenth singular value on the fixed thirteen-dimensional comparison is **0.00212019** with joint error at most **0.00133120**. Its positive lower bound, **0.00078899**, excludes any twelve-dimensional linear family closed under these ordered signed responses on this probe. Two independently stocked winding sources thus have a detectable action on separate matter. A B14 bracket or spin representation on a particular freely traveling charge still requires measuring the corresponding response on that carrier. [Derivation and finite confidence guarantee](PROOF.md) · [Checked order identity](OrderResponse.lean)

Three particular signed port readings cancel for every single mark. Their weighted sum on the two-source order record is **0.00053490** per original start, compared with **0.00053381** predicted by the exchange contact. The simultaneous error radius for this statistic is **0.00006833**. At a shorter dwell it became **0.00030303**. Halving the exchange rate at the original dwell gave **0.00012423**, while the parked comparison gave **0.00000965**. All three held quarter-phase controls are consistent with zero. The single and ordered signal sizes scale as $g$ and $g^2$ with $g=\cos^2(ht/2)\sin^2(ht/2)$. The common contrast phase cancels in these signed reads.

The complete archive has **7,976 banks** and **2,993,592,205,312 simulated starts**, sampled as aggregate multinomial counts from the computed native probabilities. It contains **2,403,554,194,167 first or second encounter failures**. On each return the probe receiver counts bright, dark and elsewhere. The three failure bins at each encounter count an empty source, two loads in the source or a returned load in another mode. Strong inference about a small signed order effect here uses a large number of independently simulated preparations. [Raw bank counts](counts.npz) · [Count-only analysis](readback.py)

On macOS arm64 the original SciPy inverse binomial assigned one trial to a neighboring unscored outcome in each of two single-mark banks. Linux x86-64 gave the other integer. Resolving near-boundary quantiles by summing binomial probabilities about their mode corrects the two original bank records while preserving every bright and dark count. All order and control banks were identical. [Portable quantiles](sampling.py)

From the repository root run

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
p=model.experiments.interval_echo_gauge.matter_class.ring_current.independent_probe
python -B -O -W error -m "$p.readback"
python -B -O -W error -m "$p.measure"
python -B -m pytest -q -W error --pyargs "$p.test_probe"
(cd Lean && lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_gauge/matter_class/ring_current/independent_probe/OrderResponse.lean)
~~~
