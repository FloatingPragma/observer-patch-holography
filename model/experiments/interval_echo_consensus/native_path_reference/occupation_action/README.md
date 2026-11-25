# Which action does an observer measure on repeated binary histories?

A local double occupies either of two adjacent ports. For the specified contact the complete history of that double has the two-state law used in the finite translation-invariant mechanics theorem. A declared quadratic action reproduces every one of its path probabilities. An observer who counts the occupied ports across independent repetitions measures a different action between fractional occupancies. Two intermediate counts discriminate them. [Apparatus and tests fixed for the comparison](PROTOCOL.md) · [Records](results.json)

## One double

Prepare a sharp double on port 2 or 3, a vacancy on port 4 and blank memories. Isolate the (2,3) contact. To read the position apply the registered echo to (3,4), publish the agreement bit and undo the echo. Agreement certifies the double is on port 3. Each branch restores that sharp position, so the same system can supply the preceding bit for the next contact. The runner checks the invariant two-column subspace of the full 99-state generator and both pointer branches. After subtracting its scalar diagonal the restricted generator has off-diagonal entries −1/2 and zero diagonal. The probability of a position change at duration t is therefore sin²(t/2). The Born instrument, isolated appointments and independent preparations are supplied by the construction. [Native generator and pointer binding](measure.py) · [Registered clock primitives](../../../../candidate_laws/interval_echo/clock.py)

A half contact from a sharp zero prepares the uniform initial law. Select the same subsequent duration for either input with cos(t)=1/3. The transition matrix has stay probability 2/3 and flip probability 1/3. This angle was chosen to realize the specified mathematical instance. The four transitions after the initial read yield 32 complete five-bit words. A separate duration of π/3 gives flip probability 1/4. Suppressing exchange in the generator preserves the two rails and their pointer. Every launch then keeps its starting position.

The positive asymmetric historical kernel also has a native implementation. Read the last bit and use the corresponding appointment for flip probabilities 107/1431 and 5/508. A contact from zero prepares the historical initial one-probability 54356/61511. The durations encode those values from the historical chain. The resulting source has the same finite path law for each depth by the sharp branch update and the Markov product. Its apparatus is the scheduled port pair. [Historical kernel and initial law](../../../../../Lean/Thermodynamics/MixingChainRealization.lean) · [All acquired five-bit banks](results.json)

## The real action supplied to the mechanics theorem

For the symmetric chain put c=log(3/2) and s=log 2. The bilinear interpolation of the four negative log transitions is

$$
L_0(x,y)=c+s(x+y-2xy).
$$

On real x and y declare the velocity-only enrichment grammar of the finite theorem, adding a y(y−1)/2, and require the constant-one modal path to be stationary at a real interior junction. The derivative of the two adjacent actions at that junction is a/2−2s. Thus a=4s. Constant zero gives the same value. For this choice

$$
L_a(x,y)=c+s(y-x)^2+F(y)-F(x),\qquad F(x)=s x(x-1).
$$

The difference telescopes on every fixed-endpoint path. The translation-invariant representative c+s(y−x)² is exactly the one in the Lean theorem. The variation gap of either constant path at one interior junction is 2s times the square of the replacement displacement. Translation gives conserved momentum 2s(y−x), zero on the realized constant mode. For the asymmetric replay the declared constant-one stationarity convention gives the original theorem's curvature 2 log(W₁₁²/(W₁₀W₀₁)) = 14.38912185. The stated grammar and stationarity convention are representation choices. All values of a agree on the binary histories. [Boundary-term and held-curvature identities in Lean](NativeOccupationAction.lean) · [Composed source and Hamilton theorem](../../../../../Lean/Variational/SourceToHamiltonianComposed.lean) · [Symmetric instance](../../../../../Lean/Variational/TranslationInvariantComposedInstance.lean) · [Original mode-extremal formula](../../../../../Lean/Variational/ModeExtremalEnrichment.lean)

## The action from occupancy readings

Prepare n independent copies with m initial ones. Each observer reads one port comparison per copy and counts the terminal ones k. The probabilities below have the full n-copy bank as their denominator. An independent microhistory specifies which labeled copies changed, whereas a count event includes every compatible labeling. Write p₀ and p₁ for the respective one-probabilities from an initial zero and an initial one. The generating polynomial of the count is exactly

$$
(1-p_0+p_0z)^{n-m}(1-p_1+p_1z)^m.
$$

The coefficient of zᵏ gives the same transition probability for every initial arrangement with m ones. The count process is therefore Markov even if earlier appointments recorded the arrangements. The contact duration and the conditional independence of the copies stay fixed. A group of size n may also be acquired by n independent relaunches of the same apparatus.

For m/n tending to x and k/n tending to y, the negative log of this coefficient divided by n tends to I(x,y). Define

$$
H_x(p)=(1-x)\log(1-p_0+p_0e^p)+x\log(1-p_1+p_1e^p),
\qquad I(x,y)=\sup_{p\in\mathbb R}\{py-H_x(p)\}.
$$

Here is an elementary convergence proof for every faithful binary kernel. With row masses s₀=1−x and s₁=x, minimize the cost of a two-by-two transition table t with those row masses and terminal one mass y,

$$
J(t)=\sum_{ij} t_{ij}\log\frac{t_{ij}}{s_i W_{ij}}.
$$

Contributions from empty rows vanish. For each integer table at finite n, the two independent binomial type probabilities lie between exp(−nJ)/(n+1)² and exp(−nJ). For a single row these bounds follow by tilting its Bernoulli probability to the observed success fraction, where that observed count is a mode among at most n+1 outcomes. There are at most n+1 compatible tables. At given margins a table is specified by one entry in a closed interval whose endpoints are maxima or minima of affine functions of the margins. The minimizing real table at the limiting margins can be approximated by clamping and rounding that entry at the finite margins. J is continuous on the closed probability simplex, as follows by expanding it into joint entropy, row entropy and a term linear in t. Compactness bounds the costs of every convergent sequence of integer minimizers from below. Thus the integer costs converge to the real minimum, and the two polynomial factors disappear in the action per copy.

Tilting both rows by the same factor exp(p) yields their one-probabilities qᵢ(p)=pᵢ exp(p)/(1−pᵢ+pᵢ exp(p)). The binary relative entropy to each tilted row splits J into py−Hₓ(p) and a nonnegative remainder. For 0<y<1 there is exactly one p with y=(1−x)q₀(p)+xq₁(p), since Hₓ″ is the weighted sum of the two positive tilted variances. The tilted table sets the remainder to zero. The endpoint limits give the same formula at y=0 and 1. This identifies a unique continuous occupancy action on the closed square. Its momentum ranges over the real line for each initial fraction, with terminal fraction confined to (0,1). For the displacement variable v=y−x the conjugate function is Hₓ(p)−px. The minimizing terminal fraction at x is (1−x)p₀+xp₁. For the symmetric experiment it equals (1+x)/3, which gives a contraction to 1/2. An arbitrary real translation of both fractions changes that typical trajectory.

### Exact coordinate symmetries of the measured action

For any faithful binary kernel, an increasing bijection T of the fraction interval that satisfies I(Tx,Ty)=I(x,y) for every x,y is the identity. It fixes zero. On the initial-zero edge the cost I(0,y)=D(y‖p₀) has its unique minimum at p₀, so T fixes p₀ as well. The derivative there changes sign according to log[y(1−p₀)/(p₀(1−y))]. The cost is strictly decreasing on the left and strictly increasing on the right. Since T preserves each side and leaves the cost invariant, it fixes every y. The symmetric kernel has in addition the reflection y↦1−y by exchanging both rail labels. Composing any decreasing symmetry with this reflection gives an increasing one, so the reflection is the only decreasing symmetry of that kernel. Each element of a one-parameter group of interval homeomorphisms is the square of the element at half the parameter, hence is increasing. Such a group acts trivially if it preserves this cost exactly. This classifies point transformations of this specific fraction coordinate with exact invariance. Symmetries up to endpoint terms require a separate condition.

## A held comparison between intermediate readings

On the initial-zero edge, I(0,y) equals the bilinear microscopic cost L₀(0,y) minus the binary entropy h(y). At finite n the exact correction is n⁻¹ log C(n,ny). Every velocity-only quadratic enrichment predicts that the correction at 1/4 is 3/4 of the correction at 1/2, independently of a. Their difference in the occupancy action tends to

$$
h(1/4)-(3/4)h(1/2)=\tfrac14\log(32/27)>0.
$$

At n=64 the two finite combinatorial terms give a difference of 0.0356751 nats per copy. The exact sign also follows from C(64,16)⁴>C(64,32)³. In 262,144 simulated banks the observer counted 10,587 occurrences of 16 ones and 652 of 32 ones. With the independently certified two-port probabilities, these yield **0.0348508 nats per copy**, with simultaneous 99% interval **[0.0331560, 0.0365844]**. The interval uses two binomial Clopper–Pearson intervals with total error at most 1%. It treats native pulse probabilities as exact properties of the specified law. Reusing one draw across 64 copies gives the same single-copy marginals while producing zero interior counts. Suppressing exchange also produced zero interior counts from the sharp-zero start. [Exact convolution, independent sampler and intervals](measure.py) · [Checks](test_measure.py)

The result supplies the declared real enrichment of PR-06 for these enrolled binary chains. The fraction observable supplies a second measured response across real intermediate values via independent-copy refinement. Its action includes the multiplicity of compatible microscopic histories. Translation and the quadratic continuation in the mechanics theorem describe a choice on the individual record alphabet. The occupancy process supplies neither physical translation symmetry for the coarse dynamics nor motion of a classical body in calibrated spacetime.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_consensus.native_path_reference.occupation_action.measure
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/native_path_reference/occupation_action
(cd Lean && lake env lean ../model/experiments/interval_echo_consensus/native_path_reference/occupation_action/NativeOccupationAction.lean)
```
