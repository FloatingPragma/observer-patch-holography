# Reusing the two heads after a complete radar cycle

The selected 80-cell entropy equilibrium changes when its heads undergo the six-unit survey. Each of the twelve directed four-station routes, for either probe charge, decreases the mean of the same left-head Pauli $X$. The exact decrease exceeds **0.0005192** in every setting. A fixed probability mixture of these 24 settings must also change the state. The result applies when the launched probe is discarded and the two heads are retained, with every timeout included. Surveying fresh independently prepared copies measures the same range means as in the finite equilibrium theorem. [Preparation and survey](../PROOF.md) · [Channel and complete counts](results.json)

## The actual post-survey state

Label the source and relay stations $i,j$. Write $U$ for one ordinary half-tick of free flight in the specified charge and left-parity sector, and $P_j$ for position at station $j$. Start the probe at $i$ with its four-dimensional coin input untouched. At each appointment the piece $P_jU$ gives a first hit; $(I-P_j)U$ continues the survivor. A hit leaves the probe at $j$, where the return leg starts on the same heads. A return hit stops the experiment. At tick 12 the experiment stops on every branch, including an outbound hit at exactly tick 12 and any surviving return attempt. The actual source has $H=h/2$ modulo a scalar, where $h$ is the sparse integer sector Hamiltonian. Each half-tick is $e^{-ih/4}$ up to a global phase. [Native generator](../native_precision.py) · [Registered first-hit clock](../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/ranging/measure.py)

Let $F_a^{ij}$ be the first-hit coin map at outward tick $a$, and $N_s^{ij}$ the no-hit coin channel through $s$ ticks after tracing out the final probe position. In this formula $F_a$ denotes the completely positive map with its one position-resolved Kraus matrix. The unconditional head map is

$$
\Phi_{ij} =
 N_{12}^{ij}+
 \sum_{a=1}^{12}
 \left(N_{12-a}^{ji}+\sum_{b=1}^{12-a}F_b^{ji}\right)\circ F_a^{ij},
 \qquad N_0=I.
$$

There are 79 disjoint terminal records before resolving the last probe position. Each tick splits into complementary projectors and each free hold is unitary, so the exact sum is trace preserving. All operations preserve left reversal parity. Expectations of the parity-blind $X\otimes I$ sum over the two unnormalized parity blocks, with the selected incoming parity weights. Position is traced only after its actual registered history. The calculation assigns the two signs separate channels before taking any lottery over them.

## An integer error bound on the complete channel

The verifier propagates all four head input columns on the actual 3,840-dimensional position and coin space. Each half-tick uses eight substeps $e^{-ih/32}$. It evaluates 22 Taylor terms using signed integer sparse products, rounding down each real and imaginary component after the division by $32k$. All terms have denominator $Q=2^{48}$. The row-sum norm of $h$ is at most 43, and $h$ is real symmetric. Write $z=43/32$ and $d=3840$.

The rounding error at one division is at most $\sqrt{2d}/Q$ in each complex column. Along the Taylor recurrence the successive errors from a given division accumulate with products $z^r k!/(k+r)!\le z^r/r!$. The sum of all term errors in one substep is consequently at most $22e^z\sqrt{2d}/Q$. The exponential remainder on a numerical input column of norm at most two is at most $2z^{23}/23!$ by the integral Taylor formula for a unitary flow. Exact substeps are isometries and survivor or hit projections are contractions. The cumulative column error over 96 substeps is bounded by

$$
b=96\left(22\cdot4\cdot88/Q+2z^{23}/23!\right).
$$

Here $e^z<4$ by $(1-z/32)^{-32}<4$. The resulting bound $b<10^{-6}$ also proves the assumed numerical column norm bound by induction from a unit input. The same estimates put every scaled Taylor power and partial sum below $9Q$ in magnitude. A sparse integer dot product is smaller than $43\cdot9Q<2^{62}$ in magnitude, including partial sums. Four input columns give single-leg operator error $\epsilon=2b$; a composed outward and return record, including an unprojected terminal probe, has error at most $\delta=2\epsilon+\epsilon^2$.

For a true terminal Kraus operator $K$ from coin to probe and coin, $\|K\|_{\rm op}\le1$. An operator approximation within $\delta$ changes $KXK^\dagger$ in trace norm by at most $(2\delta+\delta^2)\|X\|_1$. Partial trace is trace-norm contractive on Hermitian matrices. Summing the 79 records bounds the full channel error uniformly for every Hermitian input. The Kraus Gram sums and channel compositions use exact integer arithmetic. Each amplitude integer below $2Q$ is split in base $2^{24}$; each constituent dot product contains at most 960 positions and lies within signed 64-bit range. The full integer channel has denominator $Q^4$. The implementation checks that its real and imaginary coefficients have magnitude at most two. Final conversion to binary64 introduces less than $10^{-13}$ in induced Hermitian trace norm, using $\|Y\|_1\le2\|Y\|_2$ on a four-dimensional output. The uniform channel error including that conversion is below $1.670\times10^{-6}$. [Integer calculation](measure.py)

## A single head measurement distinguishes reuse

In a sham arm, prepare the selected state and read the left coin along $X$. In an independent arm, perform a complete monitored radar route and read $X$ once on the surviving heads. Both arms retain every preparation. This tests disturbance without measuring the incoming coin on the same copy. The target-state certificate places the true equilibrium within $2\times10^{-7}$ in Hilbert–Schmidt norm of its normalized Hermitian center. The corresponding trace bound is $\sqrt8$ times as large. Since $\|\Phi^\dagger(X\otimes I)-X\otimes I\|_{\rm op}\le 2$, this adds at most $1.132\times10^{-6}$ to the error in the signed change. The interval also includes $10^{-11}$ for the small final matrix sums. [Equilibrium certificate](../certificate.json)

The least disturbance at the center occurs for the negative probe on directed station indices $3\to1$, where the mean changes by $-0.00052206806$. The largest decrease in magnitude is $0.04595635$. At the true equilibrium every individual route and charge retains the strictly negative sign. The receipt reports all 24 mean changes and separately simulated binary counts for the sham and each reuse arm. Each setting has $2^{31}$ simulated starts. The two-outcome $X$ meter has plus rate $(1+\langle X\rangle)/2$. Hoeffding's inequality and a union bound give a simultaneous 99% rate radius across the 25 reported settings. The worst-case rigorous bound above is computed from the native channel and the state certificate. [Readout records](results.json) · [Regression tests](test_reuse.py)

On equal, independent reuse of all 24 settings, each sector has a unique attracting state. Its volume-matched mixture admits a correlation change at fixed area, volume and both one-head marginals with a nonzero entropy derivative. [Long-run channel and proof](../PROOF.md#what-repeated-use-of-the-radar-selects)

The same clock supplies 31 pre-cutoff bin probabilities that determine the entire incoming head state at population precision. The two independently calibrated force scores have a separate finite-count check. [Full-state readout and bounds](CLOCK_HISTORY.md)

From the repository root, replay the calculation with

    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.full_algebra.reuse.measure
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/reuse
