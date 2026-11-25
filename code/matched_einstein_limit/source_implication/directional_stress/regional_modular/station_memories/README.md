# Entropy of the four terminal loads and both head memories

At the volume-matched reused-radar equilibrium, a native preparation changes the entropy of the four terminal loads and both head memories even though it leaves the tetrahedral area, volume and all four-port sector weights fixed. The joint entropy includes the actual correlations between loads and memories. Enlarging the region to the entire closed flight also gives a nonzero signed response. [Preparation with 26 fixed readouts](../station_control.json) · [Four-port load entropy](../STATION_REGION.md)

| Entropy in nats | Derivative in the same source direction | Change from $s=0$ to $s=0.01$ |
|---|---:|---:|
| Four station loads $R$ | $[0.32465,1.36457]\times10^{-6}$ | $[-4.8150,-3.7750]\times10^{-8}$ |
| Four station loads and both memories $RH$ | $[6.9964,67.161]\times10^{-6}$ | less than $-1.6035\times10^{-4}$ |
| Both memories conditional on the four loads, $S(RH)-S(R)$ | greater than $5.63\times10^{-6}$ | less than $-1.6030\times10^{-4}$ |
| All loads and both memories | $[14.141,49.269]\times10^{-6}$ | less than $-7.4047\times10^{-4}$ |

The positive derivatives are reference-point statements. Entropy is strictly concave along these accessible preparation chords, and the quoted finite displacement is large enough to reduce the entropy below its reference value. The chord holds the twelve registered charge-pooled ranges, both four-port number weights, and the other material constraints exactly fixed. Its area and volume are the independently inferred four-station tetrahedral quantities. This addresses first-order stationarity and the finite entropy change for these specified finite algebras. A physical shrinking-ball proper area, volume and continuum generalized entropy require their own identification. [Joint certificate](joint.json) · [Closed-flight certificate](total.json)

## The joint channel

Prepare the same two-parity, two-coin head state as the 26-fixed control. Independently choose the probe charge fairly and one of ten distinct emitter ports uniformly. Propagate for four exchange units, without monitoring the mobile during this interval. Keep the four terminal load registers at 0, 1, 38 and 39 and the two original head memories. The right reversal sector is fixed by preparation, and the retained left parity is one of the physical head labels. Trace all other loads. The global charge-choice label is erased. No measurement outcome is discarded.

Within a given left parity, for charge $q$ and launch $l$, let $U_{q,l}$ be the flight isometry on the four-dimensional incoming head. The part ending at one of the four stations is a single $16\times4$ map $V_{q,l}$. Keeping all sixteen output indices preserves the position-memory coherences inside the region. If the mobile ends elsewhere at $x$, the other loads distinguish $x$ when traced. Denote the associated $4\times4$ head map by $W_{q,l,x}$. The two charges share the same terminal word when the mobile is outside, and have orthogonal number words when it is inside. For a conditional head $X$ the precise output is the direct sum of
$$
\frac1{20}\sum_l V_{-,l}X V_{-,l}^*,\qquad
\frac1{20}\sum_l V_{+,l}X V_{+,l}^*,\qquad
\frac1{20}\sum_{q,l,\,x\notin R}W_{q,l,x}X W_{q,l,x}^* .
$$
The outside block is a full four-by-four memory state. Summing the maps over one parity is trace preserving. Keeping the two orthogonal left parities gives two copies of the $16+16+4$ output algebra with their physical parity weights. [Port and head generator](../../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/ranging/measure.py) · [Integer channel contraction](channel.py)

Let $w_p$ be the weight of parity $p$, $X_p$ its normalized reference and $D_p$ its unnormalized source derivative. Both $D_p$ are traceless on this control. Writing $T_p$ for the displayed three-block channel, the signed entropy derivative is $-\sum_p\mathrm{Tr}[T_p(D_p)\log T_p(X_p)]$. The parity term $-\log w_p$ cancels exactly. An independently propagated nominal preparation gives $3.70785037\times10^{-5}$ nats per unit source. The certified interval allows the exact stationary states and the corrected physical direction to vary throughout their independent enclosures. [Entropy bound](joint.py) · [Independent sparse propagation](test_memories.py)

The full incoming head is identifiable from paired relay and return timestamps on independent copies. Applying the calibrated channel to that reconstruction also predicts the joint terminal state and its entropy. This is an indirect finite readout with counting error governed by the chosen clock inversion. The numerical comparison in this experiment evaluates population probabilities. [Single-edge tomography](../../../full_algebra/reuse/edge_clock.md)

The native integer Chebyshev recursion retains the four terminal row blocks at the original $2^{-100}$ scale. For the sum over all outside locations it also retains full output columns at $2^{-52}$. Two 26-bit limbs implement each fixed-point product without floating-point arithmetic. The extra truncation error in one propagated column is at most $6(641)\sqrt{3840}/2^{52}$. The compression for each 4-dimensional launch is within $1.2\times10^{-10}$ in operator norm, giving an induced trace-norm error below $2.5\times10^{-10}$ for the outside channel. Hermitian input matrices of trace norm one obey this bound by decomposing them into positive and negative parts. For the terminal channel the higher precision gives a bound below $10^{-18}$. The resulting integer Gram matrices retain the full quantum channel on the prepared head.

To bound the blockwise modular logarithm, the verifier checks the eigenvector Gram defect, reconstructed-state residual and logarithms of the proposed eigenvalues by integer and interval arithmetic. The reference uncertainty can be bounded at the four-dimensional input. If $X_p$ varies from its center by a traceless Hermitian matrix of Hilbert–Schmidt norm at most $r$, its operator norm is at most $\sqrt{3/4}\,r$. Putting $\epsilon=\sqrt{3/4}\,r/\lambda_{\min}(X_{p,0})$ gives $(1-\epsilon)X_{p,0}\le X_p\le(1+\epsilon)X_{p,0}$. Positivity of the channel preserves both inequalities. Operator monotonicity of the logarithm bounds its variation between $\log(1-\epsilon)I$ and $\log(1+\epsilon)I$ in every output block. Subtracting their common midpoint preserves the pairing with the traceless whole-parity output derivative. The tangent uncertainty is bounded against the measured traceless input score obtained by contracting the modular operators back through the three channels.

On the full feasible interval, eigenvalues in each parity block are at most $w_p$. The divided differences of $\log$ are therefore at least $1/w_p$ there. Along the physical chord this gives $-S''\ge\sum_p\|T_p(D_p)\|_{\rm HS}^2/w_p>3.2204$ in squared displacement units, establishing the finite decrease.

## The entire flight

For one known charge the incoming state is the head tensor a uniform classical mixture on ten orthogonal launch positions. The charge-unresolved initial state has twenty orthogonal charge-position labels. The charged sectors are distinguished by the global total load. Evolution on all positions and both memories is unitary on the enrolled sector and preserves the full spectrum. Thus the closed-flight entropy is exactly $S(\rho_{\rm head})+\log10$ for a known charge, or $S(\rho_{\rm head})+\log20$ for the fair preparation, at every dwell. Copying these independent orthogonal labels to extra records leaves the same entropy difference. A purification of the mixed head or a record of the upstream head-preparation lottery would specify different initial information. [Exact prior state and source control](../station_control.py) · [Head score and finite change](total.py)

There is a further statement for any affine extension of the output retaining the four loads. The midpoint entropy deficit $S(\sigma_0)-[S(\sigma_+)+S(\sigma_-)]/2$ equals the mean of $D(\sigma_+\Vert\sigma_0)$ and $D(\sigma_-\Vert\sigma_0)$ for opposite equally spaced endpoints. Partial trace cannot increase relative entropy. The explicit four-load midpoint lower bound $3.15266\times10^{-9}$ nats therefore also holds after retaining any additional memories from a source-independent quantum channel. This midpoint inequality alone gives no reference tangent sign. [Four-load quadrature witness](../station_region.json)

The regression verifies both archived SHA-256 hashes, recomputes the interval certificates and checks the quoted signed inequalities. An independent propagation uses 64 short Taylor steps per flight, checks the Gram matrix of all 40 launched vectors, and recovers the nominal joint entropy slope within $10^{-10}$. This check uses fresh conditional density matrices from the full port basis. [Verification](test_memories.py)

Run the replay and regression from the repository root.

    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.directional_stress.regional_modular.station_memories.joint
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/regional_modular/station_memories
