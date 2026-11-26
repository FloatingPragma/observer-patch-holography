# Does the thirty-direction spatial symbol measure energy or frequency?

Two instruments read the same directional dependence in the mean kinetic energy of a phase-patterned load. One counts antisymmetric pairs on the thirty seams. The other clocks the load's free return with every seam open. The result holds on an isolated icosahedral cell and on either tested cell of the 20-cell federation. The corresponding transverse orbital on the isolated cell has two fixed clock lines. Changing the grating transfers weight between the lines. [Energy and frequency protocol](PROTOCOL.md) · [Polarization protocol](POLARIZATION_PROTOCOL.md) · [All energy and return counts](records.json) · [All polarized receiver counts](transverse.json)

## Two independent ways to read the spatial response

The source carries one double on an otherwise empty board with blank seam memories. Equip its twelve launch ports with the registered icosahedral frame, choosing a native-seam isomorphism, and express the port coordinates $x_j$ in units of the common chart edge length. An equal-amplitude launch with phase $q n\cdot x_j$ at port $j$ has the same port populations for every gradient $q$ and direction $n$. Its amplitudes are prepared by progressively splitting a localized load along an internal path of eleven native exchanges. Private quadratic-load exposures set the relative phases. The transverse preparation uses the same method with a signed dipole amplitude envelope. The preparation program checks every resulting complex amplitude. Angles, reference phases and the chart are supplied by the apparatus. The dynamics during free flight is the candidate's native generator. [Source and clock](measure.py)

Restricted to the one-double sector, the blank two-port native operator has equal diagonal contrast and hopping $-g/2$, where $g=1$ is the symmetric-to-antisymmetric exchange gap. A separate timed contact measures transfer from a sharp rail at two short holds. Summing the unrestricted local columns over the board and applying its registered degree compensation gives $24\kappa I+(g/2)L$, with $L$ the unweighted graph Laplacian. A separately prepared uniform state of the whole board is the phase standard at energy $24\kappa$. The Ramsey comparator opens the material generator on both arms, inverts their respective preparations and recombines against that standard. Its two bright/dark quadratures give the complex return of the prepared state. All branch controls and coherent comparators are supplied.

For the seam joining $i,j$, the lower/upper sign decoder measures the symmetric and antisymmetric dimer modes. It also records the load outside the pair. The probability of the antisymmetric exit is $|\psi_i-\psi_j|^2/2$. With the bare chart rays $w_e=x_i-x_j$ of unit length, the excess energy above the unmodulated local preparation is exactly

\[
 \Delta E=\frac{g}{12}\sum_{e=1}^{30}(1-\cos(q n\cdot w_e))
         =\frac{5g}{12}\,\Lambda_1(qn).
\]

Here $\Lambda_1$ is the normalized dilated completion symbol at chart scale one. Each external edge has one zero amplitude initially, so its contribution to the mean energy is independent of the grating phase. This is why the same thirty internal terms survive on the all-open board. The identity concerns a Rayleigh quotient, with no requirement that the preparation be an eigenmode. A phase-origin shift also preserves the quotient. Homogeneous physical translations or a laboratory unit of length require separate attachments. [Exact spatial symbol](../../../../../../Lean/Screen/SeamCurrentDirichletGenerator.lean) · [Mode and support checks](test_symbol.py)

The free clock measures $s(t)=-\operatorname{Im}\langle\psi|\exp(-it(H-E_{\rm ref}))|\psi\rangle$. The estimate $[8s(h)-s(2h)]/(6h)$, $h=1/32$, uses only independent clock counts. The centered five-point derivative bound gives an energy error at most $h^4(6g)^5/30$ per preparation because the native relative generator has spectrum in $[0,6g]$. The initial uniform-cell preparation has zero relative energy on the isolated cell and $3g/24$ on either embedded owner; subtraction uses fresh zero-gradient reads on each board.

At four radians per chart edge, the seam detector gave **2.894590, 2.997730 and 3.015803** native energy increments along a port, a seam midpoint and a face centre. The port-axis free-clock increments were **2.894535** on the isolated cell and **2.895125, 2.893936** on embedded owners 0 and 7. Across all 60 native angular and radial clock comparisons the largest discrepancy with the local seam reads was **0.001029**. The simultaneous separate allowances are 0.000603 for the seam sum and 0.004351 for the baseline-subtracted free clock relative to the exact predictions. A twelve-port cosine stencil with the same quadratic and quartic moments predicts **3.223048** on the port axis. It misses the held four-radian clock results along all four directions, including the oblique one. Tripling kinetic exchange during flight with the same phase sources and detectors changes the measured energy increments by the predicted factor of three, within their larger derivative allowance.

With randomly interleaved dark, bright and half-bright standards the face-minus-port energy difference is resolved even when the receiver response depends on earlier trials. The simultaneous 99% interval is 0.0982 to 0.1431 exchange-gap units. [Randomized receiver calibration](drifting_receiver/README.md)

## Resolving the transverse clocks

For propagation along a port axis the two polar port amplitudes of a perpendicular dipole vanish. The other ten vertex projections on that axis are $\pm c$, with $c=1/(\sqrt5\,d)$ and $d=\sqrt{2-2/\sqrt5}$ the edge length on the unit sphere. Let $D$ be the unit transverse dipole and $Q$ that dipole with the lower latitude's sign reversed. These are orthogonal, independent of $q$. On the isolated cell $LD=(5-\sqrt5)D$ and $LQ=6Q$. For example, each nonpolar port has two neighbors on its latitude and two on the opposite latitude, plus a zero-amplitude pole. The same-latitude pair contributes $2\cos(2\pi/5)D_j$ to the adjacency action. The opposite pair contributes $2\cos(\pi/5)D_j$. Reversing the second latitude switches the sign of the latter term. Both pole sums vanish by pentagonal symmetry. Thus the return splits for every gradient:

\[
 \psi_q=\cos(qc)D+i\sin(qc)Q,\qquad
 \omega_D=(5-\sqrt5)/2,\quad \omega_Q=3 .
\]

Separate coherent reference arms prepare $D$ and $Q$ and read each complex transfer at zero, one sixth and one half exchange unit. The initial magnitude gives the weight; the phase increment gives the frequency. With gradient one the second line carries **0.170254** of the weight. At gradient four it carries **0.983063**. The measured frequencies were **1.381972 and 3.000001**, then **1.381975 and 3.000041**. The second line is unresolved on the unmodulated input because its template has zero overlap. Both controls use the identical holds. The one-sixth-unit pilot is unambiguous even for the tripled exchange, and resolves the integer winding of the half-unit result. [Template and clock implementation](transverse.py)

A separate test uses no spectral templates. The first three complex return lags generate the four-by-four Hermitian matrix $T_{ab}=g_{b-a}$. Each distinct sampled frequency can contribute at most one rank, since $T$ is a positive sum of outer products of phase progressions. A single oscillator, allowing both signs of frequency and a common energy shift, gives rank at most two. In the three-radian transverse launch, the two embedded owners gave smallest eigenvalues **0.002036 and 0.002036**, resolving at least four distinct line phases. On the isolated cell precisely the two transverse lines carry weight. The isolated uniform and transverse zero-gradient controls each resolve a single occupied line, as their independently known generators predict.

Each detector bank consists of $2^{34}$ simulated ideal launches compressed into bin counts. The 99% simultaneous bound for each signed mean uses Hoeffding over fewer than 5,000 detector probabilities and clock quadratures, with an extra $2/2^{25}$ for dyadic probability rounding. The quadrature error disk radius is $\sqrt2$ times that bound; each row of $T$ has at most three uncertain off-diagonals. Thus every eigenvalue's allowance is below **0.000171**. The polarized receiver uncertainty propagates the complex disks at both times through the relative phase and keeps a zero-amplitude response unresolved.

The same thirty-direction spatial symbol is measurable as the energy of this phase preparation. The clock data exclude treating its kinetic energy response as one oscillation at the square root of this symbol. The finite-aperture result leaves open whether a different extended transverse field realizes the oscillator in the completion theorem. No physical photon frequency, cofinal wave branch, SI scale or high-energy observational test follows from the measurement.

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_29_maxwell.seam_symbol
.venv/bin/python -B -W error -m "$pkg.measure"
.venv/bin/python -B -W error -m "$pkg.transverse"
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_29_maxwell/seam_symbol
```
