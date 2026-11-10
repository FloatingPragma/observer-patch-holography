# Connected finite defect-wave dynamics

## Source, state and clock

Each cell carries the complete Hilbert space of the owning
[joint defect-wave law](../PROOF.md): ten actual link words and a
`12 × 5 × 2 × 60` wave register. Its dimension is 72,000. Its unitary is
the literal `JointObserver.step`, denoted by `U`. The local link orbit,
protected boundary, apparatus and [interior constraint](../GAUSS.md) retain
their stated domains. One local constraint is imposed at vertex 5; this is
not a construction of every boundary gauge charge.

The spatial architecture is a declared cubic box, or an even periodic cube
of side at least four. It is independent of any assertion that a spherical
carrier forces three spatial dimensions. An onsite layer applies `U` in
every cell. Six edge layers follow, ordered by axis `0,1,2`, with parity
`0,1` for each axis. A layer joins `x` to `x+e_axis` when the axis coordinate
of `x` has the stated parity. Open boxes omit outgoing boundary edges.
Every layer is a matching. One complete seven-layer cycle is the clock tick.

The Hilbert space of a box is a finite tensor product. The propagator
contains gates joining neighboring factors. It cannot be replaced by a
tensor product of isolated-cell propagators.

## Gauge-covariant intercell exchange

Write a local basis vector as `|g,w>` with
`w=(port,internal,parity,address)`. Let `s=6553634` be the interior star
mask, and let `χ(w)` be minus one exactly when `port=5` and `internal≥3`.
The local constraint is

\[
\Gamma|g,w\rangle=\chi(w)|g\mathbin{\mathrm{xor}}s,w\rangle.
\]

Define the diagonal frame `F` by a minus sign exactly when
`g > g xor s`, `port=5`, and `internal≥3`. For the partner link word the
ordering reverses, so

\[
F\Gamma F=J,\qquad J|g,w\rangle=|g\mathbin{\mathrm{xor}}s,w\rangle.
\]

For neighboring cells define

\[
S_{xy}=F_xF_y\,\operatorname{SWAP}_{w_x,w_y}\,F_xF_y,
\qquad
E_{xy}=aI+bS_{xy},\quad
a=\frac{1-i}{2},\quad b=\frac{1+i}{2}.
\]

The SWAP leaves both link words at their original cells. In the frame `F`
it commutes with each pure link flip `J`. Consequently `S` and `E` commute
with `Γ_x` and `Γ_y` separately. Since `S=S†`, `S²=I`, `a²+b²=0`,
`2ab=1`, and \(\bar a b+\bar b a=0\), one has

\[
E^\dagger E=I,\qquad E^2=S.
\]

The coefficients specify an exact square-root SWAP gate. No measured mass,
coupling or velocity selects them. This gate and the cubic placement extend
the local program; neither is asserted to follow from the isolated `U`.
The onsite covariance `[U,Γ]=0` and the edge covariance prove every
interior constraint is conserved by the complete spatial circuit.

Let `Z` equal `2` on the triplet and `−3` on the doublet. Both `F` and `U`
commute with it. The edge exchange obeys

\[
[E_{xy},Z_x+Z_y]=0,\qquad
J_{x\to y}=Z_x-E_{xy}^\dagger Z_xE_{xy}
=\frac{Z_x-Z_y}{2}-\frac{i}{2}[Z_x,S_{xy}].
\]

Thus `ΔZ_x=−J_{x→y}` and `ΔZ_y=J_{x→y}`. Summing the edge currents
through one cycle gives an exact finite continuity equation. This conserved
internal code charge has no assigned electric charge unit here.

## Complete receiver experiment

`positive_local` constructs the normalized positive-constraint state from a
definite link word. It is the coherent orbit state
`(|g,w>+χ(w)|g xor s,w>)/sqrt(2)`. The complete projector instrument in the
owning package contains both signs; conditioning its output on the plus
record has probability one half. The experiment explicitly declares its
normalized physical input rather than hiding a failed projection.

Take the receiver initially in this state with `g=1`, port 6, triplet
internal level 0, blank parity and address 0. The sender has `g=1`, doublet
level 3, blank parity and address 0. Compare sender port 5 with port 6.
Execute one `E` gate followed by the receiver's literal `U`. The sender is
retained as an arbitrary quantum reference in the sparse pair interpreter;
it is not traced between gates. Sender onsite `U` commutes with the final
receiver detector and therefore does not change this marginal.

All five regional face projectors are retained. In face order
`(0,2,5,16,17)` the exact receiver probabilities are

\[
\begin{array}{c|ccccc}
\text{sender port}&0&2&5&16&17\\\hline
5&301/600&1/1800&1/1800&0&179/360\\
6&1&0&0&0&0.
\end{array}
\]

Each row sums to one. Every branch preserves the receiver constraint.
This executes a connected local gate protocol in which the incoming wave
changes the receiver's native defect motion. It is a circuit segment;
the calculation does not substitute a two-cell marginal for a full
seven-layer cubic trajectory. A separate charge preparation with sender
charge `−3` and receiver charge `2` gives receiver probabilities `(1/2,1/2)`
after exchange, while the total charge remains `−1` exactly.

## Strict causality and finite-volume comparison

An onsite gate leaves cell support fixed. In one matching layer, conjugating
an observable enlarges its support only across matched edges touching that
support. Reverse the six layers to compute its exact backward circuit cone.
After `t` cycles, every ancestor has coordinate distance at most `2t` and
Manhattan distance at most `6t` from the final support. Distances in a torus
use the periodic metric. The bounds are
independent of box size and local Hilbert-space dimension. They are circuit
bounds, not a derived Lorentz cone.

The same argument applies to the Heisenberg dual of each branch of a local
apparatus. Summing all branch weights gives the complete record law.
Attach any untouched reference algebra `R`. If two initial states agree
on the entire backward cone together with `R`, and the two boxes have
identical gates in that cone, their complete local record/reference outputs
agree exactly. No assumption of a product input is needed. There is no
shared adaptive FIFO or external controller communicating outside the
specified circuit edges. Apparatus feed-forward must itself be included in
the gate graph if it changes that graph's causal support.

This is the same partial-trace mechanism formalized in
`Lean/ObserverPatchHolography/Locality/NoSignalling.lean`.
`HistoryCausalInvariance.lean` gives the corresponding finite history
transport interface. The source-specific local gates and their support are
supplied by this construction. Those Lean modules do not certify Lorentz
covariance or a continuum limit of this circuit.

For a cube centered far from its boundary, the local operator and all its
finite-time record statistics are identical in every sufficiently large
cube. This supplies the fixed-spacing infinite-volume local dynamics by
eventual equality. It supplies no limit as spacing and tick duration vanish.

## Exact stationary background and executable band

For each triplet color `j∈{0,1,2}`, let

\[
\Omega_j=|\text{uniform 10 links}\rangle\otimes
|\text{uniform 12 ports}\rangle\otimes|j\rangle\otimes|0\rangle
\otimes|\text{uniform 60 addresses}\rangle.
\]

The address reflection fixes the uniform address. In the triplet, central
link signs act trivially, and the port-symmetric vector has no antisymmetric
component for seam repair. Each controlled defect move permutes the
uniform link state. Therefore `UΩ_j=Ω_j` and `ΓΩ_j=Ω_j`. The frame `F`
acts trivially on all three vectors. The owning
[background construction](../background/) supplies their native preparation.

An identical product of `Ω_0` over the cells is fixed by every edge gate
because `a+b=1`. Its finite-volume reductions are consistent and stationary.
One color replacement `Ω_0→Ω_1` gives a strictly invariant single-excitation
subspace. The local `U` acts as the identity there, while `E` applies the
two-site matrix `[[a,b],[b,a]]` to its spatial amplitudes. The executable
`one_charge_step` is precisely this restriction. This color excitation is
an observer-accessible internal preparation and readout, without a physical
particle assignment.

With a period-two unit cell, put `z=exp(ik)`. In one axis the literal even
and odd layers have the Bloch matrices

\[
E_0=\begin{pmatrix}a&b\\b&a\end{pmatrix},\qquad
E_1=\begin{pmatrix}a&bz^{-1}\\bz&a\end{pmatrix},\qquad
U_1=E_1E_0.
\]

The receipt checks this symbol by lifting all eight parity columns to
finite periodic cubes of side 4 and 8, including `z=i`. It also verifies
normalization, the inverse and translation by two cells. Direct algebra
gives

\[
\operatorname{tr}U_1=i(\cos k-1),\quad\det U_1=-1,\quad
\lambda_\pm=e^{-i\pi/2\pm i\omega(k)},\quad
\cos\omega(k)=\sin^2(k/2).
\]

The three-axis symbol is `U_1(k_x)⊗U_1(k_y)⊗U_1(k_z)`.
The branch containing the uniform eigenvector has phase

\[
-\frac{k_x^2+k_y^2+k_z^2}{4}+O(|k|^4).
\]

At a zone corner the phase splittings are sums of signed
`|q_axis|/sqrt(2)`. They do not give the relativistic norm
`sqrt(q_x²+q_y²+q_z²)`. Internal color multiplicity and a finite causal cone
therefore do not supply a relativistic dispersion relation for this sector.
This result concerns the exact spectator triplet sector; it does not decide
the interacting doublet sector's spectrum.

Period-two translations commute with the full periodic circuit because
they permute each matching and each identical local gate. Their finite
Fourier projectors define lattice momentum records. The polynomial
commutation interface in `Lean/QFT/PeriodicMomentum.lean` applies to these
actual translations. It does not identify a principal unitary phase with a
positive physical energy. A stationary background, a faithful reference,
and a positive-energy vacuum are separate conditions.

## Reused constructions and boundaries

The local complete-record causality argument reuses
`model/native_local_response`, with the local support proved for the
specified gates. The shared-seam flux Laplacian in the repair-law research
is a diffusion operator; its quadratic bands do not identify the unitary
clock here. The scalar Gaussian frequency hypotheses in
`Lean/SourceLaw/ContinuumTimeLimit.lean` are not supplied by this circuit.
The spatial cube and exact local continuity equation require neither an
imported Wilson action nor a Standard Model mass fit.

The result is a connected, gauge-constrained finite family with complete
local observations and exact finite-volume control. A relativistic branch,
interacting scaling limit, physical spin/statistics, measured species and
absolute energy scale require their own constructions on the chosen source.
