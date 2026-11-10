# A positive local energy law on the shared-link source

## Specified clock

At each spatial cell let `U_F,x` be the complete parity-controlled shared-link
occupation unitary in [the species construction](../species/PROOF.md).
For a neighboring pair let `S_xy` be its Gauss-dressed one-body exchange,
with both link registers held at their cells. On the full finite wave Fock
space define

\[
 A_x=I-\frac{U_{F,x}+U_{F,x}^*}{2},\qquad
 h_{xy}=\frac\pi4\,\mathrm d\Gamma(I-S_{xy}),\qquad
 H_\Lambda=\sum_{x\in\Lambda}A_x+\sum_{\langle x,y\rangle}h_{xy}.
\]

The continuous evolution is `exp(-it H_Lambda)`. These formulas, coefficient
one and coefficient `pi/4`, and dimensionless t define the architecture.
The edge choice reproduces the existing root-SWAP:
`exp[-i(pi/4)(I-S)]=(1-i)I/2+(1+i)S/2`, and exterior lifting preserves the
identity. The onsite positive functional is a declared replacement clock;
`exp(-i A_x)=U_F,x` is not asserted. The Floquet phase records of the
seven-layer circuit therefore cannot be relabelled as H energies.

## Positivity, locality and constraints

For any unitary V,

\[
 I-\frac{V+V^*}{2}=\frac12(I-V)^*(I-V),\qquad 0\le A\le2I.
\]

The local `U_F` is identity on its empty-wave sector. Hence
`0 <= A_x <= 2 P_occ,x <= 2 N_x`. The dressed S is a Hermitian involution,
so `I-S` is positive and bounded by `2I`. Exterior second quantization of
this one-body positive operator is positive and satisfies
`0 <= dGamma(I-S_xy) <= 2(N_x+N_y)`.
On a cubic graph of degree at most six, on the fixed total-number-N sector,

\[
 \boxed{\quad 0\le H_\Lambda\le(2+3\pi)N I.\quad}
\]

The upper bound is attained on N=1 by a dark mode at lattice momentum
`(pi,pi,pi)` on an even periodic cube. Thus `||H||=2+3pi` on that sector.
This bound is independent of the number of cells at fixed N. It is a bound
on energies of this declared clock; it supplies no capacity-to-SI conversion.
On the full finite occupation algebra each cell has 7200 modes. The local
norms are bounded by2 and3600pi, respectively, and the full finite H exists
without occupation truncation. Sparse N=1 and N=2 executions do not replace
the complete algebra by a truncated CAR space.

Each term acts on one cell or one neighboring pair, is even in fermion
parity, and commutes with the admitted local Gauss constraints. Disjoint
even term algebras commute. Every term preserves total N and the internal
U(3) times U(2) symmetry. Whole-cell translations, including the link
registers, permute the terms of a periodic cube, and therefore commute with
H. No full noncentral gauge action or Lorentz transformation follows from
these finite symmetries.

For finite t the continuous propagator generally has small tails beyond a
fixed circuit cone; the strict cone of a finite gate cycle is not assigned
to this clock. For a finite list of self-adjoint local terms H_r, the unitary
Duhamel identity and telescoping give

\[
 \left\|(\prod_r e^{-it H_r/n})^n-e^{-it\sum_r H_r}\right\|
 \le\frac{t^2}{2n}\sum_{r<s}\|[H_r,H_s]\|
 \le\frac{t^2}{n}\sum_{r<s}\|H_r\|\,\|H_s\|.
\]

Only intersecting supports contribute to the first sum. This is a finite
simulation error bound, not a continuum field theorem. The exact fixed-N
norm bound also controls Taylor remainders on that invariant sector.

## Reference, source energies and finite interaction

Any positive-Gauss empty-wave link state is annihilated by every term of H.
It is an exact zero-energy ground state. The choice is degenerate, includes
the specified odd regional flux in every cell, and is not identified with
a unique physical vacuum.

There are source energies independent of any fit. At carrier port2, choose
blank parity memory and address amplitudes `+1,-1` on `(0,1),(1,0)`, divided
by square root2. Its address sum vanishes, so the address reflection gives
minus one. Neither addressed seam includes port2, so seam repair acts as
identity. Port2 is outside all controlled interior edges, so their link
moves act as identity. Thus for every retained link word and every internal
component,

\[
 U\eta=-\eta,\qquad A\eta=2\eta.
\]

The positive-Gauss combination of the two link representatives preserves
this identity. The uniform-link, uniform-port/address triplet background
has `U Omega=Omega`, and therefore energy zero under A.

Two distinct internal dark modes at the same cell have exterior eigenvalue
`(-1)^2=1`; the endpoint controls are inactive at port2. Their exact onsite
energy is consequently zero, whereas one dark mode at each of two separated
cells has total onsite energy four. This is an explicit failure of an
additive free one-body energy law. It is an interaction in the specified
finite occupation Hamiltonian, without a scattering, binding-energy or
particle-lifetime identification. Hopping energy is kept separate in this
comparison.

The full two-cell occupation generator is also executed directly, including
fermionic creation/removal signs in `dGamma(I-S)`. On the constrained input
with link words `(1,3)` and occupied modes `(port5,internal3,0,0)` in the
first cell and `(port6,internal0,0,0)` in the second, it gives
`<H>=89/30+pi/2` and `<H^2>=3371/360+89pi/30+3pi^2/8`. Its generator image
has 968 nonzero exact amplitudes and preserves both constraints. These are
conditional on the normalized two-cell positive-Gauss preparation; its
complete ambient instrument retains the three failed signs as well as the
successful sign, each with probability one quarter.

The zero energy of the empty-wave ground is not a positive global gap proof.
On the product uniform-link background, one spatially uniform triplet wave
with uniform carrier ports/addresses is also annihilated by every onsite
and edge term. It is orthogonal to the empty-wave state by occupation number.

The admitted link sector has exactly one regional flux in every cell. Its
regional topological-charge operator is therefore `Q_x=I` throughout this
Hilbert space. Onsite moves relocate that flux among the cell's five faces;
the hopping term holds both link registers fixed. Consequently this law
has no intercell transport of the regional topological charge, even though
waves and local flux-position records influence one another. Its propagating
wave-energy bands are not a map of traveling protected fluxes to physical
particles. Such a map requires a link sector admitting different regional
charges and intercell charge-transfer moves, or a proved identification of
a different collective excitation with the physical particle observable.

## Exact two-wave collision in the same Hamiltonian

Retain the two dark modes at each cell with internal components `a=0,b=1`.
Their local occupation numbers are each zero or one. The actual source
wave has eigenvalue minus one on each mode, all endpoint controls are
inactive at port2, and exterior lifting gives `U_F=(-1)^(n_a+n_b)` on this
entire two-mode space. Consequently

\[
 A_x=2(n_{a,x}+n_{b,x})-4n_{a,x}n_{b,x}.
\]

The spatial hopping preserves the same dark address difference and both
link registers. On two cells, with `t_h=pi/4` and `D_x=n_a,x n_b,x`,
the full source Hamiltonian restricts exactly to

\[
 H_{\rm dark}=(2+t_h)N-4(D_L+D_R)
 -t_h\sum_{\alpha=a,b}(c^*_{L\alpha}c_{R\alpha}
                         +c^*_{R\alpha}c_{L\alpha}).
\]

This is an attractive two-site Hubbard operator in a specified invariant
source sector. `dark_collision.py` constructs the six normalized N=2
source vectors, including both positive-Gauss link pairs and all wave
address amplitudes. Applying the full `pair_energy_action` to every column
gives this six-by-six operator with zero amplitude outside the embedded
sector. The calculation therefore checks the source intertwiner as well
as its projected matrix.

The contact interaction is specific to the chosen continuous clock.
On these modes, the original local Floquet phase `(-1)^N` is an exterior
one-body phase. The product of the two original onsite ticks is identity
on N=2, verified on every embedded column. Applying `I-Re(U_F)` separately
at each cell introduces the displayed nonadditive contact. It is not
assigned to the original Floquet circuit's dark-sector dynamics.

In orbital order `(La,Lb,Ra,Rb)`, define

\[
 S=\frac{|La,Rb\rangle-|Lb,Ra\rangle}{\sqrt2},\qquad
 D_+=\frac{|La,Lb\rangle+|Ra,Rb\rangle}{\sqrt2}.
\]

Here S is an internal two-component singlet; no physical-spin interpretation
is made. The invariant block and complete six-dimensional sector spectrum are

\[
 H|_{S,D_+}=\begin{pmatrix}4+2t_h&-2t_h\\-2t_h&2t_h\end{pmatrix},
 \qquad
 \operatorname{spec}H_{\rm dark}|_{N=2}=
 \{2+2t_h\pm2\sqrt{1+t_h^2},\;2t_h,\;(4+2t_h)^{\times3}\}.
\]

The eigenvalue multiplicities follow from the exact characteristic polynomial.
Starting with S, the probability of either double-occupation outcome is

\[
 P_D(\tau)=\frac{t_h^2}{1+t_h^2}
       \sin^2(2\sqrt{1+t_h^2}\,\tau).
\]

The complete local occupation PVM on this invariant N=2 sector has outcomes
`n_L=0,1,2`, with probabilities `(P_D/2,1-P_D,P_D/2)`. Their sum is one
at every time. Higher occupation labels of the full number instrument have
zero weight on this sector. At
`tau=pi/sqrt(16+pi^2)`, the first transfer maximum, the actual source state
and complete readout are executed with probabilities

\[
 \left(\frac{\pi^2}{2(16+\pi^2)},\;
       \frac{16}{16+\pi^2},\;
       \frac{\pi^2}{2(16+\pi^2)}\right).
\]

Both local constraints and total wave number are preserved in every branch.
An initial pair with identical internal components remains one-per-cell
with zero double-occupation probability, providing a Pauli-blocked control.
The displayed collision probabilities are conditional on positive-Gauss
preparation. Its success probability is one quarter; all three other sign
records, each of probability one quarter, are retained. The receipt also
gives the occupation weights per original trial.

An additive onsite countercontrol keeps the same hopping, preparation and
detector but replaces the local positive-clock energy by `2(n_a+n_b)`.
Its S/D block has equal diagonal entries `4+2t_h`, so its double-occupation
probability is `sin^2(2t_h tau)`. At `t_h=pi/4` and `tau=1`, that probability
is one, whereas the interacting law has probability at most
`pi^2/(16+pi^2)`. The complete detector contrast is therefore at least
`16/(16+pi^2)`. The observation distinguishes the contact interaction from
independent hopping on the same finite source sector.

This is exact finite-time collision and occupation transfer in the same
declared H. It supplies no asymptotic S matrix, cross section, particle decay
width or Standard Model identification. The protected regional charge remains
one in each cell; the collision does not transport that topological charge.

## Exact joint spatial action and detector response

`model.py` consumes the exact positive-Gauss encoding of
[`cubic_readout`](../cubic_readout/PROOF.md). A basis vector stores one
occupied cell, its 7200-valued wave mode and all changed local link orbits;
the other cell links retain their prescribed background values. The
encoding is an isometry, not a frozen-background substitution.

On N=1, the sum of onsite A terms is exactly
`psi-(U_onsite psi+U_onsite* psi)/2`. The hopping term is `pi/4` times the
cubic wave-position Laplacian, with every link value held at its cell.
The complete Hamiltonian action is stored as `Apsi+pi Bpsi` with exact
Gaussian-rational coefficients. It retains source states where the defect
has moved. The same encoded preparation as the cubic circuit, localized
uniformly over eight cells with carrier mode `(port7,internal3,parity0,address0)`,
has exact moments

\[
 \begin{split}
 \langle H\rangle&=\frac{59}{30}+\frac{3\pi}{4},\\
 \langle H^2\rangle&=\frac{283}{72}+\frac{59\pi}{20}+\frac{3\pi^2}{4},\\
 \operatorname{Var}(H)&=\frac{113}{1800}+\frac{3\pi^2}{16}.
 \end{split}
\]

These moments use H; the circuit's frequency histogram is not used as an
energy histogram. Let Q record any local flux outside the initial face.
Then `Q psi=0` and exact action gives `||Q H psi||^2=1/1800`.
On N=1 take the rational bound `B=80/7 > 2+3pi`. The unitary Taylor
remainder gives

\[
 \|Q e^{-itH}\psi+itQH\psi\|\le t^2B^2/2.
\]

Because `1/43 < 1/sqrt(1800) < 1/42`, at `t=1/10000` the complete
yes/no flux measurement has probability bounded between

\[
 (t/43-t^2B^2/2)^2\quad\hbox{and}\quad(t/42+t^2B^2/2)^2.
\]

The lower endpoint is strictly positive. Both outcomes are retained. This
is an analytic finite-time probability enclosure from the executed exact
generator, not a sampled full-time histogram. The input is the normalized
positive-Gauss preparation; its ambient preparation flags are those of
the owning cubic construction.

## Exact bands and energy measurements

The dark modes are unaffected by the retained field in every cell. Their
entire one-wave subspace is invariant under H, with band

\[
 E_2(k)=2+\frac\pi2\sum_{i=1}^3(1-\cos k_i).
\]

The uniform triplet carrier mode on uniform link backgrounds has the
corresponding band `E_0(k)=(pi/2) sum_i(1-cos k_i)`. Actual side-four
plane waves in the dark band are evaluated through the full encoded
Hamiltonian action. Their energy variance is exactly zero at six displayed
three-dimensional momenta. The energy2 at zero momentum is a source band
edge in this continuous convention. It is not a measured lepton mass or a
relativistic rest-mass identification.

A complete local energy estimator follows from the owned source operations.
The two Kraus operators `(I+U_F,x)/2` and `(I-U_F,x)/2` are complete; assigning
values0 and2 gives mean `Tr(rho A_x)`. For an edge, diagonalize its dressed
one-body exchange into symmetric and antisymmetric modes and retain every
occupation count. This rotation is finite and explicit: with `P=diag(1,-i)`
on the two spatial lanes, the owned root-SWAP E satisfies
`P E P=exp(-i*pi/4) Hadamard` in the Gauss frame. The local quarter-turn
phase is a declared apparatus service; its exterior lift includes all
occupation sectors. Giving antisymmetric count n the value `pi*n/2` has
mean `Tr(rho h_xy)`. A classical recorded choice among all local terms, followed
by the selected complete instrument and inverse sampling weight, estimates
`Tr(rho H)`. This does not turn noncommuting local terms into a joint exact
energy PVM. On a finite H, controlled `exp(-itH)` with the explicit simulation
error bound supplies the Fourier energy apparatus. Resolution, simulation
error, all flags and quantum references must accompany that apparatus.

At fixed spacing, these energy laws are defined without an energy-logarithm
ambiguity. A change in the functional or coefficients is a different
architecture. Relativistic dispersion, physical Spin/current attachment,
the appropriate pole or threshold, and an independent unit comparison are
separate mathematical requirements.

## Finite theorem reuse

The exact local operator identities supply finite self-adjointness,
positivity, an invariant ground reference and translation commutation.
`Lean/QFT/PeriodicMomentum.lean` is compatible with the supplied commuting
translation and generator; it does not calculate a physical mass.
`Lean/QFT/SourceHistoryGNSDynamics.lean` uses a positive source-history
generator in its own representation and requires a separate intertwiner
before transferring its conclusion to this occupation algebra. The finite
product bound uses the same Duhamel argument as
[`spatial/ROUTE.md`](../spatial/ROUTE.md). The bound there for a free Dirac
symbol is not inherited by this interacting positive clock.
