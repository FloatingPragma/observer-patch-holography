# Native interacting field in local oscillator registers

This package defines a finite native implementation of the size-two,
one-charge quadratic field source. It retains 190 coexact oscillators and
three harmonic canonical pairs. The architecture uses a register for each
oscillator, actual twelve-control M6 words, heralded reset with retained
timeouts, and the native coin/current pointer. The large word is specified
by a finite iterator and exact resource bounds. Small complete words are
executed against the full six-level source matrices on every input column.
No large joint matrix is claimed to have been executed.

The source is the frozen quadratic construction with electric coefficient
`a=1/2` and magnetic coefficient `b=3`. Its integer incidence matrix is B,
its face-boundary matrix is A, and K=A A^T. These are architecture inputs.
No empirical electromagnetic coupling or target particle mass is fitted.
The source preparation is the coexact Gaussian vacuum, three independent
harmonic Gaussians of electric variance 1/2, and a one-particle spinor on
the eight actual root cells. Every signed route and all modes are retained.
The proof extends to arbitrary spectator reference vectors in the admitted
matter input. It supplies neither a two-fermion native route word nor a
spatial interacting limit.

## Frozen coordinates and actual source comparison

`INPUTS.json` binds read-only snapshots of the producer import and data
closure. The quadratic manifest digest is
`43fca237a676954e4a213b0888e2cdfee70bb9d8cefd8fbf113f5d2a60cbc81e`.
`geometry.py` reconstructs the integer complex, checks the frozen frame's
Gram and curvature residuals by exact integer limb products, and checks
the resulting polar-frame error inequalities. The curvature norm majorant
1504/17 is consumed from the frozen quadratic-field certificate.

Let Q and Lambda be the frozen rational frame and positive diagonal
eigenvalue approximation. The exact physical coexact frame is

    U = polar(K Q Lambda^-1).

The producer certificate gives

    ||U-Q|| <= delta,
    ||U^T K U-Lambda|| <= eta,
    ||Q|| <= sqrt(1+epsilon).

These inequalities do not identify Q as an exact eigenframe. The comparison
source uses canonical q,p coordinates in U, diagonal frequencies
`omega_j=sqrt(ab Lambda_j)`, and current coordinates Q^T J. The exact source
uses curvature C=U^T K U and current coordinates U^T J. The harmonic frame
is constructed exactly: rational null vectors of `[B; A^T]` are orthogonalized
over the rationals and divided by certified square roots of their squared
norms. The common isotropic harmonic Gaussian and p^2 Hamiltonian are
invariant under this orthogonal choice of the three coordinates.

Put `w0=sqrt(ab(min Lambda-eta))`, with an outward rational lower bound,
and `v=a/(2w0)`. The true and diagonal comparison ground states have a
positive real overlap by the real positive Gaussian convention. Their
norm difference is at most

    epsilon_vac = sqrt(190 ab eta / w0^2).

To see this, the expectation of the potential perturbation in the diagonal
vacuum has absolute value at most `190 b eta v/2`. The zero-point energies
differ by at most `190 ab eta/(4w0)`, from the square-root operator bound.
The excess energy is therefore at most `190 ab eta/(2w0)`. Divide by the
true oscillator gap w0 and use `||Omega-Omega0||^2 <= 2(1-|overlap|^2)`.
All zero-point phases are retained in the native free word.

For any comparison coherent branch with occupation at most lambda, its
q covariance is bounded by v I and its squared mean norm by 4v lambda.
Consequently

    ||q^T(C-Lambda)q psi||
      <= eta v sqrt((4lambda+190)^2+16lambda+380).

This follows by evaluating the fourth moment of a real Gaussian vector.
The free perturbation is b/2 times this quadratic form. Integrate over all
free intervals and multiply by the upper branch l1 weight. For a current
interval, the coordinate mismatch is at most delta||J|| and its generator
norm on the coherent branch is at most

    delta ||J|| sqrt(v(4lambda+1)).

Integrating every current interval gives the current term in
`geometry_transfer`. The total is the vacuum error plus the free and
current Duhamel errors. Matter operations and the unchanged harmonic and
Coulomb factors preserve the norm comparison. These bounds apply to all
branches, including interfering branches, and tensor with a reference.

## Tensor-box cutoff on all positive modes

Let `P=product_j 1_{N_j<=R}`. The free diagonal Hamiltonian commutes with P.
For a routed translation write

    L_J = sum_j gamma_j (a_j+a_j^dagger),
    gamma_j = (Q^T J)_j sqrt(a/(2omega_j)).

The retained finite generator is P L_J P. Since the modes commute, its
exponential factorizes into 190 local truncated oscillator exponentials.
This factorization avoids a joint array of dimension `(R+1)^190`.

Let L be the cumulative sum of Euclidean coefficient norms over every
current interval, including inverse and repeated kicks. The exact route
norm bounds from the interacting Gaussian producer are multiplied by
`sqrt(1+epsilon)` for Q. Every prefix coherent branch has total occupation
at most `lambda=L^2`. Let Bpath bound the number of prefix branches.
The sum of branch coefficient norms is at most sqrt(Bpath), also when the
coefficients are arbitrary reference vectors.

The complement of the box is contained in the total-occupation tail N>R.
For `R+2>lambda`, set

    p = lambda^(R+1)/(R+1)!,
    t = p (R+2)/(R+2-lambda).

Each boundary face is crossed by the lowering operator from local number
R+1 to R. Its norm is sqrt(R+1); each face's coherent probability is at
most p. Therefore

    ||P L_J (1-P) psi||
      <= sqrt(Bpath) sum_j |gamma_j| sqrt(R+1) sqrt(p)
      <= sqrt(Bpath) sqrt(190)||gamma|| sqrt(R+1) sqrt(p).

The projected exact equation has this boundary source. Integrating it and
adding the terminal tail yields

    ||psi-psi_box||^2
      <= min(4, 2 Bpath [t+190 L^2(R+1)p]).

This is a tensor-box theorem. The factor 190 is included explicitly; the
total-photon-cap flux theorem is not silently substituted. For a vacuous
cap the code reports the bound 2. The free vacuum input itself lies exactly
inside every positive cap.

## Three harmonic modes with a finite regulator

For each harmonic coordinate use an auxiliary Hermite basis with
`q=(c+c^dagger)/sqrt(2)` and `p=i(c^dagger-c)/sqrt(2)`. Its initial Gaussian
has Var(q)=Var(p)=1/2 and is exactly the auxiliary number vacuum. The
physical Hamiltonian is `a p^2/2`, which freely shears the packet. Calling
this reference state an auxiliary Hermite vacuum does not claim a
normalizable ground state for the harmonic Hamiltonian.

During any branch through free time T and cumulative current length S,
the transformed annihilator has the form `u c+v c^dagger+alpha`, with

    |u|+|v| <= 1+aT,
    |alpha| <= (1+aT) S.

These are uniform upper bounds; S sums all integer current norms. For an
integer k>=2, the number moment satisfies

    ||(N+1)^k psi|| <= M_k,
    M_k = [((1+aT)sqrt(2k)+(1+aT)S)^2+1]^k.

Conjugate each N+1 to the initial vacuum. Each application raises the
maximum polynomial degree by at most two. On degree at most 2r, the
annihilation and creation norms are at most sqrt(2r) and sqrt(2r+1).
Bounding both factors in `c'^dagger c'+1` by the displayed expression and
iterating k times proves the formula. No measured occupation is used.

For a local cap R, the norm of its tail N>=R+1 is at most
`tau=M_k/(R+2)^k`. The three-mode box terminal tail is at most sqrt(3)tau.
A q or p kick has boundary norm at most `sqrt((R+1)/2) tau`. The p^2 free
generator connects R+1 to R-1 and R+2 to R, and its boundary norm is at most
`a sqrt((R+1)(R+2)) tau/4` per mode. Thus the full harmonic error is at most

    sqrt(Bpath) [sqrt(3)tau
      + sqrt(3)S sqrt((R+1)/2)tau
      + 3 T a sqrt((R+1)(R+2))tau/4].

The code caps this state bound by 2. The finite free generator is the
compression of the full p^2, including its full diagonal at the top number.
It is not the square of the compressed p. Both its diagonal and its
number-difference-two pairs are present in the native word. For a static
holding interval the separate formula sets S=0 and Bpath=1.

## Native word and resource theorem

Each local number register has `ceil(log2(R+1))` bits. The physical numbers
are 0 through R. Every local padding number has generator zero, so its
action is defined. Tensor products containing padding also have a defined
action from the other local factors. No padding column is discarded.

For an ordered pair l<r the matrices

    X_lr = |r><l|+|l><r|,
    Y_lr = i|r><l|-i|l><r|

have exact expansions in commuting Pauli strings. Each uses 2^(k-1)
terms on a k-bit local register. Multiplying by c binary predicate
projectors gives 2^c times as many terms. `pair_word` emits every one of
these factors through the frozen rational M6 compiler. A diagonal phase
is an exact Walsh expansion, including its scalar term. Prefix Walsh sums
of a linear number energy are computed by a binary recursion in O(k),
without constructing an energy array.

The routed source diagonalizes each actual alpha axis by the producer's
Pauli basis and CNOT identities. The sign is the resulting physical rail;
the three actual root-cell bits and that rail control all 193 coefficients
of every signed route. The destination map is independently read from all
eight actual routes and checked to be the size-two coordinate flip. The
spin basis is then restored. Mass uses the actual beta=Z on spin bit zero.
The charge-sector Coulomb phase is calculated from an exact reduced
Laplacian solve. A fresh common clock wire receives `exp(-i h Z/16)` only
after the completed macro. The initial plane/clock preparation is eight
native pulses after its declared zero resets.

The local truncated translation is approximated by s repetitions of its
adjacent ladder-pair product. Only adjacent pairs fail to commute; their
commutator norms are `gamma^2 sqrt(n(n+1))`. Hence its operator error is
at most

    gamma^2 R(R+1)/(4s).

The harmonic free product contains one diagonal factor and all difference-
two pairs. If their norm sum is bounded by
`S_H=|time|a[(2R+1)+R(R+1)]/4`, its product-formula error is at most
`S_H^2/(2s)`. Rational square-root intervals price every ladder and
frequency coefficient. The whole-program certificate sums the actual
coefficient interval radii for all 48 routes and all 193 modes, together
with the local ladder and frequency rounding terms.

Each Pauli string on k bits needs at most 8k+1 native pulses. The exported
resource formulas give a finite bound of order
`macros*s*193*R^2*log R`, with fixed sixteen-route and predicate factors.
There is no exponential dependence on the number of field modes. The
Hilbert-space dimension itself is exponential and is reported as
`(Rc+1)^190 (Rh+1)^3`. The chosen s can be extremely large; polynomial
syntax is not a practical runtime claim.

Every pulse consists of twelve exact rational integrated native controls.
The frozen compiler encloses sqrt(2), sqrt(5), pi and its exact source
coefficients by outward rational intervals. A conservative universal
full-response error for the maximum declared angle is multiplied by the
entire pulse bound. Precision increases until this sum meets its budget.
The proof compares full six-level responses before coherent export, so
every first complement branch is included in the resulting channel bound.
The symbolic full iterator emits the same certified pulse records as the
small all-column replay. The replay includes both current-control values,
both field quadratures, vacuum scalar phase, local padding, noncommuting
free/current stages, and a changed-control rejection.

## Preparation and complete apparatus records

`preparation_contract` prepares every field register in number zero by the
frozen all-78-input reset contract. Repeated blocks stop at the first
heralded success, retain every fine record, and retain final timeout. The
union bound is the number of registers times the block failure bound raised
to the repetition count. This is an implemented reset schedule contract,
not a host assignment of a desired field state. The native memory load and
coherent export services remain declared architecture services, as in the
consumed rational compiler.

For arbitrary charged root sectors the physical electric observable is

    E_edge = (B^T L^+q)_edge
      + sum_j U_edge,j sqrt(omega_j/(2a)) i(a_j^dagger-a_j)
      + sum_k H_edge,k p_h,k.

The exact U frame lies in ker B and every harmonic column lies in ker B.
Consequently this declared observable has B E=q in every charge fiber.
The encoded charge label, the all-zero field preparation and this observable
define an explicit spectral-code realization of the Gauss-fiber map. They
do not require equivalence with the old finite Z16 edge representation.
Actual motion applies the nonzero transverse route factors, while static
charge preparation with a common field is a distinct declared instruction.

`edge_characteristic_word` controls `exp(i t E_edge)` on a probe. It uses
the exact rational Gauss particular phase for each source charge cell and
all 190+3 Y-quadrature oscillator factors. The implemented coexact frame Q
has a separately priced distance from U. For this observable the variance
majorant is `vp=omega_max/(2a)`; the additional coexact Weyl length is at
most `|t|sqrt((1+epsilon)vp)`, and the harmonic length increases by |t|.
`edge_transfer` includes these probe intervals in both cutoff bounds and
allows twice the branch count. The measurement-frame perturbation is
bounded by `|t|delta sqrt(2Bpath) sqrt(vp(4L_probe^2+1))`.

Native Hadamard preparation and X/Y analysis make the ideal binary laws
`(1+Re chi(t))/2` and `(1+Im chi(t))/2`. The actual pointer uses its six
native CNOT pulses and its committed C/E coin/current grammar. Both fine
instruments are executed and their effects are retained. Pointer-ready
work is obtained by blank reset and a native flip, with that error and the
probe/work reset timeouts included. State comparison, product-formula,
coefficient, native export, preparation and pointer bounds add to a complete
record total-variation bound. Failure records are never conditioned away.

This is an implemented bounded Fourier observation of an actual edge
field. A complete simultaneous sharp-bin measurement of every edge is not
asserted; finite oscillator p operators on overlapping edge combinations
need their own joint instrument. The characteristic observation and its
source transfer are sufficient to make the physical electric observable
operational, without claiming that stronger instrument.

## Common static field for a charge-energy instrument

`static_holding_contract` supplies the common field preparation, all 190
half-frequency phases, and the three harmonic p^2 words at a requested
holding time. `free_field_word(..., include_coulomb=False)` is its exact
finite iterator; the consumer supplies its own exact charge-sector phase.
For every common field state, this q-independent finite field unitary
cancels from the ideal reduced charge/clock law. Its native export errors
and reset timeouts are retained by that law. Its cutoff, spectral and
product-formula source comparisons remain relevant to full field and
edge-field records. This cancellation does not replace the nonzero
transverse disturbance in a routed moving-charge instruction.

## Scope of the executable evidence

The frozen spectral approximation contributes a nonzero error floor to
the displayed source comparison. Raising the oscillator caps does not
remove that floor. A refined spectral frame would require another frozen
integer certificate and the same priced transfer. The existing DYN13
single optical-loop periodization theorem is not used as a theorem for the
three free harmonic modes. Their finite regulator is proved above.

The executable reports the native small-word results, all actual current
coordinates, source comparison contracts, finite pulse prefix and complete
resource bounds. No reserved acceptance member is used and no official
row is changed. The construction does not infer physical seconds from
pulse counts, an observed coupling from architecture constants, a stable
interacting particle species from the prepared Gaussian, or all-31-row
validation from a field component.
