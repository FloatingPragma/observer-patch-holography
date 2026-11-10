# Preparation and complete readout for the coupled defect wave

The source is the finite unitary `U=C4 C3 C2 C1 C0 B` specified in
[the construction](PROOF.md). Its basis consists of one of the ten admissible
fixed-boundary link words, twelve carrier ports, five internal states, two
coherent parity records and sixty oriented seam addresses. Its dimension is
`d=72000`. The controlled link updates execute in the listed order after the
link-dependent wave step. One autonomous clock tick is one complete `U`.

This note constructs a faithful stationary reference from this same `U` and
states the complete preparation and measurement records. It does not inherit
the reference of a fixed-link walk. The link register and the wave exchange
information during a tick; neither is replaced by an external time-dependent
background.

## A complete preparation circuit

An owned preparation register supplies a uniform ensemble of the `d` joint
computational basis states. Equivalently, prepare
`sum_k |k>_label |k>_source / sqrt(d)` and restrict to the source algebra.
Its reduced state is `I/d`. The label is retained and identifying its value
conditions the source state. Thus the uniform ensemble is a declared finite
preparation procedure, not a claim that every physical observer has this state.

Prepare a flag in
`|a>=sqrt(2/3)|0>+sqrt(1/3)|1>`, apply `U` controlled by flag value one,
and measure the flag in the orthonormal basis
`{|a>, |b>=sqrt(1/3)|0>-sqrt(2/3)|1>}`. The two recorded branches act on the
source by

`K_+=(2I+U)/3`,

`K_-=sqrt(2)(I-U)/3`.

The flag operations and controlled `U` are declared preparation apparatus.
Their coefficients follow from this circuit and use no measured particle
parameter. For every input, including an input entangled with a reference,

`K_+* K_+ = [5I+2(U+U*)]/9`,

`K_-* K_- = [4I-2(U+U*)]/9`,

`K_+* K_+ + K_-* K_- = I`.

The star denotes the adjoint. Both outcomes are retained. The corresponding
recorded instrument is completely positive and trace preserving when its flag
output is included. The construction uses the stipulated quantum register and
instrument semantics; it is not a derivation of Born probabilities from a
classical repair law.

## Same-source stationary reference

On input `I/d`, define

`Z=(5/4)d+Re tr(U)`,

`rho_+=[(5/4)I+(U+U*)/2]/Z`.

This is the success-conditioned reduced source state. It commutes with `U`,
since both `U` and its adjoint commute with their sum. Hence
`U rho_+ U*=rho_+` under the autonomous clock. If `exp(i omega)` is an
eigenvalue of `U`, the corresponding numerator eigenvalue is
`5/4+cos(omega)`, between `1/4` and `9/4`. Therefore `rho_+` is faithful,
without assuming a gap or diagonalizing a large matrix. Its smallest eigenvalue
is at least `1/(4Z)`.

The success probability is `p_+=4Z/(9d)`. Every pure input also has success
probability at least `1/9`, because the smallest singular value of `K_+`
is at least `1/3`. These statements do not identify `-log(rho_+)` with the
autonomous clock or with a physical energy in GeV.

## Exact normalization without a dense matrix

For a fixed wave port `p`, the five controlled link updates induce a permutation
`F_p` on the ten link words. Only words fixed by `F_p` contribute to `tr(U)`.
The wave-step diagonal is independent of the central link signs: the address
reflection has diagonal `-29/30`, and the seam correction contributes a
diagonal factor `1/2` on the ten incident oriented seam addresses and `1` on
the other fifty. Summing over the two parity and five internal states gives

`sum_(internal,parity,address) (W_g)_(x,x)
 = 10(-29/30)(10/2+50) = -1595/3`.

The exact fixed-word counts for ports zero through eleven are

`[6,6,10,10,6,0,10,10,10,6,10,6]`.

Their sum is ninety. The [independent calculation](preparation.py) therefore
gives

`tr(U)=(-1595/3)90=-47850`,

`Z=42150`,

`p_+=281/1080`, `p_-=799/1080`,

`lambda_min(rho_+) >= 1/168600`.

The five rim ports each control one involutory link move, with six fixed words.
The central port controls all five moves and fixes no word. The other six ports
control none and fix all ten. This finite count fixes the preparation numbers;
no fitted mass or target spectrum enters them.

## The admitted interior gauge constraint

The [gauge construction](GAUSS.md) combines the central-vertex link-star flip
with the central code action at that vertex into an involution `Gamma`.
It commutes with this same `U`. Its plus eigenspace has dimension `36000`;
`tr(Gamma U)=0` gives `tr(U|plus)=-23925`. The same filter restricted to that
sector consequently has `Z_plus=21075`, success probability `281/1080` and
reference eigenvalue bound `1/84300`.

Starting from the ambient `I/72000` and first measuring
`{(I+Gamma)/2,(I-Gamma)/2}` adds a success factor `1/2`. Thus the joint
probability of the plus constraint and the successful reference filter is
`281/2160`. The gauge-complement branch and both filter branches are retained.
This is one declared interior central gauge constraint, with boundary records
held fixed. It is not a claim that all Standard Model Gauss constraints have
been constructed. Apparatus operating within the constrained sector must
commute with `Gamma`, or include an explicit charged reference. A detector of
regional flux position is invariant; a bare link-word detector resolves the
two frame-related words and is not an invariant detector on that sector.

## Records, apparatus and capacity

A detector with projectors `P_j` satisfying `sum_j P_j=I` has complete joint
preparation-and-detector branches `P_j U^t K_a`. Their squared norms sum to
the input norm, with every preparation failure and detector complement
included. Complete probabilities here mean normalization over all outcomes of
the specified apparatus. They do not assert that the apparatus spans every
effect on the 72000-dimensional space.

The executable starts an explicit check from the basis state
`(g,p,j,r,a)=(1,0,3,0,0)`, applies the complete filter, executes one autonomous
tick on each branch and measures the wave coordinate `(0,3,0,0)`. For success
and failure flags respectively, the two detector probabilities are
`(2267953,5832047)/14580000` and `(1132603,2107397)/7290000`.
All four joint probabilities sum to one. The preparation marginal is `5/9`
for this pure input, rather than `281/1080`, which applies to the uniform
ensemble. Keeping the input distinction prevents an ensemble reference law
from being assigned to every individual history.

The joint register can be measured in its computational basis. A coarser detector
can report the regional flux position or a wave port, retaining the complement.
Its branch changes the joint state, so subsequent propagation uses the same
`U` on that changed state. Copying the resulting orthogonal outcome label
copies a classical record; it does not copy the unmeasured joint quantum state.

If each attempt starts from a freshly prepared `I/d`, at most `k` attempts have
all-failure probability `(799/1080)^k`. Fresh preparation labels, every flag
outcome, and any exhaustion outcome must be retained within the declared
finite record budget. Applying the filter again to its failure-conditioned
state without a fresh preparation does not have this geometric law.

Stationarity applies to `rho_+` on the reduced joint source algebra after
conditioning only on success. Conditioning further on a preparation label,
detector result or detailed retained history generally changes that state.
The enlarged source-plus-archive state is not asserted to be stationary.

The preparation and readout close a finite same-unitary reference construction
for movable boundary flux coupled to a wave. They supply no physical spin,
statistics, Lorentz transformation law, mass shell, continuum limit or particle
identification. The fixed regional boundary and the controlled repair order
are explicit construction choices.

## Reproduction

```sh
python -m model.coherent_defect_wave --check
python -m pytest -q model/coherent_defect_wave/test_preparation.py
```

The package receipt includes the independent fixed-word trace calculation and
the four executed preparation-and-detector histories.
