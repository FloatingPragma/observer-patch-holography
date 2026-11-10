# Boundary charge, allowed repairs and finite spectral gaps

The central-link background in the coherent walk changes loop holonomies and
detector probabilities. Vertex rechartings cannot remove its flux. A link update
can remove its adjacent pair, so that background alone establishes no dynamical
defect protection. This note specifies an additional finite link-memory program
and proves the protection it has under an explicit move class.

## An owned finite link register

The [program](protection.py) stores one bit `ell_e` on each of the actual
icosahedron's thirty edges. It represents the central code transporter
`G_e=Z^(ell_e)`, with `Z=diag(I3,-I2)`. The corresponding face bit is

`b_f = sum_(e in boundary f) ell_e mod 2`.

A value one denotes a nontrivial central face holonomy. The field register costs
thirty bits, before the scheduler and retained event records. Its clock counts
attempted local link updates. The register, permitted moves and selected boundary
are declared architecture; their uniqueness or selection by the committed load
rule is not asserted.

Each edge borders two faces. Flipping one edge toggles exactly those two face
bits. Hence the parity of all twenty face bits is even for every link word.
The executable checks the two incidence ranks, eleven for vertex boundaries and
nineteen for face boundaries. Thus the filled spherical carrier has first
homology dimension `30-11-19=0` over the two-element field. Its graph cycles must
not be mistaken for protected noncontractible cycles of the filled sphere.

## Exact regional protection

For any face region R, multiplication of its face holonomies cancels every
interior edge twice:

`sum_(f in R) b_f = sum_(e in boundary R) ell_e mod 2`.

Fix the boundary and exterior link bits. Every permitted interior-edge flip
then preserves this regional parity. If it is odd, no finite sequence of those
repairs can reach a configuration with all region faces flat. This is genuine
sector protection **relative to the specified boundary and operation class**.
It neither protects the position of an individual defect nor protects a quantum
memory against arbitrary noise.

The actual control takes R to be the five faces meeting vertex 5. Its boundary
is `(0,1),(0,11),(1,9),(4,9),(4,11)`, and its interior edges join vertex 5 to
`0,1,4,9,11`. An initial negative link `(0,1)` creates flux on `(0,1,5)` inside
the region and `(0,1,7)` outside it. The boundary records odd interior parity.
All 32 possible interior link words retain it, with sixteen distinct odd face
syndromes. No interior word is the vacuum.

Allowing the boundary link `(0,1)` to flip back erases the entire initial pair.
The program reports this counterexample alongside the protected result. A
boundary condition must not be silently converted into an unrestricted repair
theorem. Vertex sign rechartings flip incident link bits in pairs on every
triangle and preserve all face holonomies; those rechartings alone do not move
the defects.

## Executed motion with complete records

The transport event selects an interior edge. If exactly one of its two faces
has flux, flip the edge; otherwise retain the state and emit a no-operation
record. A move transfers the defect to the neighboring face. The event is an
involution, since either direction of a move again has exactly one neighboring
defect. It preserves the number of interior defects on the one-defect class.

Enumeration gives ten reachable link states and all five possible interior
defect positions. The two link states at each position differ by retained
frame information. All fifty state/edge event pairs are checked: twenty move
and thirty retain a no-operation record. A deterministic round-robin schedule
executes thirty events before the state and scheduler return to their initial
values. Every event, including a hold, has its before/after face syndrome in
the [receipt](protection_results.json).

This is a finite source construction of movable regional flux with a conserved
boundary label. It has no particle spin or statistics, no scattering law and no
derived energy cost. In particular, the dynamic thirty-bit program is **not**
the static-link parameter of the coherent wave observer. The separate
[joint defect-wave program](../coherent_defect_wave/PROOF.md) feeds this field
register into seam transport and couples its permitted moves to wave occupation.
It supplies its own combined clock, complete readouts and
[stationary reference](../coherent_defect_wave/PREPARATION.md). That construction
does not inherit the fixed-W reference: a stationary reference of a fixed W is
not automatically stationary under a time-varying sequence of link backgrounds.

## Why a finite flux gap is not a rest mass

Let a connected graph have N vertices, M edges and one negative central link in
the doublet block. Its signed Laplacian satisfies

`x* L_signed x = sum_(edges ab) |x_a - s_ab x_b|²`.

The constant trial vector has norm squared N. Only the negative edge contributes,
with squared difference four. The variational principle therefore gives

`0 <= lambda_min(L_signed) <= 4/N`.

On the actual carrier the trial quotient is `1/3`. A triangle containing the
negative link is frustrated. A zero-energy vector would have to be parallel
across every edge and change sign around that triangle, forcing it to vanish;
connectedness then makes the signed Laplacian positive definite. Thus this
finite block has a positive minimum eigenvalue, while its full wave includes
the flat triplet's zero branch.

For the declared reflection-walk clock,

`omega_min = arccos(1-lambda_min/(2M)) <= arccos(1-2/(MN))`.

In a growing-volume family with a bounded number of such links and fixed local
spectral normalization, the trial bound tends to zero. It supplies no uniform
positive infinite-volume mass gap. A capacity-dependent conversion of ticks to
physical time changes that claim and requires its own source derivation; it
cannot be chosen solely to turn a finite-volume frequency into a desired mass.
Rest-mass identification also requires a controlled energy-momentum dispersion,
physical transformations and one energy convention on the same source.

## Reproduction

```sh
python -m model.coherent_seam_repair.protection --write
python -m model.coherent_seam_repair.protection --check
python -m pytest -q model/coherent_seam_repair/test_protection.py
```

The coefficient identities are checked on all thirty link generators and
extend by linearity to every thirty-bit state. Reachability and the complete
finite transport schedule are exhaustive on the declared region. No Monte Carlo
fit or experimental particle datum enters this result.
