# Typed finite compositions

The executable produces simulation-generated observer records on two declared
constructions. The quantum carrier has 384 states and periodic record modulus
6. The reciprocal classical carrier has two positions on 96 ports, hence
9216 joint states. The external local walk has twelve event slots per position
and a distinct generator. Equality of names, a common geometry, or a shared
count of source pulses does not identify these constructions.

## Checked admission and unchanged maps

The frozen clock helpers have mathematical preconditions: a finite Hermitian
positive density matrix of trace at most one, the declared shape, and finite
nonnegative time and rates. The number-sector semigroup and count partition
also require number-block-diagonal inputs. A conditional trajectory requires
trace one, positive native rate, a uniformization rate at least as large, and
strictly ordered interior event times with valid marks. Full-space unitary
and native maps admit cross-number coherences.

`strict_facade.py` checks these conditions before invoking a raw helper. Matrix
checks use an explicit absolute floating tolerance of 2e-12 for Hermiticity,
positive eigenvalues and trace bounds. Number-sector APIs require exactly zero
off-number-block entries. The producer's per-entry 1e-13 tolerance is not an
admitted compression domain. These tolerances classify numerical representations. They do
not change the mathematical domain to all Hermitian matrices. No input is
projected, clipped, symmetrized, or renormalized. The exact theorem concerns
positive matrices. The numerical checker does not provide an interval proof
of positivity for a matrix whose least computed eigenvalue lies within the
roundoff tolerance.

The force helper requires a finite nonnegative horizon and a finite
nonnegative 2 by 288 rate array with positive total rate, two valid source
positions, the actual seam list, and a valid target seam. Zero individual
rates and zero horizon are admitted. The all-zero-rate case lies outside the
raw exposure formula, which divides by total rate, and is explicitly refused.
The force score additionally uses the positive horizon of the frozen member.

On admitted inputs every facade method calls the identical frozen method with
identical argument objects and returns its output unchanged. This proves map
equality on the common exact input domain independently of numerical examples.
The executable also compares outputs against the independently archived
`frozen_common_clock.py` and `frozen_force_measurement.py`. Normalized states,
subnormalized reservoir states, the zero matrix, cross-number coherences on
full-space maps, tilted clocks and a marked trajectory are exercised. A positive matrix with off-block entries of 9e-14 gives a raw zero-time
compression loss whose trace norm exceeds 2e-12. The checked compressed APIs
reject it; full-space unitary, native and trajectory APIs retain their
cross-number domain. Every actual state inspected in the declared composition
and all twelve frozen clock inputs have exactly zero off-block entries. Guard
mutations use a trap which raises if a raw helper is entered, establishing
that malformed inputs are refused before raw computation. The opposing audit
sparse matrix fixtures are loaded directly; its negative-horizon and
negative-rate mutations are reproduced. The raw helpers retain their
input preconditions; only this checked facade claims checked entry points.

## Quantum composition and coarsening

Let the four reservoir preparation maps be R_a, the complete count maps be
C_0, C_1, C_ge2, and the regional detector conclusive maps be D_q. Its reset
and export failures enter separate flagged terminal spaces. Let P_g be a
complete grouped energy and global-charge projector of the actual Hamiltonian.
All regenerated local detector Kraus families, including every flagged abort
map, are compared directly with the pinned archive. The fine conclusive branch is

    P_g D_b C_j D_a C_i R_r(rho_ref) P_g.

Every map in this expression is linear. Summing a recorded label after its
operation therefore commutes with each subsequent map. The executable
regenerates all four reservoir branches, executes the complete first count
partition, executes the detector, retains all stop families, executes the
complete second count partition on the sum of detector-success states, and
executes the second detector and all energy groups. Its 348 terminal
coefficient matrices equal the sum of the frozen 576 by 348 fine coefficient
matrices within the stated numerical tolerance. This comparison checks
poststates in each spectral subspace, including two-dimensional degenerate
groups, rather than only the total probability. No success-conditioned law
is used. Failure states and their fine labels retain the pinned particle
archive's rectangular Kraus reconstruction; the coarsened run reports their
unconditional probabilities and never assigns them a C10 energy.

The coherent Gauss-code load/export and routing service is a declared
reference-preserving operation on the proper code. It does not measure fine
field or boundary records. The native detector's successful local Kraus maps
act tensor the identity on the 48 spectator coordinates. Energy readout is
a declared grouped spectral Luders instrument. The package supplies no native
spectral compiler. Source evolution is paused during detector apparatus under
the declared schedule. This convention does not identify apparatus pulse
counts with the Hamiltonian duration parameter or with physical seconds.

## Field and force domains

At modulus 6 the integer-field family is checked against the actual C10
Hamiltonian and every complete MIX Kraus operator by an explicit basis
permutation. Its lifted preparation retains all four reservoir columns and
all 384 source labels. Its complete-history comparison bound is

    min(1, 2 sum_{j >= floor((M-1)/2)-5+1} (14837/5000)^j / j!).

The modulus-6 bound equals one. Modulus 42 gives the rational upper bound in
`RESULT.json`, with outward floating value approximately 4.18541624e-6. The
executable constructs the 2688-state modulus-42 member and evolves every lifted
reservoir column, retaining its norm. The analytic certificate includes the
complete Poisson event mixture; the coherent execution is a finite witness,
not a substitute for that theorem. The useful bound changes the regulator and
uses a lifted modulus-6 preparation, not the modulus-42 ground reference. The
modulus-6 particle detector and excitation dictionary do not transfer through
a label match. A proper-code modulus-42 detector, its actual spectral readout,
and a bounded observable comparison must be supplied for that claim.

On the reciprocal classical source, the executable regenerates the 512-round
pair energy table and the 40-bit integer controller on all 22 pair positions
needed by the complete first-event instrument. Each outgoing rate is formed
from its energy difference; the native endpoint rule supplies every coin and
current outcome. The exposure score retains no-event, holding, cap and both
preparation labels. The integer controller to rate and marked-path comparison
uses the frozen exact-real certificate. Floating exponential evaluation has
no separate interval enclosure. Static smoothed three-dimensional Coulomb
refinement has no supplied trajectory clock and supplies no dynamic Maxwell
or Lorentz transfer.

## Physical claim gates

The gate binds source generator, carrier, regulator, preparation, clock,
observable, and complete outcome domain. It requires regenerated array
witnesses for the quantum join and refuses label-only calls. Mutations alter
each interface field, remove the count tail or apparatus refusal, assert a
zero field tail, or attempt a cross-construction substitution. The command
must reject each mutation. It also checks frozen dependency hashes before
and after execution.

`ALL31_GATE.json` copies the criterion wording and exact same-source missing
connection from the pinned `ALL31_READY.json` snapshot. Finite evidence is
attached only to the applicable rows. Row 1 has a finite literal charge
mechanism, and its physical interpretation is a separate gate. No physical
row is promoted by this integration. Component counts and component audits
provide no proof that all 31 physical criteria or one common physical
construction are complete.

The separate public row-1 law controls use a twelve-port native protocol and
fixed phase sampling, with 32 samples per edge and a declared absolute turn
bound of three. Native histories satisfy the tested integer-unit and regional
gates. Half-current histories preserve continuity and fail contour closure.
Random integer histories satisfy contour topology while some fail the
integer-unit or regional gate. Removing a record is detected against the
source loads; the damaged integer tape retains valid spatial contours. The
pointer reads counters, and the sampled phases are a diagnostic encoding.
These facts support a partial-at-most physical interpretation, with no native
law uniqueness award and no implied C10 instrument transfer from this separate
protocol. The code and receipt are bound in `ROW1_SCOPE.json`.


## Historical source custody

`run.py` executes the original relative-path layout in `frozen_workspace`.
All 145 reviewed input pins retain their expected bytes. The routing helper
at `openwave-push/issues/N3D_cubic_routing/n3dlib.py` is supplied by the exact
archived `consumed/n3dlib.py` bytes whose SHA-256 starts `6459c6f096c9`.
The concurrently edited live helper is preserved and receives no validation
from this historical result. No producer source contents or audit pins are
rewritten. Imported local dependencies discovered by constructor execution
are copied and hashed separately; their count is a custody statistic.

The child process disables Python site startup and environment import paths,
adds only the installed external library directories, and installs a guard
before project imports. The guard denies file access to the mutable project
outside the snapshot and rejects a deliberate attempted read of the live
routing helper. Snapshot source hashes are checked before and after the run.
`SNAPSHOT_EXECUTION.json` records the custody and isolation result. External
Python numerical libraries remain environment dependencies; the analytic
certificates do not acquire a floating-library error enclosure from isolation.
