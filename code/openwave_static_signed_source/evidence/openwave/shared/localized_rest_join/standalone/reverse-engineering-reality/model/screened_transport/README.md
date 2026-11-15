# Transient current and corrected transport of the screened code

This CP-04/#1140 packet tests moving response on the **same** canonical
screened `q=1` code used in the
[particle packet](../particle_screening/README.md) and the
[compiled active-recovery service](../screened_recovery/README.md). It
distinguishes one native attempted repair from the continuing cycle that
applies each record-selected recovery and resets the classical core. Its
[campaign](campaign.json), [exact receipt](transport_results.json),
[source](transport.py), and [independent tests](test_transport.py) bind the
complete C/E branch menu. The time parameter below is a native attempt count.
The campaign's predecessor commit records the baseline provenance; its
`source_sha256` map binds the current verification inputs, including the
recovery receipt and its producer.

## A positive same-code native current form

For both logical columns, the one-port occupation vector on the twelve
source ports is

```text
w = (1/2, 1/2, 1/4, 1/4, 1/4, 0, 1/4, 0, 0, 0, 0, 0).
```

All eight occupied-pair supports are **nonedges** of the native five-regular
port graph. On seam `(a,b)`, exactly one of the two fair coins moves a load
when one endpoint is occupied. Hence the uniform-seam current second moment
for any logical input is the exact weighted Laplacian

```text
C_code = E[(x_after-x_before)(x_after-x_before)^T]
       = (1/60) sum_(a,b seam) (w_a+w_b) (e_b-e_a)(e_b-e_a)^T.
D_code(v) = v^T C_code v/2.
```

The receipt independently reconstructs this matrix from **all 240
pair/seam/coin source columns per logical state and law**, including zero
currents, and verifies the pointwise finite continuity update
`x_after-x_before = J(e_b-e_a)`. C and E agree. The matrix has rank eleven,
trace `1/3`, and native moving-current probability `1/6`. Its exact scaled
characteristic polynomial is

```text
det(lambda I - 240 C_code)
 = lambda (lambda-4)
   (lambda^3-30lambda^2+254lambda-528)
   (lambda^3-16lambda^2+80lambda-126)
   (lambda^4-30lambda^3+310lambda^2-1270lambda+1728).
```

This is a spectrum of a **one-step current covariance**, measured in
attempts and port occupations. It is not a Hamiltonian or an energy
dispersion. The code preparation is not the stationary source measure of
the [native action-form result](../clock_defect_energy/MODE_ENERGY.md), so
the finite normalized-symmetrizer theorem does not convert this matrix to a
physical kinetic term.

The prior source eigenprofile
`v=(1,-1,-1,1,0,0,0,0,0,0,0,0)` makes the comparison sharp. On this code,
`D_code(v)=13/120` is independent of the unknown logical input, but the
static variance of `F_v(x)=sum_p v_p x_p` is `1/2` on `|0_L>` and `5/2` on
`|1_L>`. The ratios are **`13/60` and `13/300`**, respectively. For a
normalized superposition `alpha|0_L>+beta|1_L>`, the variance is
`|alpha|^2/2+5|beta|^2/2`. Therefore the earlier stationary uniform-source
ratio `D/Var=1/10` and its wave-phase interpretation do **not** transfer as
one frequency of this corrected qubit. This is an exact same-candidate
countercontrol, rather than a fitted or missing numerical mass.

## The full corrected cycle has a flat logical channel

The separate source-executed protocol checks all `90` positive C and `118`
positive E complete raw flags. Its ideal branch relation is
`R_h K_h V = sqrt(p_h) V`, with zero-probability and abort flags retained.
Summing all flags gives the identity channel on every logical density matrix
and any untouched reference. Its four logical-superoperator eigenvalues are
`(1,1,1,1)`; the per-attempt generator `Phi-I` has four zeros. The completed
source pulse execution is checked by fresh numerical diagnostics against
fixed admission ceilings, without an exact hardware fidelity certificate.
A common global phase is invisible to
this channel and cannot be called zero rest energy.

Before the explicit core-label reset, one uniform-seam attempt from core
`(0,1)` leaves it at `(0,1)` with probability `14/15`. It reaches each of
eight other recorded labels with probability `1/120`:
`(0,5)`, `(0,7)`, `(0,8)`, `(0,9)`, `(1,5)`, `(1,7)`, `(1,10)`, `(1,11)`.
This is the same for C and E and every logical input. The core changes on
`1/15` of attempts, and leaving the native edge set is included in those
excursions. The corrected service then resets the core to `(0,1)` and
archives old/new/status; at each completed attempt boundary the core and
logical state are back at the same code. With independent uniform selectors,
the chance of at least one core reset in `t` attempts is
`1-(14/15)^t`, and the expected reset count is `t/15`. These count
controller activity, not cumulative spatial displacement or decay.

Skipping the core reset leaves only `29/30` in sharp `q=1` after recovery,
so treating the pre-reset excursion as a continuing protected moving center
would change the tested protocol. An adaptive moving-code family could be
tested separately; this packet gives no autonomous coherent propagation or
momentum parameter for the current fixed-code cycle. The new protocol's
uniform-seam mean native-plus-recovery pulse counts are `37499681/60` for C
and `36870437/60` for E, with 23 or 24 fixed-format retained classical bits
per attempt including the core archive. These are the **same screened code**
and the same uniform-seam selector as the current form above. They are
compiler resource counts, not energy or power. No earlier fixed-seam code
cost enters this packet.

## Finite theorem interface and physical boundary

- The checked pointwise `J` incidence supplies the source premise for
  [RegionalContinuity](../../Lean/Screen/RegionalContinuity.lean), including
  the actual signed seam orientation. Occupation continuity is not a stress
  tensor.
- [SourceActionMeasure](../../Lean/Geometry/SourceActionMeasure.lean)
  controls normalized reversible measures on a connected response graph.
  Its earlier uniform-source application remains valid; the present code
  state is not a stationary action measure.
- [InternalEnergyInertia](../../Lean/Geometry/InternalEnergyInertia.lean)
  and [ProperTimeInternalAction](../../Lean/Geometry/ProperTimeInternalAction.lean)
  require a selected composite action and proper-length/clock premises.
  Native current variance and compiler pulse counts do not supply them.
- [WardLimitManifest](../../Lean/Dynamics/WardLimitManifest.lean) requires
  a stress realization and controlled refinement beyond finite occupation
  continuity. No such map is constructed here.

The exact result is a **positive transient source-current form together
with a flat actively corrected logical transfer**. Neither the form's
eigenvalues nor the controller's attempted-repair rate is a measured
kinetic mass. Assigning kilograms or seconds requires an independent common
observer-clock and energy dictionary. `E=mc²` remains a subsequent
same-candidate check against independently measured rest action, moving
response and conserved stress. CP-04 remains open.

```sh
python -m model.screened_transport.transport --check
python -m pytest -q -W error model/screened_transport/test_transport.py
```
