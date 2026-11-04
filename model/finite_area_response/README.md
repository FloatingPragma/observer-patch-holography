# Matter changes the finite geometry distribution

The existing reciprocal scalar model gives a controlled finite response of
geometric area to a matter disturbance. Randomizing part of the field leaves
the geometry distribution unchanged at preparation. One application of the
unchanged thermal law then changes the mean two-face area by
`-0.0006270566389 epsilon`, with an exact rational enclosure. Both geometry
branches have the same proper volume. This is a thermal two-geometry
construction; identifying a common interacting continuum and Einstein
dynamics requires further mathematical input.

## Specified geometry, matter and update law

The calculation uses the scalar branch of
[reciprocal_vacuum](../native_equilibration/README.md#reciprocal-thermal-matter-and-geometry).
Its geometry register has two values `r0=9/10`, `r1=11/10`; its field has
three oscillators at cutoff two, hence eight states per geometry. The
native angular band is six and `beta=1/sqrt(14)` is supplied. All sixteen
joint states participate. Geometry-dependent field Hamiltonians and the
thermal reservoir/readout services are specified by that model.

In a normalized local angular frame the action density is
`C(r)=diag(-1/r,1/r,r^2,r^2)`. The independent metric dictionary
`g=sqrt(-det C) C^-1` gives

```
g(r)=diag(-r^2,r^2,1/r,1/r),
normalized proper volume=1,     two-face area A(r)=2/r.
```

The geometries therefore have areas `20/9` and `20/11`. These are geometric
readouts of the declared radial/angular action family. The finite result
does not identify that family with a measured spacetime or the separate
Weyl/contact particle construction.

The update counter is one application of the specified thermal law. Its
Hamiltonian and supplied `beta` use the common model-coordinate clock.
The metric has lapse `r`, so a fixed branch has `d tau=r dt`, proper-time
generator `H/r` and proper inverse temperature `r beta`. A shared proper
temperature, a physical duration per update and laboratory calibration are
not supplied by the comparison.

Let `D=diag(2 popcount(mask)+3)` and `C_field=XXI+IXX`. Direct evaluation of
the original compressed Hamiltonians gives the exact rational matrices

```
M_g=beta H_g=a(r_g) D+c(r_g) C_field,
a(r)=r/4+1/(7r)+3r^2/28,       c(r)=-1/(7r).
```

They do not commute. They are the actual compressed oscillator action,
including its top-level diagonal, rather than the square of a truncated
coordinate operator.

For complete energy labels `i` and `j` in opposite geometries, the specified
reciprocal update has transition probability

```
T_(j,i)=(1/16) w_j/(w_i+w_j),        w_i=exp(-beta E_i).
```

All rejected and same-geometry proposals contribute to the diagonal.
The identity `w_i T_(j,i)=w_j T_(i,j)` proves stationarity of the full joint
Gibbs state. The source prescription includes energy readout, auxiliary
native proposals, a fixed-gap bath qubit for each of the 64 cross-geometry
pairs, thermal refresh and retained records. It does not execute the entire
reservoir-bank controller as one native source history. This result proves
the response of that specified complete finite law; it adds no claim of
such an execution or of autonomous reservoir preparation.

## A matter disturbance under the unchanged law

Start at the joint Gibbs state `rho`. Apply to the field the complete channel

```
R_epsilon=(1-epsilon) identity + epsilon [X -> Tr(X) I8/8],
0<=epsilon<=1.
```

It leaves the geometry marginal `p_g=Z_g/(Z0+Z1)` unchanged exactly. The
replacement has Kraus operators `|i><j|/sqrt(8)` for all 64 ordered pairs;
the mixture includes a separate identity flag. Every outcome, including
zero-probability endpoint outcomes, is retained. A possible dilation swaps
the field with three supplied maximally mixed bits and retains the displaced
input. Initialization, mixing and memory are explicit apparatus services.
No native execution of this combined intervention is asserted.

Apply one original reciprocal update, with its original Hamiltonians, bath
gaps, energy readout and transition matrix. Define

```
S=sum_(i in g0,j in g1) w_j/(w_i+w_j).
```

The reset input has populations `p_g/8`. Summing all 64 cross-geometry
transitions and retaining the rejection probabilities gives exactly

```
Delta p1=epsilon [S-64 p1]/128,
Delta <A>=epsilon (20/11-20/9)[S-64 p1]/128.             (1)
```

This equality is affine for every admitted `epsilon`, not merely an
infinitesimal fit. Exact enclosures prove

```
-0.000628 < Delta <A>/epsilon < -0.000627
```

for nonzero `epsilon`. The enclosure midpoint is approximately
`-0.0006270566389034`. The geometry-1 probability increases at rate about
`0.001551965181286`. Since both branches have identical proper volume,
the full intervention/update also preserves its expectation.

The disturbance supplies energy: the same source enclosures give
`2.697 < Delta <H>_preparation/epsilon < 2.698` in Hamiltonian units.
That transfer must be supplied by the preparation apparatus. Its division
into work and heat is not derived. Bath energy changes during the subsequent
update obey the original energy-exchange ledger. Thus the area response
does not rely on costless matter randomization or an unrecorded geometry reset.

## A separate equilibrium susceptibility

For the family `H_g(lambda)=H_g+lambda n_left`, with the same occupation
operator in both geometry blocks, cyclicity of trace in the differentiated
matrix exponential gives

```
d<A>/dlambda=-beta Cov_rho(A,n_left)
 =-beta p0 p1 (A0-A1)(<n_left>0-<n_left>1).             (2)
```

The geometry observable is scalar on each block, so it commutes with both
Hamiltonians and the perturbation. The perturbation does not commute with
either Hamiltonian; no formula `rho'=-beta O rho+...` for the full density
matrix is assumed. Explicitly,
`d exp(-beta H)/dlambda=-beta integral_0^1 exp(-(1-u)beta H) O exp(-u beta H) du`.
Multiplying by the central `A` and taking the trace reduces this integral
to `-beta Tr(A O exp(-beta H))`, which proves (2).

The certified susceptibility lies between `-0.000675` and `-0.000673`,
with midpoint approximately `-0.0006739995630072`. This calculation changes
the Hamiltonian family. Its energy readout and bath gaps must follow the
perturbed Hamiltonian under the existing family prescription. It is an
equilibrium susceptibility, separate from the fixed-law state response (1).
No native realization of the whole parameterized reservoir family is claimed.

## Exact bounds and controls

The characteristic polynomials of the rational matrices isolate all eight
levels per branch. Independent rational LDL congruences of `M_g-xI` check
the number of eigenvalues on each side of every endpoint. Every interval
contains exactly one eigenvalue and all sixteen levels are covered.
Rational 60-term exponential sums with an explicit geometric tail enclose
every thermal weight and transition rate. Arithmetic rounds outward at
denominator `10^30`; floating eigenvalues do not establish the signs.

A separate rational matrix-exponential recurrence encloses the diagonal
of `exp(3I-M_g)` without eigenvectors or spectral intervals, supplying the
occupation expectations in (2). The shift is common to both geometries and
cancels from probabilities. For either scalar or matrix series with norm
bound `B<=3`, the tail is at most
`B^61/[61! (1-B/62)]`. The source matrices, signs, positivity, inertia counts
and arithmetic bounds are replayed from their definitions.

Equal conditional Hamiltonians, constant area and `beta=0` give zero
fixed-law response. An equilibrium input without the field disturbance is
stationary. Independent tests use the complete original transition/record
law, full quantum densities and matrix-exponential Frechet derivatives.
They retain all rejected proposals and check that the geometry changes
only after the thermal update.

## Entropy and the Einstein boundary

The finite equilibrium obeys
`S(sigma)-beta Tr(sigma H)=log Z-D(sigma||rho)`.
It supplies entropy/free-energy stationarity on a specified thermal family.
It does not identify that functional with geometric generalized entropy.
The two geometry labels have one-dimensional internal blocks, so their
internal edge contribution `sum_g p_g log d_g` is zero. Since
`A0-A1=40/99`, equations `log d_g=eta A_g+c` for both labels force `eta=0`.
Those particular internal dimensions cannot supply a positive area-entropy
coefficient. The Shannon entropy of the geometry distribution is a separate,
generally nonzero term; no exclusion of all central or edge assignments
follows from this calculation.

A common interacting Einstein construction needs a microscopic
edge/area relation with derived normalization, a common stress/modular
identification and state/geometry stationarity under the required volume
variations, with controlled limiting errors. Finite responsive geometry
probabilities are compatible with that research route, but do not supply
those missing identifications.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.finite_area_response --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/finite_area_response
```

`--write` generates the source-bound receipt. Replay checks inherited source
hashes and the finite proofs above; it does not repeat the inherited native
pulse experiments or alter their receipts.
