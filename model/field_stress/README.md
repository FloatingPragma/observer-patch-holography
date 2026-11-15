# Energy continuity and entropy on the interacting field

This CP-09 / [#1145](https://github.com/FloatingPragma/observer-patch-holography/issues/1145)
and CP-13 / [#1137](https://github.com/FloatingPragma/observer-patch-holography/issues/1137)
component constructs finite energy currents and a quantum entropy first law
on the **same** [Weyl/contact action](../interacting_field/DERIVATION.md)
used by the particle calculations. It addresses missing premises without
declaring that the conditional Einstein theorem has been instantiated.

The retained open cube has 27 sites, 54 CAR modes and the complete
1,431-dimensional two-particle sector. Its action is
`H=dGamma(h)+sum_x D_x`, `h(k)=-(1/8) sigma.sin(k)`,
`D_x=n_(x,0)n_(x,1)`. The spatial lattice, action and CAR statistics are
supplied. Every time derivative below uses Hamiltonian **model time** with
the supplied convention `hbar=1`. No physical energy scale is inferred.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.field_stress --check model/field_stress/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/field_stress
```

## An actual conserved energy and its boundary flux

Assign half of every Hermitian hopping bond to each endpoint:

```
h_x = U D_x + (1/2) sum_(b incident at x) H_b,
sum_x h_x = H,
J_(x->y) = i [h_x,h_y].
```

Each `h_x` is Hermitian, hence `J_(x->y)` is Hermitian and antisymmetric in
its endpoints. Heisenberg evolution gives the operator identity

```
d h_x/dt = i[H,h_x] = -sum_y J_(x->y).
```

Summing over any site set A cancels every internal current. Its energy
changes by minus the **complete** outward boundary flux. These statements
hold for every state and untouched reference, including a propagating bound
pair, and do not use the Einstein equation. The energy assignment has a
declared ambiguity: another division of bond energy changes individual
currents. Here overlapping stars can produce currents between sites at
distance two; they must not be dropped by assuming only nearest-neighbor
energy transport. Disjoint even CAR operators commute, so greater distances
give zero. The exact construction applies to every retained finite cube;
the full operator replay here is R=1.

[energy.py](energy.py) assembles the local terms independently of the full
Hamiltonian sum, then checks the identities using sparse **Gaussian-integer
arithmetic**. Local energies have denominator 32 and currents denominator
1,024. Real and imaginary numerator matrices use signed integers; products
are bounded before evaluation to exclude overflow. All 351 unordered site
pairs and 27 continuity equations are checked. There are 81 nonzero current
pairs, with maximum range two. The 18-site half-cube has 18 nonzero crossing
current pairs. The U=0 deletion retains the kinetic terms and changes the
interaction-dependent energy-transfer channels. This is energy accounting
for the actual interaction, rather than a number current renamed as stress.

[experiment.py](experiment.py) prepares a contact/separated-pair
superposition with an existing compiled hopping word and executes the
same full field action. The initial local energy flux is nonzero. Exact
Hamiltonian evolution conserves total energy. A finite compiled palindrome
is an approximation and need not conserve that energy exactly. With
`B=6c+U=7/4`, duration t, and n subdivisions, the existing ideal word bound
is `epsilon=B^3 |t|^3/(3 n^2)`, so

```
|<H>_split - <H>_exact| <= 2 B epsilon.
```

This follows by inserting and subtracting one evolved vector in the energy
pairing; both ideal evolutions are unitary. The receipt retains both split
counts, preparation/routing resources, full-sector sparse references and
numerical source residuals. Floating contractions, sector compression and
physical hardware errors are separate from the rational ideal-word bound.
The construction concerns unitary field evolution: detector backaction and
externally driven control work require their own energy accounting.

## A first law on this action and algebra

The faithful density `tau=exp(-beta H)/Z` is constructed on the **full
finite N=2 matrix algebra**, with declared `beta=log(25/8)`. This choice
permits comparison with a two-level clock ratio; it does not prepare the
full field density from native repair or identify the whole native clock
with this Hamiltonian. Negative kinetic energies cause no difficulty for a
finite-sector Gibbs state; the state is not declared a physical vacuum.

Let delta be any Hermitian trace-zero variation. Spectral functional
calculus gives `K=-log(tau)=beta H+log(Z) I`. The derivative of von Neumann
entropy is

```
DS_tau(delta) = -Tr[delta (log(tau)+I)]
              = Tr(delta K) = beta Tr(delta H).
```

No commutation between delta and tau is required. More strongly, for any
density rho on the same algebra,

```
S(rho)-S(tau) = beta Tr[(rho-tau)H] - D(rho || tau).
```

This follows directly by expanding `D=Tr rho(log rho-log tau)` and using
unit trace. The receipt tests the noncommuting, physically specified
direction `delta=|origin contact pair><origin contact pair|-tau` and
`rho=tau+delta/64`. It retains the reference spectral bounds, the nonzero
commutator diagnostic, the nonzero entropy and energy derivatives, the
relative-entropy remainder and their numerical residuals.

For `rho_s=tau+s delta`, `0<=s<=epsilon<1`, positivity gives
`rho_s >= (1-epsilon) lambda_min(tau) I`. The logarithm derivative has
divided differences bounded by the reciprocal of this lower eigenvalue.
Consequently the entropy Hessian has absolute quadratic form at most
`||delta||_HS^2/((1-epsilon)lambda_min(tau))`. Taylor's integral remainder
therefore bounds

```
0 <= D(tau+epsilon delta || tau)
   <= epsilon^2 ||delta||_HS^2 /
      [2 (1-epsilon) lambda_min(tau)].
```

The analytic bound is conservative. Its receipt value uses numerical
eigenvalues and is not an outward-certified number. The algebraic identity
and its proof are distinct from these floating diagnostics.

The nearest existing result is
[`thermoFirstLawData_passes`](../../Lean/Thermodynamics/EinsteinPremiseLink.lean),
which proves the finite classical package on its own simplex tangent
space, conditional on its reference split. This component supplies an
explicit **noncommutative** matrix-algebra attachment to the same
interacting action used for the particles. It does not claim a general
entropy theorem or a Lean proof. Defining a formal bulk pairing as
`beta Tr(delta H)/(2*pi)` would not identify a geometric cap modular charge
or an edge entropy, so we do not use that relabeling to discharge gravity.

## Regional restriction is a substantive additional map

The [CAR restriction](../field_regions/README.md) traces the same full
reference to the origin and its +x neighbor, including fermionic reordering
signs. All local particle-number sectors 0, 1 and 2 are retained. This is
an 11-dimensional positive support. The density is not faithful on the
entire 16-dimensional four-mode Fock space, where sectors 3 and 4 have
zero weight.

On this positive support, compute `K_A=-log(tau_A)` and compare it with
the bare internal regional Hamiltonian `H_A`, which contains no crossing
bonds. Even the best affine fit `a I+b H_A+mu N_A` fails. Allowing an
additional `N_A^2` term for canonical number constraints leaves a
nonzero residual. A nonzero `[K_A,H_A]` directly excludes an affine bare
energy/number identification. These are numerical finite-model controls,
not a universal regional modular no-go. Fixed total number, boundary
couplings and the chosen region matter; no claim attributes the entire
mismatch to the contact interaction alone. Other references, corrected
regional generators and scaling limits are possible.

## The remaining Einstein premises

The actual consumer is
[`EinsteinRegisterData`](../../Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean).
Its fields require much more than conservation of a finite Hamiltonian or
the Gibbs identity. The field-to-theorem attachment here is:

| Premise | Supplied here | Open requirement |
| --- | --- | --- |
| PR-23 / G2, G4 | Same-action quantum modular first law on its declared full-sector Gibbs algebra. | Geometric cap modular flow, diamond kernel, edge-center split and a justified native reference map. |
| PR-24 / G3 | Exact local energy continuity, including interaction energy and boundary currents. | One symmetric Lorentz stress tensor and its Ward identity under the theorem's four chart steps; half-sided inclusion/null translations. |
| PR-25 / G2, G6 | No null-charge identification is inferred from the energy current. | Common stress feeding entropy and null tomography, symmetry, geometry/stress null matching. |
| PR-26 / G4 | Actual entropy tangent and a controlled finite relative-entropy remainder. | Fixed-proper-volume area variation and generalized-entropy stationarity. Gibbs stationarity alone does not supply them. |
| PR-27 / G1, G5 | The supplied finite cubic domain is stated. | Lorentzian/geodesic geometry, Bianchi identity, enough observer directions and one controlled shrinking-ball family. |
| PR-28 / G6 | No physical normalization added. | One common physical refinement tower, vacuum residue, universal coupling and an independent positive Newton scale. |

Thus CP-09 and CP-10 are open. In particular, `EinsteinRegisterData`
has not been constructed and the Einstein conclusion is not invoked.
The energy continuity identities and the same-action entropy map are
reusable concrete inputs for that further construction. No physical
comparison or integrated acceptance gate is added.
