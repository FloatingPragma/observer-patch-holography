# What the finite particle constructions identify

**The model energies, binding energy and inertial mass can be identified.**
The cube has eight exactly characterized interacting spectral energies;
its source-prepared reference and local phase kick reach each with positive
spectral weight. The moving pair
has dispersion epsilon(P)=P²/4-gamma² for the reflected generator
J_a=6/a²-H_a, so its model binding energy is gamma²
and its inertial mass is2 in those coordinates when hbar_model=1. These are derived model
quantities, not assignments to an observed species.

The [source-prepared charged defect](../ground_spectral_source/README.md)
provides a separate constructive chain: prepare the ground of the specified
full cube action, create a localized odd excitation, preserve its number
flag through the declared repair waits, and read its same-action spectral
response with complete pointer outcomes and a certified error. This identifies
a measurable excitation of that finite model. The
[identification route](ROADMAP.md) states what connects such a construction
to propagating particle sectors and measured species.

Converting them to joules or kilograms needs a common independently
identified energy unit Eunit and length unit Lunit. With
p_phys=hbar P/Lunit, the pair's inertial mass is
`M=2 hbar²/(Eunit Lunit²)` and its binding energy is `Eunit gamma²`.
The corresponding physical time unit must obey `Eunit Tunit=hbar`.
Neither the rest energy offset nor these physical units is fixed by the
bound-band theorem. A binding energy is not the total rest mass.

This audit tests specific identification maps, without fitting a mass scale
or treating simulated output as an observed particle. Owners: CP07/CP12;
physical assignments and new eligible CP11 comparisons both remain zero.

The concrete cube spectrum below uses the original H and contact unit U=1.
Decimals are displays of the exact root intervals retained in the receipt.
Weights are for the normalized centered density probe; the executed phase
target has 15/16 times these weights and an additional elastic weight1/16
at the reference energy3.

| H energy | Signed gap E-3 | Normalized probe weight |
| ---: | ---: | ---: |
| 1.816498920 | -1.183501080 | 0.002465366 |
| 1.860537501 | -1.139462499 | 0.005644624 |
| 1.955137632 | -1.044862368 | 0.005881219 |
| 2.108692846 | -0.891307154 | 0.007137664 |
| 2.167639734 | -0.832360266 | 0.004528123 |
| 3.015861346 | +0.015861346 | 0.421577939 |
| 3.030769653 | +0.030769653 | 0.415789141 |
| 3.044862368 | +0.044862368 | 0.136975924 |

These are energy eigenvalues and signed transitions, not eight rest masses.
The lowest line is not even the lowest N6 energy. Independently, the exact
one-particle identity h16²=3I and positive contact give
`-3sqrt3/8 <= E0_N6 <= -sqrt3/8`. The upper bound uses a fully
hidden-spin-polarized six-orbital Slater state, which has zero contact.
Thus all eight lines have strictly positive gaps from the actual N6
ground. The receipt gives rigorous gap intervals and separately records
an approximate sparse ground estimate of -0.366904154; its small residual
is a numerical diagnostic, not a certified enclosure. This N6 response
uses an excited reference. The separate
[full-Fock ground theorem](../interacting_vacuum/README.md) supplies a ground
across number sectors for the specified cube action; its physical-vacuum
identification requires a source and continuum argument.

| Required physical identification | Cube local density response | Moving bound pair |
| --- | --- | --- |
| Carrier and dynamics | Original open 2x2x2 Weyl/contact action; protected N6 sector under the specified complete four-block repair service | Original infinite-lattice two-body action with the declared critical coupling and global comparison J=6/a²-H |
| Number and grading | DeltaN=0 relative to N6 reference; even operator and state parity | N=2 relative to the empty Fock state; even parity |
| Electric charge and gauge representation | Number-neutral; no derived identification of N with electric charge or full physical gauge representation | Model number two; no electric-charge or physical gauge assignment |
| Spin and exchange | Local density probe has hidden spin zero; the larger eta3 sector contains hidden singlets and triplets | Antisymmetric two-fermion wavefunction, contact singlet source; no derived physical Lorentz spin or many-composite exchange law |
| Vacuum | Prepared reference R has H energy3 and is not a ground state, even within N6 | No common physical vacuum identified; energy reflection and subtraction are essential supplied conventions |
| Propagation and mass | Eight exact finite transition frequencies, not eight invariant rest masses | Controlled dispersion epsilon(P)=P²/4-gamma²; model inertial mass2 when hbar_model=1, not a physical relativistic rest mass |
| Physical clock, units and detector | Complete finite Ramsey detector source-executed; certified accelerated emulation resolves weights on the supplied eight-line support; no physical clock or unit calibration | Explicit bound-band isometry and propagator estimate; no physical detector or unit calibration |
| Interactions and empirical comparison | Contact-shifted finite energies; no observed-species comparison | One composite sector, not interacting many-composite dynamics; no observed-species comparison |

The Standard Model comparison requires physical quantum numbers and a common
energy-momentum interpretation, not only a count of spectral lines. For
example, the [PDG electron listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-electron.pdf)
specifies spin one half; the [PDG quark review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-quark-model.pdf)
identifies quarks as spin-one-half fermions and gives their additive quantum
numbers. The [PDG kinematics review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf)
uses E²-|p|²=m² in c=1 units. These references supply the target conventions;
their numerical masses are not inputs to this audit.

## Exact grading exclusion

Write Gamma=(-1)^N. The reference, phase-kicked state and every spectral
component of the cube density probe lie in N6. Their transition operators
|e_N6><R_N6| commute with Gamma and with N. `grading_certificate` checks the
entire 8008-dimensional number sector, all56 reference occupation columns,
and the local phase on all65536 occupation addresses.

Suppose a proposed dictionary is an isometry J intertwining model and physical
grading: Gamma_phys J=J Gamma. Both J|R> and J|q> then have the same parity.
An odd physical creation operator A instead anticommutes with Gamma_phys,
so A J|R> has the opposite parity and cannot equal the nonzero J|q>.
This excludes identifying this density response with a single odd fermion
under a grading-preserving dictionary. It does not exclude other model odd
excitations, a fermion-antifermion response, a composite or a different map.
Even parity alone proves neither Bose exchange nor integer Lorentz spin.
Number neutrality is likewise a model fact, not an electric-charge assignment.

## The prepared reference is not a vacuum spectral ground state

The positive greater/Wightman spectral measure of the normalized density probe
is supported at omega=z/16-1. Its exact minimal polynomial is

`(z²-16z-12)(z³-16z²-12z+64)(z³-16z²-12z+128)`.

The quadratic has one root below0 and one above16. Each cubic has one root
below0, one between0 and16, and one above16. The exact nonzero Krylov
determinant makes all eight spectral weights positive. There are therefore
**five negative and three positive** frequencies, with no zero frequency.
This statement concerns a positive spectral measure, not a retarded
commutator whose negative frequencies would occur even in a ground state.

For example, z_minus=8-2sqrt19 has strictly positive probe weight
`(1-4/sqrt19)/14`. Its normalized spectral component has energy
`5/2-sqrt19/8 < 2`, below the reference energy3. A unitary swapping this
state and R, acting identically on the orthogonal complement, preserves N
and lowers the mean energy by more than1. Thus R is not passive for the full
number-preserving field unitary algebra. This is not an executed work
extraction protocol: the implemented local phase itself preserves the mean
energy3. The source preparation and protected spectral-response results
remain valid.

Adding any scalar or any common number chemical potential to H leaves these
N6 transition differences unchanged. Global reversal of time exchanges the
two sign counts and still gives mixed signs. These changes do not make R a
ground state. A different reference or physical Hamiltonian requires a new
identification argument.

## Exact quadratic-band mass-shell exclusion

Permit any fixed affine energy map `E=alpha*epsilon+beta`, alpha nonzero,
and any invertible fixed real linear momentum map `p=L P`. Along P=r e1,
the residual `E²-c²|p|²-m²c⁴` has quartic coefficient `alpha²/16 > 0`.
The other terms have degree at most two. A polynomial vanishing on any
open momentum set would vanish identically, so the quadratic limiting band
cannot be an exact massive or massless Lorentz shell under these maps.
This conclusion is independent of a fitted mass, gamma, the scalar energy
origin, or anisotropic linear momentum units.

It does not reject a nonrelativistic limit. After matching the constant and
quadratic coefficients, let E_quad=E0+A r² and
E_rel=sqrt(E0²+2E0 A r²), with rest energy E0>0 and A>0. Their exact difference is

`A² r⁴/(E_quad+E_rel)`.

For 0<=r<=rmax, `low_momentum_gap` encloses it between
`A²r⁴/[2(E0+A rmax²)]` and `A²r⁴/(2E0)`. The discrepancy vanishes at low
momentum; no finite-precision experiment or alternative scaling is excluded
without its own domain/error analysis. The theorem does not apply to an
unspecified nonlinear momentum dictionary, whose physical validity would
have to be established separately.

## Existing comparisons remain retrospective

The retained `particle_comparison` packet assigns one measured electron mass
as the common energy scale to its dimer particle/hole doublets. Their equal
gaps give a conditional particle-antiparticle equality, but zero independently
predicted absolute masses. Its finite dictionary also lacks the physical
vacuum, charge, Lorentz spin and continuum pole joins. The exact672-case
three-lepton exclusion remains scoped to that dimer frequency menu and its
common scale/chemical-potential map; it is not applied to the new cube.

This audit reads and hashes the existing receipt and observations without
rerunning its empirical fit, modifying its exposed-target history, or adding
an eligible physical comparison. The new cube and pair candidates do not
inherit the earlier assigned units or particle names. Physical identification
requires a common vacuum/clock, actual charge and spin maps, propagating
one-particle poles, and independent dimensionless tests with errors.

Run from the repository root:

```sh
.venv/bin/python -W error -m pytest model/particle_identification -q
.venv/bin/python -W error -m model.particle_identification --check model/particle_identification/results.json
```

The receipt binds the source modules, exact tests, supplied-reference audit
and this mathematical scope. It makes no whole-OPH impossibility claim.
