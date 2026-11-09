# Vacuum-energy tests for coherent phase and attenuation

The configured coherent action admits exact energy tests without particle-data
inputs. Its free charged-CAR sea energy is phase-independent at zero momentum.
For every finite positive-weighted set of momenta, external minimization of
that free energy sets the pass amplitude to one. A concrete invariant sector
of the actual finite quantum-Higgs vertex has an exact phase-dependent energy;
external minimization sets its phase to zero and its pass amplitude to one.
The resulting coherent action has a repeated pair of eigenvalues.

Nonzero momentum in the supplied free field can favor an interior phase.
Rational bounds establish this effect and show that the ordering changes with
the kinetic scale. These results give usable selection tests. They supply no
autonomous phase register, preparation law, interacting vacuum or physical
lepton-mass prediction.

[`phase_attenuation.py`](phase_attenuation.py) uses the actual coherent
coefficients and quantum-Higgs tensor. No measured mass, calibrated phase,
calibrated attenuation, electron energy unit or target spectrum enters it.

## Existing results used as premises

| Source | Result used | Boundary |
|---|---|---|
| [Native coherent analyzer](../standard_model_source_map/KOIDE_SOURCE.md) | All128 raw records give the exclusive weights `(10,5,5,2,50)/72`; the complete analyzer supplies a coherent three-cycle action. | Its phase is a controller setting. |
| [Complete attenuator](../particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md) | The sixteen-flag instrument realizes every declared `0<=r<=1` and retains its loss outputs. | Realizing a value does not select it. |
| [Free charged field](../charged_lepton_response/RELATIVISTIC_MASS.md) | The same positive family coefficient has charged energies `sqrt(p²+|v_H|² y_j²)`, with both spins and all naive-lattice valleys retained. | The Higgs number, state, clock and free-field attachment are supplied. |
| [Finite quantum Higgs](../particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md) | The literal three-state scalar and full CAR vertex preserve the twelve gauge currents and `N_F-2N_H`. | A finite gauge-invariant vertex supplies no interacting continuum or selected vacuum. |
| [Source phase](../charged_lepton_response/SOURCE_PHASE_SELECTION.md) | The exact native ordering potential commutes with the coherent circulant; C/E reversal conjugates its phase. | Those operations do not select a nonzero unordered mass spectrum. |
| [Classical-record nonidentifiability](../particle_source_selection/PROOF.md) | Complete committed records are invariant under coherent Kraus rephasing. | The theorem permits additional accessible coherent source data. |
| [Capacity nonselection](../particle_candidate_atlas/CAPACITY_NON_SELECTION.md) | The displayed capacity closure persists under a small change of the supplied attenuation law. | Self-consistency does not select that law. |
| [Finite Koide theorem](../../Lean/ObserverPatchHolography/KoideCirculant.lean) | Positive balanced roots obey a phase-independent quotient. | A phase-independent identity cannot choose a phase. |

The computations below add energy identities to these source interfaces.
They do not replace the cited nonselection proofs by repeated searches over
their controls.

## Exact coherent moments

Let `R` be the literal regular three-cycle of the coherent analyzer, and set

\[
 a=\sqrt{10/72}=\frac{\sqrt5}{6},\qquad
 b=\sqrt{5/72}=\frac a{\sqrt2},\qquad
 C=aI+rb(zR+\bar zR^2),\quad Y=C^2/3,
 \qquad z=e^{i\theta}.                                      \tag{1}
\]

For `0<=r<=1` and `|theta|<=pi/12`, all roots are nonnegative:

\[
 c_j=a+2rb\cos(\theta+2\pi j/3),\qquad y_j=c_j^2/3.
                                                               \tag{2}
\]

The boundary at `r=1, |theta|=pi/12` has one zero root. The exact literal-matrix
calculation, using `R³=I` and `Tr R=Tr R²=0`, gives

\[
 T_1:=\operatorname{Tr}Y=\frac5{36}(1+r^2),                    \tag{3}
\]

\[
 T_2:=\operatorname{Tr}Y^2
 =\frac{25}{3888}+\frac{25r^2}{648}+\frac{25r^4}{2592}
   +\frac{25\sqrt2\,r^3}{1944}\cos(3\theta).                  \tag{4}
\]

For `r=1`, equation (4) is `(425+100 sqrt(2) cos(3theta))/7776`.
The positive-root quotient is `(1+r²)/3`; exact Koide balance therefore
corresponds to `r=1` on this family. The executable verifies (3)-(4) as
Laurent-polynomial identities and compares the full first three moments
against the owning attenuator at a complex rational unit phase.

## Free sea energy fixes no zero-momentum phase

Write `p` for the nonnegative kinetic norm `c|p_vector|`, or the corresponding
sine-lattice-symbol norm, and let `h=|v_H|>0` be the declared common Higgs
number. The charged contribution to the filled free-CAR sea is

\[
 E_{\rm sea}(p;r,\theta)
 =-2\sum_{j=0}^2\sqrt{p^2+h^2 y_j^2}.                         \tag{5}
\]

The factor two counts physical spins. The Nambu doubling is not counted a
second time. The source's massless spectators give a phase-independent term
and cancel in comparisons. Equation (5) concerns the same free charged
symbol as the cited field construction; it is not an interacting vacuum
energy obtained by replacing its quantum scalar with an expectation value.

At zero momentum, positivity of `Y` gives

\[
 E_{\rm sea}(0;r,\theta)=-2hT_1=-\frac{5h}{18}(1+r^2).
                                                               \tag{6}
\]

All admitted phases have the same energy. This conclusion uses the complete
charged CAR ground rather than a one-particle or one-family truncation.

The attenuation result extends to every momentum. For fixed `p>0`, define

\[
 f(x)=\sqrt{p^2+h^2x^4/9},\qquad
 f''(x)=\frac{2h^2x^2(3p^2+h^2x^4/9)}
 {9(p^2+h^2x^4/9)^{3/2}}\ge0.                               \tag{7}
\]

At `p=0`, `f(x)=h x²/3`, which is also convex. Each root has the form
`a+r u_j`, with `sum_j u_j=0`. Consequently
`S(r)=sum_j f(a+r u_j)` is convex and `S'(0)=0`. For the nontrivial source
cycle, `S` is strictly increasing on `r>0`. Thus (5) is minimized at `r=1`.
Summing with positive momentum multiplicities preserves the result. A
nontrivial attenuation `r<1` requires an additional term or preparation
principle if energy minimization is intended to select it. The argument does
not depend on measured masses or on any momentum cutoff value.

## Nonzero momentum supplies phase sensitivity, with scale dependence

At nonzero momentum, the square roots in (5) depend on more than `T_1`.
This produces an actual interior-phase effect, rather than a deduction from
Koide balance. Use the rational control without measured inputs

\[
 z_*=(63+16i)/65,\qquad 0<\arg z_*<\pi/12,\qquad r=h=1.
\]

Integer square-root enclosures and rational interval arithmetic prove:

| Declared kinetic norm | `E(z*)-E(1)` | `E(z*)-E(exp(i pi/12))` |
|---|---:|---:|
| `1/256` | approximately `-0.002260005521004324` | approximately `-0.00003028696865690325` |
| `1/64` | approximately `+0.005450663078193801` | approximately `-0.0004828229109534914` |

The machine certificate gives outward rational endpoints at thirty decimal
places; the signs do not depend on floating arithmetic. Since the energy is
continuous on the closed positive chamber, the first row proves existence of
an interior minimizing phase: both endpoints have strictly greater energy
than `z*`. It does not identify `z*` as the minimizer. The second row proves
that a change of kinetic scale reverses the comparison with zero phase.
No claimed global minimizer is inferred from a sampled grid.

These are exact charged-symbol controls. They supply no source-selected
momentum distribution or finite-volume boundary condition. A physical vacuum
calculation must include the entire prescribed momentum set, common cutoff,
all retained species, gauge/scalar interactions and the source's energy
reference. Selecting one momentum because its minimum agrees with a desired
phase would introduce a calibration.

## The actual finite quantum-Higgs interaction has a solvable sector

The finite source scalar has states `|0>, |1>, |2>` and creation operators
`B_h†=|h+1><0|`. Let the three families supply eighteen active CAR modes:
two spins on `L_0`, `L_1` and `e^c` per family. The literal source vertex is

\[
 W(Y)=\sum_{h,f,g}Y_{fg}B_h^\dagger
 (c_{L_h,f,\uparrow}^\dagger c_{e,g,\downarrow}^\dagger
 -c_{L_h,f,\downarrow}^\dagger c_{e,g,\uparrow}^\dagger),
 \qquad V=W+W^\dagger.                                      \tag{8}
\]

It is the same directed charged-Higgs convention as the owning source
vertex. The executable compares (8) on all192 scalar/CAR inputs for a single
family against `quantum_higgs_operators()`, including both occupied scalar
inputs. Its general-family action is defined on every occupation input.

The conserved number `Q*=N_F-2N_H` is distinct from electric charge. Its
complete `Q*=0` sector consists of the scalar/fermion vacuum `Omega`, together
with every occupied scalar and every two-fermion occupation. Its dimension is

\[
             1+2\binom{18}{2}=307.                            \tag{9}
\]

For every pair state, another application of `W` vanishes because
`B_h†|k+1>=0`. An application of `W†` either gives `Omega` or vanishes.
Thus no second pair can be created within this sector. This statement follows
from the actual finite scalar cutoff and full fermion occupancy; it would
require a different argument for an unbounded bosonic Higgs field.

Set `v=W Omega`. The distinct `(h,f,g,spin)` pair states are orthogonal, so

\[
 \langle v,v\rangle=4\operatorname{Tr}(Y^\dagger Y)=4T_2.     \tag{10}
\]

For a common declared number energy `epsilon_F N_F+epsilon_H N_H`, all
occupied-scalar/two-fermion states have the same energy
`d=2epsilon_F+epsilon_H`. The entire 307-state Hamiltonian has the exact block
form

\[
 H_{Q^*=0}=\begin{pmatrix}0&v^\dagger\\v&dI_{306}\end{pmatrix}.
                                                               \tag{11}
\]

The actual `build_join` control has default `epsilon_F=0, epsilon_H=1`, hence
`d=1`; equation (11) holds for every fixed real `d`. The two-dimensional span
of `Omega` and `v` is invariant. Its energies are

\[
 E_\pm=\frac{d\pm\sqrt{d^2+16T_2}}2.                        \tag{12}
\]

All305 orthogonal dark states are retained and have energy `d`. For the
nonzero coherent family, `E_-<min(0,d)`, so it is the ground energy of the
whole sector. The code checks every sector basis input, Hermiticity, the
bright-vector norm, `Hv=4T_2 Omega+d v`, and the full dark block, using a
complex off-diagonal family control. It also checks charge conservation on
inputs outside this sector. The other sectors of the full786432-state
scalar/CAR carrier are defined by the same vertex and are not deleted.

The source Ward identity makes both `Omega` and `v` gauge invariant. The
local quantum number terms are invariant as well. Spatial hopping, quantum
links and scalar transport need not preserve this two-dimensional span and
are excluded from (11). In particular, (12) is neither the global ground
over all `Q*` sectors nor the interacting spatial ground.

For `0<r<=1` and `|theta|<=pi/12`, the cosine coefficient in (4) is positive.
Its maximum is at `theta=0`. Every `r`-dependent term in (4) increases with
`r` on that chamber, since `cos(3theta)>=1/sqrt(2)`. Equation (12) decreases
strictly with `T_2`. External minimization therefore gives

\[
                        \theta=0,\qquad r=1.                  \tag{13}
\]

The supplied coherent family at (13) has action eigenvalue ratios
`1:1:(34+24 sqrt(2))`. Equation (13) gives no quantum-Higgs expectation value
and does not turn that ratio into a physical charged-lepton mass ratio. It
does show that this particular source-compatible vacuum-sector mechanism
cannot supply a nonzero phase or a nontrivial attenuation.

## A dynamical-control route and its required inputs

The coherent phase and attenuation in (1) label different supplied controller
programs. Their energy comparison alone does not make either variable a
dynamical observable. Even a direct sum of field Hamiltonians over a finite
phase register leaves the initial probabilities of its phase sectors
unchanged if no transition or preparation law acts on that register.

One possible construction makes the control a dynamical register. That
route requires the following specified data:

1. A finite coherent control register for the phase and attenuation, including
   its coupling to the complete analyzer and every retained loss flag.
2. A source transition, Hamiltonian mixing term or complete preparation law
   which changes or selects the register's state. A claimed ground or Gibbs
   preparation requires that state's construction and its common clock.
3. The same full field state, momentum/cutoff prescription and relative energy
   reference across register values. A scalar term `u(theta,r) I` leaves every
   fixed-control field probability unchanged and can change which control
   minimizes energy when the register is made dynamical.
4. The resulting complete observer mass readout, including the transfer from
   a quantum scalar state to any effective broken-phase coefficients.

An autonomous register is not required for every source derivation. A fixed
declared program can instead determine an effective phase directly through
its action and complete observer correlations. Architecture and software
are admissible design choices; uniqueness is not a prerequisite. A phase
chosen using a measured mass remains a calibration for that comparison,
even when stored as a program constant. Both routes require the actual
observer output and its physical observable map.

The nonzero-momentum effect gives one concrete interaction channel for such
a construction. The finite quantum-Higgs calculation supplies another
fully soluble test, whose minimum has the wrong multiplicity for a
three-distinct-eigenvalue action. These bounds constrain additional dynamics
without excluding every source extension or introducing particle-data fits.

Run the source and interval checks with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -W error -m model.particle_calibration_selection.phase_attenuation
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -W error -m pytest -q model/particle_calibration_selection/test_phase_attenuation.py
```
