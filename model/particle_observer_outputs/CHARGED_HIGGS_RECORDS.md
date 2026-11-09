# Complete charged records of the finite quantum-Higgs source

The same Hamiltonian and prepared sector ground as
[the neutral readout](HIGGS_RECORDS.md) give exact charged addition and
removal responses. The complete addition sector has 1650 states, including
1614 dark states. The complete removal sector has 36 states. A compensating
charged probe makes the experiment gauge invariant; its quantum reference
coordinates and both addition/removal outcomes are retained with their
original probabilities.

[`charged_higgs_records.py`](charged_higgs_records.py) derives these outputs
from the literal finite scalar/full-CAR vertex. The displayed spectra are
finite interacting energy differences in the same program clock. They are
not labelled as measured lepton masses, relativistic particles or
antiparticles. The family coefficient and number energies are declared
architecture data with no particle-mass inputs.

## The same operator on the adjacent sectors

Keep the actual vertex `W(Y)` and all eighteen active CAR modes. Let

\[
 H=W+W^\dagger+\epsilon_FN_F+\epsilon_HN_H
   =W+W^\dagger+dN_H+\epsilon_F Q^*,\quad
 d=2\epsilon_F+\epsilon_H,\quad Q^*=N_F-2N_H.             \tag{1}
\]

Write `G=4 Tr(Y†Y)>0`, `Delta=sqrt(d²+4G)`, `E0=(d-Delta)/2`, and
`alpha=(1+d/Delta)/2`, `beta=(1-d/Delta)/2`. The same prepared ground is

\[
 |\Psi_0\rangle=\sqrt\alpha|\Omega\rangle
              -\sqrt{\beta/G}\,W|\Omega\rangle.        \tag{2}
\]

Its domain is the full 307-state `Q*=0` sector. Adding or removing one
fermion changes `Q*` by `+1` or `-1`.

In `Q*=-1`, every state has an occupied scalar and one fermion. There are
`2*18=36` states. Creation by `W` is blocked by the finite scalar cutoff;
annihilation by `W†` needs two fermions. The entire vertex vanishes, so

\[
                     H_-=(d-\epsilon_F)I_{36}.          \tag{3}
\]

In `Q*=+1`, the scalar-vacuum block has one fermion and dimension 18. The
occupied-scalar block has three fermions and dimension
`2 binomial(18,3)=1632`. Let `W_1` be the literal source map from the first
block to the second. On all 1650 states,

\[
 H_+=\epsilon_F I+
       \begin{pmatrix}0&W_1^\dagger\\W_1&dI_{1632}\end{pmatrix}.
                                                               \tag{4}
\]

The executable constructs `W_1` by applying the original full-CAR source
vertex to each one-fermion occupation. It also replays that vertex on
all 1650 addition inputs and all 36 removal inputs. The remaining charge
sectors of the full 786432-state scalar/CAR carrier are not assigned spectra
by (3)-(4).

## An exact small Gram matrix retains the large complement

Set `w=W Omega`. Since each pair-creation term has even fermion parity,
`W c_i† Omega=c_i† w`. Thus

\[
 S_{ij}:=(W_1^\dagger W_1)_{ij}
 =G\delta_{ij}-\langle w|c_j^\dagger c_i|w\rangle
 =G\delta_{ij}-\gamma_{ij}.                             \tag{5}
\]

The actual source spin singlet gives, for each spin and weak component,

\[
 \gamma_L=YY^\dagger,\qquad
 \gamma_e=2(Y^\dagger Y)^T.                             \tag{6}
\]

The transpose in the charged-singlet family block follows from the
all-left creation-pair convention. A complex non-Hermitian family control
in the tests rejects the formula without that transpose. The trace is
`Tr gamma=8 Tr(Y†Y)=2G`, as required for a two-fermion state of norm squared
`G`.

If `s_j` are the squared singular values of `Y`, the eigenvalues of `S` are

\[
 \lambda_{L,j}=G-s_j\quad(4\text{ copies}),\qquad
 \lambda_{e,j}=G-2s_j\quad(2\text{ copies}).             \tag{7}
\]

Multiplicity in the family eigenvalues is retained. Since
`s_j<=Tr(Y†Y)=G/4`, equation (7) gives `S>=G I/2>0` even when some family
singular values vanish. Therefore `W_1` has rank 18. Its range and the
one-fermion block give 36 bright directions; its orthogonal complement
has dimension `1632-18=1614`.

For each normalized eigenvector `u` of `S`, the two-coordinate map
`(a,b) -> a(u,0)+b(0,W_1u/sqrt(lambda))` is an isometry onto the
corresponding pair of directions. Equation (4) then has matrix
`epsilon_F I+[[0,sqrt(lambda)],[sqrt(lambda),d]]`. Its energies are

\[
 E_{\lambda,\pm}=\epsilon_F+
                {d\pm\sqrt{d^2+4\lambda}\over2},
 \qquad E_{\rm dark}=\epsilon_F+d.                      \tag{8}
\]

All 1614 dark states remain in the complete spectral decomposition. Reducing
the eigenvalue calculation to `S` does not delete them from the operator.

## Charged spectral weights from the actual prepared state

For the normalized mode `u`, creation on (2) gives the two-coordinate vector

\[
 c_u^\dagger|\Psi_0\rangle
  \longleftrightarrow
  \begin{pmatrix}\sqrt\alpha\\-\sqrt{\beta\lambda/G}\end{pmatrix}.
                                                               \tag{9}
\]

It has no dark component for this particular ground input. Let

\[
 \delta_\lambda=\sqrt{d^2+4\lambda},\quad
 A_\lambda=\alpha+\beta\lambda/G,\quad
 X_\lambda=-4\lambda/\Delta-d(\alpha-\beta\lambda/G).
\]

The two addition weights are

\[
 a_{\lambda,-}=\tfrac12(A_\lambda-X_\lambda/\delta_\lambda),
 \qquad
 a_{\lambda,+}=\tfrac12(A_\lambda+X_\lambda/\delta_\lambda).
                                                               \tag{10}
\]

These are squared projections of (9) onto orthonormal eigenvectors and are
nonnegative. Their energies relative to the same prepared reference are

\[
 \omega_{\lambda,\pm}=E_{\lambda,\pm}-E_0
       =\epsilon_F+{\Delta\pm\delta_\lambda\over2}.     \tag{11}
\]

Removal gives one frequency and its original weight,

\[
 \omega_{\rm rem}=d-\epsilon_F-E_0
                  ={d+\Delta\over2}-\epsilon_F,
 \qquad r_\lambda=\beta(G-\lambda)/G.                  \tag{12}
\]

For a mode and its conjugate removal convention,
`a_(lambda,-)+a_(lambda,+)+r_lambda=1`. This is the full CAR spectral
sum, rather than a normalization of whichever pole is largest. The
addition correlator is the sum of the two terms
`a_(lambda,sign) exp(-i omega_(lambda,sign) t)`; the ordered removal
correlator is `r_lambda exp(-i omega_rem t)`. A retarded correlator uses its
appropriate opposite-time removal term. The two ordered functions are
not silently assigned the same signed-frequency convention.

The exact addition Green matrix can also be calculated without diagonalizing
the large sector. For absolute spectral parameter `z`, set
`x=z-epsilon_F`. Then

\[
 F^\dagger(z-H_+)^{-1}F
 =\left[\alpha(x-d)I+
      \left({\beta x\over G}-{2\over\Delta}\right)S\right]
      [x(x-d)I-S]^{-1},                                \tag{13}
\]

where the columns of `F` are the actual `c_i† Psi0`. The tests solve the
full 1650-state linear equation at a nonreal spectral parameter and compare
its result with (13). The poles and residues in (8)-(10) are consequences
of that same resolvent.

## A complete gauge-invariant observer experiment

A bare charged `c_i` is not a gauge-invariant pointer. Introduce a charged
probe whose addition reference transforms in the dual of the actual
18-mode one-fermion representation, and whose removal reference transforms
in the direct representation. Keep those reference coordinates coherent.
The two probe maps are

\[
 K_{\rm add}\psi={1\over\sqrt{18}}
      \sum_i c_i^\dagger\psi\otimes|i^*\rangle,
 \qquad
 K_{\rm rem}\psi={1\over\sqrt{18}}
      \sum_i c_i\psi\otimes|i\rangle.                  \tag{14}
\]

Contracting a representation with its dual makes these maps gauge
intertwiners. The probe has odd parity in each charged reference sector;
the field parity change is therefore compensated. Only the addition or
removal flag is committed at this step. Committing a weak-component basis
label would be another, gauge-dependent instrument and is not done here.

The full CAR identities give

\[
 K_{\rm add}^\dagger K_{\rm add}=I-N_F/18,
 \qquad K_{\rm rem}^\dagger K_{\rm rem}=N_F/18,
 \qquad \sum K^\dagger K=I.                            \tag{15}
\]

On all 307 admitted inputs, these effects are `I-N_H/9` and `N_H/9`.
Every component map is constructed from its literal occupation signs and
checked against (15). For the prepared ground, the original event weights
are `p_add=1-beta/9` and `p_rem=beta/9`; no branch is postselected.

Follow (14) with the complete spectral projectors of `H_+` or `H_-` on the
field while retaining the quantum probe and every spectator reference.
The energy projectors commute with the source gauge generators. A line
with multiplicity `m` has unconditional addition probability
`m a_(lambda,sign)/18`, and removal probability `m r_lambda/18`.
Their sum over all channels is one. The dark addition outcome, of rank 1614,
is retained even though its weight on the ground is zero. Every claimed
energy uses the same Hamiltonian parameter `t` as the neutral experiment.

The probe Hamiltonian is declared zero during this field-energy readout.
Its preparation and coupling are apparatus services of the chosen model.
The [complete readout contract](READOUT_CONTRACT.md) supplies the finite
pointer construction and reference-stable composition bound after the
required larger field/probe register is included. No compiled pulse tape
for (14) is claimed here. If an ambient field input lies outside the
specified charge sectors, the apparatus retains an outside-sector outcome
without assigning it one of the energies in (8).

Equation (1) makes an additional energy convention explicit. Holding `d`
fixed while changing `epsilon_F` leaves every neutral-sector output
unchanged and adds `epsilon_F Q*` to the Hamiltonian. Addition energies then
shift by `+epsilon_F`, and removal energies by `-epsilon_F`. The chosen
source specifies this number energy and the probe reference; a neutral
gap alone does not determine cross-sector offsets. This is a construction
input, rather than a requirement to prove the architecture unique.

## Concrete output and controls

The exact fixture `Y=diag(1,2,3)`, `epsilon_F=0`, `epsilon_H=1` has `G=56`
and `Delta=15`. Its removal energy is 8. Its six addition families have

| Internal source component | `s_j` | `lambda` | Addition energies relative to `E0` |
|---|---:|---:|---|
| left doublet | 1 | 55 | `(15 +- sqrt(221))/2` |
| left doublet | 4 | 52 | `(15 +- sqrt(209))/2` |
| left doublet | 9 | 47 | `(15 +- sqrt(189))/2` |
| charged singlet | 1 | 54 | `(15 +- sqrt(217))/2` |
| charged singlet | 4 | 48 | `(15 +- sqrt(193))/2` |
| charged singlet | 9 | 38 | `(15 +- sqrt(153))/2` |

The doublet rows have multiplicity 4 and the singlet rows multiplicity 2.
Their weights follow (10) and (12), and the certificate prints every
unconditional probability. The total removal probability is `7/135` and
the total addition probability `128/135`.

The certificate also evaluates the unchanged coherent analyzer at
`z=1` (phase angle zero), pass amplitude 1, `epsilon_F=0`, `epsilon_H=1`.
Its own family matrix and
`G=(425+100 sqrt(2))/1944` are used throughout. The fixture does not replace
those coefficients. No measured lepton ratio enters either calculation.

Dark response is tested on another admissible neutral input
`W(X)Omega` with `Tr(Y†X)=0`. Its complete addition instrument has a
strictly positive probability of the 1614-dimensional dark energy outcome.
Deleting that outcome because the ground does not populate it would lose
probability on an allowed gauge-invariant preparation. Zero singular family
values are also retained: at `epsilon_F=0`, their lower addition line has
zero energy and unit mode weight, with zero upper/removal weight.

The interacting readout is nonlinear in the supplied coefficient. For
`r=G-lambda`, its lower addition energy is
`epsilon_F+(Delta-sqrt(Delta²-4r))/2`; it is not the singular value of `Y`
read from a free fixed-Higgs action. This relation is derived for the stated
finite quantum scalar and prepared state. A relativistic charged-particle
interpretation requires spatial propagation, the corresponding spin and
scattering response, and a physical clock and probe dictionary on the same
construction. None is inferred from matching a finite pole to a particle
name.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -W error -m model.particle_observer_outputs.charged_higgs_records
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -W error -m pytest -q model/particle_observer_outputs/test_charged_higgs_records.py
```
