# Transverse field records, Green forces and moving mass shells

## Source and mathematical domain

The electromagnetic construction uses the local oriented cube stencil of
`model/native_wilson_maxwell/field.py`. Periodic identification of opposite
faces defines the spatial family used here. This boundary condition and the
counting metrics are declared architecture. The periodic family does not
inherit the interacting open-cube ground-state theorem. Its field and charge
normalizations are fixed once for the wave and electrostatic calculations.
The separate moving-matter calculation uses the complete native mass-channel
operator of `model/source_mass_shell`; no equivalence of the two interacting
sources is asserted.

The finite local selection theorem in `model/maxwell_local_selection/PROOF.md`
classifies the quadratic action under strict face locality, gauge symmetry
and committed rotations on its twelve-port carrier. Its energy/reversal
argument applies to any finite incidence complex with the stated positive
energy. That argument, rather than a new independent Ampere postulate, is
the shortcut used for the field evolution below. The action classification
on a twelve-port carrier does not by itself select cubic spatial gluing.

## Finite quantum fields without a full Fock simulation

Let `C` be curl, `D` gradient, `CD=0`, and let

```
W = (epsilon ||E||^2 + beta ||B||^2)/2,  c^2=beta/epsilon.
```

Faraday evolution, conservation of this positive energy for every state,
and the reversal `(E,B)->(E,-B)` give

```
dot E = c^2 C* B,       dot B = -C E.
```

For this periodic field convention, canonical coordinates are identified
as `E=-P` and `B=CQ`; thus `dot Q=P` and `dot P=-c^2 C*CQ` give the
displayed signs. The electric Gauss label is defined by `D*E=rho`. This
choice of electric orientation does not transport the interacting
open-cube current convention or alter its electrostatic quadratic energy.
Consequently `D*E` and magnetic divergence are constant. The quantum
finite-cutoff implementation requires no canonical commutator assumption
at its boundary. On total occupation at most `K`, use the actual truncated
ladder matrices and `H=sum_j omega_j a_j* a_j`. Their energy differences give
`[H,a_j]=-omega_j a_j` on every retained column, including the top boundary.
For `q_j=(a_j+a_j*)/sqrt(2 omega_j)` and
`p_j=-i sqrt(omega_j/2)(a_j-a_j*)`, therefore

```
i[H,q_j]=p_j,        i[H,p_j]=-omega_j^2 q_j.
```

The positive correction
`P_K (H+sum_j omega_j I)/2` added to the truncated quadratic oscillator
energy makes it `H+(sum_j omega_j)I/2`. Omitting this term changes its
boundary dynamics. `waves.cutoff_theorem` checks complete two-mode
cutoffs; the displayed energy-difference proof covers every finite number
of modes and every cutoff. Thus the linear field expectations follow the
cochain evolution exactly, without diagonalizing an exponentially large
field matrix. Each real electric/magnetic component is a finite Hermitian
observable; its spectral projectors define a complete measurement.

## Two polarizations and time-domain propagation

For cubic spacing `a`, the oriented forward difference has symbol
`d_j=(exp(i a k_j)-1)/a`. Curl is the exterior product by `d`; hence

```
C* C = |d|^2 I - d d*,
P_T=I-d d*/|d|^2,             rank P_T=2       (d != 0).
```

For each transverse plane field, its frequency is

```
omega_a^2 = (4 c^2/a^2) sum_j sin^2(a k_j/2).
```

With `B(0)=0`, the exact real field evolution is
`E(t)=cos(omega_a t)E(0)` and
`B(t)=-sin(omega_a t) C E(0)/omega_a`.
The executable independently evolves the real-space stencil by symmetric
kick/drift/kick updates and compares both field records to this exact
propagator. It checks both polarizations, all divergence constraints,
second-order temporal convergence, and four spatial refinements. For
`|a k_j/2|<=1`,

```
0 <= c^2 |k|^2-omega_a^2 <= c^2 a^2 sum_j k_j^4/12.
```

A fixed `c=1` is a source normalization. The theorem proves convergence to
that common speed and two transverse modes; it does not determine a
laboratory time unit or identify every mode as an observed photon.

## Electrostatic force from the same incidence operator

On the same cubic complex let `L=D*D` and `G=L^+`, with its constant mode
retained as gauge redundancy. A neutral source `rho` has the unique
minimum-norm solution `E_L=D G rho`, since

```
D* E_L=rho,     C E_L=0,
||E_L||^2=rho*G rho.
```

Any other solution is `E_L+F` with `D*F=0`, and its energy increases by
`||F||^2/2`, because `<D G rho,F>=0`. The physical source convention
multiplies this normalized electric field by the declared coupling `g`.
The energy is `g^2 rho*G rho/2`; varying `g` changes the force and does
not affect incidence identities.

The executable computes `G` from the finite lattice symbol, rather than
sampling a supplied Coulomb potential. For charges `+1,-1` at separation
`r`, it obtains `g^2[G(0)-G(r)]`. Potential differences, central-difference
radial forces, box refinement and orientation controls test the long-range
form `G(r)~1/(4 pi r)`. The continuum coefficient is therefore conditional
on the specified electric energy and charge normalization, `g^2/(4 pi)`.
No measured fine-structure constant enters the calculation, and no
prediction of that constant is claimed. The point charges are supplied
Gauss sources, not identified dynamical electrons.

## Massive wave records on the complete native one-body source

The existing mass-channel proof retains all eight lattice valleys in a
sixteen-row Bloch matrix. For supplied coefficients `(m0,m1,m2,m3)`, it gives

```
H_a(p)=(v/a) sum_j sin(a p_j) gamma_j + m0 beta0 + sum_j m_j beta_j.
```

The anticommutation identities split the operator by the commuting
involution `Q=beta0 (sum_j m_j beta_j)/|m|`. In each eight-row sector,

```
H_a^2 = [(v/a)^2 sum_j sin^2(a p_j)+(m0 +/- |m|)^2] I.
```

The second-order wave law follows from this source identity. The executable
compares its exact cosine/sine propagator to a direct matrix exponential,
checks all signed branches and multiplicities, and verifies the finite
spacing bound

```
0 <= v^2 |p|^2-(v/a)^2 sum_j sin^2(a p_j)
   <= v^2 a^2 sum_j p_j^4/3.
```

A continuum Lorentz boost preserves the limiting mass shell. Evaluating
the same boosted momentum on the finite-spacing operator exposes a nonzero
error that decreases quadratically with spacing. This is a quantified
free-sector limit; it does not assert exact finite-spacing Lorentz symmetry,
source-derived mass coefficients, a moving interacting defect, or a
laboratory particle identity. A hidden scalar energy shift breaks the
squared law and is rejected by the complete matrix comparison.

## Bound-state consumer of the computed Green field

A Dirichlet boundary replaces the periodic boundary for the atom-like
standing-wave calculation; the interior Gauss and electric-energy stencil
is identical. In the finite cube, the point-source potential is obtained
by inverting its discrete Dirichlet Poisson operator. No `1/r` array is
inserted. A supplied nonrelativistic mass and kinetic law then define

```
H = L/(2 m a^2) - (g^2/a) G_D(.,0),
m=1, g^2=4 pi.
```

The coupling convention makes the continuum Coulomb coefficient one.
It is a declared model normalization, not a measurement fit or a derived
fine-structure constant. The cube length is 32 model units. The executed
41, 63 and 95 point grids have spacings `32/42`, `1/2` and `1/3`.

Reflection in each coordinate commutes with both the kinetic matrix and
the centered Green potential. The normalized sums/differences of each
reflected site pair, with the center retained in the even sector, give
exact reducing isometries. Their tensor products reduce the 857375-site
finest grid to 110592 even-sector rows. Cubic coordinate permutations
identify the three sectors odd in one axis and even in the others. This
symmetry reduction preserves the complete source operator on each sector;
it supplies no spherical-symmetry assumption.

The numerical ground, next even-sector and lowest odd-sector energies on
the finest grid are approximately `-0.4911620284`, `-0.0766068319` and
`-0.0707661659`. Direct eigenvector residuals are below `10^-8`. These
negative finite-box levels support localized bound-state witnesses for the
specified attractive potential. Reversing its sign gives a positive
Hamiltonian and no negative level. The source is a prescribed static
charge, not a self-consistent dynamical nucleus.

The ground-to-odd spacing decreases from `0.68599` to `0.49575` to
`0.42040`; the continuum hydrogen value for these supplied units is
`3/8`. The splitting between the next even and lowest odd sectors decreases
from `0.04116` to `0.01562` to `0.00584`. The finite cubic source therefore
has bound levels and a useful refinement diagnostic. These computations
are not a precision hydrogen spectrum, and neither the supplied kinetic
law nor the measured atomic mass/coupling is derived by the field theorem.
