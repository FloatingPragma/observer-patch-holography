# Native gauge characteristics on the Dirac clock

## 1. Source and inherited trace normalization

Use the compact gauge-coordinate parent, M15 representation and full
Dirac/Fock matter of [FIELD_ROUTE](../particle_observer_outputs/FIELD_ROUTE.md).
At spatial spacing `a`, fix two positive dimensionless architecture
parameters `e0,kappa0`, with one Hamiltonian clock:

```
H = e0/(2a) sum_links C_trace
    + dGamma(h)
    - kappa0/(15a) sum_plaquettes Re Tr_M15 R(U_p).       (1)
```

The electric form is inherited from
[`native_casimirs`](../particle_field_preparation/preparation.py), whose
metric is `Tr_M15(t_r t_s)` for Hermitian generators. It is positive and
adjoint invariant. The full Fock carrier has sixty modes per spatial
cell. The ordinary Wilson choice is
`h=(beta+H_W)/a`; the normalized positive-smoothing action in
[the smooth-field proof](../particle_smooth_sign_field/PROOF.md) has
the same leading free Dirac symbol. These are explicitly distinct
Hamiltonian choices, with their own finite-cutoff bounds.

For clarity, choose a canonical internal basis. The eight color matrices
are the Gell-Mann matrices divided by two in the first three fundamental
coordinates. The three weak matrices are the Pauli matrices divided by
two in the last two coordinates. The remaining matrix is

```
Y5=diag(-1/3,-1/3,-1/3,1/2,1/2).
```

Apply the retained exterior derivative on `Lambda²(C5)+Lambda⁴(C5)`,
in the frozen ordered source basis. The twelve resulting Hermitian
matrices are `t_r`. Their full Gram matrix is

```
I_rs = Tr_M15(t_r t_s)
     = diag(2,2,2,2,2,2,2,2,2,2,2,10/3).              (2)
```

The charge convention is `q6=6Y`, `Q=T3+Y`, with the internal multiplets
written as all-left representations. The trace of `q6²` is120.
The executable reconstitutes every one of the twelve retained native
port currents from these matrices and checks rank12. It verifies the
inherited Casimir, including its charge normalization:

```
C_trace = sum_(r,s) (I^-1)_rs E_r E_s
        = C3/2 + C2/2 + E_q6²/120.                    (3)
```

Here the last expression denotes the electric Casimir on any admitted
link representation; the M15 matrices test the same quadratic form.
The representation indices `(10/3,2,2)` are the ones used in
[`coefficient_selection`](../observer_particles/coefficient_selection.py)
and `Lean/Screen/KineticFamilyCancellation.lean`. The latter proves
conditional algebra on that ray. It does not select a physical running
coupling or a family count. The carrier invariant-form theorem
`Lean/Screen/GaugeKineticInvariantForms.lean` also keeps its extra
simultaneous-invariance hypothesis explicit. Neither theorem is used
to replace the electric metric of (1).

## 2. Magnetic Hessian from ordered native links

For an oriented plaquette write its four factors as
`exp(i epsilon X1)...exp(i epsilon X4)`. Each `Xj` is a Hermitian M15
matrix; negative orientation is included in `Xj`. Set `F=sum_j Xj`.
The ordered product has coefficients

```
P0=I,  P1=iF,
P2=-(1/2)sum_j Xj²-sum_(j<k) Xj Xk.
```

Hermiticity gives `ReTr(Xj Xk)=Tr(Xj Xk)=Tr(Xk Xj)`. Therefore

```
1-ReTr(P)/15 = epsilon² Tr(F²)/30 + R3.             (4)
```

No link commutativity is assumed. With
`x=|epsilon| sum_j ||Xj||`, operator-norm product estimates give

```
|R3| <= sum_(n>=3) x^n/n! .                        (5)
```

The factor15 cancels the dimension in `|Tr A|<=15||A||`. The code
uses the inherited rational exponential-tail enclosure. This bounds
the multiplication potential on its stated weak-field coordinates;
it is not a ground-state localization estimate or a many-body
detector-probability error.

With link coordinates `theta_i^r(x)` and
`U_i(x)=exp(i sum_r theta_i^r(x)t_r)`, the first-order plaquette field is

```
F_ij^r(x) = Delta_i theta_j^r(x)-Delta_j theta_i^r(x).
```

Equation(4) gives the quadratic magnetic Hamiltonian

```
H_B^(2) = kappa0/(30a) sum_(x,i<j,r,s) I_rs F_ij^r F_ij^s.   (6)
```

The canonical cotangent symbol of the source Casimir at the identity
gives `H_E^(2)=e0/(2a) sum I^-1_rs pi_i^r pi_i^s`. This is the
quadratic tangent Hamiltonian of (1). The parent compact-group quantum
Hamiltonian has curved coordinate and interaction terms beyond that
tangent expansion. A global canonical position-momentum commutator
on a finite rotor is not asserted.

The expansion also retains a concrete nonabelian interaction. Take a
constant plaquette `exp(i epsilon X)exp(i epsilon Y)exp(-i epsilon X)
exp(-i epsilon Y)`. Its linear curvature vanishes, while direct ordered
expansion through degree4 gives

```
1-ReTr(P)/15 = -epsilon⁴ Tr([X,Y]²)/30 + O(epsilon⁵).
```

For the native weak generators `X=t_weak1,Y=t_weak2`, the coefficient
is `1/15`. Hypercharge commutes with both, giving zero for that
commutator loop. The Hamiltonian therefore contains a nonzero
nonabelian magnetic interaction; no asymptotic scattering process
or measured cross section follows from this local coefficient.

## 3. Exact lattice curl, Gauss law and common characteristic

On a periodic cubic lattice use Fourier phases `z_i=exp(i a p_i)` and
forward differences `d_i=z_i-1`. In the cyclic plaquette ordering
`(23,31,12)`, the curvature is `C(d) theta`, where `C(d)v=d cross v`.
Direct multiplication gives

```
K(d) = C(d)^dagger C(d)
     = lambda(d) I3-d d^dagger,
lambda(d) = sum_i |d_i|² = 4 sum_i sin²(a p_i/2).     (7)
```

For `d!=0`, the eigenvalues are `0,lambda,lambda`. The gauge gradient
`theta=d phi` has zero curvature. The linearized Gauss equation is
`d^dagger pi=0`; Hamilton evolution preserves it because
`d^dagger K=0`. Quotienting the longitudinal gauge coordinate and imposing
Gauss leaves two transverse oscillator responses per internal
generator. At `d=0`, all three spatial coordinates have zero quadratic
magnetic cost. Flat torus holonomies and the associated electric
zero-mode evolution are retained. They are not counted as two
nonzero-momentum polarizations.
Constant noncommuting link variations can acquire the quartic
magnetic cost above, despite their zero quadratic cost.

Hamilton's equations, with the **full** trace matrices, are

```
dot(theta) = e0 I^-1 pi/a,
dot(pi) = -kappa0 I K theta/(15a),
ddot(theta) = -e0 kappa0 K theta/(15a²).             (8)
```

Thus every gauge direction has the same quadratic dispersion

```
omega_g²(p) = c_g² [4 sum_i sin²(a p_i/2)]/a²,
c_g² = e0 kappa0/15.                               (9)
```

The inverse Gram matrix in the electric term is essential. Equal
coefficients on the canonical `Y,SU2,SU3` electric generators define
a different source. For example, that alternative has a hypercharge
to weak squared-speed ratio `5/3`. It is not the trace-metric source
in(1). Conversely, multiplying an individual trace-normalized
direction without changing the conjugate coordinate is not a basis
change. For any real invertible basis map `t'=t S`,
`I'=S^T I S`; electric momenta and link coordinates transform dually.
All cross terms must be retained and `I'^-1 I'=I` preserves(8).

In particular `q6=6Y` sends `theta_q6=theta_Y/6` and
`pi_q6=6 pi_Y`. Its trace index is120, its electric coefficient is
`e0/120`, and its magnetic index is120. The product in(8) is unchanged.
The one-direction q6 rotor with a separately chosen `E²/2` coefficient
is a distinct control; its coefficient cannot be imported as the
full field's trace metric.

The principal Dirac speed on the declared clock is one. Therefore

```
common gauge and Dirac characteristic  <=>  e0 kappa0=15.    (10)
```

This is a permitted architecture condition, without a measured input
or a uniqueness claim. For example `kappa0=1/7,e0=105` satisfies it.
The distinct member `kappa0=1/7,e0=1` has `c_g²=1/105`; for
`kappa0=1/100,e0=1` it is `1/1500`. A common change of the time unit
rescales both matter and gauge speeds and cannot remove these ratios.
Selecting the matched family changes the declared electric action.
Electric preparation-energy bounds must use its actual multiplier;
the unscaled member's state and records do not transfer by relabeling.

This condition leaves one positive action parameter. At the quadratic
level the common coefficient is proportional to the M15 trace form;
the coupling ratio obtained after canonical normalization is the
same conditional trace ray discussed above. No renormalization scale,
running coupling, particle census or comparison to measured coupling
constants is derived here.

## 4. Spatial approximation and the Dirac comparison

For the native Clifford matrices set `alpha_i=i beta gamma_i`. The
ordinary Wilson symbol, on all four spinor rows, is

```
h_W(p) = [sum_i alpha_i sin(a p_i)
          + beta sum_i(1-cos(a p_i))]/a,
E_W² = [sum_i sin²(a p_i)+(sum_i(1-cos(a p_i)))²]/a². (11)
```

The exact square follows from the native anticommutators. Both
chiralities are retained. The normalized smooth action has the same
unit leading coefficient by the calculation in its owning proof;
equality of its interacting records with Wilson records is not used.

For `|p|<=P`, elementary sine bounds give uniform estimates

```
0 <= |p|²-lambda(a p)/a² <= a² P⁴/12,
|E_W²-|p|²| <= 7 a² P⁴/12,
||h_W-alpha.p|| <= a P²/2+a² P³/6.                  (12)
```

The first follows from
`0<=x²-4sin²(x/2)<=x⁴/12`. For the second, the sine contribution is
bounded by `a² sum_i p_i⁴/3` and the Wilson square by `a²|p|⁴/4`.
The operator estimate separates the Clifford sine error from the
Wilson scalar. Hence the matched quadratic gauge and free Dirac
symbols approach the same principal cone on bounded momentum sets.
Their finite-spacing dispersions need not be equal.

For nonzero lattice frequency the gauge group velocity satisfies

```
|grad_p omega_g|²/c_g²
 = sum_i sin²(a p_i)/[4 sum_i sin²(a p_i/2)] <=1.
```

This is a quadratic dispersion bound. It does not assert an exact
sharp causal front for a continuous-time lattice Hamiltonian or
replace the full interacting Lieb-Robinson estimates.

## 5. Constant coframe compatibility

Use exactly the coframe convention of
[the geometric-response proof](../particle_geometric_response/PROOF.md):
`G=V^T V`, `v=det V>0`, `A=G/v`. Both electric and magnetic spatial
quadratic forms use `A`, while the Dirac frame is `V^-T`.
For a real continuum covector `k`, `C(k)` is the real cross-product
matrix. The matched gauge second-order principal matrix is

```
-A C(k) A C(k) = q I-k k^T G^-1,
q=k^T G^-1 k.                                    (13)
```

To verify it, use the cross-product identity
`-A C A C=det(A)[(k^T A^-1 k)I-k k^T A^-1]`.
Here `det(A)=1/v` and `A^-1=v G^-1`.
Equation(13) has the gauge kernel `span(k)` and repeated physical
eigenvalue `q`. The Dirac square is independently

```
[alpha.(V^-T k)]²=q I4.                           (14)
```

Thus the inherited coframe construction also preserves the shared
principal metric under(10). Equation(13) concerns the continuum
principal coefficients at a constant coframe; the exact complex
forward lattice curl in(7) is the flat-lattice statement. This
calculation supplies neither spacetime frame dynamics nor an
interacting Lorentz transformation law for the prepared quantum state.

## 6. Finite evidence and physical scope

The tests reconstruct all twelve native directions, check their full
trace form and inherited electric Casimir, expand ordered noncommuting
plaquettes, and compare the remainder against direct matrix
exponentials. They independently assemble the real-space curl on a
periodic `3x3x3` cube: its rank is52, with29 zero directions and
transverse eigenvalue multiplicities `3:12,6:24,9:16`.
The29-dimensional kernel comprises26 nonconstant gauge gradients
and three constant coordinates. Tests also check Hamilton-equation
normalization, charge covariance, nonorthogonal internal basis
cross terms, finite-spacing bounds and rational sheared coframes.

The calculation derives a quadratic compatibility condition from the
same native gauge action. It does not prove that a quantum ground
state is localized in this weak-field chart, that the interacting
spectral measures have relativistic poles, or that those poles are
photons, gluons, W/Z bosons or massive fermions. In particular the
twelve massless tangent directions do not identify twelve observed
massless particles. A state-sensitive treatment of interactions,
symmetry breaking, confinement and gauge-invariant readout is
required for that identification. The common-clock condition is one
necessary check for a proposed relativistic field construction.
