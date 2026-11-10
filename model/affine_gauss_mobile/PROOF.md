# Mobile compensated charge on the finite affine-Gauss carrier

## One changed Hamiltonian and its actual ground

The immutable static parent is
`model/affine_gauss_source/results.json`. The public calculation validates
its complete source pins and freshly replays that parent. A preserved
independent research receipt supplies a separate exact scientific comparison,
as stated in [PROVENANCE.md](PROVENANCE.md).
Its carrier, thirteen CAR modes, fixed compensator at vertex 1, twelve
longitudinal bits, total-photon cutoff three, counting metrics and exact
cutoff boundary term are retained. Its 30 transverse modes have finite
carrier dimension 5456. The physical source dimension is44695552; its
one-charge sector dimension is 70928. All parent Gauss constraints and
reference energy `R` are retained. This construction is distinct from the
native 119-state link Hamiltonian.

Let d_y be the compensated *annihilator*, the adjoint of the parent's
creation D_y. For each inherited oriented cone edge e=(u,v), declare the
finite Peierls matrix `W_e=exp(i A_e)` and Hamiltonian

```
H_eta=H_static+V,
V=-eta sum_(u,v) [d_v^dag d_u W_e + W_e^dag d_u^dag d_v],
eta=1/24.
```

Every A_e is the actual finite transverse coordinate of the parent.
W_e is its literal finite matrix exponential. It is unitary; no canonical
Weyl commutation or Baker-Campbell-Hausdorff identity is used. Cutoff
coordinates can have noncanonical boundary commutators. No compact
independent-link representation is inferred from the Peierls notation.
The number-conserving d_v^dag d_u changes both source positions and their
longitudinal bits. It leaves the compensator at vertex 1 with the same R.
On full kinematic inputs the individual lowering and raising operations
and their boundary zeros are retained.

All terms commute with every parent's Gauss generator and total matter
number N. On the physical space the embedding identifies d_y with the
ordinary CAR annihilator and R=N. Since W_e acts only on photons,
Cauchy-Schwarz gives for every state

```
|<d_v psi, d_u W_e psi>| <= (<n_v>+<n_u>)/2.
```

The graph degrees are 12 at the apex and 6 at every boundary vertex.
Therefore, as a physical-space quadratic-form inequality,

```
V >= -eta sum_y degree(y)n_y >= -12 eta N,
H_eta-E_vac >= H_ph+E_C+(2-12 eta)N
             = H_ph+E_C+(3/2)N.
```

The parent's blank neutral vacuum is annihilated by V. It is the unique
actual ground of H_eta on the full physical carrier. Since
`omega_min=sqrt(4-sqrt(5))<3/2` and neutral one-photon states retain their
exact energies, the actual physical gap is exactly omega_min. This is a
finite-volume gap in declared units, not an identified photon mass.
The declared blank preparation is deterministic and prepares this actual
ground. It is not a transfer of the native interacting-vacuum preparation.

The clock for every mobile experiment is `exp(-it H_eta)`. The static
parent's charged energy eigenvectors generally cease to be eigenvectors;
its exact charged energy-line and joint charged-field probabilities are
not imported. Neutral field operators preserve N=0, where H_eta equals
the parent. Thus every parent's neutral field spectral residue and complete
neutral field instrument does transfer exactly to this same interacting
finite member.

## Charge, reference and current accounting

The source electric operator is unchanged:
`E=D G L_charge-E_T`, with the parent's sign convention for E_T. On the
physical space `L_charge=n-N delta_r`. A hop u to v changes its
longitudinal flux by `D G(delta_v-delta_u)`. Consequently its electric
divergence changes by exactly `delta_v-delta_u`; the compensator remains
at r. The transverse Peierls matrix acts in the divergence-free photon
carrier. No boundary reference is moved or renamed by a charge translation.

For `T_e=d_v^dag d_u W_e`, define `J_e=i eta(T_e-T_e^dag)`. Exact finite
CAR commutators give

```
dot n=D^T J,     dot N=dot R=0,
D^T dot E=dot n
```

on the physical source. The second identity follows also by differentiating
the exact Gauss constraint with this H_eta. This statement does not require
canonical commutators for the truncated transverse coordinates. A local
untruncated Ampere law with a canonical current requires additional cutoff
control and is not asserted.

The compensated charge pump is still D_y/2 on the actual neutral ground.
It has success probability 1/4 and failure probability 3/4, with every
branch physical. Its successful state is
`psi_y=|one charge at y, matched reference and longitudinal bits> tensor |0_ph>`.
Its reference contributes one unit to all charged transition energies.

## Actual charged spectral window

Let `E_y=2+E_C(y)`, measured above the actual neutral vacuum. The static
part acts on psi_y with this energy. The hopping part maps it to mutually
orthogonal neighboring position sectors, each with a unit-norm Peierls
photon output. Hence, exactly,

```
<psi_y,(H_eta-E_vac)psi_y>=E_y,
||(H_eta-E_vac-E_y)psi_y||^2=eta^2 degree(y).
```

For the actual spectral projector on `[E_y-w,E_y+w]`, the spectral
Markov bound gives outside probability at most `eta^2 degree(y)/w^2`.
Choose w=1/4. The complete pump and energy record has probabilities

```
pump failure=3/4,
energy window=(1/4)<psi_y,P_window(H_eta-E_vac)psi_y>,
energy complement=(1/4)<psi_y,(I-P_window)psi_y>.
```

They sum to one and use the full interacting spectral projector. The
unconditional window probability is at least 1/6 for an apex charge and
5/24 for a boundary charge. These windows can overlap between source
positions. They are individual settings, not mutually exclusive
particle-species labels or exact static energy lines.

## Exact first hopping integral and all-order remainder

Let x be an actual neighbor of y, P_x its source-position projector,
`U_0(t)=exp(-it(H_static-E_vac))`, and `U(t)=exp(-it(H_eta-E_vac))`.
The exact twice-integrated Duhamel identity gives

```
U(t)psi_y=U_0(t)psi_y
 -i integral_0^t U_0(t-s)V U_0(s)psi_y ds + R(t),
||R(t)|| <= (t^2/2)||V||_(N=1) ||V psi_y||
          <= (t^2/2)(12 eta)(eta sqrt(degree(y)))
          <=24 eta^2 t^2.
```

The last inequality uses `sqrt(degree(y))<=4`. The full and comparison
propagators are unitary, so there is no exponential remainder factor.
This bounds every higher hopping insertion of the actual H_eta; the
source evolution is not replaced by a first-order Hamiltonian.

The zero-order term has no x component. Its first-order component, up
to a global phase, is

```
i eta integral_0^t exp(-i [E_C(x)-E_C(y)+H_ph]r)dr W_(xy)|0_ph>.
```

The orientation selects W_e or its adjoint; all norm bounds apply to both.
All static electric costs are between 0 and 6/31, and `H_ph<9`. Thus the
absolute spectral value delta in the integrand is less than 10. The
integral has magnitude `t |sinc(delta t/2)|`. Since
`sin z/z>=1-z^2/6`, at t=1/4 it is at least `t*71/96` on every photon
spectral component. W_(xy) is unitary, so

```
||P_x U(t)psi_y|| >= eta t*(71/96)-24 eta^2 t^2
                   =47/9216 >0.
```

The actual complete pump-and-arrival instrument has probabilities
`3/4`, `(1/4)||P_x U(t)psi_y||^2`, and
`(1/4)||(I-P_x)U(t)psi_y||^2`. It includes all transverse states at the
arrival site. Its strictly positive unconditional lower arrival bound is

```
2209/339738624 > 6.5020*10^-6.
```

Energy and arrival are separate settings; the bound does not condition on
a prior energy measurement. Source position and the retained reference
have their original meanings. No preparation trace-error subtraction is
needed because the declared preparation is the exact actual vacuum.

## Arrival with nonzero transverse occupation

The coupling changes the field. This can be proved within the finite
carrier rather than by a canonical displacement formula. For an edge e,
write `A_e=a(f_e)+a^dag(f_e)` and

```
sigma_e^2=||f_e||^2=sum_lambda (P_lambda)_(e,e)/(2 sqrt(lambda)).
```

The total-occupation cutoff is invariant under U(30). Rotating f_e into
one direction reduces its vacuum exponential exactly to the four-level
Jacobi matrix `sigma_e(a+a^dag)`. Spectator sectors have smaller cutoffs,
so on the full carrier

```
||A_e||^2=(3+sqrt(6))sigma_e^2 <11/4,
<0,A_e 0>=0,      <0,A_e^2 0>=sigma_e^2.
```

For independent variables a,b in the vacuum spectral measure of A_e,
`(a-b)^2<11`. The elementary bound
`1-cos z>=z^2/2-z^4/24` therefore implies

```
1-|<0,W_e0>|^2
 =integral [1-cos(a-b)]dmu(a)dmu(b)
 >=sigma_e^2/12 >=(P_T)_(e,e)/72 >1/108.
```

The last step uses the exact cone diagonal values 165/217 and 151/217.
This is a finite occupation statement for the actual Peierls matrix.

Let Q denote nonzero total photon occupation. Q commutes with the free
photon Hamiltonian in the first hopping integral. Its norm lower bound
therefore gains a factor `||Q W_e0||>1/11`. At t=1/64 the sinc bound is
larger than 11/12. After the same all-order remainder,

```
||P_x Q U(t)psi_y|| >= eta t/12-24 eta^2 t^2
                     =13/294912.
```

The complete detector keeps pump failure and all four commuting outcomes
of position x versus its complement and photon vacuum versus its
complement. Its unconditional arrival-with-occupation probability is at
least

```
169/347892350976 >4.8578*10^-10.
```

A twelve-state retained pointer, with unused outcomes kept, suffices for
every specified energy, arrival, joint-occupation or neutral-field setting.
It has zero reference Hamiltonian. The parent source-reference-plus-pointer
capacity count is unchanged; other control storage and work are declared
separately.

The other four records are retained with their exact finite positive
projector laws; the probabilities sum to one. Nonzero finite transverse
occupation is not identified with asymptotic emitted radiation.

## Finite spectral evaluation of the first integral

The executable also evaluates the first integral more sharply for the
actual radial edge 0 and boundary edge 12, with source y=1, at times 1/4 and
1/64. It uses the exact eight projector weights and thirty-mode spectrum
from the static receipt. Only one mode direction in each degenerate
frequency space is populated by f_e, so its free history lies in the
165-dimensional occupation space of eight mode directions with total
occupation at most three.

For `alpha=sqrt(3-sqrt(6))` and `beta=sqrt(3+sqrt(6))`, the four-level
matrix exponential is obtained by interpolating its even cosine part in
`1,Z^2` and odd sine part in `Z,Z^3`, where `Z=a+a^dag`. Expanding its
four vacuum amplitudes into the eight mode directions uses exact
multinomial weights. Each occupation contributes its exact squared first
integral `4 sin^2(delta t/2)/delta^2`, with t^2 at delta=0. Outward interval
arithmetic encloses every algebraic and trigonometric expression.

The retained probability lower bounds subtract the complete Duhamel
remainder from the first-integral norm before squaring, including the
charge-pump factor1/4. No truncated exponential of H_eta, fitted oscillator
tail or numerical eigensolver enters this calculation. The receipt keeps
both exact universal rational lower bounds and these sharper finite
interval certificates.

This is a finite interacting charged-source precursor. Its source,
metric, cutoff, reference and Peierls services are declared architecture.
The construction supplies no thermodynamic Coulomb phase, inverse-distance
continuum potential, Lorentz-covariant charged threshold, asymptotic
scattering, physical spin or observed particle dictionary. The original
native source receives no ground or propagation transfer from this member.
