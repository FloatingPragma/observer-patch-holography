# Three kinetic momentum balances on the same interacting source action

The unchanged open-box Weyl/contact action now has three exact momentum
balance columns with explicit nearest-face kinetic stresses and contact
forces. These complement the conserved energy column in
[coarse_energy_ward](../coarse_energy_ward/README.md). They close the finite
source attachment for these declared observables, while identifying exactly
why they do **not** fill the three conserved momentum columns required by
the Einstein interface.

The exact certificate uses the complete 54-mode CAR algebra of the
`3x3x3` box at `c=1/8,U=1`. It checks 81 local forced equations, 54 nonzero
kinetic face currents per column and all 27 rectangular block partitions
per column. All `2^54` Fock states, their coherences and untouched references
are covered by the operator identities; this is symbolic coverage, not an
execution of that many state amplitudes. No Gibbs reference, physical clock
or spacetime interpretation is assumed.

## Actual kinetic momentum and face stress

Let `Q_x` be the one-particle projector onto the two components at x. On an
axis-i forward nearest edge define `(D_i)_(x s,y t)=-i delta_st/2`, with
its Hermitian reverse. The source convention is

```
h0 = -c sum_i sigma_i D_i,       H0=dGamma(h0),
V=U sum_x n_x0 n_x1,            H=H0+V,
P_j=dGamma(D_j),                p_j(x)=dGamma({Q_x,D_j}/2).
```

The braces denote the **one-particle** anticommutator, before second
quantization. `dGamma(q)=sum_ab q_ab c_a^dagger c_b` includes no vacuum
subtraction. The D_i act on separate Cartesian factors, so `[D_i,D_j]=0`
and `[h0,D_j]=0`, including on the stated finite open box. Thus P_j is a
conserved free kinetic observable. Its infinite/periodic Fourier symbol
is `sin k_j`; calling it kinetic momentum does not identify it with the
logarithm of a lattice translation, physical momentum, or the generator of
translations of the finite open boundary.

For an oriented nearest edge x to y, set

```
A_xy = Q_x h0 Q_y - Q_y h0 Q_x,
Pi_xy,j = (i/2) dGamma({A_xy,D_j}).
```

This is Hermitian and antisymmetric under exchanging x,y. Although the
operator can act on sites one further step away, its *flux label* is the
nearest face. No range-two current is omitted: those terms are included in
the anticommutator. Since
`sum_y A_xy=[Q_x,h0]`, direct expansion using `[h0,D_j]=0` gives

```
i[H0,p_j(x)] + sum_y Pi_xy,j = 0.
```

The only boundary convention is the existing open action: absent spatial
edges have zero flux. Independent complete-Fock Jordan--Wigner checks on
all three axes and on a plaquette bind these matrices to the same source
Hamiltonian and its Clifford representation. Exact arithmetic uses only
Gaussian rationals with unbounded integer numerators. The symbolic evaluator
admits at most 64 sites and binary source coefficients with denominators at
most `2^16`; the algebraic identities apply to every finite rectangular box
and finite real c,U for the stated self-adjoint action.

## Interaction force and an exact all-direction obstruction

The interaction contribution is not renamed a divergence. It is
`f_j(x)=i[V,p_j(x)]`. On a forward axis-j edge define

```
F_xy = (U/2) sum_s (n_x,bar(s)-n_y,bar(s))
                    (c_xs^dagger c_ys+c_ys^dagger c_xs).
```

The opposite-spin number factors commute with each displayed hopping
operator. With the half-bond p_j convention, f_j(x) is one half the sum of
F_xy over axis-j edges incident at x. **The two endpoints have the same
force sign**, whereas the kinetic face stress has opposite signs. The
CAR commutator proves the exact local balance

```
i[H,p_j(x)] + sum_y Pi_xy,j = f_j(x).
```

The force vanishes on the vacuum and one-particle sector and when U=0. It
is quartic in Majorana operators on the complete algebra. For distinct
edges, its Majorana monomials are orthogonal in normalized Hilbert--Schmidt
inner product. Each full edge force has squared norm `U^2/8` in that inner
product. Consequently, if E_i counts axis-i spatial edges,

```
Tr[(sum_x f_i(x))(sum_y f_j(y))]/dim(Fock)
    = delta_ij U^2 E_i/8.
```

At the retained 27-site instance this Gram matrix is exactly `(9/4) I_3`.
Every total force has operator norm at least 3/2, and every nonzero real
linear combination has nonzero total force. Hence none of these three
interacting kinetic momenta, or their nonzero linear combinations, can
satisfy a closed-boundary conservative balance. Adding any antisymmetric
internal face flux cannot remove a nonzero summed force. The calculation
does not assert that a physical external force acts on the simulator: it
identifies the finite lattice/contact term for this supplied momentum
observable. Restriction to a smooth sector or another limit would require
its own quantitative force estimate.

## Complete number-preserving quadratic conservation class

The obstruction is stronger than testing a guessed current. Let
`Q=dGamma(q)` be any number-preserving quadratic observable and assume
`[H,Q]=0` on the complete Fock algebra at one fixed `U>0`.

1. Restrict to the one-particle sector. V vanishes there, so `[h0,q]=0`.
   The CAR lift identity then gives `[H0,Q]=0` on *every* sector.
   Therefore `[V,Q]=0` separately; no cancellation with the kinetic term
   has been assumed away.
2. For x different from y and arbitrary s,t, use the two-particle input
   `c_x,bar(s)^dagger c_yt^dagger|0>` and output doublon at x. Their V
   energies are zero and U. The selected matrix element of Q is exactly
   `+q_(xs,yt)` or `-q_(xs,yt)`: only that single hop reaches the output.
   Thus `[V,Q]=0` forces every intersite coefficient to vanish.
   Conversely any site-local bilinear commutes with the local doublon
   projector, so this condition is also sufficient for contact commutation.
3. For `c>0`, a connected nearest axis-i bond imposes
   `q_x sigma_i = sigma_i q_y`, or `q_y=sigma_i q_x sigma_i`.
   Pauli conjugations along different axes commute, even though the
   Pauli matrices themselves anticommute. Therefore the entire onsite
   family is the transported root matrix
   `q_x=R_x q_0 R_x^dagger`, with
   `R_x=sigma_1^x1 sigma_2^x2 sigma_3^x3`.

Every Hermitian root 2x2 matrix is allowed. The conserved homogeneous
number-preserving quadratic space has **real dimension four**: total number
and the three existing hidden-spin charges. An additive scalar identity is
allowed separately. The certificate isolates all 2,808 offsite coefficients
on the 27-site box by exact occupation-bit matrix elements and verifies the
four transported generators. An independent small-cube linear constraint
rank check gives dimension four. This proves the absence of a conserved
quadratic observable retaining an intersite kinetic momentum component,
at the fixed positive interaction; no assumption of commutation for every
U is needed. The theorem does not exclude higher-body dressings, anomalous
pairing observables, nonlocal nonquadratic translation charges or a different
physical action/reference. It also does not call hidden spin Lorentz spin.

## Finite time cells and coarse-graining retain the force

Heisenberg evolution uses the full interacting H. Integrating the exact
operator balance over a cell gives

```
p_j(x,t1)-p_j(x,t0) + sum_y integral_t0^t1 Pi_xy,j(t) dt
    = integral_t0^t1 f_j(x,t) dt.
```

Any normalized state/reference expectation obeys the same equation. With
lower-face indexing, it supplies `ddiv(T)_j=integrated_force_j`, rather
than the zero right side in the Einstein premise. These are finite cell
quantities; no pointwise symmetric stress tensor or metric volume factor
is inferred. The previous energy column remains conserved, with the same
H and arbitrary input state. The four columns have not thereby acquired
Lorentz symmetry or the symmetry `T_ij=T_ji`.

Spatial blocks sum the momentum densities and forces and cancel their
internal kinetic faces. The divergence commutes with pushforward under
any finite block map, and successive block maps compose as previously
proved. Merging time cells telescopes endpoints and adds both stress and
force integrals. `Lean/QFT/LatticeMomentumBalance.lean` proves the generic
forced block equation, cancellation of total antisymmetric divergence,
the total momentum change and exclusion of absorbing nonzero total force
into a closed internal divergence, without admissions or new axioms. The
concrete source and quadratic classification are separately justified by
the exact CAR certificate and the sector argument above.

An independent complete 64-state three-site spectral computation checks the
finite-cell equation at duration 3/4. It evaluates both operator integrals
through a sinc formula that retains degenerate frequencies. Deleting the
force or replacing its integral by an endpoint sample fails by a nonzero
operator residual. These floating controls are not the all-time proof.

## Source approximation and finite sampling

The same locality and complete-field source-splitting bounds used by the
energy column apply here. Conservatively, `||p_j(x)||<=1`,
`||Pi_xy,j||<=2c`, and `||f_j(x)||<=U`. The supports fit respectively
radius/size pairs `(1,3)`, `(2,6)` and `(1,3)`. For the stress norm,
`||A_xy||_1=2c` and `||D_j||<=1`, so the one-particle anticommutator has
trace norm at most 2c; `||dGamma(q)||<=||q||_1`. The force norm follows
from the explicit two-spin edge expression and half-bond allocation.

For any observable of norm bound B and support radius r with at most s
sites, let L be the existing normalized locality tail and epsilon the
patch unitary error, uniform for `0<=u<=t`. Its endpoint error is bounded
by `B min(2,L+2 epsilon)`. Its integral error is at most t times that
bound. Since `||dO/dt||<=2 B s(U+6c)`, replacing the integral by m midpoint
samples adds `t^2 B s(U+6c)/(2m)`. Thus a momentum-cell residual is bounded
by two momentum-endpoint errors, six kinetic-face integral errors and one
force-integral error. Supplied instantaneous readout errors add twice the
endpoint error and t times each of the six face errors and force error;
none is certified without a readout model.

The receipt evaluates the exact rational budget at `R=18,t=1/8,n=10^7`
source split steps and `m=10^12` midpoint samples on an analytic `37^3`
patch, distinct from the 27-site symbolic instance. This is a resource
bound, not source execution or measurement of those counts. It is uniform
in the larger ambient volume, not in the size of an arbitrarily extensive
coarse block. Local input restriction, initialization, measurement accuracy,
floating/hardware errors and controller/export-failure energy remain
separate. The inherited tail domain requires `24ct<=32` and `2<R<=4096`.

```sh
python -m model.lattice_momentum_balance --check model/lattice_momentum_balance/results.json
python -m pytest -q -W error model/lattice_momentum_balance
cd Lean && lake env lean QFT/LatticeMomentumBalance.lean
```
