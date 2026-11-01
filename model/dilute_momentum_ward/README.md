# Collision-suppressed local momentum Ward cells

The same interacting periodic Weyl/contact action now has a controlled
**local** momentum-conservation defect on a specified initial-state domain.
The bound uses the small contact weight of a dilute momentum band; it
improves the previous full interaction-norm estimate. Exact face currents
and spatial coarse-graining retain the contact force. The result is an
approximate Ward theorem with an explicit error, not an exactly conserved
local stress or a completed Einstein premise.

No projection, filter or effective Hamiltonian is inserted during evolution.
The band projector below describes the allowed initial states, including
their arbitrary untouched reference correlations. Preparing these states
and measuring the currents remain separate source tasks.

## Same local columns, now on the periodic source

Retain the action and graded boundary convention of
[periodic_momentum](../periodic_momentum/README.md),
`H=H0+U D`, `D=sum_x n_x,up n_x,down`, with U>=0. Active periods are at
least three; a side of length one is inactive. Let t_j be the one-particle
positive shift and `q_j=(t_j-t_j^*)/(2i)`. Its symbol in the retained
`exp(-ikx)` convention is sin(k_j). Thus q_j is the **negative** of the
central difference in the earlier open-box momentum attachment.

With Q_x the one-particle site projector, define

```
s_j(x) = dGamma({Q_x,q_j}/2),
A_xy   = Q_x h0 Q_y - Q_y h0 Q_x,
Pi_xy,j = (i/2) dGamma({A_xy,q_j}),
F_j(x) = i[U D,s_j(x)].
```

The braces are one-particle anticommutators before the CAR lift.
Cartesian shifts commute with h0. Therefore the unchanged local CAR proof
gives, including wrapping bonds,

```
i[H,s_j(x)] + sum_y Pi_xy,j = F_j(x).
```

The nearest-face stress is antisymmetric in x,y. For an oriented positive
axis-j edge (x,y), the contribution to **each** endpoint's force is

```
-(U/4) sum_s (n_x,bar(s)-n_y,bar(s))
                  (c_xs^* c_ys+c_ys^* c_xs).
```

The same sign at both endpoints matters: this term is not silently made
an internal conservative flux. The exact certificate on `3x3x3`, c=1/8,
U=1 checks all 81 local equations, 81 nonzero kinetic faces per column,
and all 27 rectangular block partitions per column on the complete
54-mode Clifford algebra. This covers all `2^54` Fock states symbolically;
it does not execute their amplitudes. The evaluator admits the unchanged
nonnegative rational source coefficients and at most 64 sites.

## A band has small contact weight, with all scattering outputs retained

Let B be any set of spatial Fourier modes on a finite torus of volume V.
Retain both spin components at every allowed momentum. Write
`f=|B|/V` and let Omega be the Fock projector with every occupied momentum
in B and `N<=Nmax`. Define `d=floor(Nmax/2)`.

For the band annihilator b_xs, translation invariance of the Fourier
projector gives `||b_xs||=sqrt(f)`. The projected contact operator is

```
Omega D Omega = sum_x b_x,up^* b_x,down^* b_x,down b_x,up
                    <= f N_up,  and <= f N_down.
```

This operator preserves N_up,N_down. Decomposing into those integer
sectors gives `Omega D Omega<=f d`. Also `0<=D<=d` on N<=Nmax, hence
`D^2<=d D`. Consequently

```
||U D Omega|| <= U d sqrt(f).
```

This is the interaction on **all outgoing states**, not just
`||Omega U D Omega||<=U f d`. Confusing those two norms would lose the
leading outside-band scattering. H0 preserves Omega because its spin
matrix acts at fixed momentum. Duhamel comparison with the free evolution
therefore gives the all-time finite bound

```
||(I-Omega) exp(-itH) Omega|| <= min(1, |t| U d sqrt(f)).
```

No invariance under the interacting H is assumed. The bound is uniform in
V at fixed f,Nmax,U,t; how f changes with a family of boxes or bands must
be specified.

In N=2 there is a sharper exact diagnostic. At total momentum K, the
normalized contact vector d_K has squared band projection

```
lambda_K = |B intersection (K-B)|/V.
```

The contact has rank one in each fiber, so
`||D Omega||^2=max_K lambda_K` and
`||(I-Omega)D Omega||^2=max_K lambda_K(1-lambda_K)`.
For B=-B, the first maximum is f at K=0. If f<=1/2, the second is f(1-f).
These are exact principal-angle identities, not cutoff perturbation
theory. On the full 28-dimensional N=2 sector of the existing four-site
source, the spin-complete band `{0,1}` gives lambda values
`1/4,1/2,1/4,0` in its four total-momentum sectors; all sectors and all
outside-band outputs are retained.

## Local and coarse integrated Ward bounds

Each full edge force is bounded by U: the two opposite-spin number
differences commute with their hopping terms, each of norm at most one.
The half allocation over the two incident axis-j edges gives
`||F_j(x)||<=U`. Normal ordering the full edge expression produces eight
quartic terms. Each has two annihilators; on Omega their norm is at most
`sqrt(f)sqrt(f)=f`, while the two creation factors have norm at most one.
The same allocation therefore gives

```
||F_j(x) Omega|| <= min(U,4U f).
```

For Nmax<=1 the force is exactly zero. In the remaining cases, decompose
the evolved state with Omega and use the leakage bound. If preparation
is at time zero and `0<=t0<=t1`, define the **actual** local Ward residual

```
R_x[t0,t1] = s_j(x,t1)-s_j(x,t0)
                  + sum_y integral_(t0)^(t1) Pi_xy,j(t) dt
            = integral_(t0)^(t1) F_j(x,t) dt.
```

Unitary invariance of norm and integration yield

```
||R_x[t0,t1] Omega|| <= min{
  U(t1-t0),
  4U f(t1-t0) + (U^2 d sqrt(f)/2)(t1^2-t0^2)
}.
```

The initial-band clock cannot be restarted at t0: the state there need
not lie in Omega. `budget(..., start_time=t0)` keeps this preparation
origin. Merging adjacent cells telescopes the density endpoints and
adds the force integrals; the displayed un-clipped polynomial bounds
also telescope correctly. No point-sampled current is substituted for
the time integral.

For a spatial block C, sum the actual local equations. Antisymmetry
cancels every internal face, and the norm defect is bounded by the sum
of its constituent site budgets (at most |C| times the uniform bound).
Successive spatial pushforwards compose. `Lean/QFT/DiluteMomentumWard.lean`
checks the finite-sum error transport, the coarse forced Ward bound and
merged error bound. The state-domain and CAR estimates are supplied by
the proof above and the exact Python certificate, not encoded as new
Lean axioms.

These operator bounds imply the same absolute expectation error for
every normalized state in Omega and arbitrary untouched references.
They are uniform for a fixed-size block, not for an extensive block
whose volume grows without an error budget. Source splitting, current
readout and numerical/hardware errors would have to be added separately.

## A sharp total-force formula and a better local momentum approximation

In N=2 at total K the nearest-neighbor total momentum S_j has symbol

```
sin(k_j)+sin(K_j-k_j) = 2 sin(K_j/2) cos(k_j-K_j/2).
```

The uniform contact vector has mean zero and variance
`2 sin(K_j/2)^2` for this multiplier on every active period at least
three. A rank-one projector commutator has norm equal to the square
root of that variance. Hence, exactly on the complete K fiber,

```
||i[U D,S_j]|| = sqrt(2) U |sin(K_j/2)|.
```

The nonzero commutator has rank two. Exact full N=2 matrices check its
cubic polynomial and squared trace in all four L=4 momentum sectors;
the norm squares are `0,1,2,1` in U-squared units. At K=0 this **total**
force vanishes while the local forces need not vanish.

Now assume initially `|k_j|<=kappa` for each particle and `2kappa<pi`.
The conserved principal total P_j=K_j then stays in `|K_j|<=2kappa`,
even for scattered relative momenta outside B. On that entire total-K
subspace, `||P_j-S_j||<=4kappa`; initially the sine remainder improves
this to `kappa^3/3`. Consequently

```
||(P_j-S_j) exp(-itH) Omega||
    <= kappa^3/3 + 4kappa min(1, |t| U sqrt(f)).
```

The spectral restriction here is two-particle total momentum, preserved
by the actual H, rather than an incorrectly preserved individual cutoff.
For f=O(kappa^3), fixed U,t and N=2, even division of this error by kappa
tends to zero. This is a controlled local-sine approximation in a dilute
long-wavelength domain, not a derivation of physical momentum units.

## Quantified examples and the critical-scaling boundary

The finite counting example uses a `64^3` torus and B=`{-1,0,1}^3`:
f=27/262144, 125 possible contact total-momentum sectors, Nmax=2,
t=1/8, U=1 and `kappa=pi/32<11/112`. Exact outward rational bounds give

* outside-band leakage below 0.00127, compared with the old 0.125 bound;
* one-site Ward defect below 0.000132;
* eight-site block defect below 0.00105;
* rescaled local momentum error below 0.0083.

This is an exact mode-count and analytic resource example, not execution
of a `2^(2*64^3)` state. The separate eight-site N=2 diagnostic executes
the original finite Hamiltonian, its local densities, kinetic stresses
and force integrals. At t=1/2, a zero-momentum pair has local unforced
defect about 0.0437282, while the complete forced identity agrees to
below 1e-12. The exact initial collision and escape norm squares are
1/8 and 7/64. Floating time evolution illustrates the all-time proof;
it does not establish the bounds.

The source coupling is visible in every formula. In the distinct
critically tuned band-edge family, f proportional to a^3 and U
proportional to a^-2 make `U sqrt(f)` proportional to a^-1/2. The present
leakage estimate then gives no vanishing bound. This failure of the
estimate is not a disproof of the moving bound-pair limit or of a
restricted renormalized theory. The exact force norm likewise does not
automatically bound a selected bound-state expectation from below.

Remaining obligations include finite-density or critically scaled local
interacting stress, a symmetric Lorentz tensor, physical time and units,
geometric modular data, entropy/area stationarity and an Einstein-register
inhabitant. None follows from taking the present dilute errors to zero.

```sh
.venv/bin/python -m model.dilute_momentum_ward --check model/dilute_momentum_ward/results.json
.venv/bin/python -m pytest -q -W error model/dilute_momentum_ward
cd Lean && lake env lean QFT/DiluteMomentumWard.lean
```
