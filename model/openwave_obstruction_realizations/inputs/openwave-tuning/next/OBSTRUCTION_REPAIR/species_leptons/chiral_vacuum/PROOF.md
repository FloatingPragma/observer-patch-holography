# The full compact-loop vacuum and its boundary response

## Source and charge accounting

This is a one-loop, one-physical-spin control of the frozen wall
Hamiltonian. There are S=2 or4 complete two-component Dirac slices,
2S complex CAR modes, and one full `L2(U1)` loop coordinate theta.
Every occupation sector N=0 through2S is retained. The physical
Hamiltonian is

```
H=E^2/4+dGamma[A+C exp(i theta)+C* exp(−i theta)],
E=−i d/dtheta, h=1/8,
A=(I_S tensor beta+H_fifth+H_boundary)/h,
C=−I_S tensor(beta+i alpha_z)/(2h).
```

Here beta=Z, alpha_z=X and the exact fifth and boundary terms are those
of `chiral_service`. The two bare masses are zero and the actual432
electron coefficient. Their source vacua are computed separately. There
is no three-family gauge vacuum or three-dimensional momentum space in
this control.

On a one-vertex loop the incoming and outgoing electric divergences
cancel. A static reference charge −N cancels the matter number N.
Its energy is zero by declaration. The direct sum of all reference
sectors forms one Gauss-invariant source. A charged fermion insertion
also changes the reference charge, making the complete measured probe
neutral. Reference energy, locality and physical external-charge
interpretation cannot be inferred from this finite construction.

The unitary one-particle matrix `U=reflection_s tensor Y` obeys
`U A^T U*=−A` and `U C^T U*=−C`. Since A and C are traceless, its
canonical particle-hole transformation maps N to2S−N and leaves the
Hamiltonian and loop flux unchanged. The executable checks these
identities and retains every paired number sector. Particle and hole
boundary residues are computed in their actual distinct CAR bases.

## Vacuum and complete boundary propagator

For each S, mass and flux cutoff K=4,8,16, every number-sector ground
energy is compared. The lowest sector is N=S in the executed controls.
The complete eigensystems also retain the empty, full, bulk and
boundary occupations. No positive-band or three-particle reduction is
used. In the half-filled ground state Omega, the boundary field is
`b_a=sum_j E_ja* c_j`, with the frozen isometry
`E=|0>P+ +|S−1>P−`. Reference translations accompany b and b*.

For every eigenstate in the N+1 and N−1 sectors, respectively, define

```
v_n,a=<n,N+1|b_a*|Omega>,
w_n,a=<n,N−1|b_a|Omega>.
```

The output retains all excitation energies and both2x2 residue matrices,
including zero weights. Their sum is exactly the boundary CAR matrix I.
At zero frequency the retarded matrix is

```
G_ab(0)=−sum_n v_n,a* v_n,b/(E_n,+−E0)
         +sum_n w_n,b* w_n,a/(E_n,−−E0).
```

The coefficient of beta in `−G(0)^−1` is a specified finite boundary
chirality-mixing diagnostic. It is distinct from the lowest charged
transition energy and from a Lorentzian continuum pole mass. Its field
normalization is fixed by E, and its complete spectral sum is checked.

The bare-mass-zero S2, K16 control gives
`<cos theta>=−0.9360412511988625`, a charged reference transition near
12.27572654517, and a beta coefficient near31.92965874060.
The nonzero input electron coefficient changes the transition by only
about3.67e−5. Thus this finite dynamical vacuum does not preserve the
light pole of its flat-background boundary service.

For S4, m=0, K16, the corresponding values are
`<cos theta>=−0.9531267170862867`, transition9.944208317540458,
and beta coefficient127.14751931671904. The infinite-flux charged
transition interval is `[9.94410796568,9.94430854118]` under the stated
arithmetic model. Since cosine is bounded below by−1, the measured
negative mean implies at least95.3 percent probability in the
`cos(theta)<=0` sector, where W>=1. The scalar localization hypothesis
therefore fails in a directly measured part of the actual vacuum.

The source has no magnetic pinning potential or normal-ordering
counterterm. Its filled-band energy favours the cutoff side of the sole
spatial holonomy. At theta=pi the scalar Wilson shift is W=2, outside
the |W|<1 localization regime of the scalar wall bound. At zero
frequency the static boundary Schur mass is exactly `W^S/h`, giving
32 forS2 and128 forS4. These numbers explain the direction of the
dynamical boundary response without fitting it. A ring with several
spatial cells retains several momenta, so its full filled-band energy
and Gauss-law electric terms require a separate calculation. This
one-cell failure is not a theorem against domain-wall fermions.

## Infinite-flux ground enclosures

Write the finite matrix as H_K and let A_N,C_N be the literal CAR lifts
in a number sector. Put n=min(N,2S−N). The zero-mass A_N has lower edge

```
a_N=−16 sum_(j=1..n) cos[j pi/(2S+1)].
```

The boundary input changes that edge by at most |m|. The flux-changing
operator has norm
`c_N=8 sqrt(floor((n+1)^2/4))`: each singly occupied slice is a two-state
raising operator, and the maximum total-spin block gives this norm.
The discarded flux subspace therefore has the rigorous lower bound

```
d_N=(K+1)^2/4+a_N−|m|−2 c_N.
```

If E_K is the finite ground energy, the two-block quadratic estimate is

```
E_infinite >= (E_K+d_N−sqrt((E_K−d_N)^2+4 c_N^2))/2.
```

This coarse bound is retained even when it does not separate sectors.
For a trial lower endpoint x<d_N, the full Schur complement is bounded
below by `H_K−x−delta P_edge`, with
`delta=c_N^2/(d_N−x)` and P_edge the two retained boundary-flux blocks.
The full finite eigensystem supplies a low subspace P with J states.
Let p bound `||P_edge P||^2`. Its orthogonal complement is bounded by
the first excluded finite eigenvalue. A second two-block estimate gives

```
lambda_min(H_K−delta P_edge)
 >= lower_root(E0−delta p, E_J−delta, delta sqrt(p)).
```

When this is at least x, x is a full infinite-flux ground lower bound.
The declared numerical enclosure tolerance is1e−4 source-energy units;
it selects no physical parameter. The actual boundary probabilities and
all failed coarse/enclosure checks are retained. At S2, m=0, K16,
combining the enclosed neutral and charged sectors gives the interval
`[12.27562654255,12.27582654639]` after outward decimal rounding.
It is separated from zero by the bound, rather than merely by cutoff
agreement. Reference-sector selection also uses these lower and upper
ground bounds.

## Arithmetic scope

The finite complete eigensystem is checked by its full matrix residual
and Gram matrix. A conservative standard-IEEE envelope uses at most16d
real operations per complex dot of length d and a1e−10 assembly margin.
If `||Q*Q−I||<=eta<1` and `||H Q−Q D||<=r`, the diagonal eigenvalues
have error at most `r/sqrt(1−eta)`. Polar orthogonalization gives the
operator comparison used for the low-subspace edge probability, with
all stated roundoff terms retained. The receipt records this arithmetic
model; it is not a claim that the linked BLAS library performs interval
arithmetic.

The independent rational review imports no source module. Machin's
arctangent formula bounds pi by alternating rational series; alternating
cosine series and integer-square-root bounds then check the outward
values of a_N,c_N,d_N. It checks every number sector and complete
particle/hole record count, and independently sums the CAR weights.
The analytic Schur inequalities apply to the infinite rotor, while
matrix residuals and the displayed decimals have the stated numerical
arithmetic scope.

## Physical boundary

The source-controlled scalar-shift suppression theorem survives as a
conditional theorem. Its localization hypothesis is not satisfied by
the tested compact-loop vacuum. The failure includes actual gauge
electric dynamics, all CAR occupations, both particle and hole probes,
the static charge reference and the inherited clock. It does not solve
the muon ratio, establish an observed electron mass, or transfer the
Euclidean Ginsparg-Wilson Ward identity to this Hamiltonian. A larger
spatial source and its common charged-family vacuum are separate
operator and state constructions.
