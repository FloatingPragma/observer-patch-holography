# The full matrix spectrum of the original grounded cube

The entire eight-site down-addition response of the original cube has the
same 78 energy lines as the previously resolved origin probe. Every matrix
residue has rank four. The observed cyclic subspace therefore has dimension
312 inside the full 560-state charged sector. Its response, inverse response
and exact memory kernel are all affine functions of the signed cube
adjacency. This gives an exact finite spatial structure, explicit static
coefficients and a stronger controlled frequency disk.

The action and ground are those of
[the retained source response](../retarded_matter_response/README.md):
`c=1/8,U=1`, original zero chemical potential, ground
`E_G=-sqrt(3)/4`, and all eight columns `V_x=c_x,down^dagger G`.
The source restriction `H16=16H_(N_up=4,N_down=1)` keeps all 560 states.
No projected approximate dynamics or fitted mass is introduced.

## Exact matrix resolvent from a finite source recurrence

The previous source polynomial `q` is monic of degree 78 over `Q(sqrt(3))`.
The new calculation proves `q(H16)V=0` for **all eight columns**, using
unbounded integer arithmetic. It also evaluates every complete matrix
moment through order 77 and verifies

```
1296 V^dagger H16^j V = d_j I + a_j A,  0<=j<78,
A=A^T,  A^2=3I,  tr(A)=0.
```

The factor 1296 is the denominator of the exact source vector Gram matrix;
`V` itself is an isometry. The polynomial recurrence fixes all subsequent
moments. Equivalently, polynomial division of `(x-H16)^-1 V`, with its
remainder killed by `q(H16)V=0`, gives the exact rational identity

```
G(z) = 16 [p_I(x) I+p_A(x) A] / [1296 q(x)],
x = 16z-4sqrt(3).
```

The coefficients of `p_I,p_A` are finite convolutions of `q` and the exact
moments. This proves the identity wherever the resolvent exists, not just
as a large-frequency approximation. The independent complete 560-state
numerical resolvent is a diagnostic cross-check, not an admission rule.

## Full matrix residues and the observed cyclic space

The calculation further verifies exactly

```
q divides p_I^2-3p_A^2.
```

At each of the 78 distinct real roots `x_i` already certified for the
origin probe, let

```
alpha_i=p_I(x_i)/(1296q'(x_i))>0,
beta_i =p_A(x_i)/(1296q'(x_i)).
```

The polynomial identity gives `3 beta_i^2=alpha_i^2`. Strict outward
integer interval arithmetic determines the sign `s_i` of every beta.
Consequently the full residue is

```
R_i=alpha_i I+beta_i A=2alpha_i P_(s_i),
P_+/-=(I +/- A/sqrt(3))/2.
```

Both projectors have rank four: `A` is symmetric, `A^2=3I`, and its trace
is zero. Every `R_i` is therefore positive with rank four. There are 40
lines in the plus band and 38 in the minus band. Each band's weights
`2alpha_i` sum to one, because `sum R_i=V^dagger V=I`.

For a spectral projection `E_i` of the full source, the rank of `E_i V`
equals the rank of its Gram matrix `V^dagger E_i V=R_i`. Distinct energy
spaces are orthogonal, so

```
C=span{H16^j V u:j>=0,u in C^8},  dim(C)=sum_i rank(R_i)=312.
```

This is a reducing subspace of the self-adjoint source. The other 248
dimensions are orthogonal to all these probes and their time translates.
Their eigenvalues and degeneracies have not been resolved. No claim of
completing the full 560-state spectrum follows.

Every site has the same local measure because `A` has zero diagonal.
The 78 energies remain excitation energies of this finite source, rather
than identified physical species or relativistic rest masses.

## Exact inverse response and a larger controlled disk

Write `G=F I+J A`. For `Im z>0`, both band functions
`g_s(z)=sum_(i:s_i=s) 2alpha_i/(z-E_i)` have strictly negative imaginary
part. Thus neither vanishes, and

```
G^-1=(F I-J A)/(F^2-3J^2).
Sigma=zI-(A/16+I/2)-G^-1.
```

The full finite inverse response and memory therefore contain only
onsite and the original signed nearest-neighbor matrix entries. This
special cube identity does not assert spatial locality on larger boxes
or identify a Lorentzian principal symbol.
There are only two adjacency eigenvalues on this cube. Closure of that
two-dimensional matrix algebra does not supply a spatial derivative limit,
propagating three-dimensional dispersion or a common limiting light speed.

The certified common support lies strictly between the rational constants

```
g=299/10000=0.0299,    b=37/20=1.85.
```

These bounds apply to `K=(H-E_G)|C`. Let `Q_C` be the complement of the
observed columns inside `C`. The compressed operator `B_C=Q_C K Q_C` has
the same lower and upper bounds. The exact coupling `W=KV-VM1` lies in
`C`, so the full Schur memory equals `W^dagger(z-B_C)^-1 W`; the
orthogonal unobserved space contributes zero. This is a reducing-space
identity, not deletion of leakage.

Since `W^dagger W=I/4`, the memory Taylor series on `|z|<=r<g`, truncated
after order `p`, obeys

```
||R_p|| <= r^(p+1)/[4g^(p+1)(g-r)].
```

The response itself has four times this error bound. The controlled disk
is larger than the earlier conservative `7/1024` full-sector gap disk.
The stronger bound is not applied to the remaining 248 source dimensions,
an excited reference or an infinite-volume ground.

## Exact static values

At `z=0`, evaluate the rational functions at `x=-4sqrt(3)` to obtain
`C=-G(0)=c0 I+c1 A` exactly over `Q(sqrt(3))`. The retained receipt gives
all four rational coefficients of `c0,c1`. It also gives the exact inverse
and static memory using

```
C^-1=(c0 I-c1 A)/(c0^2-3c1^2),
Sigma(0)=C^-1-(A/16+I/2).
```

These exact values satisfy the prior certified approximations
`c0≈12.9058381883`, `c1≈-4.9026412320`. A floating linear solve is no
longer required to calculate these static coefficients. The physical
meaning remains an exact finite source response, not an inferred metric.

## A complete finite-gap particle-dictionary test

The full support also permits a stronger identification test than assigning
a few attractive-looking lines. Consider the explicit proposal that any
visible positive single-charge ground-probe gap is a charged-lepton rest
energy, with one common positive unit factor and no additive energy shift.
All 78 down-channel gaps and the positive up-channel line lie in the
certified interval used here. Its upper-to-lower ratio is below 62.

The already retained [PDG 2025 comparison baseline](../particle_species_detector/physical_targets.json)
has a muon/electron ratio above 206 and tau/electron ratio above 3477,
even throughout the declared boxes consisting of each reported central
value plus or minus its reported uncertainty. Those boxes are benchmark
inputs, not exact guarantees about physical confidence coverage. No choice
of two visible lines and positive unit factor can reproduce either pair.

The mismatch has a quantitative bound. If the source ratio is at most `R`,
the target ratio at least `Q>R`, and both relative errors are at most
`epsilon<1`, then

```
(1-epsilon)/(1+epsilon) <= R/Q,
epsilon >= (Q-R)/(Q+R).
```

Thus at least one member of a proposed electron/muon pair has relative
error exceeding 54 percent; the corresponding electron/tau bound exceeds
96 percent. This is an exact rational lower bound over all assignments
and common positive scales, not a fitted calibration. The target data
never enter the source spectrum.

The tested finite-gap-as-rest-energy dictionary is rejected. This does
not prove that those gaps are relativistic rest masses, or exclude other
probes, unobserved sectors, composite particles or a continuum refinement.
No whole-model falsification or eligible prospective physical test is
inferred. The remaining identification must derive the appropriate
relativistic sector and its actual detector dictionary.

Run `python -W error -m pytest -q model/matrix_spectral_response` and
`python -W error -m model.matrix_spectral_response --check` from the
repository environment. All earlier source proofs and their boundaries
remain separately retained.
