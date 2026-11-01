# Original-action retarded matter response

This package computes a state-dependent interacting response of the unchanged
three-dimensional cube action. It supplies the full eight-site addition/removal
resolvent, an exact nonzero memory kernel, and controlled finite-frequency
expansions. These are matter-response inputs to the existing geometry theorem
stack. A common effective metric or an Einstein equation is not inferred from a
static density shift.

## Source and reference

The source is the open `2 x 2 x 2` cube with two CAR modes per site,
`c=1/8`, `U=1`, and all 65,536 Fock states. There is no chemical potential or
new mass term. The existing onsite unitary frame turns its Weyl hopping into
the real spin-independent signed adjacency `h=A/16`, where `A^2=3 I`.
All calculations below use this hidden spin frame; the physical response is
its sitewise unitary conjugate. Hidden spin is not assigned a physical Lorentz
spin or electric charge.

The retained `interacting_vacuum` theorem certifies the full global ground
energy `E_G=-sqrt(3)/4`, its five-dimensional spin-two ground space, and a
strict gap `delta=7/1024` to every different particle-number sector.
We choose its explicit normalized highest-weight Slater vector `G`: the four
negative-band hidden-up orbitals are occupied and the down orbitals are empty.
This is a source-defined reference, not an executed native preparation.
The second stationary reference is

```
A_state = sqrt(2) c†_(0,up) G,     E_A = -3 sqrt(3)/16.
```

Both references are exact eigenvectors of the original action, with zero
contact energy. The second reference is excited; it receives no ground-gap
claim. `response.py` restricts the actual complete source matrix to invariant
hidden-number sectors, retaining every state in those sectors. The calculation
does not replace source evolution by the observed eight-dimensional space.

## Exact retarded dilation and memory

For a normalized stationary reference `psi` of energy `E`, define, for
`Im z > 0`,

```
G_xy(z) = -i integral_0^infinity e^(izt) <{c_x(t),c_y†}> dt
        = <c_x [z-(H-E)]^-1 c_y†>
          + <c_y† [z+(H-E)]^-1 c_x>.
```

The addition vectors are `c_y† psi`; the removal vectors live in the
**conjugate** removal Hilbert space and are `conj(c_y psi)`. Thus

```
L = (H_(N+1)-E) directsum -conj(H_(N-1)-E),
V_y = c_y† psi directsum conj(c_y psi),
V†V = I,                         G(z)=V†(z-L)^-1 V.
```

The conjugate removal space gives the correct matrix indices for a general
complex reference. Here all arrays are real in the inherited frame. CAR gives
the isometry exactly, including addition and removal, with no survival
normalization. On `P=VV†`, `Q=I-P`, let

```
M1 = V† L V,       W = Q L V = L V - V M1,
B  = Q L Q restricted to Ran Q.
```

Block elimination of `z-L` gives the exact identity

```
G(z) = [z I - M1 - Sigma(z)]^-1,
Sigma(z) = W†(z-B)^-1 W,          W†W = M2-M1^2.
```

Both inverses exist for `Im z>0`: `L` and `B` are self-adjoint, and
`-Im G(z)=Im(z) V†[(z-L)^-1]†(z-L)^-1 V` is positive definite.
Similarly `-Im Sigma(z)>=0`. Eliminating the unobserved states leaves their
exact frequency-dependent memory; dropping `Sigma` would change the response.

For both references, the **up channel** obeys `LV=Vh` exactly, so
`G_up(z)=(z-h)^-1` and `W=0`. This is a complete cyclic-space identity, stronger
than agreement of a few moments. The ground up dilation has 56 addition and
56 removal states; the charged up dilation has 28 and 70 respectively.

For the **ground down channel**, removal vanishes and the addition isometry
has eight columns in the full 560-state `N_up=4,N_down=1` block. Exact source
arithmetic gives

```
M1 = h + I/2,                    W†W = I/4.
```

In particular, every nonzero observed vector couples to eliminated states;
the coupling `W` has norm `1/2`. An independent operator derivation is
`W_x=(n_(x,up)-1/2)c†_(x,down)G`. Distinct down locations are orthogonal,
and each up occupation has mean `1/2` and variance `1/4`.

The charged down channel retains all 448 states with `N_up=5,N_down=1`.
In site order
`000,001,010,011,100,101,110,111`, its exact results are

```
M1 = h + diag(1, 2/3, 2/3, 1/2, 2/3, 1/2, 1/2, 1/2),
W†W = diag(0, 2/9, 2/9, 1/4, 2/9, 1/4, 1/4, 1/4).
```

The change in `M1` is exactly the added orbital density
`diag(1/2,1/6,1/6,0,1/6,0,0,0)`. It is a computed response to a change in
matter state at fixed microscopic graph and law. Removing the contact term
in the explicit `U=0` comparator changes the ground down response to
`(z-h)^-1` and removes the memory entirely. It does not make `G` the full
free-Fock ground.

## Complete local charged spectrum

The original physical origin probe `v=c†_(0,down)G` has unit norm.
`spectrum.py` resolves **all 78 distinct positive lines** in this probe's
spectral measure; the previously unresolved weight `4/5` is now accounted
for. This resolves its cyclic response, not the entire 560-state block or
its degeneracies. There is no new executed preparation or detector history.

The retained exact monic polynomial `p` has degree 78 over `Q(sqrt(3))`.
Direct unbounded-integer arithmetic verifies `p(16H)v=0`. A nonzero
78-by-78 Hankel determinant modulo 97, with `sqrt(3)` mapped to 10,
proves the matching characteristic-zero cyclic-rank lower bound.
Hence the local cyclic dimension is exactly 78. Modular reconstruction
only proposes `p`; the exact annihilation and rank checks admit it.
The exact first 78 moments determine its resolvent numerator `q`, with

```
<v,(z-16H)^-1 v> = q(z)/(1296 p(z)).
```

Here `q` uses the unnormalized integer state numerator of squared norm
1296. A floating eigensolver only proposes brackets. Exact 1024-bit outward
dyadic interval arithmetic verifies 78 disjoint strict sign changes of
`p`, exhausting its degree, then refines every bracket to 384 bits and
encloses every normalized residue `q(x)/(1296 p'(x))` strictly above zero.
The gap is `(x+4 sqrt(3))/16`. The receipt publishes all gap and weight
intervals rounded outward to 24 decimal places, with total weight enclosed
within `1 +/- 4*10^-23`.

The smallest gap is approximately `0.029950938342053226`, with weight
`0.3130113772155167`; the largest is approximately `1.8360744654423855`,
with weight `2.711810625844545e-10`. These displayed decimals are summaries
of the retained rational intervals. Exact polynomial substitution proves
that the lines `sqrt(3)/16` and `1+sqrt(3)/16` have total weight exactly
`1/10` each, excluding coincident residual contributions. The other 76
lines sum to `4/5`; unresolved local probe weight is zero. These positive
finite charged gaps are neither infinite-volume poles nor relativistic
rest masses.

## Controlled finite-frequency domain

For the interacting ground down channel, `K=H_(4,1)-E_G` satisfies
`delta I < K < 2 I`. The lower bound uses the retained full global-ground
certificate, not the ground-vector eigenidentity alone. For the upper bound,
`D<=1` on this block, the kinetic norm is at most `5 lambda`, and
`-E_G=4 lambda`, where `lambda=sqrt(3)/16`; hence
`K<=1+9 lambda<2`. Compression gives the same bounds on `B`.

Consequently the low-frequency series converges in operator norm:

```
Sigma(z) = -sum_(n>=0) z^n C_n,       C_n=W† B^(-n-1) W >= 0.
```

For every `|z|<=r<delta`, truncation after order `p` has error

```
||R_p|| <= r^(p+1) / [4 delta^(p+1) (delta-r)].
```

This follows by a norm-convergent geometric series for `B^-1`; it also gives
`0<=C_(n+1)<=C_n/delta`. The same Taylor expansion of `G` itself has four
times this remainder bound. These are actual finite gapped frequency
expansions. The coefficients are exact functions of the supplied source
matrices; this does not assert they are spatially local.

The package also computes all matrix moments `M0,...,M4` exactly over
`Q(sqrt(3))`. For `|z|>=R>2`,

```
G(z) = sum_(n=0)^p M_n/z^(n+1) + remainder,
||remainder|| <= 2^(p+1) / [R^(p+1)(R-2)].
```

`bounds.py` evaluates these bounds using rational arithmetic and rejects
the closed gap boundary, the large-frequency boundary and inexact inputs.
The finite gap is not transferred to the excited reference or infinite volume.

The static coefficient is also computed, with an exact error certificate.
Write `C=V†K^-1V=-G_down(0)`. In operator norm,

```
||C - (12.9058381883 I - 4.9026412320 A)|| < 10^-7,
||Sigma(0) - (-0.3633623064 I - 0.0105943755 A)|| < 10^-6.
```

`static.py` uses a floating-point solve only to propose an integer matrix
`X/2^40`. It encloses every entry of the exact residual `KX-V` over
`Q(sqrt(3))`, using an integer square-root bracket for `sqrt(3)`, and bounds
its Frobenius norm. Since `K>delta`, the inverse-vector error is at most
`||KX-V||/delta`. Compression by the exact isometry then certifies the
displayed rational approximation to `C`. Because `K<2`, `C>I/2`, so
`||C^-1||<2`; a second exact residual, using `A^2=3I`, certifies a fixed
rational inverse and `Sigma(0)=C^-1-M1`. The identities are not admitted by
the numerical solve. A corrupted proposal fails. We do not infer exact
equality to an affine function of `A` from these enclosures.

## Interface to the existing geometry stack

`thermodynamic_ground` supplies a translation- and gauge-invariant ground
subsequence for the same original action on growing boxes, and analytic
local retarded responses uniformly on compact subsets of `Im z>0`.
It does not compute the limiting state or provide a real-axis limit. The
present cube is an exact finite response calculation; it is not declared a
convergent approximation to a particular infinite ground without an
additional state-convergence estimate.

The relevant existing interface is
`model/relative_observer_geometry/intrinsic_geometry.py:wave_density_metric`.
It reconstructs the metric from a complete local four-dimensional action
density `C=mu diag(-w,a)` with positive `mu,w` and positive definite spatial
`a`, through `g=sqrt(-det C) C^-1`. The existing
`model/observer_particles/geometric_matter_join.py:source_tetrad` then
constructs a coframe and proves common scalar, Weyl, Dirac and Maxwell
principal-symbol relations for its supplied comparison density.

The new output is the actual source response from which a candidate effective
principal symbol must be extracted, including its memory and state dependence.
The next join requires controlled spatial and low-energy limits of the
**inverse** response, its normalization and a common coefficient across the
relevant matter channels. A Hartree diagonal has no spatial derivatives and
does not determine that principal symbol or its conformal normalization.
The current finite gapped disk is not a Lorentz cone. No array from the cube
is silently substituted for the distinct screen-source coefficients in
`source_tetrad`.

Fixed microscopic geometry is compatible with such emergent response
geometry. A new fundamental metric register is not required by this route.
After the common response/action-density join, the actual modular, area,
stress, coarse-graining and null-balance hypotheses of `EinsteinRegisterData`
still require proof on the same family. None of G1--G6 is marked complete here.

## Replay

Run `python -W error -m pytest model/retarded_matter_response` and
`python -W error -m model.retarded_matter_response --check` from the repository
environment. The receipt binds the earlier global-ground proof and receipt;
its costly all-number Gram proposal is not rerun as part of this package.
New moment, isometry and memory identities are replayed exactly. Independent
floating-point resolvent tests are diagnostics, not admission criteria.
