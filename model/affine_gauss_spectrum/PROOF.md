# The charged ground of the actual finite interacting source

## Source and decomposition

Keep the exact Hamiltonian, total-photon cutoff three, physical Gauss
embedding, charge unit, compensator at vertex1 and clock of
`model/affine_gauss_mobile`. The actual neutral ground is unchanged.
Work in its full one-charge sector of dimension `13*5456=70928`.
Let P select photon vacuum and Q=I−P. The source Hamiltonian above the
neutral vacuum is the exact block matrix

```
H = [ A   R† ],       A=P H P, R=Q H P, C=Q H Q.
    [ R   C  ]
```

P has dimension13. The reference energy remains one on this entire sector;
no charged static eigenline is transferred after adding hopping.

For edge e let `sigma_e²=<0|A_e²|0>`. Total occupation cutoff is invariant
under orthogonal changes of the thirty real normal modes. The vacuum
Peierls vector `exp(i A_e)|0>` therefore has the exact four-level
oscillator amplitudes. Put `alpha=sqrt(3−sqrt6)`,
`beta=sqrt(3+sqrt6)`, and

```
b_e = [cos(sigma_e beta)−cos(sigma_e alpha)]/(2 sqrt6),
k_e = cos(sigma_e alpha)+(sqrt6−2)b_e,
c_2(e) = sqrt2 b_e.
```

Here k_e is the vacuum overlap and c_2(e) the amplitude in the normalized
two-photon state of the e direction. The exact projector polynomials of
the parent show there are only two edge values: radial r and boundary b.
The executable verifies this using every diagonal of every required power
of the actual curl matrix. Outward intervals prove `0<k_r,k_b<1` and
`c_2(r),c_2(b)<0`. Thus

```
A_yy=2+E_C(y),      A_xy=−eta k_e  if x,y share edge e,
eta=1/24.
```

The two units in the diagonal are the source and retained-reference
energies. The graph is connected, and every off-diagonal edge entry is
negative. Perron–Frobenius therefore makes A's ground unique, with a
strictly positive real normalized eigenvector u.

## Bounds for the entire complement

The cone adjacency has characteristic polynomial

```
(x²−5x−12)(x²−5)^3(x+1)^5.
```

The executable checks this identity exactly. Blockwise norm comparison
with its nonnegative adjacency gives the full one-charge hopping bound
`||V||<=v=(5+sqrt73)/48`. Consequently

```
C >= c I,       c=2+sqrt(4−sqrt5)−v >3.
```

This includes every photon occupation allowed by the source. In R, the
edge norm is `eta w_e`, where `w_e=sqrt(1−k_e²)`; it is the norm of the
actual Peierls output after its vacuum component is removed. Block norm
comparison and the radial/boundary equitable partition give

```
||R|| <= r = eta[5w_b+sqrt(25w_b²+48w_r²)]/2.
```

The quotient's positive eigenvector extends to a positive vector on all
thirteen vertices. Perron–Frobenius makes its eigenvalue the spectral
radius of that nonnegative weighted adjacency, so this is a full operator
bound rather than a bound on a selected trial.

Write a_j for the ordered eigenvalues of A. For lambda<c the exact Schur
complement is `A−lambda−R†(C−lambda)^−1 R`, bounded below by
`A−lambda−r²/(c−lambda)`. Inertia and the min–max principle imply

```
f(a_j,c,r) <= E_j(H) <= a_j,      j=1,...,13,
f(a,c,r)=[a+c−sqrt((c−a)²+4r²)]/2,
E_14(H)>=c.
```

The upper bounds use the full13-dimensional variational subspace P.
The lower bound for the fourteenth level uses its codimension13
complement Q. Every a_j is below c, so precisely thirteen eigenvalues
occur below the certified complement threshold. Interval LDL congruence
and inertia count each a_j, including multiplicities, between rational
endpoints. Every pivot must exclude zero; otherwise the program refuses
to certify an inertia. No approximate eigenvector is needed for these
enclosures. Monotonicity of f in a,c and its decrease in r propagate the
outward endpoints to the actual H.

## Strict energy lowering by the actual transverse field

Let `|B>=sum_(j=1)^30 |2_j>/sqrt30`. For every real unit normal-mode
direction f, `<B|2_f>=1/sqrt30`. The one-photon and three-photon pieces
have zero overlap with B. In each arrival position x,

```
<x,B|R u>=−eta sum_(y~x) u_y c_2(xy)/sqrt30.
```

All terms on the right have the same sign. If
`d=min(|c_2(r)|,|c_2(b)|)>0`, positivity of u and minimum graph degree6
give

```
||R u||² >= (eta² d²/30) sum_x(sum_(y~x)u_y)²
          >= eta² d²/5.
```

This is an actual matrix element, not a generic assertion that coupling
usually shifts energies. The full one-charge Hamiltonian has norm below12:
its static energy is at most `2+6/31+3 omega_max`, with `omega_max<3`,
and its hopping norm is v<0.3. In the variational span of u and
`R u/||R u||`, the lower eigenvalue is strictly below a_1. Since a_1>0
and `||R u||<=r<1`, its reduction is at least `||R u||²/13`. Hence

```
E_1(H) <= a_1−eta² d²/65 < a_1.
```

The bound certifies transverse-field dressing in this finite Hamiltonian.
It does not use an infinite-oscillator displacement formula.

The computed interval is approximately
`[1.86354168535,1.87398449270]`, whereas a_1 is approximately
`1.87398476077`. The next actual level is greater than `2.00219970691`.
Thus the charged ground is unique with a strictly positive certified
charged-sector gap. All displayed decimal values are explanatory;
`results.json` retains rational outward endpoints.

## Complete charged-ground preparation and readout

Prepare the parent's actual neutral vacuum. The compensated creator in
the normalized spatial mode u is `D_u=sum_y u_y D_y`. The physical CAR
embedding gives `||D_u||=1`. Its joint-even complete pump `D_u/2` has
success probability1/4 and failure3/4. The success state is the exact
zero-photon compression ground `|u,0>`, whose actual-H mean energy is a_1.

Let L_1 be the certified actual ground lower bound, L_2 its second-level
lower bound and U_A the compression ground upper bound. The spectral
measure of this state gives

```
<u,0|P_ground(H)|u,0> >= (L_2−U_A)/(L_2−L_1) = F_lower >0.
```

An energy window enclosing E_1 but below L_2 isolates that ground in the
prepared N=1 sector. The complete records are

```
pump failure:      3/4,
ground window:     <u,0|P_ground(H)|u,0>/4,
energy complement: [1−<u,0|P_ground(H)|u,0>]/4.
```

They sum to one. Conditional on the selected window, the state is the
unique actual charged ground. Every measurement commutes with the Gauss
constraints; no compensating charge is moved or erased. The trial-mode
preparation and finite spectral detector are explicit finite matrix-functional
calculus services, including the declared Peierls exponentials;
an efficient native circuit is not proved.

For its actual ground vector, the Q block eigen-equation also gives

```
Pr(nonzero transverse occupation)
 <= r²/[(c−U_A)²+r²].
```

This probability is also strictly positive. Write delta for the proved
energy reduction below a_1 and q for the nonvacuum probability. The block
form gives `E_1>=a_1+(c−a_1)q−2r sqrt(q(1−q))`. Together with
`E_1<=a_1−delta` and `c>a_1`, this implies
`q>=delta²/(4r²)>0`. The actual charged ground therefore contains a
nonzero transverse-field contribution.

These statements refer to the same source, state preparation, reference
and evolution clock. Subtracting the known reference energy1 reports the
source contribution to the transition energy. The compensator stays at
one vertex, so the construction still describes a finite boundary problem.
Its thirteen spatial levels are not thirteen particle species. Neither a
translation-invariant mass shell, Lorentz spin, GeV mass nor Standard Model
identity follows from these bounds.
