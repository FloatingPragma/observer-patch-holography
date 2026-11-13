# Full compact-group source and a finite-volume vacuum gap

## Carrier and positive Hamiltonian

Take the open cubic graph with vertices `{0,...,L-1}^3`. Each oriented edge
carries `L2(SU(2),dU)` with normalized Haar measure. The physical Hilbert
space is the range of the orthogonal local gauge average over all vertices.
This retains every Peter--Weyl representation, including arbitrarily large
spin. There are no matter factors, electromagnetic photons, or delta-function
link preparations. Let `C_l=-Delta_l` have eigenvalues `j(j+1)`. Set

```
E=alpha sum_l C_l,    V=beta sum_p (1-W_p),
W_p=Re Tr(U_p)/2,     H=E+V,
alpha=g^2/(2a),       beta=2/(g^2 a).
```

Both coefficients are positive. Each `W_p` is bounded between minus one and
one. Thus `V` is bounded and nonnegative, and `H` is the unique self-adjoint
bounded perturbation of the positive compact-manifold Laplacian. Its domain
is the electric domain, and finite Peter--Weyl polynomials form a form core.
Gauge averaging commutes with the electric operator and plaquette
multiplication. The physical range reduces the source.

The compact connected product group has a strictly positive heat kernel at
every positive time. The bounded-potential Feynman--Kac formula retains
strict positivity. Restriction by gauge averaging preserves it. Compact
resolvent and positivity improvement imply a unique, strictly positive
physical vacuum and discrete spectrum at each fixed finite graph. These
facts give no graph-independent gap estimate.

## The electric physical gap

An invariant spin network with any nontrivial edge cannot have a vertex
incident on exactly one nontrivial edge: an isolated nontrivial irreducible
representation has no invariant vector. Its nontrivial support therefore
contains a cycle. The cubic graph has no cycle shorter than four edges.
Each nontrivial SU(2) edge costs at least `3/4`, hence every nonconstant
invariant electric state has energy at least `3alpha`. A fundamental
plaquette character attains this value. This proves the complete electric
gap; it does not replace the interacting vacuum by the electric vacuum.

Let `P` denote the number of elementary plaquettes and let `1` be the
normalized constant wavefunction. Haar orthogonality gives

```
<W_p>=0,    <W_p^2>=1/4,    <W_p W_q>=0  (p != q),
mu=<1,H1>=beta P,    r^2=||(H-mu)1||^2=beta^2 P/4.
```

For distinct elementary faces an edge occurs in one face alone, so its Haar
integral kills the mixed term. The graph producer verifies the complete face
census and this distinct-edge condition. Positivity of `V` and min--max give
`E_1(H)>=3alpha`, while the constant trial gives `E_0(H)<=mu`. Thus

```
gap(H) >= d=3alpha-beta P,
E_0(H) >= max(0,mu-r^2/d)                    when d>0.
```

The second bound is Temple's inequality with the independently established
excited-energy floor. It concerns the full untruncated source.

## Local finite-energy noncollapse

Write the normalized true ground as `c 1+v`, with `v` orthogonal to constants.
The compressed block on that complement is at least `3alpha`. Its ground
equation gives `||v||<=|c| r/d`, so

```
||v||^2 <= r^2/(d^2+r^2).
```

For `L=2,g^2=8,a=1`, the graph has 12 links and six plaquettes. Here
`alpha=4,beta=1/4,mu=3/2,r^2=3/32,d=21/2`, and the last probability is
strictly below `1/1024`. The distance of the ground projector from the
constant projector is consequently below `delta=1/32`. A bounded observable
of norm at most one changes its expectation by at most `2delta`. Applying
this to `W_p,W_p^2` gives

```
Var_ground(W_p) >= 1/4-2delta-4delta^2 = 47/256.
```

The group normalization is fixed by `C chi_(1/2)=3chi_(1/2)/4`.
For one edge, the squared gradient of its normalized fundamental trace is
`(1-W_p^2)/4`; all four edges contribute. Plaquette multiplication commutes
with every magnetic term, so the operator identity on the form core is

```
one half [W_p,[H,W_p]] = alpha(1-W_p^2).
```

The centered vacuum vector therefore has energy form at most `alpha=4`.
Its spectral measure is nonzero and lies above the gap, proving

```
21/2 <= inf Spec(H-E0)|vacuum-perp <= 1024/47.
```

This instantiates a finite-volume version of the local-correlation
noncollapse criterion in `particle_yang_mills_gap_contract`. No equality of
prelimit and continuum spectra is assumed. The local curvature vector need
not generate the entire vacuum complement to give the upper bound; the
lower bound uses the complete min--max argument independently.

## Domain of the bounds

At fixed coupling the displayed lower bound becomes inconclusive as the
plaquette count grows. At fixed graph all energies and both bounds scale as
`1/a`. The local energy bound also diverges as `a` decreases. These are
executed controls, so neither uniform infinite-volume decay nor finite
continuum noncollapse follows from this construction. No renormalized
four-dimensional correlation limit, physical glueball GeV mass, or result
for every compact simple group is supplied. The Yang--Mills existence and
mass-gap problem has those additional requirements; see the primary
[Jaffe--Witten statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
