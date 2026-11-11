# Complete scalar-site response of a localized band

## One spatial source and one unit interaction

Use the literal compact-pair Hamiltonian in
[`compact_defect_pair.spectral.source`](../compact_defect_pair/spectral.py),
with radius three, texture amplitude 64 and half-separation one. Its complete
686-row matrix has a four-dimensional band in `[-r,r]`, `r=7/240`, and every
complementary eigenvalue has absolute value at least 62. The
[parent proof](../compact_defect_pair/PROOF.md) supplies the complete-band,
localization and exterior-volume estimates. Its separate exact certificate
excludes zero for this source. The overlap calculation diagonalizes all 686
rows and keeps the entire band, including both signs and any degeneracy.

Let `W` be an orthonormal frame for that spectral subspace, `Q` any complete
orthonormal complement, and `C=D_stagger i sigma2`. The parent identities give

```
hW=WA, C conjugate(W)=WS,
S=W† C conjugate(W), S^T=-S, S†S=I4.
```

For site projector `P_x` define every local scalar intervention by
`Delta_x=P_x C`. These are antisymmetric creation-pair kernels. The same
spatial factor can multiply a unit symmetric native internal Higgs tensor
`D`, giving `Delta_x tensor conjugate(D)`. Tensoring retains the native
current intertwiner of the parent source. This calculation uses the spatial
factor with unit primitive coefficient; no family matrix is supplied. It
does not transfer the local gauge-link or physical Gauss-state results of
a different Higgs Hamiltonian.

The complete scalar-site kernel and density matrix are

```
B_x=W† P_x W, Y_x=W† P_x C conjugate(W)=B_x S.
```

Every one of the 343 sites is included in `site_kernel.npz`, together with
`W`, `S`, `B_x` and `Y_x`. This is a basis for independent scalar interventions,
not a selected scalar vacuum. The matrix `Y_x` is antisymmetric and has rank
at most two. A scalar profile `f` gives `Y[f]=sum_x f_x Y_x`. Its singular
values occur in equal pairs. Four spatial modes therefore supply at most
two different response strengths in this construction.

For any unitary change `W -> WU`, the kernel transforms as
`Y_x -> U† Y_x conjugate(U)`. Singular values, ranks, complete leakage and
response spectra are frame independent. Selecting an individual eigenvector
inside a degenerate band is unnecessary. The full normalized band metric is
identity, so no nonorthogonal spatial coordinate norm creates the split.

## Every complementary mode

The local omitted amplitude is `L_x=Q† P_x C conjugate(W)`. Its complete
cross-site Gram has the exact identity

```
L_x† L_y = S† (delta_xy B_x-B_x B_y) S.
```

The stored `B_x,S` therefore determine every cross-site leakage amplitude
Gram without storing a redundant 343-by-343 array of four-row matrices.
For any real `f`, put `B[f]=sum_x f_x B_x`. Then

```
L[f]† L[f]=S† [sum_x f_x^2 B_x-B[f]^2] S >= 0.
```

The implementation independently forms `f C conjugate(W)-W Y[f]` on all
686 rows and checks this identity. No complementary mode is discarded by
assuming a strictly invariant local Higgs insertion.

The uniform profile gives `sum_x B_x=I`, `Y[1]=S` and `L[1]=0` exactly.
Its four singular values equal one. A uniform primitive scalar interaction
on this carrier therefore supplies no spatial hierarchy. Localized profiles
produce a different measurable response on the same spatial band.

## Geometric intervention controls

The complete local kernel is calculated before these diagnostic contractions.
The profiles use only the declared core coordinates. A unit scalar insertion
at the positive core gives rank two. The indicator of the 27-site cube
`||x-(0,0,1)||_infinity<=1` has numerical singular values

```
(0.999999981, 0.999999981, 0.0000610314301, 0.0000610314301).
```

Its largest/smallest strength ratio is about 16,385, and the complete
complement leakage norm is about 0.00781202. The small strength measures
the spatial tail of the other localized channel in the specified patch.
These are numerical source diagnostics at the stated tolerance, not
outward interval certificates. No empirical mass selects the patch or
any source coefficient. The complete site kernel permits every alternative
scalar profile to be evaluated on the same domain.

The weak strength has a direct hopping interpretation. The negative core has
exactly one nearest neighbor inside the positive-core patch, at `(0,0,0)`.
The literal hopping column into that site has Gram `I2/4`; the exterior
onsite inverse contributes `1/M` in amplitude. Its first-hop squared overlap
is therefore `1/(4M^2)`. The full source calculation at the declared controls
`M=32,64,128` retains every spatial row and compares this leading term with
the actual weak strength. This is a source expansion diagnostic, with no
outward remainder certificate. Since `Y[f]=B[f]S` and `S` is unitary, a
nonnegative profile has strengths equal to the eigenvalues of `B[f]`, not
their square roots. The overlap vertex involves two localized wavefunctions.

Dividing a profile by its Euclidean norm divides every strength and the
leakage by that common norm; the strength ratio is unchanged. This operation
is an explicit coordinate normalization. A physical scalar kinetic metric
or pole residue is not supplied. In particular, the patch divided by
`sqrt(27)` is not asserted to be a canonically normalized Higgs mode.

## Resolvent residue and the finite intervention

The parent's four core-site coordinates are not themselves orthonormal
physical localized eigenstates. In the complete core/exterior split write
`h=[[0,B†],[B,D]]`. At each actual band energy `E`,

```
psi=(u,-(D-E)^-1 B u),
Z(E)=I+B†(D-E)^-2 B,
||psi||^2=u† Z(E) u.
```

The calculation reconstructs every complete normalized eigenstate from this
resolvent and checks `u†Z(E)u=1`. The resulting `W` supplies the normalized
spatial response above. This spectral residue is a finite resolvent metric;
its identification with a relativistic field residue is a separate premise.

For an additional real amplitude `t`, retain the complete Nambu matrix

```
mathcal B(t)=diag(h,-conjugate(h))
            +t [[0,fC],[(fC)†,0]].
```

Its eight-row band compression is
`[[A,tY[f]],[tY[f]†,-conjugate(A)]]`. Let `ell=||L[f]||`,
`eta=|t| max_x |f_x|`, and assume `d=62-r-2eta>0`. The complementary
resolvent at any low-band energy has norm at most `1/d`. The exact Schur
self-energy consequently has norm at most `t^2 ell^2/d`. Its associated
residue metric is

```
Z_t(E)=I+t^2 V_PQ (E-B_QQ(t))^-2 V_QP,
0 <= Z_t(E)-I <= (t^2 ell^2/d^2) I.
```

The full calculation diagonalizes all 1,372 Nambu rows for the uniform and
patch profiles divided by their Euclidean norms, at `t=1/8`. It compares
all eight low eigenvalues with the band compression. The numerical `ell`
is used only for a diagnostic evaluation of these exact formulas. The
Schur comparison is not promoted to an outward certificate by a floating
residual. Both full and compressed matrices retain the nonzero spatial
matrix `A`; setting it to zero changes the source.

The source establishes a concrete way for same-current localized modes to
have very different responses to one scalar intervention. A scalar Hamiltonian
and its selected collective mode determine which contraction `Y[f]` is
physically relevant. Three generations, a physical Higgs kinetic residue,
relativistic propagation and measured Yukawa values do not follow from this
finite two-strength example.
