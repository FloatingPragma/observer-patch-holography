# Magnetic spin response of the gauged native kinetic family

## Declared source and retained modes

The sixteen native kinetic matrices `gamma_j` and four mass matrices are
those of `model/native_mass_channels/algebra.py`. They retain all eight
folded microscopic valley labels. The scalar mass matrix `beta0`
anticommutes with every `gamma_j`. A finite gauged spatial extension is
specified by

```
H = sum_j gamma_j tensor Pi_j + m beta0 tensor I,
Pi_j = i(T_j-T_j*)/(2a),
(T_j psi)(x)=exp[i q a A_j(x)] psi(x+a e_j).
```

The finite spatial calculation has open boundaries in x and y and
`Pi_z=0`. No native internal row is discarded. The uniform external
field is `A_y=Bx,A_x=0`. This Peierls gauge extension is a declared
architecture; equivalence to the primitive interacting repair source
under an external electromagnetic field is not proved here.

The exact link current is obtained by differentiating this Hamiltonian
with respect to the link potential. It is proportional to `gamma_j`,
with coefficient `-q/2` multiplying the Hermitian link sum. Its derivative
is checked independently by finite differences. Simultaneously changing
the link gauge and the probe phase preserves the response; changing only
the links is rejected by a controlled comparison.

## The curvature coefficient follows from the source square

The Clifford identities give the exact finite operator identity

```
H^2 = m^2 I + I tensor (Pi_x^2+Pi_y^2)
      + Sigma_z tensor i[Pi_x,Pi_y],
Sigma_z=-i gamma_x gamma_y,        Sigma_z^2=I.
```

`Sigma_z` commutes with `beta0`. Each of the four joint eigenspaces of
`beta0` and `Sigma_z` has rank four. Hence the complete square reduces to

```
K_+ = Pi_x^2+Pi_y^2+i[Pi_x,Pi_y],
K_- = Pi_x^2+Pi_y^2-i[Pi_x,Pi_y].
```

The code verifies the complete 144-row unreduced Hamiltonian on a small
box, both signed energy branches, and all multiplicities against this
reduction. Neither a Pauli term nor a gyromagnetic factor is inserted.
The spin-curvature term is the cross term of the actual gauged source
operator. Dropping it fails the complete square identity.

For smooth fields and low model momentum,
`i[Pi_x,Pi_y] -> qB I`. The positive-branch energy has the expansion

```
sqrt(m^2+K) = m + K/(2m) + O(K^2/m^3).
```

With spin `S_z=Sigma_z/2`, the magnetic term is
`q B Sigma_z/(2m)=q B S_z/m`, up to the stated electric-charge orientation.
Thus the magnitude of the tree-level coefficient defined by
`H_mag=g q B S_z/(2m)` is exactly `g=2` in this limit. The factor two
follows from the native Clifford algebra and the spin normalization;
it is not fitted to a measured magnetic moment. This argument supplies
no anomalous magnetic moment.

## Controlled finite response and coherent phase

The executed source has `m=64`, `q=1`, `B=1/512` and a fixed Gaussian
probe width of four model length units, multiplied by a cosine taper
that vanishes on the box boundary. Spacings `1,1/2,1/4` use grids
33,65,129 on the same spatial interval. The largest unreduced
spin/valley/spatial carrier has 266256 one-body rows. Exact internal
reduction leaves two 16641-row operators.

Since `||Pi_x||,||Pi_y||<=1/a`, the complete source bound is
`0<=K_+/-<=4/a^2`. For `x=K/m^2<=1`, the sixth-order square-root
polynomial has a uniform remainder at most

```
R = m (33/2048) [4/(m^2 a^2)]^7.
```

The executable evaluates the polynomial as a sparse matrix, subtracting
only the common displayed rest coefficient. Both branches retain a
uniform analytic truncation error `R`; their difference has error at most
`2R`. Floating-point evaluation and eigensolver diagnostics are separate
numerical checks, not interval-certified rounding bounds.
The operational response coefficient is

```
g_eff = 2m [<E_+>-<E_->]/(qB),
|g_eff-g_exact| <= 4m R/|qB|.
```

The same two source branches evolve a coherent spin superposition for
ten source-time units. Their overlap gives the relative phase and
visibility. Its phase rate agrees with the energy response. The probes
are finite wave packets, not exact stationary localized eigenstates;
this is an external-field spin response with a finite observation time.
The energy absolute-value branches are those of the source spectrum,
or equivalently its positive block after the exact massive Clifford
block diagonalization. The spin pointer is specified in that branch
frame; a microscopic spatial spin identification is a separate claim.

The values approach two under refinement, with finite-spacing corrections
and the relativistic remainder retained. Scalar-current replacement,
zero link flux and omitted curvature controls each remove or corrupt
the predicted response. The construction supplies no stationary
self-sourced magnetic dipole, physical electron identification, measured
Bohr magneton or radiative correction. The mass coefficient, external
field, gauge extension and clock normalization are explicit inputs.
