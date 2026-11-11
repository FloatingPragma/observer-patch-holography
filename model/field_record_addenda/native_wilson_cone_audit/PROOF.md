# Matching the free relativistic dispersion coefficient

This is an audit and a separate proposed architecture, not a modification
of the frozen native Wilson–Maxwell cube. Its matter rest frequency m
is the same selected native excitation. No measured mass or speed is
inserted. The source matrices alpha_j and beta are read directly from
`model/native_wilson_maxwell/current.py`.

## Increasing volume does not remove the current mismatch

On the translation-invariant bulk symbol, the existing hopping matrices
and onsite term give

`h(p)=sum alpha_j sin(a p_j)/a + beta[m+sum_j(1-cos(a p_j))/a]`.

Put `w_j=1-cos(a p_j)`, W=sum w_j, and let the lattice Maxwell
frequency be `Omega_gamma^2=2W/a^2`. Clifford anticommutation and
`sin^2(a p_j)=2w_j-w_j^2` give the exact identity

`E_F^2-m^2=(1+ma)Omega_gamma^2+2 sum_(i<j)w_i w_j/a^2`.     (1)

The frozen cube chooses a=1/m. Along each coordinate axis, (1) says
`E_F^2-m^2=2 Omega_gamma^2` at every momentum. In particular the
massive energy curvature and the photon frequency have different
small-momentum coefficients. One common time or length conversion
does not change this ratio. This is a free-regulator statement. It
does not assume a pole or a dispersion law for the interacting
many-body ground and cannot by itself exclude all possible interacting
renormalization effects.

The coefficient is not the group velocity at zero momentum, which
vanishes for a massive branch, or a proved maximum signal speed. It
is the coefficient required by `E^2=m^2+c^2|p|^2` when the same
physical momentum and clock are used for matter and the Maxwell field.

The usual spacing refinement is one valid repair route: at fixed m,
let a tend to zero while increasing the physical volume as separately
specified. Then ma tends to zero. This changes the source family;
increasing the number of sites at a=1/m is not that refinement.

## An explicit finite regulator change

There is also a simple alternative that removes the ma cross term
without a second mass input. Define, on the same four components,

`Gamma=alpha_1 alpha_2 alpha_3 beta`.

The exact matrices give Hermitian Gamma with square I, anticommuting
with every alpha_j and beta. It commutes with the three physical
spin-rotation generators. Change only the Wilson regulator direction:

`h_perp(p)=m beta+sum_j alpha_j sin(a p_j)/a+Gamma W/a`.    (2)

This is an explicit architecture choice. Its squared spectrum is

`E_perp^2=m^2+sum_j sin^2(a p_j)/a^2+W^2/a^2`
`       =m^2+Omega_gamma^2+2 sum_(i<j)w_i w_j/a^2`.        (3)

Along every coordinate axis the massive and photon relations now agree
exactly. In all directions their small-momentum quadratic coefficient
is one. This does not make a finite cubic lattice exactly Lorentz
covariant: the last term in (3) is a directional cutoff correction.
Using `0<=w_j<=(a p_j)^2/2` and
`0<=p_j^2-4 sin^2(a p_j/2)/a^2<=a^2 p_j^4/12`,

`-a^2 |p|^4/12 <= E_perp^2-(m^2+|p|^2) <= a^2 |p|^4/6`.

For m>0 this gives
`|E_perp-sqrt(m^2+|p|^2)|<=a^2 |p|^4/(12m)`.
At m=0, the only zero in the Brillouin zone is at p=0: W is a sum
of nonnegative terms, vanishing only where each cosine is one.
The usual corner doublers have W>0 and are lifted. The four-component
field is still a Dirac field, not an unpaired chiral Weyl species.

## Local gauge coupling, spin and the exact finite background bound

The new directed hopping matrices are `B_j=-(Gamma+i alpha_j)/2`.
Each remains a rank-two partial isometry. On the complete two-site
one-body carrier the bond eigenvalues are -1/a twice, zero four times,
and +1/a twice. Its full Fock norm is 2/a. Thus the local charge,
Peierls phase coupling and the current-algebra norm estimates have the
same structural bounds as before; their actual ground contractions
and preparation receipts must be recalculated for the new member.

For arbitrary Hermitian finite graph links, including a finite gauge
background, write the one-body operator as `h=m beta+A`. Every
hopping/onsite term in A contains alpha_j or Gamma and hence
anticommutes with beta. Therefore `h^2=m^2+A^2>=m^2`. This statement
concerns that one-body operator. Adding a dynamical field Hamiltonian
and filling a many-body sea does not inherit this gap by notation.

The original charge-conjugation antiunitary C0 K, `C0=i beta alpha_2`,
has `C0 alpha_j* C0^dagger=alpha_j`,
`C0 beta* C0^dagger=-beta`, and `C0 Gamma* C0^dagger=-Gamma`.
It preserves the required particle/antiparticle energy transformation
when momenta and gauge phases are conjugated. Ordinary parity at
nonzero spacing changes the sign of the new regulator; exact ordinary
parity is therefore not claimed. The charge representation must be
common to these four components. This is a vectorlike charged-field
construction; it does not gauge the Standard Model's distinct weak
left/right representations without a Higgs construction.

## What this does and does not import from the literature

Rotating the Wilson and mass terms relative to one another is related
to twisted-Wilson improvement. Frezzotti and Rossi discuss improved
lattice observables at maximal twist in
[hep-lat/0306014](https://arxiv.org/abs/hep-lat/0306014).
Aoki and Baer explain why interacting improvement depends on the
definition of the twist in
[hep-lat/0409006](https://arxiv.org/abs/hep-lat/0409006).
These are useful construction precedents, not proofs that this OPH
Hamiltonian has an interacting relativistic continuum or renormalized
mass. Equations (1)-(3), the finite matrix bounds and the stated
symmetries are proved here from the actual matrices, independently.

To use (2) in the complete native source one must recalculate its
Foldy-Wouthuysen code map, completed full ground, physical preparation,
compensated currents and energy/spatial detector law. A common free
quadratic cone is a check on that new construction, not a substitute
for those interacting calculations or for measured GeV masses.
