# A charged band does not determine its electromagnetic current

The native Dirac code supplies a free positive-energy dispersion and an
explicit spinor readout. Coupling its scalar square-root dispersion by
Peierls phases is a distinct action choice from coupling that readout's
Dirac current. They have different finite matrix elements.

For a mass m>0 and Omega(p)=sqrt(m^2+|p|^2), let U(p) be the normalized
positive-frequency frame of `model/native_block_dirac`. Compare incoming
p=(0,0,-k) with outgoing p'=(0,0,k), and take polarization epsilon=(1,0,0).
The momentum transfer is q=(0,0,2k), transverse to epsilon. Direct
four-component multiplication gives

```
U(p')^dagger U(p)=(m/Omega) I2,
U(p')^dagger alpha_x U(p)=i(k/Omega) sigma_y,
U(p')^dagger alpha_z U(p)=0.
```

In contrast, differentiating the positive scalar square root by the
Sylvester equation `Omega deltaOmega+deltaOmega Omega=delta(P^2)` gives
the scalar minimal-coupling vertex

```
<p'|deltaOmega|p> proportional to
       (p'+p) dot epsilon/[Omega(p')+Omega(p)]=0.
```

At k=3m/4, the density is 4I2/5 and the Dirac transverse spin-flip
matrix element has magnitude 3/5. The actual native selected m sets the
unit; no measured mass is used. These values persist at fixed physical
momentum in the continuum comparison. This is not a discretization
correction. The zero longitudinal current is compatible with the
zero transferred energy; the continuity equation does not determine the
nonzero transverse current from the charge density alone.

The result excludes identifying scalar square-root Peierls motion with
the established free Dirac readout's electromagnetic current. It does
not exclude a charged scalar-band model on its own action. Equal charge,
statistics and free energy do not fix an interaction vertex.

## A local canonical matter interface

A concrete permitted regulator on the same sixteen-occupation native
code is

```
h_W(p)=sum_j alpha_j sin(a p_j)/a
       +beta [m+sum_j(1-cos(a p_j))/a],
Omega_W(p)^2=sum_j sin(a p_j)^2/a^2
       +[m+sum_j(1-cos(a p_j))/a]^2 >= m^2.
```

The corresponding positive native kinetic matrix Omega_W-m is
nonnegative. Its momentum-dependent Bogoliubov map supplies all four
canonical site fields Psi. Its zero-energy native occupation state is
the filled negative-band state in these field coordinates. The local
charge is rho_x=Psi_x^dagger Psi_x-2, with integer spectrum -2 through2;
the total charge is N_plus-N_minus. Local rho_x differs from the
original onsite native particle-minus-antiparticle number and contains
pair terms in those coordinates. The exact control gives a nonzero
U(p)^dagger V(0) density matrix and U(0)^dagger alpha_x V(0)=sigma_x.
Dropping those terms changes both the charge distribution and current.

A Wilson bond has block -(beta+i alpha_j)/(2a). Its full two-site
one-body spectrum is +/-1/a, each with multiplicity2, and four zeros.
The full Fock norm is 2/a, the sum of its positive one-body eigenvalues.
For genuinely local unit-norm link operators, six incident bonds of
two-cell support give the weighted interaction bound
24 lambda_spatial/a. This is the precise matter input to the existing
finite causal observer estimate. A finite-mode Coulomb-gauge photon
representation with a total occupation cutoff does not automatically
give such local tensor link operators; its support and cutoff error
need their own proof. The filled-band state, complete Gauss reference,
photon carrier and coupled Hamiltonian are additional data before this
interface is an interacting source.

No native Hamiltonian or frozen free package is changed by this audit.
An interacting extension must specify its actual ground, or label an
explicit finite nonstationary preparation, and keep every source and
photon complement. These current identities alone supply neither a
Maxwell source nor interacting relativistic scattering.
