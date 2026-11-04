# The complete one-particle spectral instrument

Consider the original action on a periodic cube of side eight,

\[
H=d\Gamma(h)+U\sum_x n_{x0}n_{x1},\qquad
h(k)=-\frac18\sum_{j=1}^3\sin(k_j)\sigma_j,
\quad k_j=\frac\pi4\ell_j,\quad \ell_j\in\{0,\ldots,7\}.
\]

The contact operator annihilates every one-particle vector, for any U.
Number is conserved by both hopping and contact. Hence this entire sector
is exactly invariant under the original action. No particle-number
postselection, shifted chemical potential or replacement vacuum is used in
this restriction. In particular its lowest energy is negative,
`-sqrt(12)/16`; the empty Fock state is not a positive-energy ground reference.

## Exhaustive spectrum

The eight values of `2 sin(pi*ell/4)` lie in `Q(sqrt(2))` and have squared
values zero, two and four with multiplicities two, four and two. Let
`n=2 sum_j sin(k_j)^2`. The Pauli relations imply

\[
(16h(k))^2=2nI_2,\qquad \operatorname{tr}h(k)=0.
\]

There are `d_n=8 binomial(6,n)` momenta in shell n, since their generating
polynomial is `(2+4z+2z^2)^3=8(1+z)^6`. For `n>0` the two energies and their
momentum-block projectors are

\[
E_{n,\pm}=\pm\frac{\sqrt{2n}}{16},\qquad
P_\pm(k)=\frac12\left(I_2\pm\frac{16h(k)}{\sqrt{2n}}\right).
\]

The square identity proves idempotence, orthogonality and completeness;
the trace proves that each projector has rank one. At n=0 the whole spin
block has energy zero, and both component projectors are retained. Exact
arithmetic checks every distinct block and all 512 momentum labels.

| n | Energy squared | Multiplicity of each signed energy | Maximum one-site probability |
|---:|---:|---:|---:|
| 0 | 0 | 16 total | 1/64 |
| 1 | 1/128 | 48 | 3/64 |
| 2 | 1/64 | 120 | 15/128 |
| 3 | 3/128 | 160 | 5/32 |
| 4 | 1/32 | 120 | 15/128 |
| 5 | 5/128 | 48 | 3/64 |
| 6 | 3/64 | 8 | 1/128 |

The multiplicities add to 1,024. There is no discarded momentum sector or
additional one-particle band hidden among unlabelled outputs. The zero
momenta are the eight vectors with coordinates zero or pi. They are the
nodes of the signed Weyl bands. Finite nonzero shell energies are not
relativistic rest masses or independent particle species.

## Complete detector and source implementation

Use the original [native Fourier network](../particle_species_detector/README.md),
which retains all 1,024 physical modes. Write `16h=d.sigma`,
`e=|d|=sqrt(2n)`. When `d_x^2+d_y^2>0`, apply the two-mode matrix

\[
R(k)=\frac1{\sqrt{2e(e+d_z)}}
\begin{pmatrix}e+d_z&d_x-id_y\\-d_x-id_y&e+d_z\end{pmatrix}.
\]

It has determinant one, is unitary, and satisfies
`R h R*=diag(e,-e)/16`. For a positive z-axis direction use the identity;
for a negative z-axis direction use `[[0,1],[-1,0]]`; for d=0 use the
identity and label both outputs zero. These charts include every momentum
without division by zero. The phases of the two zero-mode components fix
an apparatus convention, not a preferred basis of the zero eigenspace.

After the Fourier network, each pair `(k,0),(k,1)` is adjacent in the
retained mode ordering. Apply its fermionic lift `Gamma(R)`, then copy all
occupations to fresh blank pointers. On N=1 this is a projective instrument
with 1,024 effects. Their energy labels agree with the original H exactly
in the ideal map. Vacuum, every multiple-occupation string, and all native
source/export failure branches are also outcomes of the complete device.
On N>1 these occupation strings do not generally label the interacting H.

For each nontrivial rotation, put
`alpha=arg(d_y+i d_x)` and `t=acos(d_z/e)/2`. Its chronological word is

\[
D^*\ ;\ \exp[-it(c_0^*c_1+c_1^*c_0)]\ ;\ D,
\qquad D=\exp(i\alpha n_0).
\]

Thus the one-particle upper off-diagonal is
`-i exp(i alpha) sin(t)`. The two native commuting XX and YY Pauli words
implement the middle operation, and both scalar and Z phase terms are
retained in D. The vacuum and double-occupation entries of this determinant
one Fock lift are one. All four columns of each of 125 distinct local blocks
are independently executed through the existing compiler.

Inverse trigonometric functions at fixed 80-digit precision only propose
literal binary64 words. `native.py` uses the integer dyadic interval kernel to
enclose their sine and cosine, the exact radical target entries and every
native port-control error. The full Fock operator error is bounded by the
entrywise error plus fixed-angle and source-control errors. Summing these
bounds over 492 nontrivial rotations is valid on arbitrary Fock input and
untouched references. The ideal rotations are Gaussian; native residuals
are retained in the full source output rather than projected onto N=1.

The added rotations cost 12,792 native pulses. Including Fourier routing and
all 1,024 occupation pointers, detection costs 32,925,688 formal pulses and
retains 121,573,504 first-export abort slots. The certified complete
isometry error is below `4.155e-9`. For the same supplied input, contractivity
of the subsequent readout gives unhalved output L1 error at most twice that
bound. These counts exclude state preparation and field evolution. The
1,024-column ideal detector and local four-column native blocks are
executed; no `2^1024`-dimensional native array is executed.

## Finite stationary localization

Let `P_E` be the full spectral projector at one of the energies above.
Translation invariance makes its on-site block independent of x. For a
positive shell the exact sum of each component of d over its momenta is
zero, since momentum inversion preserves the shell and reverses d. Therefore

\[
P_E(x,x)=\frac{\dim P_E}{1024} I_2.
\]

The same formula holds for the sixteen-dimensional zero eigenspace by
direct counting. If a normalized state lies in this eigenspace, its
occupation probability at x is at most `dim(P_E)/1024`: the nonzero
eigenvalues of `P_E P_x P_E` and `P_x P_E P_x` agree. The bound is attained
by `P_E|x,s>/sqrt(dim(P_E)/1024)` for either component s.

Any stationary one-particle density matrix commutes with H, hence is a
convex combination of density matrices supported in its distinct energy
eigenspaces. Its one-site probability is consequently at most `5/32`.
For a region C of m sites, summing the site bounds gives

\[
\operatorname{tr}(\rho N_C)\leq\min(1,5m/32).
\]

At least `1/16` of the occupation lies outside every six-site region. This
is a necessary support bound; no claim is made that seven sites suffice.
It applies to stationarity under this original N=1 Hamiltonian. Prepared
wave packets, moving cores, defects in other number sectors and changes to
the repair or memory law are outside this finite exclusion.

## Infinite lattice

On `Z^3` with arbitrary `c>0`, the unchanged one-particle symbol is
`h(k)=-c sum_j sin(k_j) sigma_j`. For any real E,

\[
\det(h(k)-EI_2)=E^2-c^2\sum_j\sin^2(k_j).
\]

For compactly supported wavefunctions the Fourier components are Laurent
polynomials in `z_j=exp(ik_j)`. Multiplication by the adjugate gives the
determinant times each component equal to zero. Its coefficient at `z_1^2`
is `c^2/4`, so it is a nonzero Laurent polynomial at every E, including
E=0. The Laurent polynomial ring is an integral domain; both wavefunction
components must vanish. `laurent_certificate()` checks the symbolic
determinant, adjugate identity and nonzero coefficient.

There is also no nonzero square-summable stationary one-particle eigenvector.
For fixed `(k_y,k_z)`, the zero-determinant equation has at most four
solutions for k_x in the Brillouin circle. Fubini's theorem makes its full
zero set a set of three-dimensional measure zero. Away from it the Fourier
matrix is invertible. Plancherel's theorem then forces any square-integrable
eigenvector to vanish. This assertion concerns point spectrum, not the
continuous Weyl spectrum or the finite torus spectrum. At c=0, H is zero
on N=1 and one-site stationary states are counterexamples; that parameter
is explicitly excluded.

The exclusion also covers every nonzero positive trace-class stationary
one-particle density operator on the infinite lattice. Such an operator is
compact and has a strictly positive eigenvalue with a nonzero
finite-dimensional eigenspace. Commutation with the bounded self-adjoint
one-particle H makes that eigenspace H-invariant. The restriction of H to
it has an eigenvector, contradicting the point-spectrum result. This does
not exclude translation-invariant states of the infinite many-body algebra,
which need not be trace-class densities in the one-particle representation.

This sector supplies no stationary localized core to promote to a protected
particle. The result does not exclude interacting many-body bound states,
localized resonances or a different common source realization. It fixes
the domain of an exhaustive search and its obstruction without assigning
physical charge, spin, units or observed-species names.
