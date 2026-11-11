# Complete scalar modes and native filled-sea backreaction

## Source, scalar cap and reference

Use the complete radius-three, amplitude-64, half-separation-one source in
[compact_defect_pair](../compact_defect_pair/PROOF.md). Its 343 sites and 686
spin coordinates include the entire exterior. The native matter carrier is
`F(C^686 tensor R15)`, with 10,290 fermion modes and no inserted family factor.
Each of the three native Yukawa tensors has primitive coefficient one. The
local pair has the same fixed phase `D_x=(-1)^(x1+x2+x3)` as the compact
particle-hole operator `C=D i sigma2`. All source parameters are declared
architecture inputs; empirical masses do not enter them.

The scalar carrier has a single global occupation cap, rather than a tensor
product of per-site caps. Its states are the vacuum, every site of a particle
doublet, and every site of an antiparticle doublet. The dimension is
`1+4*343=1373`. Define

```
phi_xk=|0><particle_xk|+|antiparticle_xk><0|,
K=I+L_graph+64^2 diag(|n_x|^2),
Omega=positive_sqrt(K),
H_H=diag(0,Omega tensor I2,conjugate(Omega) tensor I2).
```

The graph Laplacian uses every internal nearest-neighbor edge of the open
cube. The texture amplitude is zero at the two cores and one elsewhere.
`K>=I`; its square root is the unique positive one. This is a declared
finite scalar architecture, with no global canonical commutator assertion
at the occupation cap. The free law is
`H0=dGamma(h tensor I15)+H_H`. The interacting law adds every literal native
site term `H_Y`. Its clock is `exp(-it(H0+H_Y))`, with hbar one.

The scalar-vacuum response below uses the scalar blank under H_H. The
backreaction reference is the Slater determinant that fills all 343 negative
spatial eigenmodes with all fifteen native internal labels, tensored with
the scalar blank. It has 5,145 fermions. Neither reference is identified
with the ground of H0+H_Y.

## Complete scalar spectrum and canonical response

The archive retains every orthonormal scalar eigenvector `u_j`, every
frequency `omega_j=sqrt(lambda_j)`, K, Omega and the coordinate metric.
There are 343 modes for each of the four scalar branches. The real symmetric
K has nonpositive connected off-diagonal entries, so its lowest eigenvector
is unique and strictly positive by the Perron theorem. Core-exchange
reflection makes that mode even. The implementation resolves the two bound
roots through the complete exterior Schur complement to avoid mixing their
close numerical eigenvalues. It also retains all 341 exterior modes.

For either doublet component, the retarded scalar matrix on the blank is

```
D_phi(z)=(z-Omega)^(-1)-(z+Omega)^(-1)
       =2 Omega (z^2-K)^(-1),             Im z>0.
```

This follows directly from the two nonzero matrix units in each phi. Its
zeroth spectral moment is zero and its first moment is 2 Omega. Thus the
canonical vacuum coordinate and its response are

```
q=(2 Omega)^(-1/2) phi,
D_q(z)=(z^2-K)^(-1).
```

A primitive coupling to phi_j therefore has coefficient
`sqrt(2 omega_j)` times its unit-ladder mode overlap when written as a
coupling to q_j. Dividing the overlap by that factor would invert the
specified change of coordinate. This normalization is a statement about
the free blank reference. The interacting spectral moment must be evaluated
in the interacting state before assigning its own pole normalization.

The lowest stiffness eigenvalue is approximately 6.99829101533, with
frequency 2.64542832361 and equal core amplitudes 0.707106633672. Its four
unit-ladder spatial strengths are, in pairs,

```
(0.706890883991, 0.706890883991,
 0.706804630458, 0.706804630458).
```

Their ratio is 1.00012203306. Canonical vacuum strengths are approximately
1.62598067048 and 1.62578227129, each twice. The complete exterior leakage
norm of the ladder insertion is 0.01460788368. These are numerical
diagnostics, not outward interval values. The large ratio of the separately
declared positive-core patch in
[defect_mode_overlap](../defect_mode_overlap/PROOF.md) is not the response
of this source-selected free scalar mode. Switching the scalar defect term
off gives K=I+L_graph and the constant lowest profile `1/sqrt(343)`, with
four equal spatial strengths. That is a distinct declared control.

## Exact global singlet and complete native coefficients

Let T_A be any of the twelve actual native current matrices. Filling all
fifteen internal labels of a fixed spatial mode transforms as the determinant
representation; its infinitesimal charge is `Tr T_A`. The code checks that
all twelve traces vanish. A product of these filled blocks and the scalar
blank is therefore a GLOBAL native gauge singlet. The actual native tensor
Ward identities imply `[Q_A,H_Y]=0`, so its first interaction image is also
a singlet. H0 is identity on native internal indices and has the same
symmetry. This argument supplies no quantum links or local Gauss law.

A retained idle conjugate copy of the sea serves as reference. Its odd
fermion parity compensates the odd source parity, giving even combined
parity; H_Y preserves both parities. The reference is unchanged by every
instrument in this package. It is part of the preparation and each output,
including failures, rather than a scalar factor removed by conditioning.

For either scalar component, the Hilbert-Schmidt Gram of the three literal
symmetric native pair tensors, ordered up, down, charged lepton, is exactly
`diag(6,6,2)`. These coefficients are derived from the actual fifteen-label
representation, without a lepton-only projection. Cross-channel terms
vanish by this Gram. The complete quantum CAR control separately evaluates
an eight-mode Fock space and all five states of a two-site, one-component
scalar cap, checking both particle and antiparticle pair kernels and their
fermion energy moment.

## Complete first interaction image

Let P be the negative spectral projector of the full 686-row h, and let
P_xy be its two-by-two spin block. The phase-weighted pair Gram is

```
G_xy=D_x D_y det(P_xy).
```

Wick's theorem gives this determinant from both spin contractions. The
identity `h C=-C conjugate(h)` and the absence of zero energies imply
`P+C conjugate(P) C^dagger=I`. Hence `Tr P_xx=1`, creation and annihilation
pair Grams are conjugates, and G is real for this source. The code retains
its complete matrix and rejects an imaginary part beyond tolerance.

The literal up channel multiplies phi, and the down and charged-lepton
channels multiply phi-dagger. On the scalar blank, the unnormalized reduced
scalar matrices for each component are consequently

```
R_particle=8G+6 conjugate(G),
R_antiparticle=6G+8 conjugate(G).
```

Both components of each branch occur. In particular,
`c^2=||H_Y |sea,0>||^2=56 Tr G`, and the complete scalar frequency record is

```
p_j=2 [u_j^T R_particle u_j+u_j^T R_antiparticle u_j]/c^2,
sum_(j=0)^342 p_j=1.
```

The modal projector sums over every component and both particle branches,
so it is gauge invariant. The component entries in the receipt are
unnormalized reduced-state weights, not a claim that an uncompensated
charged component is a gauge-invariant classical record.

The executable obtains `c^2=29.01797916421`. The lowest even scalar mode has
weight 0.000411882689, and the odd bound mode has weight 0.964506586904.
Every remaining mode, including its small weight, is retained. The largest
eigenvalue of G is approximately 0.499786435425. Its normalized profile has
opposite core values very close to plus/minus `1/sqrt(2)` and four equal
spatial strengths approximately 0.70684782963. The profile characterizes
the reduced first image. It is not a dressed Higgs pole or an interacting
vacuum expectation.

For a complete preparation instrument, set `B=28*343`. At each site the
literal scalar vacuum/occupied interaction is a self-adjoint star with
four outgoing CAR blocks A_b. The absolute coefficient sum of each block
is exactly fourteen, hence `||A_b||<=14`,
`||H_Y,x||=||column_b A_b||<=sqrt(4*14^2)=28`, and the triangle inequality
gives `||H_Y||<=B` on the full Fock/scalar carrier. Define

```
K_success=H_Y/B,
K_failure=positive_sqrt(I-H_Y^2/B^2).
```

The effects sum to identity. After success, measure every invariant scalar
mode projector. The original mode probabilities are `c^2 p_j/B^2`; the
failure probability is `1-c^2/B^2`. There are 344 classical outcomes, with
both scalar orientations and all doublet coordinates retained. The
success probability is small and is not replaced by one. This construction
is a finite spectral instrument definition, with no claim of an efficient
laboratory circuit.

## Full-law variational bound

Write `|s>=|sea,0>` and `|v>=H_Y|s>/c`. They are orthonormal because their
scalar occupations differ. H_Y only links the blank and occupied scalar
sectors, so its diagonal expectation on either vector is zero. In this
entire two-dimensional subspace, the exact Rayleigh-Ritz matrix relative
to E_sea is

```
[[0,c],[c,d]],
d=<s|H_Y(H0-E_sea)H_Y|s>/c^2.
```

The scalar part of the numerator is `56 Tr(Omega G)`. If
`Q=P h`, its fermion part is
`-56 sum_x Tr(adj(P_xx) Q_xx)`, retaining both particle and hole excitations.
The direct small Fock control checks this energy identity independently.
The complete source gives d approximately 9.29053640750 and the smaller
Ritz eigenvalue approximately -2.46785318204. This is a numerical upper
bound estimate for the true singlet ground energy relative to the sea,
not an exact ground-state reconstruction or an outward certificate for
that decimal value.

A separate exact margin establishes a strict inequality without relying
on those floating values. Every run regenerates the parent's complete
686-eigenvalue rational certificate. It verifies `|h|>=delta=1/1000`; the
parent hopping estimate gives `||h||<=66`. On a normalized core spin vector
u, `<u,h u>=0`. Splitting into positive and negative spectral parts gives

```
a I2<=P_xx<=(1-a) I2,     a=delta/(66+delta)=1/66001.
```

Together with `Tr P_xx=1`, this implies
`det P_xx>=a(1-a)` at each of the two cores and thus
`c^2>=112a(1-a)=7392000/4356132001`. The two fermion excitations cost at
most 132 and `||Omega||<=sqrt(4109)<65`; hence `d<=197`. The singlet trial
vector `|s>-H_Y|s>/197` therefore proves

```
E_ground(H0+H_Y) <= E_sea - 1456224000/169057134218809
                 < E_sea - 0.0000086137.
```

The ground here is the lowest energy in the same global native singlet
domain. The strict bound applies to the full untruncated fermion carrier
and all scalar states. Neither the lowest free scalar profile nor the
filled free sea can be substituted for that interacting ground.

## Counting, artifacts and scope

The complete compact fermion band has dimension four. Tensoring with the
single native fifteen-label module gives sixty band CAR coordinates. For
any one real scalar profile, its four-by-four antisymmetric overlap matrix
has at most two distinct singular strengths, each repeated twice. This
linear-algebra count is neither two physical generations nor three physical
generations. The collection of scalar mode tensors need not share a common
family diagonalization, and the local interaction leaks into the complete
fermion exterior.

`results.json` retains every scalar frequency, every modal probability,
original success and failure weights, the exact rational margin, scalar
normalization, independent controls and transitive source hashes.
`scalar_modes.npz` retains every scalar mode, both charge-component Gram
matrices, G, its full spectrum and its dominant profile. The checker
regenerates the source, certificate and distributions, then validates
both artifacts. Mutation tests reject truncated modes, missing branches,
forged Grams, postselection and physical generation claims.

The results describe finite simulation-generated response records on the
declared source and preparation. No observed Yukawa hierarchy, exactly
three physical generations, interacting Higgs pole, relativistic mass
shell or continuum gauge theory is inferred from them.
