# One calibrated quadratic action and its complete finite records

## Construction and input ancestry

The source is the calibrated native 45-component all-left chiral action of
`model/common_source/calibrated_matter.py`, its optional Weinberg kernel, and
canonical quadratic gauge and Higgs fields. The source has two distinct
momenta p=(0,0,1) and -p in GeV. Natural units have c=hbar=1, the energy unit
is one GeV, and programmed time is inverse GeV. This convention is an
empirical calibration rather than a reconstructed native clock.

The nine Yukawa values, Higgs background v, three gauge couplings and scalar
quartic are empirical inputs. The charged parameters use the frozen pure
MSbar calibration at mu=91.1876 GeV. Its CKM matrix is an additional calibrated
boundary prescription. The optional Weinberg extension uses normal ordering,
zero lightest mass, oscillation mass differences, a calibrated PMNS matrix
and zero Majorana phases. Its low-energy-to-MZ effective-theory transport is
unproved. No native repair or lepton compiler coefficient supplies these
Yukawa inputs. The scalar potential is the explicit tree comparison
V=-lambda_SM v² H†H+lambda_SM(H†H)². Its stationarity at the calibrated loop
v is an extra prescription, distinct from the matched loop effective potential.

The program checks every inherited calibrated-action binding, reconstructs
its matrices and validates 72 exact charged current/Higgs identities and 12
exact Weinberg identities. It freezes all owning code, tests, proof and data
through the repository import closure before and after calculation. Separate
composite records authenticate their own archived source bindings.

## Complete fermionic action

Let K be the reconstructed symmetric 45-coordinate kernel and W its exact
unitary family and neutral mass basis. The program verifies W^T K W directly.
It consists of 21 off-diagonal Dirac 2-by-2 blocks, with color counted, and
three neutral diagonal Majorana blocks. The rank is 44 because the first
neutrino mass is the declared zero input. No extra sterile field occurs in
this 45-coordinate action.

At the two distinct momenta the actual native CAR source has 180 annihilation
operators. `kinetic.momentum_pair_source` builds its h and antisymmetric Delta;
the full Nambu matrix is [[h,Delta],[Delta†,-h^T]], dimension 360. The complete
numerical diagonalization is checked against all positive and negative
frequencies, not selected lines. Its positive subspace is a canonical
Bogoliubov frame: V†V=I and V†tau V*=0. The actual Nambu electric charge
commutes with the source, and its complete positive-subspace spectrum gives
the listed particle and antiparticle charges. The exact mass basis resolves
neutrino masses whose p=1 frequency differences binary64 cannot resolve.

The arbitrary numeric eigensolver frame is only a subspace diagnostic.
The named field frame is constructed separately. Transform by
T=diag(I4 tensor W,(I4 tensor W)*) to the real mass-block basis. For every
coordinate j with kinetic diagonal +1, the unique mass-paired coordinate
has kinetic diagonal -1. If R is the off-diagonal part of this mass-basis
Nambu matrix, its normalized positive eigenvector is

    [e_j+R e_j/(omega+1)] / sqrt(1+m²/(omega+1)²).

Its native vector is T times this vector. Diagonal charge, momentum and spin
generators assign every named mode before any degenerate numerical
diagonalization. For neutrinos the spin is converted to helicity relative
to the actual signed momentum. The formula retains m/(omega+1) for tiny
masses and avoids subtracting numerically indistinguishable frequencies.
The program checks the complete named frame's eigenvector equations and
both CAR identities against the actual native matrix.

The canonical two-component block [[p,m],[m,-p]] squares to (p²+m²)I;
canonical spatial rotations replace p² by |p|². Thus every listed mode obeys
omega²-|p|²=m². The test also diagonalizes the actual full native symbol at
p=(1,2,2), retaining every band. Each normalized canonical fermion mode has
unit spectral weight. The one-momentum p=0 Nambu source has dimension 180,
rank 176 and nullity four; it is not doubled at a self-conjugate momentum.

After vacuum-energy subtraction H_F=sum_j omega_j a_j†a_j on the complete
2^180-dimensional Fock space. `car_action` implements its signed occupation
grammar for every one of these states. Its energy projectors are exactly
P_E=sum_{n:sum_j omega_j n_j=E}|n><n|, combining any degeneracies. A source
particle cap is absent. The declared state is its normalized Bogoliubov vacuum.
This vacuum is a color singlet by uniqueness and SU(3) invariance. Its electric
charge is zero: at zero pairing each internal coordinate contributes one
occupied spin at each momentum, giving 2 Tr(Q45)=0, and scaling every mass
continuously from zero leaves the positive p=1 gap open and preserves charge.

There are 42 probe groups: nine charged species times particle/antiparticle
times two momenta, plus three neutral mass fields times two momenta. Color
and the two spin or helicity components remain coherent within each group.
For M=180 the creation and removal maps are

    K_g,+ = sum_(j in g) a_j† tensor |dual j> / sqrt(M),
    K_g,- = sum_(j in g) a_j  tensor |j>      / sqrt(M).

Orthogonal reference states give sum_g(K_g,+†K_g,++K_g,-†K_g,-)=I by the CAR
on every input state. Vacuum preparation gives 42 positive creation records
with original probabilities dim(g)/180 and 42 zero removal records. Every
other energy projector has zero weight and is retained by an explicit outside
projector. This is a compact exact description of the full finite PVM, not a
postselection. The reference dimension is 361 including its blank state.

## Gauge and scalar kinetic normalization

The gauge kinetic term is -F_a²/4 and the doublet kinetic term is |D H|².
The program takes the actual native doublet matrices t8,...,t11, forms
v_i=g_i t_i(0,v/sqrt(2)), and computes M²_ij=2 Re(v_i†v_j). Exact unitary
diagonalization gives W+,W-,Z,photon. The adjoint electric generator has
eigenvalues +1,-1,0,0 on those vectors. The eight color generators leave
the scalar background fixed, so their quadratic mass coefficients vanish.

The complete four-real-coordinate scalar Hessian has rank one and radial
eigenvalue 2 lambda_SM v². The three null directions are Goldstone coordinates.
The physical coefficients are mW=gW v/2, mZ=sqrt(gW²+gY²)v/2,
mH=sqrt(2 lambda_SM)v, and zero for the photon and gluons. The radial scalar
has kinetic term (partial h)²/2. Eliminating the nondynamical Proca time
component gives longitudinal raw kinetic coefficient m²/(|p|²+m²) and
potential m². The canonical coordinate qL=m A_L/omega therefore has
frequency omega=sqrt(|p|²+m²), equal to both transverse polarizations.
Massless gauge fields contribute two transverse oscillators and no time-like
or longitudinal physical mode. These are free quadratic constraints, not a
local interacting Gauss construction.

The physical inventory is 28 bosonic modes per momentum: W+ and W- each
three, Z three, Higgs one, photon two, and eight gluons each two. The source
therefore has 56 oscillators at the two momenta. The four distinct frequencies
are the W, Z, Higgs and massless values; their multiplicities are 12,6,2,36.

## Finite boson source and complete instrument

Use the global total-occupation cap N<=2, invariant under all mode rotations
that preserve the source charges and frequencies. The Hilbert dimension is
1+56+56*57/2=1653. A separate cap on each color component would lack this
invariance. On this declared domain set a_j=P b_j P and H_B=sum omega_j n_j.
The actual sparse creators contain their bosonic sqrt(n+1) coefficients.
Their commutation relation changes at the top cap; the vacuum and one-particle
field response agree with the uncapped quadratic theory. For
q_j=(a_j+a_j†)/sqrt(2 omega_j), the vacuum positive-frequency spectral weight
is 1/(2 omega_j), and the retarded normalization 2 omega_j times that weight
is one. No zero-frequency massless oscillator is assumed because |p|=1.

Literal multiplication of all 56 creators gives S=sum_j{a_j,a_j†}=56,58,2
on total-occupation sectors zero, one and two. The 12 coherent boson groups
have creation/removal maps of the fermionic form with denominator sqrt(58),
plus K_empty=sqrt(I-S/58). Thus the complete effect is I on all 1653 input
columns. The no-operation squared coefficients are 1/29,0,28/29 respectively.
The quantum reference has dimension 113 including its blank state.

The program enumerates every occupation state and all 15 energy projectors:
vacuum, four single-frequency sectors and ten unordered two-frequency sectors.
Their ranks sum to 1653. It retains all 25 pointer outcomes crossed with all
15 energy outcomes, giving 375 records. Vacuum preparation gives twelve
positive creation outcomes, with probability dim(g)/58, and a no-operation
outcome of probability 1/29. The other 362 joint records have exact zero
probability. Tests reconstruct every creation and removal energy probability
directly from the actual sparse source operators and energy projectors.

## References, clock and neutrino flavour records

The reference carries the conjugate color/electric representation for a
creation, and the original representation for a removal. Fermionic reference
states have odd parity, so the joint maps have even parity. Color components
remain quantum reference indices rather than classically measured color
labels. These maps preserve the global color and electric symmetries of the
frozen Higgs background. They do not assert an unbroken electroweak symmetry
of that selected background or an interacting local Gauss implementation.

The blank reference has energy zero and each nonblank reference state costs
one GeV. The clock is exp[-it(H_source+H_reference)]. Source energies in the
receipt subtract this declared reference cost from the joint transition.
The maps are declared finite Stinespring probes whose completeness guarantees
a unitary extension on the same finite source and apparatus; external work
allows particle creation. An autonomous energy-conserving detector law is
not derived. The product source state is the two quadratic vacua.

A separate neutral probe uses twelve mass-mode annihilators, two helicities
at two momenta for three fields. Define f_alpha,nu=sum_i U_alpha,i a_i,h=-
and f_alpha,antinu=sum_i conjugate(U_alpha,i) a_i,h=+. Exact PMNS unitarity
gives canonical analyzers. Creation and removal divided by sqrt(12) form a
complete instrument. The receipt retains every mass eigenvalue with original
creation probability |U_alpha,i|²/12, all removal zeros, and a zero outside
projector for each pointer. The 96 records include the antitau-neutrino
analyzer. These are opposite-helicity flavour analyzers of three Majorana
fields, not six independent fields or sterile particles. Their weights are
normalized mass-mode overlaps, without a weak production-rate claim.

## Particle dictionary and boundaries

The 17 elementary rows group the nine charged fermions, three neutral mass
fields and five boson types, counting the W pair and eight gluon components
within their types. Their values are computed from this action's matrices.
Every row reports its mass convention, charge/antiparticle multiplicity,
kinetic normalization and missing physical attachment. In particular,
running quark coefficients are not asymptotic quark poles, tree bosonic
coefficients are not measured complex poles, and the charged lepton compiler
uses a distinct source and scale. The receipt marks all physical native mass
predictions as null. Recovering a calibrated input is no independent prediction.

Separate comparison rows link the quenched pion triplet, nucleon isospin
partners and rho anchor. Their source, empirical tuning, alternative scale
and missing physical transfer are retained, without merging them into the
elementary quadratic Fock space. Additional sterile architectures, compiler
mirror valleys, neutral composites and scalar or polarization branches are
classified separately. No completeness claim for every composite or every
possible source architecture follows from the 17-type elementary inventory.
