# NATIVE QUANTUM LINK: AN EXACT SPECTRAL PROJECTOR WITH NONINTEGER FOCK LIFT

Scope

This example uses the source's literal 60-state link and its e^c raising
operator. Its two-site scalar onsite energies are 0 and 64, with hopping
1/2. It tests the proposed one-particle-to-Fock inference in the native
operator algebra. It does not identify the spectral projector of the
343-site compact texture, prove the failure of every alternative quantum
link construction, or compute a physical particle mass.

Construction

Let |0> be the link vacuum, |e> its e^c link state, and

 U=|e><0|, P_e=|e><e|, P_0=|0><0|,
 H=[[0,U/2],[U^dagger/2,64I]].

This is a 120-dimensional one-particle matter/link Hamiltonian. Only
|x,e> and |y,0> mix. The coupled block is [[0,1/2],[1/2,64]], with
eigenvalues (64 plus or minus sqrt(4097))/2. The other 59 x states have
energy 0; the other 59 y states have energy 64.

Set R=sqrt(4097) and t=(1-64/R)/2. The spectral projector onto the 60
lowest eigenvalues is exactly

 P_xx=I-t P_e,
 P_yy=t P_0,
 P_xy=-U/(2R), P_yx=P_xy^dagger.

The identities P^2=P, P^dagger=P and [P,H]=0 hold over the algebraic
field Q(sqrt(4097)). The low-to-high spectral gap is 64. In particular,
this is a genuine isolated spectral projector rather than an arbitrary
choice of invariant subspace.

Fock lift

For two canonical matter annihilators c_x,c_y, define

 N(P)=sum_ij c_i^dagger P_ij c_j.

The complete Fock/link dimension is 4 times 60 = 240. In the matter
vacuum sector N(P)=0. In the one-particle sector it equals P. In the
two-occupied sector it equals

 P_xx+P_yy=I+t(P_0-P_e).

Consequently its exact spectrum is

 0, multiplicity 120;
 1, multiplicity 118;
 1+t, multiplicity 1;
 1-t, multiplicity 1.

The value t is strictly between 0 and 1. The last two values are not
integers. Numerically t is approximately 0.000061024 for orientation;
no floating approximation is used to prove nonintegrality.

Gauge invariance

The e^c matter representation has primitive charge 6 and is a singlet
under the source's nonabelian factors. Its link has opposite charge at
x and matching charge at y. The exact native current matrices obey

 [L_x(T),U]=-T_e U,
 [L_y(T),U]=+T_e U.

These identities prove the one-particle Gauss invariance of H and P and
the full Fock/link Gauss invariance of N(P).

Place one static dual e^c reference at each endpoint. The state with
both matter sites occupied and a blank link has zero total charge for
each endpoint and every native current. It is therefore a Gauss singlet.
Its N(P) eigenvalue is 1+t. Gauge invariance does not restore integer
occupation for this operator-valued spectral projector.

The executable uses all 12 native generators and verifies 72 exact
endpoint commutator and singlet identities, as well as the spectral
projector and Fock-sector identities. This supplies a concrete native
link witness for the algebraic boundary described in
QUANTUM_LINK_PROOF.md. It leaves open a different many-body construction
that explicitly establishes canonical fermion coordinates.

Reproduction

 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_gauge_band.native_link_projector_witness

The JSON receipt binds the source, link-construction dependency, exact
spectra, and scope. The program performs exact symbolic identities.
