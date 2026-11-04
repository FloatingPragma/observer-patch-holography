# Observable charged poles of finite interacting cores

The Hamiltonian spectral measures of two complete finite sources have
isolated charged poles with a positive, quantitatively controlled residue.
One source is the supplied linear texture of
[the defect-core proof](../topological_defect_spectrum/CORE_PROOF.md).
The other is the original massless Weyl hopping on an odd open cube, with
no imposed mass field, constructed in
[the odd-box proof](../native_odd_box_core/PROOF.md).
Both retain the original centered onsite contact and every Fock state.
Their backgrounds and zero modes differ; their energy levels are not
identified across the two sources.

## Exact source inputs and spectral subspaces

The three-site-per-axis cubes have 54 one-body modes, 26 negative levels,
two zero modes and 26 positive levels. Filling the negative levels and
subtracting only their fixed energy gives

~~~
H0 >= Delta (1-P0),       rank(P0)=4.
W=U sum_x (n_xup-1/2)(n_xdown-1/2),   ||W||=w=27U/4.
~~~

The four zero-energy Fock states have particle numbers 26,27,27,28.
The parents prove their exact contact compression, full-complement
spectral bounds, particle-hole symmetry and an appropriate Kramers
symmetry. The linear texture has Delta>6/5 and a real zero orbital
proportional to (1,2,1) in each coordinate in its displayed phase frame.
The unmodified source at speed one has Delta>2/3 and zero orbital
proportional to (1,0,1) in each coordinate. Both have inverse participation
I=1/8. These rational lower bounds, not floating gaps, enter the estimates.

Write d for either proved gap lower bound. Suppose

~~~
w<d/4,    epsilon=w^2/(d-2w),    U I/2>2epsilon.       (1)
~~~

The parents give the exact interacting ground projector G in N=27,
of rank two, and one lowest charged state in each of N=26 and N=28,
with projectors C_- and C_+. Their excitation energies are equal:

~~~
g in [U I/2-2epsilon, U I/2+2epsilon],    g>0.         (2)
~~~

Every other state in these charged sectors has excitation energy at
least d-2w. Let G0,C_-0,C_+0 denote the corresponding noninteracting
projectors, and put rho=G/2 and rho0=G0/2. No charged sector or virtual
state is removed from the Hamiltonian.

## A bound on the actual interacting projectors

Set t=w/(d-2w). In any of the three low sectors an eigenvector psi of
H0+W has eigenvalue lambda in [-w,w]. If P is its unperturbed zero
projection and Q=1-P within that particle-number sector, then

~~~
Q psi = -[Q(H0+W-lambda)Q]^-1 Q W P psi,
||Q psi|| <= t ||P psi|| <= t ||psi||.                (3)
~~~

The inverse exists since its operator is at least d-2w. The interacting
ground is an exact degenerate doublet, so (3) holds for every vector in
its two-dimensional space, without summing separate eigenvector errors.
The charged spaces are one-dimensional. Equal finite ranks and the
principal-angle identity for orthogonal projections give

~~~
||G-G0|| <= t,     ||C_+-C_+0||, ||C_--C_-0|| <= t,
||rho-rho0||_1 <= 2t.                                (4)
~~~

To obtain the projection identity, decompose the pair into its
two-dimensional principal-angle blocks; both relevant norms are the
sine of the same angle. The trace-norm estimate uses
rank(G-G0)<=4 and division by the ground multiplicity two.
These are full-Fock estimates. The number 2^54 is not a four-state
execution claim.

## Positive pole measures and their residues

For a one-body smearing f, let c(f)=sum_i conjugate(f_i)c_i.
The exact CAR gives ||c(f)||=||f||. The positive addition and removal
measures are

~~~
mu_+(f,A)=Tr rho c(f) 1_A(H-E_ground) c(f)^dagger,
mu_-(f,A)=Tr rho c(f)^dagger 1_A(H-E_ground) c(f).
~~~

Only N=28 or N=26 contributes to the respective target. The atoms at g
have one-body residue matrices

~~~
(M_+)_ij=Tr rho c_i C_+ c_j^dagger,
(M_-)_ij=Tr rho c_j^dagger C_- c_i.                   (5)
~~~

Both quadratic forms f^dagger M f equal the displayed nonnegative pole
weights. Let Z be the complete rank-two one-body zero-mode projector.
Wick contraction of the filled sea and the equal mixture of the singly
occupied zero modes gives

~~~
M_+0=M_-0=Z/2.                                       (6)
~~~

Changing rho costs at most 2t in a normalized quadratic form; changing
C costs at most t. On the entire 54-dimensional one-body space,

~~~
||M_+-Z/2||, ||M_--Z/2|| <= 3t.                      (7)
~~~

Each matrix is the Gram matrix of transitions from two ground vectors
to one charged vector, so its rank is at most two. If 3t<1/2, (7)
makes its compression to ran(Z) positive definite, so the rank is
exactly two. For both retained cubes, (1) gives t<1/216 because
U I/2=w/108 and epsilon=w t. Thus 3t<1/2 follows from the proved
ground-ordering condition and is no extra numerical premise.
The source antiunitary leaves the ground mixture and each
charged projector invariant and squares to -1 on the one-body carrier.
It commutes with these Hermitian residue matrices. Their nonzero
eigenvalues therefore occur in Kramers pairs. This is a degeneracy of
a matrix measure, not an assignment of two observed particle species.

Choose zero-mode smearings f1,f2 in the source's Kramers basis. The
antiunitary preserves their span, so the compressed residues are
z_+ I_2 and z_- I_2. Particle-hole conjugation exchanges addition and
removal and preserves that span, giving z_+=z_-=z. The same symmetries
and CAR give total addition and removal weight 1/2 for each normalized
zero-mode smearing. Consequently

~~~
1/2-3t <= z <= 1/2,
mu_+(f_sigma,{g})=mu_-(f_sigma,{g})=z,
mu_+(f_sigma,[d-2w,infinity)),
mu_-(f_sigma,[d-2w,infinity)) <= 3t.                  (8)
~~~

The remaining positive weight is retained above its proved threshold.
For every real Hamiltonian time s this implies

~~~
Tr rho {c(f_sigma,s),c(f_sigma)^dagger}
       =2z cos(g s)+R(s),       |R(s)|<=6t.           (9)
~~~

For a native site and spin, (7) bounds the pole weight by
rho_zero(x)/2 +/-3t. The textured cube's center has rho_zero=8/27,
giving unperturbed weight 4/27. Each of the unmodified cube's eight
corners has rho_zero=1/8, giving weight 1/16. At U=1/65536 both lower
bounds are strictly positive by rational arithmetic. The unmodified
center has zero unperturbed weight and is retained as a dark-detector
control. The full residue's two eigenvalues and its compression z are
different quantities; no equality between them is assumed.

## Number, energy and complete repair records

Both sources satisfy [H,N]=0 and [N,c(f)^dagger]=c(f)^dagger.
The two pole measures carry relative CAR charges +1 and -1 and flip
fermion parity. This fixes neither electromagnetic units nor a
Lorentz-spin or exchange-statistics interpretation.

For an oriented source bond define

~~~
J_xy=i(c_x^dagger h_xy c_y-c_y^dagger h_yx c_x).
~~~

The CAR and the onsite mass/contact terms give
d n_x/ds=i[H,n_x]=-sum_y J_xy. Energy transport includes the interaction.
Assign each bond half to each endpoint:

~~~
e_x=c_x^dagger B_x c_x+U(n_xup-1/2)(n_xdown-1/2)
    +(1/2)sum_y(c_x^dagger h_xy c_y+c_y^dagger h_yx c_x),
J^E_xy=i[e_x,e_y].
~~~

Distribute the fixed filled-sea reference subtraction as scalar constants
among the e_x; it changes no current. Then H=sum_x e_x,
J^E_xy=-J^E_yx and
d e_x/ds=-sum_y J^E_xy. Every overlapping pair is retained, including
distance-two supports and interaction terms. A regional sum leaves
exactly the crossing energy currents. This supplies finite Hamiltonian
energy continuity, without identifying a Lorentz stress tensor.

The complete native C/E instrument has 32 coin/raw slots per seam.
For each branch K_r on x,y its recorded current j_r obeys

~~~
N_x K_r-K_r N_x=-j_r K_r,
N_y K_r-K_r N_y=+j_r K_r,     [N,K_r]=0.              (10)
~~~

The collision-hold and padded slots are included. Equation (10) follows
on every occupation input from the actual branch map and extends to
arbitrary quantum inputs and untouched references. In a conditional
branch it compares the output to the selected input subensemble;
Bayesian measurement updating is not an additional transport current.
Summing all branches gives unconditional charge balance. Repair can
change Hamiltonian energy; its exact mean increment is
Tr rho[sum_r K_r^dagger H K_r-H]. Hamiltonian energy conservation is
not transferred to that increment.

## Actual event clocks and the limit of spectral protection

Use the admitted independent Poisson seam scheduler, each of the 54
cube seams at rate Gamma, interleaved with the Hamiltonian. At duration
T, no event has probability p=exp(-54 Gamma T). The output retains
every event time, coin/raw record and quantum state. Its empty-record
block is exactly p times the Hamiltonian output. Its trace distance
from the Hamiltonian output with a sure empty archive is therefore
exactly 1-p, including any untouched reference. In this section trace
distance means (1/2)||rho-sigma||_1; equation (4) states the unhalved
trace-norm estimate.

Starting in rho0 instead of rho gives by (4) the comparison bound
1-p+p t for the uncorrelated ground input. Marginal ground closeness
does not remove arbitrary initial correlations with a spectator.
Conditioning on the empty archive gives the Hamiltonian law, but its
success probability p and all other outputs remain in the protocol.
The receipt stores the exact rational value of -log(p), so very small
probabilities are not rounded to zero.

There is a concrete positive-rate continuing window. Set
T=1/g_lower and Gamma=g_lower/54000, using the positive lower bound in
(2). Every seam operates at a strictly positive fixed rate, and
54 Gamma T=1/1000. The complete-output distance is
1-exp(-1/1000)<1/1000 for the exact ground input, or at most
1/1000+t for rho0. These are bounds on the whole archive, without
postselection. This rate is an explicitly chosen scheduler parameter,
not a derived universal rate or a reinterpretation of the unit-rate
control. At unit seam rate the same duration has its stated much smaller
no-event probability. Neither result executes a spectral estimator.

This derives the finite no-event comparison and exposes its cost.
It does not prove repair preservation of a pole, ground or inversion
index. The admitted repair moves a hole next to a singly occupied site
and holds a doubly occupied site. The Hamiltonian's particle-hole
symmetry is not imposed on the recorded instrument.

## Physical and operational domain

The two-point functions are even CAR quantities on the actual finite
sources. The proof gives positive charged spectral poles, detector
weights, exact currents and the complete-archive comparison. It does
not execute a long-time spectral estimator or a 2^54 matrix. A compiled
mode decoder alone is not such an estimator. Preparation, services and
clock resources must be specified by any implementing protocol.

The textured source retains its supplied mass field. The unmodified
source has boundary standing-wave zero modes; a finite-box interaction
gap is not a thermodynamic particle mass. Neither branch identifies a
measured species, calibrates physical units, proves an interacting
relativistic continuum, or supplies Einstein's stress/entropy premises.
