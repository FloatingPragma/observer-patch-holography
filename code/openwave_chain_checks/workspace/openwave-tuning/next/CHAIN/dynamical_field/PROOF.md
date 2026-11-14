# Interacting integer-field cutoff in the C10 clock

The finite member at modulus six is the actual 384-dimensional Hamiltonian
in `C10/joint_source.py`. The executable regenerates that member, the global
interacting reference, and all four charged-reservoir outcomes. It compares
the Hamiltonian and every complete MIX Kraus matrix entry with the source.
Both maximum differences are zero in the executed double-precision arrays.
The source Hamiltonian hash agrees with the common-clock contract.
The complete density and each reservoir branch are compared directly
with the common clock's archived input, alongside its duration and rates.

## One declared family and preparation

There are six CAR orbitals, six half-edges, and one cycle. Let
`G=(1,1,1,0,0,0)`, `B=3-popcount(mask)`, and let `D` have a positive entry at
each half-edge port and a negative entry at its transit. The infinite
comparison basis consists of every pair `(mask,F)` satisfying

    mask in {0,...,63},  F in Z^6,  n(mask)+D F+B delta_0=G.

This is a mathematical comparison carrier, not a physical infinite register.
Use the source's exact spanning tree with chord 3. Assigning integer `k=F_3`
and eliminating the five tree edges gives a unique record `F(mask,k)`.
Consequently this carrier is `C^64 tensor l2(Z)`.

For every integer modulus M>=6, the finite member has all 64M labels and
reduces these records modulo M. Its electric energy uses the original
componentwise representatives `r_M(x)=(x+floor(M/2)) mod M-floor(M/2)`.
Its mass, half-hopping phases, electric coefficient .37, magnetic
coefficient .61, and complete C/E/MIX instruments have the original laws.
These are declared architectural choices. Increasing M changes register
capacity; it changes no coefficient or time unit.

The isometry `J_M` sends each finite label to its integer tree solution
with `k=r_M(F_3)`. It need not send every component of F to a centered
representative. This distinction is essential for the electric energy.
The existing spatial readout uses a different tree with chord 5. On all
384 labels the two integer charts differ by a multiple-six cycle, and
their divergence and residues agree exactly. Thus their spatial face
degrees and the augmented boundary degree `-B=N-3` agree. Their electric
energies are not identified. No real-valued minimization is used.

Let `rho_6` be the actual global interacting ground reference followed by
the complete fair-0/1 charged-reservoir instrument, with its four labels
retained. The preparation of member M is its isometric lift using the
same `(mask,k)` labels. The infinite preparation is `J_6 rho_6 J_6*`,
with the same retained reference. This is a supplied preparation family,
not the ground of member M. Every one of the 384 original labels obeys
`max_h |F_h|<=5`. The proof uses this exact integer support bound, rather
than discarding small numerical amplitudes. The numerical ground has
multiplicity one, so four subnormalized vectors execute all reservoir
branches; no outcome is selected or renormalized.

## The interacting limit operator

On the finite integer-record core define

    H_infinity = H_mass + (kappa/2) sum_h F_h^2
                 + g(2I-W-W*) + H_hop,

where W adds `(1,-1,1,-1,1,-1)` to F, and H_hop is the original CAR
half-hop sum with the matching integer record increment. Every term
preserves the displayed integer Gauss constraint and the existing B.
Define `D0=H_mass+(kappa/2) sum F_h^2+2gI` and
`V=H_hop-g(W+W*)`. Each row of H_hop has at most six entries of absolute
value one. The same holds for columns, so Schur's inequality gives

    ||H_hop||<=6,   ||V||<=beta=6+2g=361/50.

D0 is a real diagonal self-adjoint operator on its maximal domain.
The bounded perturbation theorem defines H_infinity on `D(D0)`. Its
clock is exactly `U(t)=exp(-it H_infinity)`. Every finite member uses
`exp(-it H_M)` in that same source unit. Finite-particle oscillator,
quadratic Maxwell, relativistic, or physical-second identifications
are absent from this operator definition.

This field interacts with matter. The actual M6 commutator between field
and matter energies is nonzero. Its Frobenius norm is approximately
21.4472, and the executed reservoir preparation exchanges field energy
under the full autonomous Hamiltonian. No prescribed expectation current
or postprocessed Green field supplies this dynamics.

## Wrap bound for passive evolution

Put `P_r=1_(max_h |F_h|<=r)` and `L_M=floor((M-1)/2)`. One application
of V changes each integer half-record by at most one. Diagonal free
evolution changes no support. Every Dyson term of order n starting
inside P_r therefore lies inside P_(r+n).

If `r+n<=L_M`, no component wraps and no cyclic identification occurs.
The free electric phase, CAR phase, and every target label are then
identical in the finite and infinite interaction-picture words. This
uses all six record components, not only the cycle coordinate.

Let `d=L_M-r+1`, `x=beta |t|`. Each omitted Dyson series has terms
bounded by `x^n/n!`. Therefore, for any normalized vector in P_r,
including an arbitrary retained quantum reference,

    ||J_M U_M(t) psi-U_infinity(t) J_M psi||
       <= min(2, 2 sum_(n>=d) x^n/n!).                 (1)

For d<=0 the right side is defined as two. The half-trace distance
and every complete terminal outcome law have error at most the same
quantity, capped at one. No assertion of small full-carrier
operator-norm commutator error occurs.

For an arbitrary initial vector with weight p outside P_r, splitting it
into its two orthogonal components gives the sufficient norm allowance
`delta_good+2 sqrt(p)`. Mixed states follow by purification. The displayed
small-tail probabilities of the actual preparation are double-precision
diagnostics. No theorem in this package needs them as certified bounds.

## All native events in the common marked clock

The common-clock law has native attempts at rate `lambda=4/5`, duration
`T=37/100`, and independent idle opportunities at rate `4/5`. Each
native attempt executes the full MIX instrument, including zero-current
outcomes. Its isometry has norm one and moves every half-record by at
most one. Its definition extends to the integer carrier by the exact
same targets, CAR phases, and marks. Completeness follows from the
same injective fixed-current branches and fair law/coin weights.

Condition on N attempt times and retain every quantum instrument result.
Expand only the coherent intervals in V. Summing the time-simplex
volumes over all allocations to the intervals gives `T^n/n!`. Through
order `n+N<=L_M-r` all complete words intertwine. Thus (1) holds with
`d=L_M-r-N+1`, independently of the times. The probability of N is
`exp(-lambda T)(lambda T)^N/N!`, including every no-current attempt.
All count classes are retained. Convexity and the binomial convolution
of the two factorial series give

    TV(full marked law_M, full marked law_infinity)
      <= min(1, 2 exp(-lambda T)
                     sum_(j>=L_M-r+1) [(beta+lambda)T]^j/j!)
      <= min(1, 2 sum_(j>=L_M-r+1) y^j/j!),           (2)
    y=(beta+lambda)T=14837/5000.

This compares the complete marked histories and their retained terminal
source/reference states, after embedding the finite state. No Poisson
tail is dropped. Independent idle marks have identical laws, do not
change support, and preserve total coherent duration T; they do not
increase y. The proof covers arbitrary terminal measurements on the
common embedded carrier. A further native terminal MIX isometry has
its own one-step record increment, so replace r by r+1 in (2).

For `y<d+1`, the explicit rational majorant is

    2 sum_(j>=d) y^j/j!
       <= 2 y^d / [d! (1-y/(d+1))].                 (3)

The receipt stores the exact rational upper bound and its outward
floating display. At M=42, r=5 and d=16, (2)-(3) give
`TV<4.185417e-6`. Including an additional terminal MIX event gives
the separately recorded bound. For fixed T, lambda and r, d tends
to infinity with M, so (3) proves convergence of complete interacting
record laws. This limit preserves the original clock and preparation
ancestry. It is a field-register cutoff limit on a fixed triangle.

## Executed controls and mathematical boundaries

The executable retains all 64 occupation outcomes and all 30 terminal
MIX labels at M=6,12,18,30,42. It also executes a complete event at
time .17 followed by .20 coherent time at M=6 and M=42, retaining all
120 combinations of native and reservoir marks. The entire Poisson
sum in (2) is proved analytically; these executions do not numerically
enumerate its infinitely many possible counts.

The passive integer comparison uses a Hermitian Dirichlet truncation
`|k|<=26`, dimension 3392. The initial chart has `|k|<=3`. Terms through
order 23 do not encounter its boundary, so the analogous factorial
tail at d=24 bounds its error to the infinite comparison. The deleted
couplings belong only to this diagnostic approximation, and its output
is never normalized after evolution. The declared finite sources are
periodic and retain all their wrap transitions. Floating propagation
errors are diagnostics, not interval certificates. The analytic error
bounds use exact rational arithmetic and do not depend on those errors.

At M6 the integer chart differs from periodic electric energy by as
much as 13.32. The actual preparation's diagnostic weight on nonzero
energy defects is about .000127989, and its mean defect is about
.000907615. The passive vector discrepancy at T=.37 is about .0243425.
The complete-clock support theorem yields only TV<=1 at M6. It does
not promote the much smaller numerical passive discrepancy to a
certified full-clock estimate. The explicit bounded repair is the
larger record modulus with the lifted same preparation; M42 supplies
the nonvacuous bound above without clipping or postselection.

The finite centered coordinate obeys a wrap commutator defect of norm
M against a unit shift. The integer comparison instead obeys
`[F_h,W]=delta_h W` on the finite-support core. Hence this is a
prepared-record convergence theorem, not uniform operator-norm
convergence of the finite coordinate algebra. Removing electric energy
changes the actual M6 evolution by an executed vector distance about
.293685 and is a rejecting mutation.

The existing native Wilson-Maxwell refinement uses completed oscillator
energy, literal finite displacement links, a cube with five transverse
modes, and its own preparation. C10 has a periodic integer electric
record, one graph cycle, and different hopping and preparation laws.
The source-aware bounded-Dyson method transfers; those operator and
spatial conclusions do not. The scalar static three-dimensional field
control also has a different spatial family. This result supplies no
three-dimensional Maxwell, continuum-force, Coulomb-minimization,
particle-spectrum, or physical-clock conclusion. The exact coherent
services retain the common clock's declared, uncompiled status.
