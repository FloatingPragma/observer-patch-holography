# Original spinor-channel records on the selected source

## 1. Same source and actual probe representations

Fix the [primary contract](../particle_common_contract/PROOF.md), including
its interacting Wilson action, actual Gauss and product link cutoffs,
selected polynomial density rho, idle quantum reference and source clock.
At level n write a=2^(-n), L=2/a², V=L³. The density is stationary and
invariant under the actual finite binary-cube action and literal source
charge conjugation. No trial vector replaces this input density.

The [local composite construction](../particle_primary_composite_spin/PROOF.md)
computes the complete six-dimensional weak-doublet-pair times charge-six
singlet kernel. Both positive-beta spin rows are included. In the actual
native mask chart its six active source modes are (6,10,11,21,25,26),
ordered as weak1, weak2, singlet in spin row0, then the same entries in
spin row1. All other CAR modes remain present as spectators.

Let F_(lambda,m) be the normalized creation polynomials obtained by the
actual Dirac spin Casimir projectors: eigenvalue3/4 gives the two rows
of G1 and eigenvalue15/4 gives the four rows of H. The code uses the full
native Gauss-kernel columns, not a guessed assignment of elementary fields.
Every polynomial has degree three and native occupation increment

    Delta=(0,2,0,0,1), in native charge order (-4,-3,1,2,6).

It is invariant under all twelve native gauge generators. Total gauge
charge is zero. Its actual binary-cube transformation is

    U(g) F_(x,lambda,m) U(g)*
      = sum_k D_lambda(g)_(k,m) F_(g x,lambda,k).              (1)

The code restricts all48 exterior matrices and verifies the character
inner product one for each irreducible coefficient representation.
These are creation/probe multiplet labels. The original density is not
assumed a cube-singlet vector, and these labels are not the final state's
total spin or a Lorentz spin assignment.

## 2. Exact local effect and complete apparatus

On all64 occupations of the six active positive-beta modes define

    S_lambda=sum_m F_(lambda,m)* F_(lambda,m),
    N=sum_(i=0)^5 n_i,             d_G1=2, d_H=4.

Exact full-CAR matrices give

    0 <= S_lambda <= d_lambda I,
    S_lambda >= d_lambda (I-N/2),
    S_G1+S_H <= 6I.                                             (2)

The receipt includes the exact spectra of S_lambda and of
S_lambda-d_lambda I+d_lambda N/2. All eigenvalues of the latter are
nonnegative. On the empty active modes S_lambda=d_lambda; on one
active occupation S_lambda=d_lambda/2. At N>=2 the lower bound in
(2) is nonpositive, so positivity suffices. This also gives a short
proof of the lower inequality once the zero- and one-occupation blocks
are computed. The upper bounds are checked on every occupation block.
Each individual normalized F is a contraction, though their sum of
success effects needs the factor six.

For a retained uniform address x use the seven-outcome instrument

    K_(x,lambda,m) = F_(x,lambda,m)/sqrt(6V),
    K_(x,fail) = [I-(S_(x,G1)+S_(x,H))/6]^(1/2)/sqrt(V).        (3)

All six success rows and the failure are retained at every address.
The failure is the positive square root in (3); it is not the projector
failure of a single Slater monomial. Equation (2) proves completeness.

An idle memory supplies the compensating native increment -Delta for
success, and zero for failure. The success memory subspace has six
orthogonal components. Its combined parity with the odd source operator
is even. The failure preserves occupations and parity. The isometry
specified by (3) can be completed to a unitary on finite physical source
and memory blocks; preparation of that memory is a declared apparatus
resource. The gauge-singlet input source remains rho before the probe.
The intervention changes its conditional and unconditional outputs, which
are retained together with memory, pointer and reference correlations.
No autonomous or efficient implementation of this algebraic service is
asserted.

## 3. Original selected-state probability

The [selected odd response](../particle_primary_selected_odd_response/PROOF.md)
proves on this same source that charge conjugation sends positive-beta
occupations in each native block to negative-beta holes. Hence

    <N_L^+ + N_e^+> = <D_L+D_e>/2.

This is an expectation identity for the actual invariant density, not
site independence or an empty-state approximation. Averaging (2) over
all sites gives the original success probability in each probe channel:

    p_lambda = (1/(6V)) sum_(x,m) Tr(rho F_(x,lambda,m)* F_(x,lambda,m))
       >= (d_lambda/6) w,
    w = 1-<D_L+D_e>/(4V).                                    (4)

The inherited exact operator bounds give w>0.84402 at member zero and
w>59/400 throughout the source family. Therefore member zero has
p_G1>0.28134 and p_H>0.56268. These are disjoint events in one complete
apparatus, with no division by their success probability.

## 4. Exact twelve-mode residual certificate

The onsite commutator is [H_on,F]=9F/a. The electric and magnetic terms
commute with each matter-only gauge-invariant F. For any fixed link
coordinates, let T be the one-particle hopping matrix and let a_i, i=0..5,
be the six active orthonormal positive-beta modes at x. Every T a_i has
support at neighboring sites and is therefore orthogonal to all original
a_j. Distinct neighbors are orthogonal. On each incident bond the actual
Wilson and native representation matrices satisfy

    P_+ S* S P_+ = P_+ S S* P_+ = P_+/2,
    R(g)* R(g)=I.

These identities apply jointly to both spin rows and every active weak
or singlet internal component. Consequently

    <T a_i,T a_j> = z_x delta_(i,j)/(2a²),                     (5)

where z_x is the actual number of neighbors. The six normalized hopping
images b_i are orthonormal, including across the weak and singlet blocks,
and orthogonal to the six a_i. They are twelve canonical CAR modes.
Their occupations are unrestricted; hopping images are not assumed empty
or positive-beta modes of their destination sites.

If a monomial F contains three creation factors, let R(F) be the sum
obtained by replacing each one, in its existing order, by its b partner.
On the full link-function carrier,

    [H,F]-9F/a = sqrt(z_x/(2a²)) R(F).

Put G_lambda=sum_m R(F_(lambda,m))* R(F_(lambda,m)). Exact arithmetic
on all4096 twelve-mode occupations certifies

    G_G1 <= (19/3)I,           G_H <= 12I.                    (6)

The certificate first extracts each polynomial's common radical factor;
the resulting Gram numerators are integer matrices, with denominators
six and three respectively. Their exact nonzero-connectivity blocks have
size at most32. Positive rational LDL factorizations of
(19/3)I-G_G1 and12I-G_H certify every block. For G1 there are1922 blocks
and minimum positive numerator pivot7906850/17765029. For H there are1858
blocks, one zero block, and minimum positive numerator pivot6. The code
checks L D L^T equals each matrix and every nonzero pivot is positive.
No numerical eigenvalue is used to certify (6).

The tempting bound G_lambda<=3d_lambda I is false for G1 on this full-CAR
domain. In masks (224,336,392,644,770,1092,1281,2114,2177), the integer
vector (-55,-32,32,-28,8,28,-8,40,-40) has exact Rayleigh quotient
378070/59814>6. This countercontrol concerns the uncompressed twelve-mode
operator, not an eigenvalue in the source's physical Gauss sector.

The physical cutoff P commutes with each complete gauge-invariant F.
Thus the compressed residual is P R(F) P, and

    sum_m (P R_m P)* (P R_m P) <= P G_lambda P.

Pointwise uniform bounds pass to the link direct integral and to this
compression. No finite truncation is assumed to preserve intermediate
hopping-image subspaces.

## 5. Spin-resolved source-energy windows

Use the original stationary transfer measures of the complete instrument:

    mu_lambda(domega)=(1/(6V)) sum_(x,m,i,f)
       rho_i |(F_(x,lambda,m))_(f,i)|² delta_(E_f-E_i)(domega).

Equations (4)-(6) give

    mu_lambda(R)=p_lambda,
    integral (omega-9/a)² dmu_lambda <= C_lambda zbar/(12a²),
    C_G1=19/3, C_H=12,      zbar=2 B/V <=6,                  (7)

where B is the actual number of interior bonds. Member zero has zbar=3.
Chebyshev's inequality therefore gives a positive original measure on
[1/(2a),35/(2a)] in each channel. Uniform lower bounds before clock error
are

    G1: 59/1200 - 38/867,
    H:  59/600  - 24/289.                                    (8)

Use the inherited two-energy-clock construction, with the same H and
energy zero, resolution1/(8a), and each retained tail and clock budget
1/16384. Its complete probability error is at most1/4096. The reported
joint events, including success, native increment and probe-irrep label,
therefore obey

    Pr(G1, reported transfer in [1/(4a),71/(4a)]) > 1/200,
    Pr(H,  reported transfer in [1/(4a),71/(4a)]) > 3/200.      (9)

The clock error is charged once to each stated event. All failure,
out-of-band and other spin-component records remain in the output.
Initial energy recording preserves the stationary source marginal; its
reference disturbance is retained in the complete quantum output.

At member zero the narrower exact window[3,15] has centered moment
bounds19/12 and3. The reported window[11/4,61/4] satisfies

    Pr(G1, reported narrow window)>121/512,
    Pr(H,  reported narrow window)>61/128.                    (10)

The receipt supplies stronger exact rational values. The original
transfer-to-excitation theorem also applies with initial cooling allowance
1/16 and energy enlargement one; its actual bounds are retained separately.
These are source energy differences, not GeV masses or isolated particle
lines. The uniform windows in (9) have energy scale1/a and do not prove
fixed-energy survival as a tends to zero.

## 6. Conserved native-current readout without initial charge measurement

Each actual polynomial obeys [N_r,F]=Delta_r F. The failure in (3)
commutes with all five native numbers. They commute with H as well.
The [coherent transfer-meter identity](../particle_primary_charged_readout/PROOF.md)
therefore applies with five radix-three pointers. Write S for the cyclic
three-state shift and W=product_r controlled(S_r^(N_r)). Then

    W (K tensor I) W* = K tensor product_r S_r^(delta_r).     (11)

The two couplings cancel the unknown initial numbers coherently. The
admitted transfer alphabet is0,Delta,-Delta. Their residues modulo3 are
distinct, so the complete243-state pointer exactly records the native
increment. All inaccessible records have zero probability, not discarded
probability. Additional unproved transfers could alias these residues.
Equation (11) factors the pointer label from every full source/reference
Kraus history and leaves its energy law unchanged. There is no initial
absolute-number measurement, and no electromagnetic charge is inferred
from these gauge-neutral conserved-number increments.

Literal source charge conjugation maps the positive creation multiplets
to the negative-beta hole multiplets with increment -Delta. It preserves
H and rho and intertwines their proper cube representations. Conjugating
the complete apparatus, including the compensating memory and meters,
gives identical original energy-record probabilities. The idle reference
is not acted on; its conjugate transformation belongs to the canonical
comparison of the two complete outputs.

## 7. Component laws, spatial analyzers and transport scope

For any energy event, form the original coefficient Gram matrix by
summing the same local address on its bra and ket and averaging x
uniformly. The actual cube simultaneously permutes those addresses and
acts by D_G1 or D_H on the components. Invariance of rho and of both
energy effects makes this Gram matrix commute with that coefficient
representation. The irreps are inequivalent and occur once. Finite Schur
orthogonality therefore kills cross-irrep blocks and makes each diagonal
block scalar. Each retained component has exactly1/d_lambda of its
channel's original energy-event weight. This conclusion uses the uniform
address average and source symmetry; it does not require the input state
to be a group singlet or equate a probe label with final total spin.
Measuring component rows chooses an apparatus basis. Under a rotation,
that basis rotates coherently; non-Abelian component labels are not
relabelled as if every rotation were a classical permutation.

The success-address register can instead remain coherent until a complete
finite Fourier or Walsh unitary is applied. Unitary mixing of the V Kraus
rows for each fixed component preserves their total effect and complete
channel measure. The member-zero code checks the full eight-row Walsh
matrix exactly. Thus every aggregate lower bound above also holds after
summing all spatial analyzer rows. At least one spatial row is populated;
a prescribed row, a momentum eigenvalue and a dispersion curve are not
selected by this aggregate statement. A classical which-address record
made before this rotation would define a different apparatus.

On the separately prepared anchor, the
[actual six-state transport theorem](../particle_primary_composite_transport/PROOF.md)
gives neighboring third moment -3I_6/4 and vanishing lower neighbor-transfer
moments. This identifies the same polynomials and spin coefficient spaces
in that preparation. It does not transfer those moments or an isolated
branch to the original selected-state measures used here. A controlled
relativistic branch and its physical mass dictionary require spectral
information beyond (7)-(10).

The finite complete-instrument, retained-reference, group-character and
positive-energy-window results are reused on their supplied premises.
No new Lean formalization of the native matrices, relativistic
spin-statistics theorem or physical Standard Model assignment is claimed.
