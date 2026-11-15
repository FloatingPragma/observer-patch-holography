# Independent full-criterion Maxwell review

Reviewer: Codex worker `kg_complete`, distinct from the ROW29 producer DYN21.
The review owns only this directory and makes no producer or official-row edits.

The candidate is `native-local12-dirac-noncompact-quadratic-v4-proposed`,
ROW29 manifest `993744d14c2da8d66a93dc76422f81799074febad9914ae33bb60bf455cf00b6`.
Its size-three native member is frozen by manifest
`90178885fcaaaba405b950a6f132fd47d519d21a130a4f982ec09423db129ec2`.
The criterion is exactly **“Maxwell recovered; transverse waves at c”**,
under criterion hash `0c5d2a8e4212d05f5cfde6f163ca3e3c9cf32d09e97f8c23c347e5306daf980f`.

The review assesses the whole criterion, not just a dispersion coefficient.
The final signed decision and the results it requires are in `REVIEW.json`.
One review does not itself authorize an official OpenWave status change or
stand in for the required second producer-independent review.

## Construction and source boundary

This is the admitted **free noncompact quadratic field** on the actual LOCAL12
oriented face complex. The positive coefficients a=1/2 and b=3, spectral-code
preparation, source time step h, h/8 native clock phase, and terminal scheduler
are declared construction inputs. The coefficients were selected to match the
admitted Dirac convention's speed one. Their provenance is explicit; this is
neither a parameter-free prediction of c nor an empirical electromagnetic
coupling measurement. The Maxwell curl, field normalization and transverse
polarizations are derived from the finite face geometry; a continuum Maxwell
stencil is not substituted for that geometry.

**This pass is source-specific.** The scalar Klein–Gordon construction reviewed
elsewhere has speed 1/sqrt(48) in its own convention. Nothing in this review
identifies its dynamics, clock conversion or mass coin with this v4 field or
the routed Dirac source. A shared clock attachment primitive does not prove
equality of the attached source dynamics or their physical speed. The requested
unification remains a separate transfer problem.

The observed domain is free vacuum radiative dynamics and field first moments,
including arbitrary finite sums of bounded nonzero Fourier momenta. It does not
include interacting charged radiation, a compact-cosine transfer, arbitrary
continuous-field spectral measurements, or an ultraviolet limit. None is an
additional requirement in this literal row criterion. The supplied architecture
is not described as assumption-free or unique.

## Independent oriented geometry and field recovery

`review.py` rebuilds A(k), B(k) and their translated integer entries from the
actual size-three cell complex. It constructs the pseudoinverse through an
independent rank factorization and the harmonic projector through a nullspace
basis. It checks B(k)A(k)=0 as an exact Laurent polynomial, not at selected
momenta. The zero-symbol ranks are 22 and 11; the inverse squared Frobenius
norm is 3383/240 and the face squared Frobenius norm is 88.

The signed macro-face columns independently give V. All three exact identities
P D_i Q = W cross(e_i) V^T hold with W^T W=(2/3)I and V^T V=I. The cross-product
polynomial has exactly two transverse positive singular values for every
nonzero direction. In-cell connectivity forces any left-null vector of B(k)
to be constant; a nontrivial intercell phase removes that constant. Therefore
rank B(k)=12 away from zero, rank A(k)<=24, and the 22 separated optical plus
two positive low modes exhaust the positive spectrum. There is no hidden
third longitudinal propagating mode. Harmonic zero modes remain separately.

For A_phys=W^T theta/h, E_phys=−a W^T E/h² and
B_phys=V^T A(k)*theta/h², the microscopic Hamilton equations give
E_phys=−d_t A_phys exactly. On the low sector,
A(k)*=i V cross(k) W^T+O(k²), so the signs give B=curl A,
d_t B=−curl E and d_t E=ab(2/3)curl B in the controlled limit.
The last coefficient is exactly one. Both divergence constraints follow from
cross(k)k=0 and the low transverse domain. This checks the oriented E/B
dictionary, not just the squared wave operator.

The translation list independently bounds the linear Frobenius coefficient
by sqrt(53+22sqrt(2))<11 and the quadratic coefficient by sqrt(117)/2<17/2.
Reconstruction of the absolute curvature convolution gives positive weights
w with M w <=(1504/17)w, establishing the uniform curvature bound used in the
source certificate at every periodic refinement.

## Continuum error and field normalization

The block-rotation proof is valid: left and right skew eliminators each have
norm at most 4 eta, eta=11|k|. The first-order action on the perturbation costs
8 eta², while the second derivative of the two-sided unitary conjugation is
at most (8 eta)²(10+eta). Its integral remainder is half that quantity.
Adding the translated-face Taylor error gives the producer's
D=8 eta²+32 eta²(10+eta)+(17/2)|k|².

The low Hermitian-dilation cluster is separated by at least (4/5)|k|. At every
reported finite refinement this review separately checks the high-cluster
gap, in addition to all advertised admission inequalities. Projection onto
the comparison complement bounds the unwanted component by
D/((4/5)|k|−D). Two columns, polar alignment and extraction of the paired
left/right halves cost at most the stated factor four. Restoring rotations
is covered by 32 eta. Thus zeta and dw=(5/4)D/h are genuine outward bounds,
including split polarizations whose frequency matrix need not be diagonal.

The Fourier normalization h^(3/2), canonical oscillator amplitudes and W,V
normalizations give A amplitude sqrt((1/3)/(2|p|)) and E/B amplitude
sqrt((1/3)|p|/2). The conservative relative dictionary bound
64(zeta+dw/|p|) controls orientation, frequency square roots and inverse square
roots together. It tends to zero on every fixed bounded nonzero momentum set;
weighted finite sums therefore give controlled radiative wave packets. The
fixed native precision floor is retained rather than asserted to vanish.

Every rational bound at L=2^34,2^36,2^38 is recomputed independently. These
enormous lattices were **not allocated or run**. They are finite members covered
by the symbol theorem, whose observed terminal spectral instruments were run.
The complete-record TV bounds are approximately 5.70437e−5, 1.42841e−5 and
3.59419e−6. Field-mean errors add twice the TV bound and the dictionary bound.

## Preparation, cap, spectators and phase reference

The code |vac>,|y>,|z> maps to 00,01,10; 11 is retained padding. A finite native
reset and Hadamard prepare vacuum plus one selected photon. This is the
declared spectral encoding, not an inferred encoder from an old edge-register
state. The cap is imposed on the full free Hamiltonian, retaining its visible
zero-point scalar. It commutes with that Hamiltonian. Consequently compression
preserves the linear Heisenberg equations even though it does not preserve
canonical commutators. The review proves [H,a_j]=−sum_l Omega_jl a_l for a
general symbolic Hermitian two-mode Omega and rejects squaring capped field
matrices in place of compressing H. The cap CCR defect is explicitly two.

The autonomous spectator reduction is exact for trace-class joint states,
including correlations with an external reference, subnormalization and each
local instrument branch. Cyclicity of the partial trace removes U_s before
any comparison. Fresh correlated numerical tests verify that identity. The
finite-box coexact oscillator vacuum and three normalized harmonic Gaussian
references are admissible; no harmonic free-particle ground state is assumed.
Spectator vacuum phases cancel by this theorem, not by unbounded vacuum-energy
subtraction. Spectator field fluctuations remain. Reported radiative first
moments equal full-field means only when spectator first moments vanish.

The comparison diag(p,2p,2p,0) preserves the observed zero-point phase. Its
physical-code Hamiltonian distance is at most 2dw: one dw from half the trace
and one from the one-photon block. Padding has the same declared action. An
independent native all-column diagonal test retains this scalar without phase
alignment. Deleting it changes a coherent vacuum/padding reference; it cannot
be silently deleted even though some vacuum/one-photon sign laws are unchanged.

## Native records, clock and opposing controls

The native executor is implemented here from the raw six-state response
matrices and control coefficients. It does not call the producer or compiler
executor. Each pulse applies the success compression and accumulates its
orthogonal export complement after all preceding successes. The summed Gram
identity is checked on **all eight input columns** for all 32 radiative words.
The six-pulse native pointer is also reconstructed with its literal coin and
current transition records. Every tensor-product fine effect is contracted:
125 fine C outcomes and 512 fine E outcomes for each setting, without dividing
by success. Both fine laws reproduce the stored eight-outcome records.

The exact primitive angle/control certificate is reconsumed, and the bound
adds three blank resets, three READY resets and three pointers. The clock
certificate already includes one pointer; the additional two are counted
once. Source approximation adds 2Tdw. Reset timeouts and first export failures
are outcomes, not discarded runs. Floating Gram residuals are diagnostics;
the outward bounds supply the mathematical instrument certificate.

The stored eight-outcome arrays are the initialized source/pointer component
law, before the reset-success weights. For a specified native input the
complete instrument composes the rank-one reset maps, multiplies the later
success records by their actual success weights, and retains each first reset
timeout. Appending reset-failure numbers to the unweighted component array
would be wrong. The uniform lower reset-success effects bound this complete
composition. The review's explicit first-pointer-failure effects also telescope
with the full successful tensor effect to the identity.

The projected bit X/Y operator equals the specified compressed field matrix.
The other polarization has zero first moment and a fair sign law. Analysis
can populate padding, which remains in the output record. This realizes
field first moments, not an unsupported continuous-field distribution.

The independent replay checks both field and clock contrasts. Terminal
coalescing of free evolution preserves the exact source and clock group law;
the emitted word has a newly priced apparatus bound. No per-macro failure
tape or intermediate native observations are inferred from it. The split
clock control agrees, and doubling the actual clock tick is detected while
the field word is held fixed. Doubling b is separately detected. Reversing
an actual face entry breaks the complex. Native source phase deletion and
naive capped-Hamiltonian controls also reject the corresponding false claims.

## Full size-three native member

`size3_check.py` independently recomputes the complete 646-mode frame Gram and
residual certificate using **four centered 15-bit limbs**, distinct from the
producer's three 18-bit limbs. Exact integer GEMM has absolute dot-product
sums below 2^53; Python integers recombine all terms. Both squared rational
Frobenius values match exactly. The resulting cap-one Hamiltonian bound,
including its 646-mode zero-point scalar, is 1.152700483032225e−6.
No eigenvectors are reselected or silently rotated inside degeneracies.

The whole 1024-target, 25606-pulse diagonal certificate is consumed. Fresh raw
native execution covers vacuum, first and last photon, and two padding columns.
A complete prepared field/clock word is executed again and compared across
all 2048 outcomes of both pointer laws. The other settings use the same frozen
diagonal and the local X/Y analysis already independently checked above.
The visible full vacuum phase is retained. The size-three member supplies a
literal finite source execution alongside the analytically controlled family.

The initial auxiliary size-three review attempt completed the independent
frame proof, then stopped on a harness return-arity error before native replay.
The error was corrected without producer changes; the attempt log is retained
and the complete corrected command rerun. It is not a scientific rejection.

No extra interacting-pole, empirical-alpha or uniqueness condition is imposed.
No scalar/Dirac source equality, inter-row common-speed unification, absolute
vacuum-energy measurement, or official OpenWave update follows from this review.
