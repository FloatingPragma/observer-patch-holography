# Exact identity tests for neutral source composites

The source, polynomial density, idle conjugate reference and clock are the
[primary Wilson/Gauss member](../particle_common_contract/PROOF.md). The
calculation uses its complete sixty-mode cell CAR algebra. The
[composite construction](../particle_primary_composite_spin/PROOF.md)
supplies the actual gauge kernels and their rotation matrices. No continuum
particle dictionary, baryon or lepton charge, or mass parameter is inserted.
The receipt binds member zero. The algebraic charge and meter identities
hold on every finite member with these block-preserving source terms.

## 1. Actual conserved currents and conditional constituent names

In native charge order \((-4,-3,1,2,6)\), let \(N_r\) count all four Dirac
rows of block \(r\), over the source. Each native link matrix preserves
these five blocks. Onsite Dirac matrices act within each block; electric and
plaquette terms contain no CAR factors. The link cutoff and Gauss
projection commute with every \(N_r\). Consequently, on the full physical
carrier,

\[
 [H,N_r]=[\rho,N_r]=0.                                      \tag{1}
\]

The second identity follows from the actual selected polynomial
\(\rho=f(H)/\operatorname{Tr}f(H)\), including its physical projection.
It does not require a unique ground state. The twelve exact local gauge
generators annihilate each column of the computed singlet kernels. The
three creation/hole polynomials and their conjugates have increments

| Family | \(\Delta N\) | CAR parity | Local source rotation coefficients |
|---|---|---|---|
| weak pair and charge-six row | \((0,2,0,0,1)\) | odd | \(G_1\oplus H\) |
| color triple | \((1,0,0,2,0)\) | odd | \(G_1\oplus H\) |
| even determinant-hole | \((0,0,0,3,-1)\) | even | \(T_1\oplus E\oplus T_2\) |

Every monomial obeys \([N_r,F]=\Delta N_rF\), checked on the actual
native-mode chart. Opposite increments give the charge-conjugate spaces.
The charge-six, weak and color currents of these composite operators all
vanish. Gauge-neutral total charge alone would not establish the twelve
generator identities; the calculation checks them separately.

Under the *conditional* all-left constituent naming
\((u^c,L,Q,d^c,e^c)\), the source charges

\[
 B_*={-N_{-4}+N_1-N_2\over3},\qquad L_*=N_{-3}-N_6       \tag{2}
\]

assign \((B_*,L_*)=(0,1),(-1,0),(-1,1)\), respectively.
Their opposites reverse these values. Equation (1) proves conservation of
these particular linear combinations. It does not identify them with
measured baryon and lepton numbers. The weak candidate has the conditional
bookkeeping of \(LL e^c\); the color candidate has that of
\(u^c(d^c)^2\). These motivate a neutral-lepton-composite and a neutral
antibaryon-type comparison, subject to an actual field/current dictionary.

There is an exact direct-map test. A current-intertwining map from a gauge
singlet kernel into a nontrivial irreducible weak doublet or color triplet
has zero image: each target generator must annihilate that image, whereas
the target quadratic Casimir is a positive multiple of the identity.
Thus these literal local singlets cannot be directly identified with an
unbroken active weak doublet or a color triplet while preserving those
currents. This excludes no Higgs-dressed dictionary, multipole response,
or emergent interaction with a separately constructed field. The even
candidate cannot represent a single fermion through a CAR-parity-preserving
map. Source cubic labels are not Lorentz spin measurements.

## 2. Normal and anomalous scalar response selection

Set \(\alpha_t(F)=e^{itH}Fe^{-itH}\). Equation (1) and trace cyclicity
give, for any angles \(\theta_r\),

\[
 \begin{aligned}
 \operatorname{Tr}\rho F_i^\dagger\alpha_t(F_j)
 &=e^{i\theta\cdot(\Delta_j-\Delta_i)}
       \operatorname{Tr}\rho F_i^\dagger\alpha_t(F_j),\\
 \operatorname{Tr}\rho F_i\alpha_t(F_j)
 &=e^{i\theta\cdot(\Delta_i+\Delta_j)}
       \operatorname{Tr}\rho F_i\alpha_t(F_j).
 \end{aligned}                                             \tag{3}
\]

Normal kernels therefore vanish between distinct native increments.
Anomalous kernels require opposite increments. These are necessary
conditions, not assertions that an allowed correlator is nonzero. The six
signed candidate increments are pairwise distinct, so their normal
cross-family scalar kernels vanish. Operators within the same charge
family require rotation information too.

For complete uniform spatial-orbit sums, the source cube action and
\(\rho\) are invariant. Their normal kernels are intertwiners. Exact
binary-cube character contraction gives

\[
 \begin{aligned}
 G_1^*\otimes G_1&=A_1\oplus T_1,\\
 G_1^*\otimes H&=E\oplus T_1\oplus T_2,\\
 H^*\otimes H&=A_1\oplus A_2\oplus E\oplus2T_1\oplus2T_2.
 \end{aligned}                                             \tag{4}
\]

The scalar \(G_1/H\) cross-kernel is zero for this orbital construction;
a vector \(T_1\) apparatus can connect them. Arbitrary localized probes
also carry orbital transformation data, so (4) does not assert that every
localized cross-kernel vanishes. The positive odd uniform multiplets have
inversion parity g and the opposite multiplets u. Their charge-allowed
anomalous *uniform scalar* kernel is therefore also killed by inversion.
The even multiplets have parity u on both sides, so this extra exclusion
does not apply. No physical Majorana-mass theorem follows.

Equation (3) is a scalar source-marginal statement. The reference-valued
partial trace need not vanish. For example the purification
\((|00\rangle+|11\rangle)/\sqrt2\) of \(I/2\), with source raising
operator \(|1\rangle\langle0|\), has zero scalar expectation but a
nonzero reference-valued contraction. A compensating memory or a reference
intervention can carry the missing charge. The record algebra and
preparation must specify those interventions explicitly.

## 3. A complete seven-state charge-difference apparatus

Use the actual conserved integer operator

\[
 S=N_{-4}+N_2+N_6.                                        \tag{5}
\]

Its three positive candidate shifts are \(1,3,2\); the opposite
candidates have negative shifts. A cell has 28 contributing modes, so
\(0\le S\le28V\). A blank cyclic pointer of radix seven separates the
specified alphabet \(\{-3,-2,-1,0,1,2,3\}\). No smaller pointer with
ordinary cyclic addition separates these seven values. This is not a
readout of an arbitrary total occupation or arbitrary collection of
excitations.

Let \(F_{x,r,j}\) be each of the forty computed local singlet columns,
including both signs, at each of V cells. Its anchor norm is one. If its
polynomial contains \(t_{r,j}\) monomials, orthonormal coefficient columns
and the norm-one CAR monomials give
\(\|F_{x,r,j}\|\le\sqrt{t_{r,j}}\). Define the complete instrument

\[
 K_{x,r,j}={F_{x,r,j}\over\sqrt{40Vt_{r,j}}},\qquad
 K_f=\left(I-\sum_{x,r,j}K_{x,r,j}^\dagger K_{x,r,j}\right)^{1/2}.
                                                               \tag{6}
\]

Each success effect is at most \(I/(40V)\), so the square root exists.
Every effect commutes with all five occupations; functional calculus
therefore gives \([S,K_f]=0\). The input is the original selected density;
the outcome maps are explicitly \(K_b\rho K_b^\dagger\), and their sum
can change the output source density.
An isometric dilation retains the complete branch label and a finite
memory carrying the opposite native increment and parity. Since total
charges and total parity match on every branch, the isometry extends to
a unitary separately in each finite total-charge/parity block. These are
declared apparatus resources, with all failure and reference records kept.
There are \(40V+1\) original labels, or 321 at member zero.

Let X be the cyclic pointer shift and \(P_s\) the actual spectral
projectors of S. The two unitary couplings are

\[
 W=\sum_sP_s\otimes X^s,\qquad
 W(K_b\otimes I)W^\dagger=K_b\otimes X^{d_b},             \tag{7}
\]

where \(d_b\) is the branch increment and \(d_f=0\).
The apparatus applies \(W^\dagger\), then (6), then W, without a pointer
read or source wait inserted between these operations. Equation (7)
follows directly from \(P_{s+d_b}K_bP_s=K_bP_s\). It is an operator
identity on the full source, so unknown initial charges, coherent
superpositions and arbitrary entanglement with the idle reference obey it.
Measuring the first pointer would define a different apparatus.

On a blank pointer, only record \(d_b\bmod7\) occurs. The fine-history
maps are exactly \(K_b\rho K_b^\dagger\) with that deterministic
additional label. The formal refined alphabet has \(7(40V+1)\) labels;
zero-probability labels are retained. Both couplings need \(56V\)
controlled native-mode additions. The complete finite source-energy
instruments can be appended on both sides; their factors \(f_i(H)\)
commute with S and the refined history is
\(1_{r=d_b\bmod7}f_j(H)K_bf_i(H)\). This preserves every original
source/reference energy-history probability without success normalization.

The source-marginal faithfulness theorem gives \(\rho\ge cP_G\),
where member zero has \(c=2^{-582407133}\). The filled-beta, constant-link
anchor \(\Omega\) is physical and
\(\|F_{x,r,j}\Omega\|=1\). Hence every success row obeys

\[
 \operatorname{Tr}\rho K_{x,r,j}^\dagger K_{x,r,j}
     \ge {c\over40Vt_{r,j}}>0.                            \tag{8}
\]

This is original unconditioned probability, with the same reference kept.
It supplies no uniform physical residue or energy-line identification.
No operator domination of the joint pure source/reference state is used.

The conserved-charge finite theorem
[`chargeExpectation_preserved_of_dual_fixed`](../../Lean/Dynamics/ProtectedCharge.lean)
applies to the actual source-time record kernel in a joint occupation
eigenbasis: \(K_t(i,j)=|\langle j|e^{-itH}|i\rangle|^2\).
Equation (1) makes nonzero entries preserve each charge; unitarity gives
row normalization, so the pullback fixes every \(N_r\). The finite
push/pull maps satisfy the theorem's explicit dual-pairing premise.
This invocation concerns the time evolution, not charge-changing (6).
The quantum cancellation (7) is proved directly and is stronger than
that classical marginal theorem. For the final pointer read, its seven
orthogonal projectors instantiate
[`ProjectivePartition.kraus_complete`](../../Lean/EventAlgebra/PartitionPinchingCP.lean).
Completeness of the source instrument follows independently from (6).

## 4. A nonzero multi-excitation alias and the identity boundary

Use the two orthonormal \(G_1\) columns extracted by the actual source
spin-Casimir projector. Their two cubic creation polynomials have a
nonzero product on \(\Omega\): the exact CAR calculation gives one
filled-six-mode vector with coefficient \(-1\) and norm squared one.
Its native increment is \((0,4,0,0,2)\), twice the weak trion increment.
Its S increment is 2, equal to the single even determinant-hole's pointer
record. The full five-charge increments differ. This is a literal source
alias, so the seven-state decoder only identifies its stipulated single
intervention alphabet. Neither a final family census nor a measurement of
how many physical particles exist follows from it.

The six-creation product is gauge invariant and nonzero. Adding it plus
its adjoint to the autonomous H would violate (1), because its charge
increment is nonzero. A compensating external memory is allowed by (6).
This establishes a source interaction-selection constraint, with no
identification of that term as a physical Majorana mass or decay.

The calculation therefore distinguishes the three source composite
families operationally and rejects specific current/parity-preserving
direct maps. Physical confirmation requires a common spatial/clock
readout with original finite-energy weight, relativistic transformation
and interaction data. Conditional constituent names do not determine a
neutrino, neutron, mass, lifetime, or GeV conversion.
