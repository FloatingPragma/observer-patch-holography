# A native spatial loop field on the colour string

This is a new architectural attachment to the exact finite SU(3) colour-chain source frozen at e222d339. It is not a claim that the parent's uniform conjugate-angle density already contained a texture. Its spatial graph has one 16-node circumferential ring around each actual chain cut. The added local Ramsey memories store a colour-frame field on those rings, coherently locked to the cut flux. A conjugated Hamiltonian preserves this field support and every parent matter observation. Source parameters and representation weights remain declared inputs. No uniqueness, spontaneous formation or common electromagnetic-source transfer is claimed.

## The group, the physical unit and the representation trace

Let z=exp(2 pi i/3), and define

    G = (SU(3) x R) / <(z I, 2 pi)>.

This is isomorphic to U(3) by [S,theta] -> exp(-i theta/3)S. The kernel of this map from the cover is exactly the displayed infinite cyclic subgroup. SU(3) is simply connected, so pi1(G)=Z. The path g(theta)=exp(i theta C), C=diag(1/3,1/3,-2/3), together with theta from0 to2pi, closes in G and generates its fundamental group. Under the U(3) isomorphism its determinant is exp(-i theta), with degree-1.

A fundamental colour state of physical electric charge q represents [S,theta] by exp(i q theta)S precisely when exp(2pi i q)z=1. Both q=2/3 and q=-1/3 satisfy this. A colourless state represents the same parameter by exp(i q theta), requiring integer q. The actual fifteen-coordinate carrier fixes the colourless charged coordinate (3,4) to q=1, and the colourful coordinates (c,3),(c,4) to2/3,-1/3. `carrier.py` checks those exact matrices, every Q+C integer and every3Q integer. The declared hypercharge weights are inputs, as in the parent source; the group condition alone would not select the integer parts of the charges.

On a cut containing p quarks to its left, the actual colour representation is wedge^p C^3. It is fundamental for p=1, antifundamental for p=2 and singlet for p=0,3. Let n=2 N_up,left-N_down,left=3 Q_left. The physical U(1) generator on a fixed cut sector is therefore

    Q_cut=(n/3) I_d,
    d=3 for p=1,2; d=1 for p=0,3.

The corresponding traceless colour generators are C_1=diag(1/3,1/3,-2/3), C_2=diag(-1/3,-1/3,2/3), and C_0=C_3=0. Since n+p=0 mod3, the closed representation loop

    W_p,n(theta)=exp[i theta((n/3) I_d+C_p)]

has integer component exponents. Examples are (1,1,0) for a single up quark, (0,0,-1) for a single down, (0,0,1) for up+down and(-1,-1,0) for two downs. The charged two-up sector additionally has(1,1,2). The largest absolute exponent is2 throughout the declared finite domain.

Its integer determinant degree and its normalized trace winding obey

    deg(det W)=(1/(2pi i)) integral tr(W* dW)=d n/3,
    (1/(2pi i d)) integral tr(W* dW)=n/3=Q_left.

The factor1/d is the normalized trace in the actual colour representation, not a change to the electron unit. The colourless electron has d=1, Q=1 and W(theta)=exp(i theta), with winding1 in exactly the same theta parameter and the same native phase reader. For a nontrivial cut d=3, the determinant degree is n, while the trace winding is n/3. No fractional element of pi1(U(3)) is asserted: the integer winding of the full matrix loop and its electric partial winding are explicitly different quantities. The fractional electric contribution closes with its colour contribution, rather than being an unclosed scalar curve renamed a winding.

For any smooth single-valued spatial colour frame R(theta), W changes to R W R*. The determinant is unchanged pointwise. Differentiating this conjugation adds a commutator term whose trace vanishes, so the normalized trace integral is also unchanged. This holds for theta-dependent frames and arbitrary component permutations. The native field records eigenphase components in a declared diagonal frame; their product and sum of degrees are invariant. No physical expectation of a preferred colour orientation in the singlet is assumed. One can regard the measured field as the conjugacy class of W; its determinant is an actual frame-independent observable.

## Actual electric registers and an actual spatial field

The parent three-quark Gauss isometry and its SU(3) link image are retained. Add three binary memory bits per edge storing n+2 (n lies in[-2,4] on the two admitted weak sectors). The word7 remains an invalid-flux outcome. These are physical cut-electric registers in thirds of the same electron unit, with electric coefficient zero. The physical constraint is n_e=3 sum_(x<=e) Q_x. In the neutral primary sector both external fluxes are zero; in the charged comparison sector the right boundary flux is3, matching its total charge1. The extra U(1) field is an exact bounded-support Gauss dilation on this fixed three-quark domain, not a claim about a truncated rotor outside it.

Every cut also has16 spatially labelled ring vertices theta_j=2pi j/16, each with three Ramsey memories. In the singlet cut only the first component is active; the other two are explicitly constant padding fields. If w_r is an integer exponent of W, the actual local complex field is

    f_(e,j,r)=<X_(e,j,r)+iY_(e,j,r)>=exp(i w_r theta_j).

It is prepared by a native Hadamard followed by exp(-i w_r theta_j Z/2), controlled by the full source code. The reference coherence makes the phase of each matrix component observable. A bare Choi ray would lose a common matrix phase and would not suffice for the electron comparison; that is not the proposed instrument.

The new field is a dual-electric colour-frame field on a spatial ring. It is not the conjugate-angle wavefunction of an electric eigenstate. The old uniform compact-angle density remains uniform. The new local coherences vary across actual distinct memory sites when the cut carries a fractional flux. This addition is disclosed architecture, not a consequence retrospectively attributed to the old source.

The analytic interpolation along each ring is part of the declared field definition: W(theta)=exp(i theta(qI+C_p)). Integer exponents bounded by2 prove exact closure and exclude aliasing. The16 samples alone would not certify arbitrary continuous functions. Here each component's nearest-neighbour phase increment has magnitude at most pi/4; the determinant increment is at most pi/2. Thus principal-argument sums on the ring equal the continuous degrees, with positive branch margins. Native phase errors are controlled below these margins; no winding is inferred from an unrestricted interpolation.

## Hamiltonian attachment and preserved source dynamics

Let C be the complete finite unitary in `native_extension.py`. It first XORs the exact n+2 bits into the cut registers, then prepares every controlled Ramsey field. It is defined on every source padding column and every complementary initial register state. Let |0_A> be the added blank memories and define V psi=C(psi tensor|0_A>). Orthogonality of source codes makes V an isometry even when different conditional field states overlap.

Define the full finite attached Hamiltonian by

    H_tilde=C(H_parent tensor I_A)C*.

On im(V), H_tilde V=V H_parent and exp(-it H_tilde)V=V exp(-it H_parent). The orthogonal complement is explicitly defined by the same unitary conjugation; it is not removed or assigned an implicit successful probability. Parent padding energy16 is retained. For every parent operator O, its physical lift V O V* has the same matrix elements on im(V). The source occupation law, colour string energy and source time are preserved exactly. The added field changes whenever a source hop changes a cut flux; it is not only a terminal label written after an unrelated evolution.

On a rightward hop of charge q, n_e changes by-3q. This is exactly the Gauss-preserving link shift; spin exchange terms change no cut charge. All 1,678 source off-diagonal pairs and 2,432 source-code/cut configurations across L=2..5 and both weak sectors are checked. The native free word is C*, followed by the complete parent source word, followed by C. The initial C is also actually generated. Native pulses may leave the physical image during compilation; the comparison concerns the complete declared endpoint channel, as in the parent bridge, with every export retained. There is no claim that the compilation pulse count is physical source time.

Since n_e+p_e=0 mod3, n_e is nonintegral in electron units exactly when p_e=1 or2. On the physical image these are precisely the cuts between the first and last quark, with colour Casimir4/3. Thus the partial winding is confined to the actual colour string, not to an arbitrary normalization convention. Local charge is the difference of neighbouring partial windings, including the explicit boundary fluxes. The direct-field readers for an up-left and a down-left separated preparation therefore see2/3 and-1/3 on their first cut, with the same electron control.

## A persistent separated string at every source time

The primary member is L=5,kappa=2,m=1,t=1/8,U=1,J=1/4, neutral weak count1. Let P be the colocated subspace. The exact parent block estimates are

    PHP <=75/16,  ||PH(1-P)||^2<=3/32.

For the first maximal-span basis codeword, the exact energy mean and variance are mu=707/48 and sigma^2=3/64. For the two additional specified codewords (0,9,19) and(1,8,19), with individual quarks at cells0,2,4, the same mean holds and sigma^2=1/16. The latter distinguish an up-left and down-left first endpoint; these labels follow the parent's actual cell/spin/weak indexing.

Split the exact spectrum at energy10. Chebyshev applied to(H-mu)^2 gives initial low-energy weight at most108/51529 for the first preparation and144/51529 for either individual-endpoint preparation. These weights are conserved. The subspace Sylvester bound between PHP and the high spectral subspace gives

    ||P Pi_[10,infinity)(H)||^2 <=(3/32)/(10-75/16)^2=24/7225.

For every time, the norm of the colocated component is at most the sum of the square roots of these two bounds. Using(a+b)^2<=2a^2+2b^2 proves

    Pr(colocated) <=2(108/51529+24/7225)<0.011

for the primary preparation, and less than0.013 for the individual-endpoint preparations. Hence a separated string persists with probability above0.987 at every source time in all three preparations. Finite native errors add to the complete law, including failure outcomes counted as loss of a certified string. This is not extrapolation of a finite-time grid.

For kappa=2 the local colour coefficient vanishes. In a fixed span-s subspace the electric energy is exactly(8/3)s, the hopping contribution lies in[0,3/2], and the local spin contribution lies in[3/16,15/16]. Therefore the lowest compressed energy at span s+1 exceeds the highest possible bound on the lowest energy at span s by at least8/3-3/2-3/4=5/12. The supplied Hamiltonian has positive string tension. The source has exactly three quarks and no pair-creation channel; all cut strings stay locked to their endpoints.

## Complete native observed channel and error budget

`receipt.program` specifies conclusive resets of every source, flux, field and clock memory; the exact native maximal-span codeword; C; C* exp(-itH) C; the same clock exp(-itZ/16); and native C/E pointers for all source bits, every flux bit, all192 field components on L5, and the clock. Fixed X orY analyzers give two apparatus settings. All source padding, invalid flux, source/flux mismatch, first-export failure, reset timeout, fair-coin fine label and terminal bit is retained. No successful subset is silently normalized.

For an ideal final source code i, the field state is the tensor product of its actual Ramsey states and the flux word is exactly n(i)+2. The full terminal pattern/flux/field law is therefore the source code weight multiplied by the conditional product of the field's binary quadrature laws and the clock law. This conditional product comes from V, not an independence assumption about the unobserved quantum source. With an external reference, the same formula is a CP map retaining each conditional reference block. The exact fine C/E maps refine these binary projectors with their literal native coins and endpoint records. First failures are K_failure,j K_success,j-1 ... K_success,1 with no survival division.

The source matrix entries and all source Lie-product errors are those of the finite parent. The attached field uses only exact rational multiples of pi; it introduces no eigensolver or algebraic field cutoff. Its Walsh diagonal forms and XORs define every false-control and padding column. Source evolution uses n=ceil(t rate), rate>=A^2 T/(2 tolerance), so the full product error is at most tolerance=10^-18 on0<=t<=T. The frozen resource window isT=4000, with the same construction available for any positive finite T.

Unlike the parent's narrower angle assumption, new attachment Walsh angles are bounded by |q_pi|<=2. `uniform_control_error` explicitly doubles both the pi magnitude and pi-interval-width contributions. Every other angle has q_alpha=0, and the rational clock angle is at mostT/16. The native full-six-level Duhamel enclosure covers all twelve rational controls, not merely proper-code amplitudes. Reset calls include every large auxiliary memory and both work resets per pointer. Their conditional timeout union bound requires no independence. The completed-reset comparison and successful fine-record decoration are exactly the parent construction: no unit-mass success law has timeout probability appended to it.

The exact count includes all three attachment calls, two native pulses per source bit flip, all field/clock analyzers, and every pointer. The primary L5 receipt has227-bit controls, approximately5.02e35 native pulses and complete half-diamond/TV error about1.003e-18. Exact final integers and fractions are authoritative in BUDGET.json. This enormous finite word is not claimed executed.

## Winding error margins and conditioning

For a specified sharp-code preparation, the exact field has unit-modulus coherences. A native complete-channel TV error delta changes each X orY mean by at most2delta. Thus each reconstructed complex field changes by at most2sqrt(2)delta, and its phase changes by at most4delta for sufficiently small delta. Comparing neighbours adds at most8delta to each component increment; the determinant adds at most24delta. These are far below the exact component3pi/4 and determinant pi/2 branch margins. The observed nonzero polygonal field is homotopic to the declared analytic field and has the same winding.

At a later time the complete source-code record can be correlated with the field. The ideal conditional field is exact for every code i, but a uniform joint error does not justify an unqualified conditional error for a vanishing-probability branch. If the native code probability is at least eta=10^-6, a conservative conditional TV error is2delta/(eta-2delta). This still preserves every branch margin. Other code records remain in the complete law; no winding-confidence statement is made on an arbitrarily rare conditioned record. Their total native probability is at most D eta, with D<=256. In the primary separated experiment, the all-time string probability minus this amount and the complete error is still above0.988. At the specified up/down initial preparations the selected source code has probability near1, so the fractional endpoint winding has no rare-branch denominator.

These statements concern exact population laws of the native finite programs. No finite-shot estimator or statistical confidence is claimed executed. A single random field bit does not determine a winding. The complete joint binary law and its refinements are specified, so a separate sampling protocol could be budgeted without changing the source.

## Failing controls and scope

Removing C_p while keeping the same fractional electric charge leaves endpoint discrepancy sqrt3 (both quadratures of this and the noncongruent/flat controls are executed in CONTROL_NATIVE_RESULT): exp(2pi i q)I is not closed. Its16-node polygon has an integer scalar degree that differs from q. Changing n while holding p fixed to violate n+p=0 mod3 likewise fails closure. Replacing the actual attachment by the flat-field native word leaves charges, fluxes and parent dynamics unchanged but gives zero winding on every ring; every fractional cut then fails the charge/winding identity. These controls distinguish the added physical texture from merely reading n and printing n/3.

An omitted field-update control is executed on all columns of a small dressed source slice and differs from the target. The colourless integer electron has a closed rank-one loop of winding1 under the same native reader. A normalization-only colourless fractional charge cannot furnish this quotient representation: its rank-one loop fails closure. Removing the string term removes the positive supplied electric tension; the phase construction alone does not bind an unbound triple. The literal criterion still relies on both the actual parent string dynamics and the new native loop field.

Executed:5,132 actual six-level pulses (5,024 primary blocks plus108 explicit failing-control reader pulses), including three complete flux XORs, all-sector loop components, a conjugated source slice on all64 input columns including the ancilla complement, and all native C/E fine reader laws for component weights-1,0,1,2 in both quadratures. Every first-export complement is computed. Floating residuals are diagnostics distinct from exact rational error certificates. The full ring word is not executed. The attachment is supplied architecture, not a spontaneously formed soliton, a three-dimensional Yang-Mills solution or a transfer to the calibrated Maxwell source. Two producer-independent full reviews remain necessary before a row award.
