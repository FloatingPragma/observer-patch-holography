# Rotations, exchange and the proposed spin-dependent transport

The source `native-local12-current-v0-proposed` has an ordinary geometric rotation representation. The proposed Dirac source below is a replacement of its matter dynamics. Its Clifford matrices, root-port code, mass and gauge couplings are architectural inputs. They are not selected by the original twelve-port process.

## The existing matter source

The routed cubic graph admits the twelve even signed cyclic coordinate permutations. They permute the actual icosahedral ports, cubic cells and direction events. The event permutation at a port is determined by the image of its actual neighbour. These permutations preserve the holding/moving event bit, the admitted-event subspace, the padding and the uniform-event reflection. Consequently they commute with `W=S R`. `spin.py` checks the 144 group products and the complete 2048-dimensional source permutation.

On the four-dimensional acoustic span, take the frozen continuum basis `(e,sqrt(48)g_x,sqrt(48)g_y,sqrt(48)g_z)`. The actual geometric action is `diag(1,O)`. Its infinitesimal angular momentum is zero on `e` and the ordinary vector representation on the other three components. Thus `J²=diag(0,2,2,2)` and a full `2π` rotation is the identity. The propagating vector component is longitudinal, and its helicity is zero. The scalar mass coin has `R=diag(1,-1,-1,-1)` on this span and commutes with these rotations. It does not change the spin conclusion.

The same source's velocity matrices are `|e><g_a|+|g_a><e|`. Their rescaled squares are rank-two projections, not the identity. They are not Dirac matrices. This is an exact obstruction to interpreting the existing scalar/vector acoustic coordinates as a Dirac spinor by relabelling them.

The supplied CAR lift gives a minus sign when two occupied modes are swapped. The controlled native swap experiment records that relative sign, and a double swap records its square. An exterior lift of an ordinary geometric representation remains single-valued: `Gamma(I)=I`. Exchange statistics therefore do not turn the above geometric rotation into half spin. The native geometric generator programs and exchange probe records are in `ROTATION_PROGRAM.json` and `ROTATION_EXCHANGE.json`.

## A common spin-dependent replacement

The proposed source `native-local12-dirac-routed-v3-proposed` uses four co-located modes at port zero of each cubic cell. In the native compiler's little-endian tensor convention choose

`alpha_x=XX`, `alpha_y=XY`, `alpha_z=XZ`, `beta=ZI`, and `Sigma_a=I sigma_a`.

These matrices obey the four-generator Clifford relations exactly. `S_a=Sigma_a/2` has `S²=3/4`. The spin lift `U(O)` satisfies `U alpha(p) U*=alpha(O p)` and commutes with `beta`. This ties the proposed spin transformation to spatial transport. A spectator spin lift on the old scalar source would not satisfy this spin-dependent transport requirement.

The matter macro is, in chronological order, the on-site mass word `exp(-i h m beta)` followed by the three conditional axis shifts. Each shift diagonalizes `alpha_a`, translates its positive and negative eigen-rails in opposite directions, then restores the spin basis. The diagonalizers are finite words in the retained native control generators, and their complete four-column execution and export failures are recorded. On the chosen spin lift, the controlled native `2π` and `4π` spin words have respectively relative phases `-1` and `+1` in both C and E pointer laws.

These native spin words are the spin part of a proposed physical rotation. They are not a completed finite apparatus that rotates the carrier, reference field and detector. A finite full spatial intervention also needs the coordinate permutation and its Gauss dressing. In addition, the ordered xyz lattice macro does not have exact rotational covariance. Its controlled continuum generator does. No finite geometric half-spin measurement of the old source is claimed.

For a Fourier convention with a positive shift represented by `exp(-i k)`, the free symbol is

`U_h(p)=exp(-i h p_z alpha_z) exp(-i h p_y alpha_y) exp(-i h p_x alpha_x) exp(-i h m beta)`.

For `H=alpha.p+m beta`, `H²=(|p|²+m²)I`. The eigenvalues are `±sqrt(|p|²+m²)`, each twice repeated. For `a=(m,p_x,p_y,p_z)`, the unitary Lie-product estimate and telescoping give

`||U_h(p)^n-exp(-i n h H)|| <= min(2,n h² sum_{i<j}|a_i a_j|)`.

The same operator bound holds after tensoring an arbitrary reference; a complete terminal measurement has total variation at most its minimum with one. On a bounded momentum set, its supremum gives the direct-integral free limit. This statement has no interacting-field, arbitrary-momentum or single-species extension. At the eight Brillouin corners for zero mass the full symbol is `±I`; all eight zero/π nodes are retained in the receipt. No Wilson term, projection of doublers or particle-species identification is hidden in the construction.

The bound concerns the explicitly stated homogeneous translation family with that symbol. The full routed native program executed and compiled here is its finite size-two member. A native implementation and preparation contract for every spatial refinement is not inferred from the symbol calculation.

The common clock is `t=n tau0 h` after the complete macro. Its free continuum speed is `1/tau0`, whereas the scalar walk has `1/(sqrt(48) tau0)`. These are distinct dynamics and must not share a speed label without changing the contract. DYN08's harmonic field can match this proposed speed with eight field substeps, `alpha_field=1/16` and `beta_field=3/8`: `8 sqrt(2 alpha_field beta_field/3)=1`. Its uniform curvature norm bound 158 gives the harmonic stability product `237/64<4`. This schedule is a declared architecture. It is not a proof that the finite quantum Z16 field has a canonical Maxwell limit.

The encoded one-particle native compiler uses data wires 0 through 1156 and scratch wires 1157 through 1161. Its clock and additional probe must use fresh wires starting at 1162. The older field source's wire allocation cannot be reused unchanged. A full Fock encoding has a different allocation and requires its own contract.

## The same Gauss register and field intervention

Each positive axis route is an actual path from root port zero to its chosen outgoing port, across an actual external seam, and from the incoming port to the next root. The negative route reverses the corresponding positive path. For the incidence convention `B_e=+1` at the lower endpoint and `-1` at the higher endpoint, its signed edge vector `j` satisfies

`B j=delta_start-delta_end`.

The source shifts the existing 288-edge Z16 field by `q j` while moving a charge `q`. Therefore `n+B F` is unchanged modulo 16, including every intracell leg. The internal spin basis words leave the charge at the root and require no field shift. Omitting a nonzero intracell current is a rejecting control.

The one-particle source is a unitary on the complete finite matter and modular-field space: each spin basis change is unitary, and each conditional position/field translation is a bijection. Every finite field character of the same translations gives the path phase `exp(-i q theta.j)`. Under `theta -> theta+B^T chi`, it transforms by `exp(i q chi_end) ... exp(-i q chi_start)`. The receipt checks this covariance on all 32 root-spinor columns for admitted discrete characters. Arbitrary continuous characters of Z16 are not used as physical finite-field states.

The optional fermionic extension is explicit: in the finite field-angle decomposition take the exterior lift of this unitary at each character, and transform back to field labels. This defines `Gamma(U_D)_F` with the same current convention. It retains the supplied CAR statistics. The one-particle native words alone do not constitute an execution or error certificate for every Fock sector.

The electric and plaquette kicks are precisely the Z16 field operations of DYN08, now acting on four-spinor coefficients. `dirac_proposal.py` executes complete sparse two-macro experiments for both charge signs, retaining every generated field outcome. Diagnostic experiments activate one declared face; all other face couplings are zero. They verify norm, modular Gauss, field-dependent energy and spin records. They do not stand in for the homogeneous full-face, eight-substep schedule.

At size two, positive and negative shifts reach the same neighbouring cell. A full xyz macro therefore alternates the cell deterministically. Its electric-energy intervention does not by itself establish a matter-response signal. `cycle_probe.py` supplies a stricter finite test. Let `D=j_plus-j_minus` be the difference between the actual opposite-axis routes from one root. It is a closed primitive Z16 cycle. The two-root field code is `F=t D` at the first root and `F=j_plus+t D` at the second. This code is invariant under that entire axis step and satisfies the same Gauss equation. Its coordinate updates are `t->t-1` on the first root's negative rail, `t->t+1` on the second root's positive rail, and identity on the other rails.

For the source preparation with cycle amplitudes `exp(i theta t)/4`, an axis step sends a root spinor through `P_plus+exp(i theta)P_minus`. Thus the characters `theta=0` and `theta=pi` take the same input `|00>` to `|00>` and `|11>`, respectively. The native compressed circuit executes this exact restriction and reads its spin bit through the complete C/E pointer laws. This is a field-dependent matter readout. It is not a magnetic-moment measurement or a proof that the full-cycle character obeys an integer-field cutoff bound. Native preparation of the full 288-edge embedding and reset remain separate from supplying the normalized compressed source state.

`dirac_native.py` emits a complete factorized encoded one-particle native macro. Each route-current block is checked on every clean and dirty local input, with scratch return. The original cell address is retained until all route field increments are complete, and then its axis coordinate is flipped. This compiles the endpoint channel. It does not make every low-level gate prefix a physically local, Gauss-preserving transit. An implementation using actual transit modes, its initialization/reset and a full accumulated physical-error certificate remain separate obligations. No numerical pulse residual is substituted for that certificate.

The audited positive axis routes use respectively 3, 4 and 5 actual seams. The negative routes reverse the corresponding positive paths. Every route is simple, and each edge occurs in at most two of the three axis support families in the finite source. Therefore one-particle geometric transport over `n` macros has per-edge prefix displacement bounded by `2n`; `3n` is also a valid looser bound. The common cyclic-to-integer field theorem can use this actual motion bound, together with its initial field support and the sum of all magnetic coefficients. A bare bound of one shift per whole xyz macro is not justified. For `Q` particles, `2nQ` is a conservative bound before cancellations.

## Magnetic response and its boundary

If the same routed gauge source has a controlled canonical covariant-derivative limit, set `Pi=-i grad-q A`, with `[Pi_i,Pi_j]=i q epsilon_ijk B_k`. The Clifford relations then give

`(alpha.Pi+beta m)²=Pi²+m²-q Sigma.B`.

For positive energy, weak field and momentum small compared with nonzero `m`, the leading expansion is `m+(Pi²-q Sigma.B)/(2m)`. Its intrinsic moment is `q S/m`, corresponding to `g=2`. This is a conditional algebraic prediction of the explicitly chosen Dirac coupling. It is not an extracted magnetic moment of the old walk or a demonstrated interacting Z16 continuum limit. The finite oscillator control keeps the exact cutoff commutator defect and verifies the Pauli identity with that defect; it does not discard the highest oscillator state.

The rejecting spectator construction has a nontrivial supplied SU2 lift and `2π=-I`, but its transport and magnetic Hamiltonian are proportional to the identity on the spin index. Its two spin sectors have identical energy, so its intrinsic spin response is zero. A phase test or CAR sign alone cannot establish `g=2`.

The remaining common obligations are the full physical rotation apparatus and reference dressing, the homogeneous coupled native initialization/error contract, an interacting quantum field limit with explicit cutoff and angle bounds, control of the additional lattice species, and comparison with the actual OpenWave validation criteria. No official row is changed by these artifacts.
