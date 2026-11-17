# Retaining the actual current-image block

This calculation keeps the frozen full625 auxiliary source unchanged. It strengthens spectral bounds by retaining more of its actual states. It does not replace the current or the source by an effective interaction.

Let R be the complete 20-dimensional old-source low SU4 multiplicity frame at the stated physical spin and isospin. Let A=A_L+A_R create the adjoint/vector auxiliary together with its quark current, without the coefficient g. The retained frame is `[R, A R]`. Every one-auxiliary current output is included. The frame has Gram matrix `diag(G,D)`, where D=R†A†AR and G=R†R.

The formulas below are written in orthonormal coordinates; the code inserts each required inverse Gram matrix. Write H for the old quark Hamiltonian, K for its signed hopping part, D for ΣG_xai², and Δ=3κ.

A intertwines each local quark colour charge with the corresponding total quark-plus-auxiliary colour charge. Therefore the electric term intertwines exactly. Mass commutes with A. The only nontrivial commutator is C=[K,A]. The proven double-commutator identity gives

    Z = A†C = [H,D]/2 − K,
    F = A†H_Q1 A = ΔD + {H,D}/2 − K.

The retained H block is

    H40 = [ H    gD ]
          [ gD   F  ].

For the exact H² block, H applied to a zero-auxiliary retained vector lies entirely in this 40-dimensional frame. Hence its leakage covariance has exactly zero first20 rows and columns. On a one-auxiliary current vector,

    H_Q1 A = A(H+Δ)+C.

Thus, setting J=H+Δ,

    W = A†H_Q1²A = J D J + J Z + Z† J + C†C.

The zero-auxiliary image contributes g²D². The two-auxiliary image is **2g A_L A_R R**, so its complete norm contributes **4g² D_L D_R**. The actual lower-right H² block is consequently

    g² D² + W + 4g² D_L D_R.

Here D² is formed from the full CAR current images, not from a compressed square. The exact leakage covariance is

    V = g² D² + W + 4g² D_L D_R
        − g² D G⁻¹ D − F D⁻¹ F.

All these are finite source identities. No resolvent denominator has been replaced inside H40 or its H² covariance.

## Literal CAR commutator outputs

For a single spin-colour generator g_ai, define

    B_ai = c_L† g_ai c_R − c_R† g_ai c_L.

The two cell commutators have opposite signs and distinct auxiliary outputs, so

    C†C = 2t² Σ_ai B_ai†B_ai.

`moments.py` constructs every one of the24 B_ai maps directly on the complete CAR occupation carrier. Imaginary generator phases cancel only within each exact norm; the eighth colour generator retains its extra factor1/3. All output words and output norms are retained, including colour/spin complements. The integer normalizations are recorded per generator.

The two local D images are constructed separately using the frozen exact local CAR identity. Their cross Gram is D_L D_R, and their summed image exactly reproduces the previously frozen full D image. This supplies the entire two-auxiliary term without dropping outcomes.

## Direct full-carrier control

`full_carrier_control.py` independently applies the Hamiltonian on the entire Nq=1 carrier, dimension24×625=15,000, without a Gauss projection. It applies quark hopping/mass, total electric charge squared, both auxiliary gaps, and all48 cell/current transition operators directly. None of the compressed H/H² formulas enters that action. At κ=1/2,1,2, direct H and H² compression agrees with the derived blocks to below9×10⁻¹⁴. This is a control of the operator identities; it is not the Nq=9 baryon spectrum.

## Full-source bounds and proof parameters

The complement contains all other zero-auxiliary physical S,I irreps and Q states orthogonal to A R. Their diagonal lower bounds are e_high and q from the frozen spin-refinement packet. Their mutual coupling norm is at most a=g√d_max. The uniform two-by-two bound gives a complement lower bound α. Exact positivity of

    H40 − E G40 − V40/(α−E)

proves a full-source lower bound E when E<α.

A less uniform bound separates the leakage into the two complement sectors. Its exact squared images are

    V_P = g²(D²−D G⁻¹D),
    V_Q = V−V_P.

For any positive η with

    h_P=e_high−E−η>0,
    h_Q=q−E−a²/η>0,

the complement is bounded below by diag(h_P,h_Q), by the operator Cauchy inequality. Therefore exact positivity of

    H40−E G40−V_P,40/h_P−V_Q,40/h_Q

also proves a lower bound for the actual full source. The rational allocation parameter η changes only a proof inequality. It does not change H, its inherited couplings, a physical unit or any observable target. Every accepted endpoint is checked by exact integer inertia; floating computations only propose endpoints and proof allocations.

The compressed H40 ground upper bound is variational for the actual source and is independently enclosed by exact generalized inertia. A finite-current trial space is not identified with a full-source eigenstate. Any unsuccessful bound is retained, and previously certified stronger bounds remain valid.
