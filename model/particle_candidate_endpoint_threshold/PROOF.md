# Separated compensated charged responses

## Domain and result

Fix member zero of the [primary contract](../particle_common_contract/PROOF.md), including its sixty-mode-per-cell Wilson action, eight-cell cube, physical Gauss carrier, selected polynomial state, purification reference and source clock. The apparatus below changes none of these objects. Its energy units are source clock units. A complete measurement with a chosen shortest path length one, two or three has a strictly positive recorded, dressed energy-transfer weight exceeding respectively **4/25, 2/25, 7/250** in the intervals

\[
 [55/32,1337/32],\quad[239/32,2161/32],\quad[263/32,3145/32].
\]

Every row retains its opposite-charge endpoint and Wilson flux. These finite compensated pair responses supply no single-endpoint mass, asymptotic charged threshold, physical species assignment or GeV conversion.

## 1. Paths and complete apparatus

For each starting cube vertex choose an ordered list of distinct axes and flip these bits. Lengths \(l=1,2,3\) give \(N_l=24,48,48\) oriented paths. Each endpoint occurs \(N_l/8\) times and each undirected link occurs \(N_l l/12\) times. Each path touches respectively two, three or four plaquettes. The squared coordinate endpoint separation is \(l\). These statements are enumerated independently in the receipt and tests. They concern shortest cube paths, with no extrapolation in volume.

Let \(L_e=\chi_6(U_e)P_{0,e}\), using the orientation of the path, and \(L_p=\prod_{e\in p}L_e\). Each character has Casimir \(3/10\). The cutoff is per link: each factor has degree one, admitted by member zero. Multiplying different links does not require a larger per-link word cutoff.

Use the positive-beta spin modes at \(x\), negative-beta modes at \(y\), and the four orthonormal spin coefficient matrices \(\tau_a/\sqrt2\), where \(\tau_0=I\) and the other three are Pauli matrices. The coherent rows are

\[
 A_{p,a}=L_p\sum_{s,t}(\tau_a)_{st}c^\dagger_{x,+,s}c_{y,-,t}/\sqrt2.
\]

The full occupation calculation in [the pair channels proof](../particle_primary_pair_channels/PROOF.md) gives \(\|A_{p,a}\|\le1\) and

\[
 I-A_{p,a}^\dagger A_{p,a}\le D_{xy}/2+\sum_{e\in p}(I-P_{0,e}).
\]

Divide each row by \(\sqrt{4N_l}\) and include the positive square-root failure. Changing from these four coherent spin rows to the four raw spin entries is a unitary outcome transformation, preserving the summed complete spectral measure. The dressed raw rows use the existing local odd field \(D_{y,t}\):

\[
 C_{p,st}=L_p c^\dagger_{x,+,s}D_{y,t}/(3\sqrt{N_l}).
\]

The exact full-CAR four-row effect is at most identity, as proved in [the charged dressing proof](../particle_candidate_charged/PROOF.md). Multiplication by the link partial isometry preserves this inequality. Each length setting has \(4N_l+1\) complete outcomes, including its failure. A uniformly recorded choice among all three settings has 483 outcomes before temporal and energy pointers; each joint setting weight has an additional factor \(1/3\). The displayed lower bounds are for a fixed apparatus setting.

## 2. Gauss matching and compensation

The product of character shifts has one unit of oriented electric flux on each path link. The divergence cancels at every interior vertex. Endpoint matter increments are \(+6,-6\), matching the endpoint divergence. The complete row commutes with all Gauss generators and with global native occupations. It is even under CAR parity; its regional factors are odd spinors. The exact finite chain control checks these identities for all three lengths and every local occupation in its two-mode control. The full local dressing and spin identities are inherited from the complete Fock calculations, rather than inferred from this small chain control.

The source has zero total charge on its closed physical carrier. Compression of a bare operator with nonzero total charge to that carrier is zero. A regional charged factor can be described between boundary charge sectors, while the physical probe retains its compensator. Neither operation identifies the complete even probe with a single odd particle.

A coherent local increment meter can be placed immediately before and after the instantaneous row inside each controlled temporal conjugation. For each path row, the difference is the same endpoint vector \((+6,-6)\) independently of the temporal node; the failure has zero increment. Thus the meter record factors from the temporal superposition and preserves its interference. This measures the intervention increment. It does not measure the final absolute local charge after subsequent source evolution or energy readout. The existing [complete charge meter](../particle_primary_charged_readout/PROOF.md) supplies the finite pointer and retains its quantum reference.

## 3. Positive weight in the original selected state

The averaged row deficit obeys

\[
 1-w_l\le\langle D_6/16+(5l/18)\sum_e C_e\rangle_\rho
 \le l\langle D_6/16+(5/18)\sum_e C_e\rangle_\rho.
\]

The first coefficient follows from endpoint uniformity; the second uses the character Casimir and link uniformity. The last inequality uses positivity and \(l\ge1\). The authenticated [charged-string allocation](../particle_primary_charged_string/PROOF.md) gives \(w_1>993/1250\). Consequently

\[
 w_l>1-l(1-993/1250)=(0.7944,0.5888,0.3832).
\]

All weights refer to the original selected state. Constant-link success is part of the apparatus and is not a substituted preparation or conditioning event.

## 4. Full-Hamiltonian energy moment

The onsite and electric terms satisfy the exact commutator identity

\[
 [H_{\rm onsite}+H_E,A_{p,a}]=\lambda_l A_{p,a},
 \qquad \lambda_l=6+63l/4.
\]

This is a center for the apparatus energy-transfer moment. It is not a minimal energy, interparticle potential or rest mass.

For the selected path links, let \(P\) project all those links to the constant sector and \(Q\) project them all to the charge-six character. The row maps \(P\) to \(Q\). A single selected-link hopping term maps the corresponding constant or character sector orthogonally away from that sector: both diagonal compressions vanish. For distinct selected links, the outputs of \(T_e A\) are orthogonal since exactly one of their path links leaves its character sector. The inputs of \(A T_e\) are orthogonal for the corresponding constant sectors. The two aggregate blocks have mutually orthogonal input and output supports. The full-Fock Haar bound \(\|T_eP\|,\|T_eQ\|\le\sqrt{40}\), supplied by the charged-string allocation, therefore yields

\[
 \left\|\sum_{e\in p}[T_e,A]\right\|\le\sqrt{40l}.
\]

The endpoint off-path hopping contributes at most two. At either endpoint the two remaining Wilson image columns are an isometry. They are orthogonal to the opposite original endpoint modes because the path is shortest; at length one their common link is selected, and at greater lengths the endpoints are not neighbors. The corresponding two-raising-pair Fock operator with coefficient \(\tau_a/\sqrt2\) has norm at most one. Sum the two endpoint contributions. Interior matter hopping off the path commutes with the row. Overlap between the two image spaces does not invalidate the triangle bound. The executable checks the literal Wilson columns for every path, endpoint and spin mode.

Each incident magnetic plaquette has norm at most \(1/7\), giving commutator norm at most \(2/7\). Thus

\[
 \|[H,A]-\lambda_l A\|\le\sqrt{40l}+2+2(l+1)/7,
\]

whose square is strictly below \(v_l=(80,140,199)\). The native product word and Gauss projections commute with the row; compressing the untruncated identity preserves the norm bound. No hopping or plaquette term is removed from the source.

Because \(\rho\) commutes with the source Hamiltonian, its complete two-energy spectral measure satisfies

\[
 \mu_l\{|\omega-\lambda_l|\le r_l\}>w_l-v_l/r_l^2,
 \qquad r_l=(20,30,45).
\]

These ideal intervals are \([7/4,167/4]\), \([15/2,135/2]\), \([33/4,393/4]\). They overlap; no separation of physical energy levels is inferred from their differing centers.

## 5. Dressed transfer and retained records

The exact full-CAR dressing inequality is

\[
 \sum_t(D_t-\sqrt2 h_t)^\dagger(D_t-\sqrt2 h_t)
 \le(3/2)(4-N_L).
\]

Uniform endpoint occurrence and \(L_p^\dagger L_p\le I\) give, for every path length, the same complete-history error

\[
 \operatorname{Tr}\rho\sum(C-rA)^\dagger(C-rA)\le\epsilon=\delta/48,
 \quad r=2\sqrt2/3,\quad\epsilon<(129/400)^2.
\]

Here \(\delta\) is the authenticated original weak-plus-singlet defect-count bound. The factor 48 follows by summing the two outside positive spin modes (at most two), averaging the eight endpoints and using the original negative weak-hole count at most \(\delta/2\). It does not grow with path length. The triangle inequality on complete purified energy histories gives dressed weight at least \((\sqrt{(8/9)p_l}-\sqrt\epsilon)^2\), with \(p_l=w_l-v_l/r_l^2\).

Use the same normalized triangular filter of width \(T=1/128\) as [charged propagation](../particle_candidate_charged_propagation/PROOF.md). Its probability retention throughout an interval of positive upper endpoint \(u\) is at least \((1-(Tu)^2/24)^4\). The complete midpoint and Cayley compiler is independent of the number of probe outcomes: the column identity is \(\sum K_j^\dagger K_j=I\). Its event error is at most \(6/65536\). Two complete source energy clocks of resolution \(1/64\), tail allowance \(1/32768\) per clock and implementation allowance \(1/32768\) per clock give event loss at most \(1/8192\) and expand the transfer interval by \(1/32\) at either end. All analyzer rows, probe failures, clocks and reference records are retained. Exact rational arithmetic proves the lower bounds stated above. No success-normalized probabilities enter them.

## 6. What the separation records determine

The construction establishes positive dressed responses with independently labelled separation, flux and compensating endpoint on the same source. The three-length menu is finite. It supplies no large-separation plateau, charged-sector ground threshold, populated zero-momentum pole, dispersion relation or relativistic spin dictionary.

Even a complete finite pair-energy table cannot determine a single endpoint mass through \(E_{\rm pair}(R)=2m+V(R)\) unless the interaction/compensation convention is fixed: \(m\mapsto m+d\), \(V\mapsto V-2d\) preserves every pair datum. This algebraic ambiguity does not assert that arbitrary decompositions are physically realized. It identifies the missing information in an inference from pair data alone. A defensible endpoint rest mass requires a charged-sector or asymptotically controlled compensator construction, a joint source energy/momentum threshold, the corresponding relativistic field dictionary, and one common physical unit calibration. Dividing the neutral pair response by two supplies none of these premises.
