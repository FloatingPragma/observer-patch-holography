# Actual shared-vacuum auxiliary source

The extension has the same24 quark CAR modes and mass/hopping coefficients as the frozen `../SPIN_FLAVOUR_MANIFEST.json` source. Each cell now carries an independent25-state auxiliary: a scalar vacuum plus a colour-adjoint, physical-spin-vector excitation. The full auxiliary carrier has625 states, including simultaneous excitations at both cells. The full finite carrier has10,485,760,000 states; the proof uses exact operator identities and bounds, not a claim that this full matrix was diagonalized.

On the total colour-Gauss singlet space, declare

`H = m(NL-NR) + H_hop + kappa(Q_quark,L+Q_aux,L)^2`

`    +3 kappa(n_aux,L+n_aux,R) +g sum_x,a,i (|ai><0|+|0><ai|)_x G_xai`,

`H_hop=-(1/8) sum_alpha,c(c_L† c_R+h.c.)`,

`G_xai=sum_flavour c_x†(T_a tensor sigma_i/2)c_x`,

with `tr(T_a T_b)=delta_ab/2`, `g=1/8`, the inherited kappas1/2,1,2, and the actual432-analyzer light coefficient m. The auxiliary colour and spin act through the adjoint and vector generators already frozen in the local interface. Total angular momentum includes auxiliary spin. Isospin acts on quarks; the auxiliaries are weak neutral. This is an explicit new construction, not an inherited Dirac or chromomagnetic field theorem. No physical target sets a coefficient.

The total Gauss constraint includes quarks and auxiliaries at both cells. With zero exterior flux, the electric energy is the left total colour Casimir. Its use is explicit throughout; retaining the quark-only electric term after adding a charged auxiliary would define a different source.

## Exact local algebra and denominator images

`local_algebra.py` proves on all4096 one-cell CAR states

`24 G² =4N(12-N)-12I²-6C_colour-4S²`.

It also reconstructs the complete colour/spin/isospin character, with4096 dimensions, and verifies every branch of all92 local colour singlets. The N6 singlet representation50 contains `(S,I)=(1,0)` with G²=17/3 and `(0,1)` with G²=5, along with its other complete branches. The largest G² on the full local carrier is17/3. Exact two-cell norms of `G_L²+G_R²` at N=0,3,6,9,12 are0,4,8,29/3,34/3, followed by their conjugates. These norms include the quark-colour complements needed by the auxiliaries.

For strong-coupling reduction, let `P_el` be both auxiliaries in vacuum with quarks colour singlet separately at each cell. Its complete dimension is92²=8464. Define `H0` to be the electric and auxiliary-gap terms. The complement has gap at least4kappa/3. Actual creation and removal images of every92 local source column have colour Casimir4/3; all24 tests pass. All24 independent mixed-current images have quark-adjoint Casimir3. The exact adjoint-dual Ward identity shows that the current plus the auxiliary created at the same cell has total colour zero. Thus a hopping image has H0 energy4kappa/3 and an auxiliary-current image, at either cell, has energy3kappa. These statements act on the complete P_el domain, not a selected three-quark state.

All144 local projected bilinears satisfy

`P_el c_(alpha,c)† c_(beta,d) P_el=delta_cd F_alpha,beta/3`.

Together with the CAR identity this gives

`P_el K_dimensionless² P_el=N-(2/3)sum_alpha,beta F_L(alpha,beta)F_R(beta,alpha)`.

The second-order operator is consequently

`H_eff2=m(NL-NR)-(3t²/4kappa)[N-(2/3)sum F_L F_R]`

`       -(g²/3kappa)(G_L²+G_R²)`.

No spin penalty has been fitted or inserted. The local current and actual excitation denominators generate its spin dependence.

## Complete effective calculation and its limit

`effective_second_order.py` includes all8464 states using204 exact spin/isospin highest-weight blocks, of maximum dimension14. It verifies highest-weight kernels, exact dimension reconstruction and612 member/block inertia certificates. All local occupancies and baryon-number sectors remain. The effective vacuum has N12,S=I=0; its lowest B=+/-1 states have S=I=1/2. The effective B2 ground still has competing `(S,I)=(1,0)` and `(0,1)` sectors. Competition chooses local occupancies3+3, so the on-cell N6 preference does not by itself split that nuclear pair.

These are properties of H_eff2. The full-source remainder is essential. Exact CAR Howe identities bound the complete P_el hopping-square norm by20, attained in the N12 SU4 singlet, and the P_el-to-complement coupling squared by47/96. The perturbation norm is bounded by

`v=12 sqrt(t²+m²)+2g sqrt(17/3)`.

At an actual low energy E, resolvent comparison gives

`||F(E)-H_eff2|| <= (47/96)(v+|E|)/[gamma(gamma-v-|E|)]`,

where `gamma=4kappa/3`, only when the denominator is positive. The actual full-source energy interval below supplies |E|. This bound fails at kappa1/2 and1; at kappa2 it exists but is1.016959 source energy units. It is far too broad to identify a physical mass or transfer the effective baryon spin. No effective mass is promoted to a full-source or physical result.

## Nonperturbative full625-carrier bounds

Use a different projection `P_aux`: both auxiliaries in vacuum, retaining the entire frozen interacting quark source, including nonlocal colour flux. Its complement Q_aux includes all one- and two-auxiliary states. Q_aux H_g Q_aux connects those two excitation sectors; it is retained and bounded.

At fixed quark number N, let a_N² bound the P_aux-to-Q_aux coupling and b_N² bound the one-to-two-auxiliary coupling. Their exact values come from the complete local character and the identity above. The free quark one-body energies are +/-omega, with `omega=sqrt(t²+m²)` and12 levels of each sign. Thus the free fixed-N ground is at least `-min(N,24-N)omega`.

At N12, one auxiliary forces the quarks to transform as a colour adjoint. The filled negative band is a unique colour singlet, so this sector requires at least one particle-hole pair and has free energy at least -10omega. Two auxiliaries can form a colour singlet, so their conservative free bound remains -12omega. This exclusion is used only where applicable.

Write d1 and d2 for the resulting one- and two-auxiliary diagonal bounds, including their gaps3kappa and6kappa. Then

`q_N=(d1+d2-sqrt((d1-d2)²+4b_N²))/2`

bounds the full Q_aux block below. If p_N is the certified old-source ground lower bound,

`E_N >= (p_N+q_N-sqrt((p_N-q_N)²+4a_N²))/2`.

All radicals are enclosed with exact rational arithmetic. No scalar3kappa denominator replaces the full Q_aux resolvent.

A literal two-dimensional dressed vacuum trial improves the upper bound. Its first vector is an exact integer-coordinate trial in the frozen35-dimensional SU4-singlet frame. Its second is the actual current-created one-auxiliary vector. The SU4-invariant matrix elements obey

`G_L²+G_R²=32/5-(2/45)M²+E6/180`,

where `M=NL-NR` and `E6=6 C_left²` in the old source. Exact CAR Howe identities and a transitive Clifford conjugation graph on the15 traceless Pauli generators verify the symmetry-average constants. The actual kinetic double commutator is `sum[G,[G,H_hop]]=2H_hop`; the electric term intertwines with the total colour including the auxiliary. Therefore the auxiliary trial diagonal is exactly

`3kappa+(<G² H_old>-<H_hop>)/<G²>`.

The integer trial coordinates, exact norm, all rational matrix elements and a rational upper enclosure of the two-dimensional Rayleigh eigenvalue are retained. A floating eigensolver selects a trial vector only; its bound is then exact for the declared vector on the full625 source.

Every baryon-number sector is bounded. The resulting full-source vacuum is in N12 at all three kappas. At kappa2, the upper bound of the S=I=0 trial is -0.142598506, below the lower bound -0.142339201 for every nonzero-spin or nonzero-isospin competitor. Thus its full-source vacuum has S=I=0. That spin conclusion is not established at kappa1/2 or1.

| kappa | Conditional actual lowest B=+/-1 gap bounds, GeV | Full vacuum N | Full vacuum spin/isospin |
| --- | --- | --- | --- |
| 1/2 | [0.362735,3.187658] | 12 | Unresolved |
| 1 | [0.429103,1.535633] | 12 | Unresolved |
| 2 | [0.241161,0.678238] | 12 | S=I=0 |

The conversion is the earlier conditional common-energy map; a common spatial Dirac clock is still absent. The kappa2 lowest-baryon interpretation cannot reach the proton mass under that map. The bounds do not identify the full-source baryon spin or exclude an arbitrary higher excited line. Full-source B2 binding remains unresolved; its interval and every N-sector threshold are retained in ALL_THRESHOLDS.json.

This packet establishes an actual shared-vacuum extension, exact generated interactions, a complete effective calculation, and nonperturbative bounds on the full source. It does not solve the physical baryon masses, neutron-proton gap, nuclear levels/halos, relativistic current or native execution requirements. The source remains isospin symmetric. Further spin-resolved Feshbach or direct-sector work must use these full-source states and complements rather than treating H_eff2 eigenvectors as exact states.
