# Actual shared-vacuum auxiliary source: spin and mass refinement

This packet keeps the frozen `shared_auxiliary` Hamiltonian and all its source parameters. It strengthens bounds on that actual source. It introduces no fitted coupling, no change of electron coefficient and no replacement of the full source by the second-order effective Hamiltonian.

The full carrier is the 24-quark-mode CAR space times both complete 25-dimensional auxiliary registers. The physical subspace obeys global **total** colour Gauss law, including the auxiliary charges. The Hamiltonian is

\[
H=m(N_L-N_R)-\tfrac18K_{\rm hop}
 +\kappa(Q_{q,L}+Q_{a,L})^2+3\kappa(n_{a,L}+n_{a,R})
 +\tfrac18\sum_{x,a,i}(|ai\rangle\langle0|+|0\rangle\langle ai|)_xG_{xai}.
\]

Here `K` in the code already includes the negative hopping sign. The mass coefficient is the unchanged rational 432-analyzer coefficient. The declared coupling list remains κ=1/2,1,2. All 625 auxiliary outcomes remain in the source, including both-excited states.

## Complete spin counts and coupling norms

`branching.py` uses exact Gelfand–Tsetlin characters and Weyl differences. It independently sums the complete colour-singlet quark spaces: 14,112 states at N=9 and 24,696 at N=12. Their largest physical spin/isospin highest-weight blocks have dimensions 152 and 274. The SU4-to-spin/isospin decomposition includes every old-source irrep.

Using the frozen full local CAR character and the exact local mixed-current identity, it computes the maximum of

\[
D=\sum_{xai}G_{xai}^2
\]

on each colour-singlet physical S,I sector. Thus the exact P₀-to-Q₀ squared coupling norm is at most g²||D||. The free one-auxiliary Q sector requires quark colour 8 and an auxiliary spin 1. The two-auxiliary sector retains quark colours 1,8,10,10-bar,27 and auxiliary spins 0,1,2. Colour tensor products, spin/isospin triangles and band occupations are enumerated exactly. An absent sector is explicitly recorded; it is not silently treated as a zero-energy state.

## A stronger bound on the actual Q Hamiltonian

Within each fixed auxiliary-number Q block, project onto zero electric energy. Both local **total** colours must then be singlets. An auxiliary-vacuum cell has colour-singlet quarks; an excited auxiliary cell has colour-adjoint quarks coupled to its auxiliary adjoint. These statements determine all allowed local occupancies and flavour representations.

Between local total-colour singlets, Schur invariance gives

\[
P c^\dagger_{\alpha c}c_{\beta d}P=\delta_{cd} F_{\alpha\beta}/3.
\]

It applies equally to the quark-adjoint/auxiliary-adjoint singlet. Hopping changes local triality, so PKP=0, and

\[
PK^2P=N-\frac13(C_{U4,\rm total}-C_{U4,L}-C_{U4,R}).
\]

The local flavour Casimir is `7N-N²/3-6 n_aux`. `q_colour_bounds.py` decomposes every local U4 tensor product by exact characters, verifies its full dimension and takes the exact maximum hopping-square eigenvalue compatible with total S,I. The complement has electric energy at least 4κ/3. Combining that gap with the already bounded free Hamiltonian and the exact hopping image yields a rigorous two-by-two lower bound for each fixed auxiliary-number block. Combining the one- and two-auxiliary blocks retains their off-diagonal coupling. This bounds the actual Q Hamiltonian; it does not use an effective eigenvalue.

## Actual baryon states and their entire current images

`baryon_states.py` embeds all 20 old-source multiplicity states into actual CAR occupation bases. For Δ, the physical S=I=3/2 frame is the SU4 highest frame. The physical S=I=1/2 frame is obtained by the literal CAR intertwiner

\[
(2E_{30}+E_{31}E_{32})/\sqrt{18}.
\]

Exact Gram, colour-raising, physical spin/isospin-raising and all three Hamiltonian component identities are replayed. The complete colour-Cartan/Sz/Iz carrier dimensions are 4,408 for Δ and 14,376 for the nucleon. The integer mixed-current operator is applied on those complete carriers. Its images outside the original SU4 representation remain present.

Consequently the covariant matrices `D=R†DR` and `D2=R†D²R` are computed separately. In particular **D2 is not replaced by D G⁻¹D**. Their difference

\[
V=D2-DG^{-1}D
\]

is the exact Gram matrix of the current image outside the complete low SU4 block.

Literal integer-coordinate source vectors produce actual full-source two-component trial states `(ψ, Aψ)`, with `Aψ` containing every one-auxiliary current outcome. Their compressed auxiliary diagonal is

\[
3\kappa+\frac{\langle D H_{old}\rangle-\langle H_{hop}\rangle}{\langle D\rangle}.
\]

This uses the frozen exact identity Σ[G,[G,Hhop]]=2Hhop. All trial moments and coordinates are retained. These states are variational inputs, not asserted full-source eigenstates.

## Matrix Feshbach positivity on the full source

For fixed S,I, split the zero-auxiliary quark space into the complete 20-dimensional Sym³ SU4 multiplicity block L and every other SU4 representation H. The old Hamiltonian has no L–H cross term. Let `e_high` be the certified lower energy of all those other irreps, `q` the rigorous actual-Q lower bound, and `d_max` the full-sector ||D|| bound.

At a candidate E<q, set

\[
c=\frac{g^2}{q-E},\qquad h=e_{high}-E-cd_{max}.
\]

If h>0 and

\[
H_L-EG-cD-\frac{c^2}{h}(D2-DG^{-1}D)>0,
\]

then the **entire actual full625 source in that S,I sector** has energy above E. The first Schur bound retains Q through its resolvent lower bound; the second retains all other zero-auxiliary SU4 representations through their exact current covariance. Exact integer inertia proves strict positivity in dimension 20. Floating eigenvalues merely propose rational endpoints. Every accepted endpoint and its inertia count are saved.

All remaining physical S,I sectors use the rigorous scalar bounds. Q-only spin sectors have positive energy by the frozen uniform Q bound, and cannot undercut the negative baryon trials.

## Full-source particle–hole map

`particle_hole.py` verifies the finite generator identities for

`c_L -> c_R†`, `c_R -> -c_L†`, together with a cell swap and a signed orthogonal map on both complete auxiliary registers. For colour generators define `-T_a^T=o_a T_a`; for spin define `-S_i^T=r_i S_i`. The auxiliary excited state |ai> receives sign `-o_a r_i`.

The quark mixed current and auxiliary transition have the same sign, so their product is invariant. Total colour transforms covariantly; on global Gauss zero the left and right electric Casimirs agree. Both hopping and staggered mass are invariant. N maps to 24−N and B to −B, while total S² and I² are preserved. These generator identities lift to the entire finite CAR algebra and all auxiliary outcomes. Thus the B=−1 bounds transfer to B=+1 for this actual source, rather than by relabelling a selected old state.

## Certified consequences and remaining failures

The original source-wide N=12 vacuum-sector certificate remains valid. The stronger bounds now certify its spin/isospin zero at κ=1 and κ=2. At κ=1/2 the nearest nonzero-spin bound −0.596775654 still lies below the spin-zero trial −0.596119666, so its spin remains unresolved.

The **actual full-source** lowest |B|=1 spin and isospin are 1/2 at κ=1 and κ=2:

| κ | competing Δ lower bound | nucleon trial upper bound | conclusion |
|---|---:|---:|---|
| 1/2 | −0.453572841361 | −0.432543617884 | spin not certified |
| 1 | −0.213659134693 | −0.213995114981 | nucleon spin certified |
| 2 | −0.099834773690 | −0.100762808832 | nucleon spin certified |

Every other S,I sector is also excluded at κ=1 and κ=2. This removes the earlier spin-3/2-only and exact SU4 nucleon/Δ-degeneracy obstructions for these finite sources. It does not derive an observed nucleon mass or the physical Δ splitting.

Subtracting the same interacting vacuum gives conditional |B|=1 mass intervals:

| κ | conditional mass interval, GeV |
|---|---:|
| 1/2 | [1.059059714, 2.195481157] |
| 1 | [0.702948908, 1.054976747] |
| 2 | [0.331083321, 0.477614664] |

The endpoints κ=1/2 and κ=2 exclude a lowest-baryon/proton identification under this unit map. κ=1 has only a broad enclosing interval; containment is not a successful prediction, and κ=1 is not selected by fitting the target. The 432 electron anchor is still only a **conditional** common energy conversion: no common spatial Dirac clock or interacting pole identification has been proved.

The source remains isospin symmetric, so neutron/proton degeneracy persists. Full-source B2 binding, mobility, nuclei/halos, a physical string tension, native full-state preparation/readout and the Dirac current/clock join remain open. `SUMMARY.json` explicitly records `physical_rows_solved:false`.
