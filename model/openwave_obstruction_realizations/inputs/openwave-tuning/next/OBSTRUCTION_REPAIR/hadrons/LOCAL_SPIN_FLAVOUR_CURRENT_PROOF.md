# Local spin and flavour current mechanism

This result is a finite source construction. It is not a physical mass prediction, a completed two-cell vacuum calculation, or a transfer from a Dirac or chromomagnetic field source.

The one-cell carrier has two physical spins, two weak flavours and three colours: 12 CAR modes, hence 4096 full Fock states. The three-quark sector has dimension220. Its complete colour-singlet subspace has dimension20. Actual CAR spin and isospin currents split it into `(S,I)=(3/2,3/2)` of dimension16 and `(1/2,1/2)` of dimension4. The script computes both highest-weight kernels and checks that their full multiplets exhaust the colour-singlet sector.

With `tr(T_a T_b)=delta_ab/2`, define

`G_ai = sum_flavour c† (T_a tensor sigma_i/2) c`.

The complete three-quark colour-singlet operator identity is

`sum_ai G_ai² = 9/2 - (2/3) S²`.

The script builds this operator from literal CAR bilinears using the colour Fierz identity, and checks the identity on every colour-singlet basis vector. Its eigenvalues are4 on the spin/isospin doublet and2 on the quartet. A positive current square alone would therefore favour the wrong spin; the sign must come from an actual interaction, not a desired ordering.

One explicit extension adds an auxiliary scalar vacuum and a colour-adjoint, physical-spin-vector excitation, of total dimension25. Declare

`H_aux = Delta P_exc + g sum_ai (|ai><0|+|0><ai|) G_ai`,

with `Delta=3 kappa` and `g=t=1/8`. The gap uses the adjoint quadratic Casimir3 and the inherited electric coefficient; using the inherited hopping coefficient as this new current coupling is a declared architectural choice. Neither choice is forced by OPH or derived from the previous two-cell Hamiltonian. Neither is selected using a target splitting.

The mixed current maps a colour singlet to a quark adjoint. The actual colour Casimir has a40-dimensional adjoint zero-weight eigenspace in the three-quark carrier, hence20 adjoint copies. Coupling these to the auxiliary adjoint gives60 physical colour singlets after including its three spin components. Together with the20 zero-excitation states there are80 physical states. The coupling has rank20 because its Gram operator is `G²`, whose eigenvalues are strictly positive. Thus there are20 coupled two-level pairs and40 retained dark states at energyDelta.

For each eigenvalue lambda of `G²`, the exact pair eigenvalues are

`E_±(lambda) = [Delta ± sqrt(Delta²+4 g² lambda)]/2`.

The lower spin1/2 state has lambda4 and is lowered more than the spin3/2 state with lambda2. This is an exact finite coupling calculation, not a truncated perturbation formula or an inserted spin penalty. Under the conditional existing432 common-energy map, the resulting quartet-minus-doublet shifts are0.169633,0.087368 and0.044022 GeV at the three inherited kappas1/2,1,2. No observed splitting is claimed.

The zero-quark vector with the auxiliary in its vacuum remains unchanged. That fact does **not** establish that it remains the source-wide ground state: negative baryon shifts can change the vacuum. The interaction has not been joined to the two-cell source or re-grounded there. Its actual preparation, native execution, all observer records, physical field transfer and baryon mass comparison remain open.
