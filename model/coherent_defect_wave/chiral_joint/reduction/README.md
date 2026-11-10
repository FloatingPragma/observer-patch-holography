# Exact sectors for a two-cell charged response

This enumeration identifies a small computational target within the declared
finite construction. It does not change its frozen Hamiltonian receipt.
The control has two cells, a single z-directed quantum link, both spin
components of the native e^c and L module, and a five-state quantum Higgs
memory at each cell. Set both local charged-lepton Yukawa couplings to one
and the quark Yukawa couplings to zero. Empty quark modes and unexcited
colored links form reducing spectators. Each fixed boundary-charge sector
starts with `2^12 * 6 * 5^2 = 614400` active source states.

```sh
python -m model.coherent_defect_wave.chiral_joint.reduction.sector_counts --summary
python -m model.coherent_defect_wave.chiral_joint.reduction.sector_counts
python -m pytest -q model/coherent_defect_wave/chiral_joint/reduction
```

The second command emits every block as JSON. All diagonal charge weights
are read from the actual native, scalar and quantum-link current matrices.
The fixed charge minus one is a compensating boundary apparatus at the
right cell, rather than an additional particle species. Its retained fixed
value contributes no additional dimension within that sector.

| Right boundary charge | Diagonal Gauss support | Number of (L, 2Sz) blocks | Largest block |
|---|---:|---:|---:|
| 0 | 4140 | 33 | 642 at L=2, 2Sz=0 |
| -1 | 2496 | 25 | 354 at L=1, 2Sz=+1 or -1 |

Here `L=N_L-N_ec` is the conserved lepton-number label of this control.
Its one-axis kinetic term and spin-singlet Yukawa terms preserve the source
spin projection `2Sz`. These labels are algebraic conserved quantities of
the stated finite control; they do not establish relativistic particle spin.

The displayed dimensions impose local hypercharge and weak T3 constraints
only. They **are not the dimensions of the full physical Gauss kernels**.
The next operation is the joint kernel of the actual local weak raising
operators, mapping the zero-weight domain into the adjacent weight sectors.
Compressing their output back to zero weight would make them vanish and
would falsely count nonsinglets as physical. On a zero-weight vector in a
finite unitary SU(2) representation, vanishing of the raising operator also
implies vanishing of the lowering operator, by their commutator and norm
identity. The exact full weak-Gauss projection is not computed here.

After that projection, construct and check the complete Hamiltonian on
every physical conserved block. Its full ground projector or its physical
Gibbs state must precede a claimed stationary charged energy response.
The vacuum cyclic block alone need not contain the ground. A compensated
e^c probe, its same-Higgs conversion channel, and complete energy and site
records then test a finite charged excitation using the same source and
clock. The Hamiltonian, ground, charged response and mass are not computed
by this enumeration.

The open gauge graph supplies no proved conserved lattice momentum.
Fourier smearing may specify a detector but does not establish a momentum
eigenstate. A translation test and controlled boundary, volume and cutoff
comparisons are required before interpreting an energy response as a
relativistic dispersion or a physical rest mass.
