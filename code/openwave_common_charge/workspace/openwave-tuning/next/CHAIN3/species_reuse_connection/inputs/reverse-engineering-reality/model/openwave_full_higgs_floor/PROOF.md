# A charged-sector floor from all native Higgs channels

## Declared action and inherited source

The source uses the full fifteen-coordinate native matter module, two spin
coordinates and three analyzer characters. Its homogeneous periodic Weyl
kinetic term and unit neutral-Higgs background are those of the free
[charged branch construction](../openwave_lepton_mass_shell/PROOF.md).
All three native invariant Higgs tensors, up, down and charged lepton,
are present. Their normalization and signs are taken from the explicit
same-current intertwiner in `native_matter_attachment.py`.

The architecture sets Y_up=Y_down=Y_lepton=diag(e_k), where e_k are the
same positive analyzer effects. This equality is a declared extension of
the charged-lepton-only action. It is not forced by the gauge module;
the native invariant tensor lines also admit independently chosen family
matrices. The phase2/9, neutral background norm and electron calibration
retain their input ancestry. No quantum Higgs condensate or interacting
ground of a different source is transferred.

## Exact native kernel

Let L be the weak-doublet projector obtained from the actual weak Casimir.
For each native symmetric invariant tensor C_a at the neutral component,
put D_a=LC_a(1−L), with its original sign. The complete directed sum D
obeys

```
D^2=0, rank(D)=7,
P=D*D+DD*=P_charged, rank(P)=14,
Q(1−P)=0, |Q_ii|≤1.                              (1)
```

These are exact identities on all fifteen native coordinates. The seven
pairs consist of three up colours, three down colours and one charged
lepton. Their supports are disjoint and every nonzero coefficient is±1.
The single remaining internal coordinate is electrically neutral.
For K=(D+D^T) tensor diag(e_k),

```
K* K=P_charged tensor diag(e_k^2),
Q_45^T K+K Q_45=0.                               (2)
```

The spin-epsilon CAR creation kernel uses the conjugate of this all-left
symmetric matrix. Exact three-axis Clifford anticommutation in the full
180-row Nambu symbol gives charged excitation energies
sqrt(sum sin(ap_i)^2/a^2+e_k^2) and neutral energies|sin(ap)|/a.
Each character has seven charged Dirac channels per valley. All three
neutral family coefficients and all eight sine valleys are retained.
The side-four execution includes all64 momenta and11520 Nambu rows;
its48 zero-energy physical CAR modes are all electrically neutral.
The stated reference fills negative bands and leaves neutral zeros vacant,
retaining its2^48 degeneracy.

## Complete charged-sector lower bound

Let m=min_k e_k>0. Diagonalize the same charge-preserving quadratic source
in quasiparticle coordinates. Every charged occupation has nonnegative
excitation energy E_j≥m and charge |q_j|≤1. Every zero-energy occupation
has q_j=0. The chosen reference has total electric charge zero, since each
native pairing ground is a superposition of vacuum and a zero-charge
opposite-leg pair. On each simultaneous energy/charge occupation state,

```
H−E_ground=Σ_j E_j n_j ≥ m Σ_j |q_j| n_j
                          ≥ m |Σ_j q_j n_j|=m|Q|. (3)
```

Spectral calculus makes(3) an operator inequality on the complete Fock
space, including arbitrary references. A charged-lepton rest mode attains
the charge-one lower bound m. A three-quark colour determinant assembled
from same-family rest quasiparticles has energy3m and is not a zero-cost
charge-one counterexample. It is a free multiparticle state, without a
bound-baryon interpretation. Turning off either native quark channel
restores six charged internal zero modes; these failing controls verify
that a label change alone cannot establish(3).

## Finite implementation and complete records

The same address analyzer can act on each of the seven disjoint native
pairs, with all324000 enlarged CAR modes retained per site. Ancillary
number phases remove unwanted flag terms. Parity of the complete
weak-doublet leg reverses every pair term, and spatial inversion reverses
homogeneous kinetics. There are42s disjoint four-mode factors, or168s
active CAR coordinates. The existing exact fourth-order word and the
[factor bound](../openwave_factorized_compiler/PROOF.md) give complete
unhalved trace error at most

```
min(2,42s R^5|t|^5/(10n^4)), R=(7/3)(4/a+1), n≥R|t|.
```

This is a controlled mean-action implementation, with every ancillary
outcome and external reference retained. Its global control services
are declared; no autonomous implementation or local resource optimality
is inferred.

For each actual pair and character, a fixed native leg-creation
contraction applied to its pairing ground has success probability1/2.
An opposite idle charge reference compensates its charge. All four
occupation outcomes and the failed creation are recorded under the same
source evolution. The second spin orientation is related by the exact
spin-epsilon CAR intertwiner. Independent four-mode full-Fock checks
retain both spin orientations and prove the charge-one floor directly.

One electron energy anchor gives the same conditional e/μ/τ branch values.
The equal up/down/lepton family coefficients also give equal quark-family
spectra in this extension. These are model outputs that fail to reproduce
the unequal observed quark masses. The charged-sector floor supplies no
flavour selection, weak-decay lifetime, removal of spatial valleys or
physical interacting particle dictionary.
