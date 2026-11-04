# Radial defect excitation catalogue

This catalogue calculates excitations of the supplied radial lattice model.
It associates a protected zero-mode doublet and two actual interacting
charged states with a defect fixed at the origin. Observed particle species
and absolute physical masses are unassigned. The rows describe distinct
source domains; they do not count five observed species.

The finite source has 27 sites, two fermion components per site,
`a=v=1`, radial mass parameter `m=64`, and centered contact `U=1/65536`.
Its interacting reference is the actual `N=27` ground doublet, averaged
with equal weights. The complete Fock complement is retained in the proof;
the `2^54` matrix is not numerically diagonalized.

| Excitation | Energy in model units | Multiplicity | Position and momentum |
| --- | --- | --- | --- |
| Protected one-body zero doublet | 0 | 2 one-body modes | Localized at the supplied defect origin |
| Interacting removal state, `N=26` | [0.0000076231015466583, 0.0000076297367371032] above ground | 1 state; pole residue rank 2 | Central spectral weight at least 99.959000283954%; every Fourier probe has positive weight |
| Interacting addition state, `N=28` | [0.0000076231015466583, 0.0000076297367371032] above ground | 1 state; pole residue rank 2 | Same localization bound; signed pole opposite to removal |
| Negative exterior one-body branch | `-sqrt(p²+64²)` in the continuum | 8 rows per wave number | Exterior wave packets, with explicit residual bound |
| Positive exterior one-body branch | `+sqrt(p²+64²)` in the continuum | 8 rows per wave number | Exterior wave packets, with explicit residual bound |

The charged states have equal positive excitation gaps `g`; their retarded
poles are `-g` and `+g`. Every other charged excitation is at least
`8126437/131072 = 61.999794006347656` above the ground.
The two zero modes generate the low-energy interacting core; they are not
additional particle species beside its charged states. The index protects
the one-body zero sector under its stated symmetries, not the value of `g`
under arbitrary changes of the interaction.

## Defect association and localization

The [bounded radial index proof](../../model/bounded_defect_index/PROOF.md)
gives Fredholm index two for the inversion-even to inversion-odd block.
The grading is spatial inversion. Compact self-adjoint perturbations that
anticommute with this grading preserve the index. This is not an
identification of texture winding with electric charge. The stronger
[radial gap theorem](../../model/radial_defect_gap/PROOF.md) proves exactly
two zero modes in its strong-mass domain. The
[interacting core proof](../../model/radial_interacting_core/PROOF.md)
derives the charged states from this same core, including the virtual
complement. Their charges `-1,+1` are conserved source fermion-number
differences; physical electric charge is unassigned.

For either actual low-pole residue `M`, define the spatial spectral-probe
distribution `P_x=Tr(M_xx)/Tr(M)`. This measures where that pole can be
excited by a local fermion probe. It is not a particle trajectory.
The theorem gives `||M-Z/2||<=3t`, where `Z` is the rank-two one-body zero
projector and `t=27/16252874`. Both `M` and `Z` have rank two, so eigenvalue
perturbation gives `1-6t<=Tr(M)<=1+6t`. At the origin,
`Tr(M_00)>=alpha-6t`, where `alpha>=7688/7691`. Consequently

```
P_origin >= (alpha_lower-6t)/(1+6t)
         = 62475424685/62501049938.
sum_x |x|² P_x <= 3a²(1-P_origin) <= 0.0012299914813633 a².
```

For the 27 orthogonal Fourier probes at `p_j=2*pi*k_j/(3a)`,
`k_j in {-1,0,1}`, put `r=alpha_lower-101/1000`.
The quaternionic zero frame satisfies `W_0=sqrt(alpha) I2` and
`W_x*W_x=rho_x I2`. Cauchy-Schwarz bounds the sum of the other 26 frames:
`s_min(F_p*W)>=(sqrt(alpha)-sqrt(26(1-alpha)))/sqrt(27)>=r/sqrt(27)`.
The exact inequality `26(1-alpha_lower)<(101/1000)^2` justifies the
rational choice of `r`. Thus each diagonal low-pole block has eigenvalues
at least `r²/54-3t`, and its normalized two-component probability is

```
P_p >= (r²/27-6t)/(1+6t) >= 0.029897135012027.
```

The same singular-value argument for two different Fourier labels proves
nonzero off-diagonal momentum blocks. The full radial Hamiltonian has no
translation symmetry: translating it moves its unique zero onsite mass
block. Each Fourier probe sees the same localized gap `g`; its probe label
does not define a moving-defect dispersion or a conserved momentum.

## Bulk rest energies and units

The [exterior packet construction](../../model/source_mass_shell/PROOF.md)
uses the same infinite radial one-body operator. Its finite-spacing energy
is `E=+/-sqrt((v/a)^2 sum_j sin(a*p_j)^2+m²)`. Compact normalized packets
at `(0,0,R)` with `L` cells per axis have residual at most
`3*pi*v/[a(L+1)]+4*sqrt(3)*m*L/R`. Their continuum limit has
`E²-v²p²=m²`. All sixteen folded rows are retained, corresponding to four
Dirac copies in the pure-spin bulk family. These multiplicities have no
observed-species assignment. The bulk rest energy 64 is supplied through
the source mass parameter. It is distinct from the finite interacting gap.
The fixed-mass continuum refinement leaves the finite interacting core's
strong-mass proof domain; no interacting relativistic pole follows from it.

With the explicit identification `v=1` with physical light speed,
`length_unit=299792458*T` metres, where `T` is seconds per Hamiltonian time
unit. The [unit conversion](../../model/source_mass_shell/CALIBRATION_PROOF.md)
gives `energy_unit=hbar/T` and the charged energy interval
`[5.0176165870614E-30, 5.0219839487489E-30]/T GeV`. A positive clock gap identified with cesium in the same
generator would determine `T`. No such input is supplied here. The
catalogue therefore leaves absolute GeV masses and the localized defect's
rest energy null.

## Reproduction

Run from the repository root:

```sh
.venv/bin/python code/particles/derive_radial_excitation_catalogue.py --write
.venv/bin/python code/particles/derive_radial_excitation_catalogue.py --check
.venv/bin/python -m pytest -q code/particles/test_radial_excitation_catalogue.py
```

The [exact JSON catalogue](runs/catalogue/radial_excitation_catalogue.json)
retains rational endpoints, source hashes, all five entries, and their
domain restrictions. Its checker recomputes the catalogue and verifies the
parent source manifests. The owning source packages provide independent
proof controls and cold replays.
