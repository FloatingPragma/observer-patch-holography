# Complete neutral quantum-Higgs spectroscopy

## Same source and gauge-neutral operator

The Hamiltonian, neutral Gauss ground and clock are exactly those of
[`periodic_spectrum`](../periodic_spectrum/PROOF.md): two cells, two oriented
six-state quantum links, two five-state quantum Higgs memories and all twelve
active lepton CAR modes. The native quark spectators are empty in this
reducing `y_u=y_d=0` control. The declared coefficients are
`y_e=omega=kappa=1`, and the oriented hopping magnitude is `1/4`.

At site one and spin zero, use the actual operator

\[
 J=\phi_0 c_{0,11}-\phi_1c_{0,10},\qquad
 \phi_k=|0\rangle\langle1+k|+|3+k\rangle\langle0|.
\]

All twelve native gauge-current commutators vanish, as proved in
[`higgs_composites`](../higgs_composites/PROOF.md). The operator is fermion odd
and noncanonical. Its neutral gauge character identifies a source sector;
it does not identify an observed neutrino.

The unique neutral source ground has `L=2`, `2Sz=0` and translation involution
`+1`. Applying `J` reaches `L=1,2Sz=-1`; applying its adjoint reaches
`L=3,2Sz=+1`. Both targets have zero boundary charge. Each complete physical
sector has dimension 298 and splits into two 149-dimensional translation
sectors with eigenvalues `-i` and `+i`. The source ground block has dimension
251. All transfer maps are exact matrices in the inherited physical Gram
metrics, and their raw scalar/CAR images are checked against the full Gauss
embeddings.

## A complete finite instrument

On the complete two-active-mode CAR space times the scalar memory,
`A=J^dagger J` has eigenvalues `0,1,2` with multiplicities `14,4,2`.
The same spectrum holds for `JJ^dagger`. This 20-state calculation retains
all five scalar alternatives. The identities also hold on the complete
native thirty-mode Fock space, because the independent normal-ordered CAR
polynomials retain every spectator occupation.

For either setting `J_kind=J` or `J^dagger`, define

\[
 K_{\rm success}=J_{\rm kind}/2,\qquad
 K_{\rm failure}=P_0+\frac{\sqrt3}{2}P_1+
                     \frac{1}{\sqrt2}P_2,
\]

where `P_r` is the exact spectral polynomial of `A_kind` at eigenvalue `r`.
The executable checks projector orthogonality and completeness, positive
failure eigenvalues, and

\[
 K_{\rm success}^\dagger K_{\rm success}
  +K_{\rm failure}^\dagger K_{\rm failure}=I.
\]

Both settings are selected with probability one half, with that choice
retained in the records. The full instrument therefore has four branches,
and each successful energy/phase record has weight
`|<n|J_kind|0>|^2/8`. The failure state is the unnormalized positive state
`F rho0 F/2`; its energy and phase are not measured in this protocol.
The finite apparatus retains the complementary fermion parity, lepton
number and spin. It carries no gauge charge. A contraction has a finite
unitary dilation, so this specifies a complete finite apparatus. That
controller is a declared intervention, rather than part of the idle source
Hamiltonian. On the complete 40-state matter/scalar/reservoir memory, the
code constructs `U=[[F,-K^dagger],[K,F_reverse]]` and verifies both unitary
products exactly. It also checks conservation of total source-plus-reservoir
fermion parity and lepton number; the selected spin component has the
corresponding conserved spin projection. This dilation retains every
occupation and scalar alternative.

The failure operator is constructed before translation restriction. Its
image has a nonzero component of opposite translation sign. Computing the
square root after compression to the 251-dimensional ground block would
omit that component. The exact projected completeness check uses the full
unresolved source physical sector for the failure map.

## Complete spectral probabilities and errors

Let `T_kind` denote the raw composite transfer into the unresolved target.
For its translation involution `S`, the exact maps are constructed from
`(I+sign*S)T_kind/2`, using the kernel basis for the chosen eigenspace.
Their Gram effects sum to `T_kind^dagger G_target T_kind`. Each complete
energy projector, including zero-weight alternatives and degenerate
subspaces, contributes one probability interval in the receipt.

The normalized ground approximation and all spectral-projector bounds are
inherited from the same Hamiltonian's complete exact-residual certificates.
If `delta_ground` bounds the ground trace distance and `delta_P` bounds the
projector error, the conservative joint-probability error
`(delta_ground+delta_P)/2` applies because the setting and success
normalization contribute `1/8` and `||J||<=2`. The sharper exact bound is
`||J||=sqrt(2)`. All reported rational probability intervals are rounded
outward at denominator `2^27`.

The ground preparation itself succeeds with probability `1/3424`. Its
failure `3423/3424` is retained. The receipt reports both the law conditional
on that recorded success and the original experiment law obtained by
multiplying all subsequent branches by `1/3424`. No successful spectral
branch is renormalized. Both complete probability enclosures contain one.

The first energy line is populated in all four setting/translation branches.
Its source-unit transfer interval is
`[2681604835/4294967296,2681605091/4294967296]`. Each original joint weight
conditional on the recorded ground preparation lies in
`[1206821/33554432,4827285/134217728]`, approximately `0.03596607`.
The factor of one half for the setting and one quarter from the success
operator are included. Both positive lower endpoints are verified from the
same source's complete spectral bounds.

## Exact local-parity selection

The conserved operators

\[
 C_x=(-1)^{N_{\rm matter,x}}
      \prod_{e\ {\rm incident\ to}\ x}(-1)^{N_{\rm link,e}}
\]

commute with the complete source Hamiltonian and all gauge currents. A local
`J_x` changes only the star character at its own cell. The unique ground is
an eigenstate of each `C_x`. Consequently

\[
 \langle0|J_x(t)J_y^\dagger|0\rangle=0\quad(x\ne y)
\]

at every source time. This selection rule holds even though `J` is a gauge
singlet.

In an odd total-fermion sector, `C_0 C_1=-I`. Translation exchanges the two
stars, so `C_0` anticommutes with the translation involution. The executable
checks both `[C_0,H]=0` and this anticommutator in each full 298-dimensional
target. It follows that the `+i` and `-i` spectra coincide. Since a local
composite applied to the ground has a definite star character, its energy
weights in the two translation sectors coincide as well. The phase pairing
is an exact internal selection rule, rather than a relativistic dispersion
relation. Retained link dressings or a separately declared Hamiltonian can
have different transport questions; no exclusion of those constructions is
asserted here.

## Scope and custody

This result is a complete gauge-neutral composite measurement on one finite
interacting source. It supplies no moving neutrino state, observed flavor
assignment, physical rest mass or conversion to GeV. All apparatus and
Hamiltonian parameters have construction ancestry; no empirical mass datum
is used.

The receipt pins package code and proofs, transitive local Python imports,
the consumed Higgs-composite and periodic-spectrum receipts, and every input
pinned by those receipts. Their custody is validated before generation or
replay. `--check` independently recomputes the operator identities and
complete record law without writing the receipt.

The projector bound uses the inherited rational residual and Gram formula.
For an isolated contiguous energy cluster, ordered interval endpoints and
containment of every approximate energy imply that its closest outside
interval is an adjacent interval. The implementation verifies these premises
on the complete spectrum, then uses the nearest left upper or right lower
endpoint. This gives the identical separation to a full-complement scan;
representative singleton, degenerate and endpoint clusters are checked
against that scan. The optimization changes no probability bound.
