# Actual preparation and finite readout custody

The source is the original full native eight-cell h, its retained
preparation partner and the full thirty-mode static CAR reference with
energy4N. The existing physical two-copy preparation supplies the actual
conditional state and its trace-distance error. No additional ground
projector is used. The outer affine field and replicated matrix-unit
couplings are not part of this readout Hamiltonian.

The producer validates the frozen closures of
`model/native_charge_resolution`, receipt
`f22c9a266a43de8aaae776b1293894081c3dc20b7d8e5678e30ef6af993dd19e`,
and `model/affine_gauss_charged_clock`, receipt
`7c8f0214f1412e38bdc410c607d0a8f8765df323c3a33a778da923f3892e8613`.
Only source-independent phase-tail and rational-unitary identities are
reused from the latter. No affine energy, state or record probability is
transferred to the native source.

The complete native-resolution producer and the physical preparation
producer are freshly replayed. The code also replays every actual
hopping path and all59 PW119 blank-coordinate controls, then rebuilds
the exact centroid, conditioned preparation and finite-pointer bounds.
The full native source-plus-pointer matrix is not numerically executed.
Native bin admission, controlled finite matrix functions and state
synthesis remain declared apparatus services, with controller and gate
synthesis workspace separate from logical storage.

`research_reference.json` is the immutable research receipt from
`repair-law/work/CODEX_JOIN_QUANTUM/native_charge_centroid`, SHA256
`412146a491d580c14280b77739a4d997d04e24b021c1b8fd3d09982786c7faa4`.
The public replay compares every scientific field with it, excluding
only schema and path-custody metadata. Portable imports and command
entrypoints preserve the complete calculation.
