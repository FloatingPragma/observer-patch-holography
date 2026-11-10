# Native source custody and exact replay

The native constituent h and its reference/Gauss/G1 inputs are pinned through
`model/native_abelian_field_join`, receipt
`f6aadcd79eec889de0c3243925c950e482bd456ea5c22fcab6a8406978472327`.
This readout uses the original native h and its static charged reference;
the outer affine coupling is not part of this readout Hamiltonian.
No native Hamiltonian coefficient, frozen spectral selection or outer
coupling is changed by this package. It derives nested energy readouts
on the complete native h and an explicitly different refined selector.

The producer validates all inherited source pins, freshly replays the
native join, rebuilds the actual local CAR target polynomial and positive
weights, isolates all nine charged roots and the neutral ground root,
and evaluates the full-source probability and neutral-energy inequalities.
It requires no numerical eigensolver. It does not enumerate or truncate
the full native charged spectrum.

`research_reference.json` is the immutable research receipt from
`repair-law/work/CODEX_JOIN_QUANTUM/native_charge_resolution`, SHA256
`9a88dc18a69a9d12eb7058d1bd62b42f015f7e0531147d3a885e82559ab6a086`.
The public replay compares every scientific field with this receipt,
excluding only schema and path-custody metadata. Portable imports and
command entrypoints do not alter the source calculation or theorem.
