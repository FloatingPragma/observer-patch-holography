# A quartic error bound for quantum hopping

The two-cell-per-axis source has an exact sublattice parity. On its
atomic ground band this makes the canonical effective Hamiltonian an
even function of hopping strength. The [proof](PROOF.md) gives a
controlled fourth-order error bound while retaining every charged
intermediate state. The [receipt](results.json) rebuilds the atomic
gap, full hopping norm and actual source parity.

The bound applies to the declared small-hopping interval. Longer
spectral lattices contain same-sublattice hops and require their own
symmetry argument.

```sh
python -m model.coherent_defect_wave.chiral_joint.atomic_even_remainder --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/atomic_even_remainder
```

The stable exact APIs are
`bounds.uniform_quartic_remainder(gap,tau,coupling)` and
`bounds.choose_quartic_coupling(gap,tau,second_order_gap)`.
The latter requires a separately proved second-order spectral gap.
