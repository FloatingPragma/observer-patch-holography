# Modular local-coarse C/E instrument prototype

This separately named candidate module implements a complete connected one-load Gauss
sector with Z6 half-seam registers. It leaves the hard-wall, fully monitored
`JointRecords` model unchanged. It uses the existing triangle development
configuration. No reserved acceptance case is executed.

The coherent Hamiltonian is the lifted port/transit generator H_m. Register
shifts wrap modulo 6; there is no wall rejection. Gauss eigenvalues are exact
residue classes modulo 6. Integer-valued representatives of those residues do
not have a fixed integer Gauss vector. The modular calculation supplies no
isometry to the complete integer model and no general long-time approximation
bound. The triangle sector contains 36 states; it is the complete connected
Gauss sector, not the full tensor product of all unconstrained registers.

For a drawn seam (a,b), the classical native rule determines current j and
pre-repair contrast d for each port position and fair coin. The quantum lift
coherently groups sources with the same local C record `(coin,j)` or local E
record `(coin,d)` into one Kraus operator. C/E laws have equal probability.
Native arithmetic and classical record content alone do not uniquely determine
this grouping as a quantum instrument; it is a supplied local-coarse extension.
No source position or record value is added to the outcome ledger.

Let T_ab be the completed transit translation which moves a load a to b and
shifts the two half-seam registers by (+1,-1), modulo 6. Conditional on record
law C, the coin-1 operators are `(I-P_b)/sqrt(2)` and
`T_ba P_b/sqrt(2)`; coin-0 operators are `(I-P_a)/sqrt(2)` and
`T_ab P_a/sqrt(2)`. Conditional on E, each coin resolves P_a, P_b and
P_elsewhere with the same completed translations for moving outcomes.
The code includes an additional factor `1/sqrt(2)` for the law choice.
Each moving partial translation is injective; records are shifted and are
not measured. Therefore sum K^dagger K=I, and the channel is completely
positive and trace preserving.

Transit loads are included in "elsewhere" and left untouched. No transit
indicator is emitted. In particular a coherent superposition of distinct
transit states is unchanged by a local commit; a fully monitored control
destroys that coherence. A null C outcome preserves coherence between the
nonmoving endpoint and all other positions, including their record fibers.
E resolves both endpoints separately. These are operationally different
channels with the same one-step classical position transition.

`commit` accepts a specified seam or an explicit uniform seam mixture. It
returns every possible local record outcome, including zero-probability
outcomes. It has no boundary or global-source outcome. The supplied transit
policy does not refine the recorded instrument.

The continuous scheduling API implements the exact finite generator

```text
L(rho) = -i [H_m,rho] + r sum_e (C_e(rho)-rho).
```

Here r is the supplied rate of each seam; the total attempt rate is M*r.
Adding seams does not renormalize the rate assigned to an existing seam.
`evolve_per_seam` exponentiates this finite generator by sparse Krylov action.
The separate `step` method reproduces the declared one-event Bernoulli split
for comparison. Setting `r=p/(M*t)` identifies their first-order rates; it
does not assert equality of their finite-duration channels or a fast-commit
theorem. No clock-convergence sweep or many-load generator is studied.

The executed triangle uses mass 0.23, time 0.37, split probability 0.29 and
eight steps. The local-coarse and fully monitored channels are compared in
the same modular sector to isolate the instrument effect. The cutoff-1
hard-wall fully monitored model is retained as a distinct control, which
combines a regulator change with a monitoring change. A single continuous
interval uses duration 0.37 and r=0.29/(3*0.37). This is development evidence,
without an official cell or physical species/charge-unit identification.

Run the prototype and tests:

```sh
python -m oph_openwave.modular_native --output /tmp/C4-modular.json
python -m pytest -q tests/test_modular_native.py
```

The independent K01/K02 audit is complementary evidence. Finite numerical
position-marginal agreement with an integer calculation is not equality of
record spaces or proof of all observables. This prototype does not adopt the
uncertified numerical cutoff-bound values or execute the reserved bridge
acceptance. A live construction change requires its own explicit state,
register, instrument and rate declaration.

The component-rate control prepares one load on an edge. Adding an unoccupied
disconnected edge leaves the connected prepared sector and coherent generator
unchanged. The remote seam channel is the identity on this sector, so its
contribution r(C_remote-I) vanishes. The reduced occupied-component state
therefore agrees at the same supplied per-seam rate. This is a one-load
component test, without a many-load generator claim.

For an arbitrary external reference R, completeness implies
`Tr_system sum_i (K_i tensor I_R) rho (K_i^dagger tensor I_R) = Tr_system rho`.
The tests retain an entangled reference qubit and verify this identity, trace
and positivity. This proves complete positivity with the displayed Kraus
operators; it does not identify a physical reference apparatus.

Portable implementation: [modular_native.py](../oph_openwave/modular_native.py).
Regression controls: [test_modular_native.py](../tests/test_modular_native.py).
The class is separately named `LocalCoarseModular`; it does not replace
`JointRecords` or select a record energy. State tuples are (position, residue
tuple); Gauss values are obtained through `gauss_residue`, not integer Gauss.
`commit` returns coarse local record probabilities and does not expose global
source/record readout. `step` and `evolve_per_seam` name distinct schedules.
