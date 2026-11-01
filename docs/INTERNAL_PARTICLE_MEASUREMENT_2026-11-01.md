# An internal observer measures the protected excitation response

This report records the first one-time native readout. The
[subsequent extension](RESOLVED_SPECTRUM_AND_RELATIVISTIC_SECTOR_2026-11-01.md)
now resolves all eight source-supported weights through certified long-time
emulation and separately constructs a massless relativistic limit. The
prospective statements below describe this earlier packet's scope.

The [native readout construction](../model/protected_particle_readout/README.md)
connects the protected cube response to an observer's pointer records in
the specified computed universe. The connection concerns the local probe's
eight-dimensional cyclic response and its elastic reference component.
One probe does not identify every state in the larger protected sector.

The simulator's source encoding, initialized blank state, field law,
repair operations, pointer and repetition service are explicit inputs.
Their use is consistent with reading the simulator specification. The
construction does not import an external particle name, physical mass
scale or laboratory observation.

This continues the five-issue work in
[the previous report](THREE_DIMENSIONAL_PROGRESS_2026-11-01.md), directly
advancing [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139)
and [#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143).

## What the observer does

1. Prepare the same signed eta-pair reference R on the sixteen field modes,
   with number six and Hamiltonian energy three.
2. Prepare a seventeenth pointer and apply the local number phase
   F=exp(i pi n_origin/2) controlled by that pointer.
3. Evolve both arms under the original cube Hamiltonian for one model-time
   unit. Two half-time segments each use 32 symmetric splitting steps.
   Apply all five declared C repairs between the segments and all five E
   repairs after them.
4. Apply the inverse controlled phase and one of two pointer analyzer
   settings, then export the pointer through its complete six-outcome map.

The full field/pointer success state has 131,072 amplitudes. The execution
does not project it into a number or spectral sector. Local native words
are composed before being applied to this full state; their complete
first-failure maps are retained and reconstructible. All source aborts,
native repair terminations and final pointer outcomes enter the probability
accounting. No success state is renormalized.

## Which quantity those records measure

For the ideal model protocol, the two pointer settings give

```
C(t) = 1/16 + (15/16) sum_j w_j exp[-i (E_j-3)t],
p_zero(phi,t) = [1 + Re(exp(i phi) C(t))]/2.
```

The eight energies and positive weights are the
[same exact spectral values](../model/particle_identification/README.md)
derived from the original Hamiltonian. Their relation to pointer records
is now a specified source/detector map. It can be checked independently
through the field matrix, the spectral measure and the native circuit.

At t=1, the ideal probabilities are approximately 0.9942494663 for phi=0
and 0.5019797169 for phi=pi/2. The full-Fock splitting contribution to each
probability-error bound is exactly 629/1048576. Literal source-control
errors are added separately; together the bound remains below 0.000600.
This is an exact mathematical source/evolution bound. Floating execution
has a separate numerical admission.

The complete native execution and the reproducible simulated records give:

| Pointer phase | Native pointer-zero probability | Zero records out of 16,777,216 | Observed frequency |
| --- | ---: | ---: | ---: |
| 0 | 0.994249701120 | 16,680,454 | 0.994232535362 |
| pi/2 | 0.501979889450 | 8,423,380 | 0.502072572708 |

The counts are two seeded binomial samples of the unconditional event
“final pointer zero”; the complement includes every other complete output.
They are simulated records, not laboratory data or independent evidence
for the simulator law. With the specified independent-reset service,
Hoeffding's bound gives a simultaneous sampling allowance of 1/2048 with
failure probability below one percent. Both recorded frequencies lie in
the certified intervals enlarged by this allowance and the separately
reported host-numerical admission. The computed complete-output error
against the original-H comparator is about 2.27e-6, below the theorem's bound.

The fixed-number two-arm protocol measures energy differences E_j-3.
Adding a scalar to H, or a common number chemical potential, leaves its
pointer probabilities unchanged. The absolute energies use the supplied
Hamiltonian's chosen origin; they are not an independent measurement of
absolute rest energy. The reference is excited, and the response is
neutral under model particle number.

## Identifying the interaction and resolving the lines

Removing the contact interaction gives a useful exact control. Bipartite
occupation parity reverses the free hopping Hamiltonian and fixes both R
and FR. Since the free Hamiltonian annihilates R, its corresponding C(t)
is real for every time. Its second pointer probability is therefore
exactly one half. This is an exact ideal comparator; an additional native
zero-contact run is not executed. The interacting imaginary-quadrature fringe is large
enough to distinguish these two specified models with the stated source
and sampling errors.

Two pointer settings at one time do not reconstruct eight unknown lines.
A separate exact resolution theorem uses both settings at model times
0, 39, 78, ..., 273. For the eight source-specified frequencies, the
sampling matrix has determinant magnitude greater than 1400 and inverse
infinity norm below 2744/525. It therefore determines their weights
uniquely, with an explicit error-amplification bound.

One sufficient schedule uses 77,309,411,328 independent shots across those
sixteen settings, alongside a specified systematic-error bound at every
setting. That long-time experiment and its accuracy requirements are
**not executed here**. Finite-precision matrix construction and inversion
would also need their own error allocation. The theorem is about the
given eight-frequency support; it does not identify an arbitrary spectrum
from finitely many samples.

The new result is an operational identification of a finite, protected
neutral response inside the model. Controlled spatial refinement,
propagating particle poles, a physical vacuum and a relativistic
interacting field remain separate obligations. The moving pair's model
inertial mass and binding energy from the previous pass are not silently
assigned to these cube modes.
