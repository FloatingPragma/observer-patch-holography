# Calibrated CAR propagation on a native-seam ring

This finite spatial extension uses all twelve native ports, with the supplied
Hamiltonian cycle `(0,1,5,4,2,3,9,8,6,7,10,11)`. Every ring edge is checked
against the actual native seam graph. Ring position and orientation are
design choices. Native balancing remains a separate instrument.

At unit spacing and speed the Hermitian difference operator has
`D[r,r+1]=-i/2` and eigenvalue `sin(k)`. The action is
`h=-D tensor sigma_z tensor I45`, with the same local calibrated pairing
`Delta=I12 tensor Delta90`. Its exact algebraic mass/kinetic anticommutation
gives `B(k)^2=B(0)^2+sin(k)^2 I`. The numerical check retains all twelve
momenta and all 2160 Nambu directions. Both zero-momentum and pi doublers
remain, including the 24 zero Nambu directions of the massless neutrinos.
Three tensor-product copies of D and the three Pauli factors give the same
identity with the sum of three squared sines and eight nodes. That larger
construction is an algebraic observation, not a source execution here.

The full carrier has 1080 CAR modes and all `2^1080` occupation states.
Same-port collisions are retained: hopping out of the earlier hard-core
two-particle sector can create them. Likewise, kinetic evolution need not
preserve the static onsite positive-quasiparticle number. The calibrated
pairing graph decomposes the action into 24 two-mode, six six-mode and six
one-mode internal components. Each complete spatial component has at most
72 CAR modes, so the existing source Gaussian service executes it without
allocating an exponentially large Fock vector. This is an exact invariant
decomposition of the specified quadratic action, not a mode selection.

All three native current values drive the entire action for
`tau_eff=tau*(1+j/4)`, agreeing with the common-clock convention. The executed
law is the declared symmetric product of local terms. Complete four-column
Fock source gates and all routing FSWAP signs determine its Majorana map;
all local scalar phases and resources are retained in the source catalog.
The packet separately reports the numerical discrepancy from the exponential
of the summed Hamiltonian. It does not call a finite product exactly equal
to that exponential. Its static dispersion refers to the sum action.
Every local source matrix is compared with an independent full four-state
CAR Hamiltonian, retaining scalar phases that the Majorana map alone cannot
see. The local numerical admission composes to a full-Fock operator and
probability bound. A separate conservative Taylor bound compares the ideal
symmetric product with the static sum; it is explicitly allowed to reach
the uninformative unitary bound 2. Neither bound upgrades the floating local
checks to outward-rounded arithmetic.

A faithful Gaussian preparation gives each work mode occupation 1/4 and
retains a CAR purifier for each mode. A second preparation raises one mode
at native port zero to 3/4. Both preparations are actually source-executed.
Their propagated occupation difference displays spatial transport. A bond
current is read by an actual source rotation followed by complete neutral
occupation pointers; both outcomes and the full conditional covariances,
including every purifier, are retained. Current is the difference between
two separate complete pointer experiments. There is no postselection.

The complete native balancing extension outside the old two-particle domain
is a separate supplied instrument. This package supplies the field action
and its source propagation on the full carrier for composition with that
instrument; it does not invent that extension by applying an exterior power
to a nonunitary Kraus map. Static per-site Bogoliubov coordinates connect the
bare fields used here to the earlier colored internal chart.

This is a one-dimensional finite witness with calibrated matter inputs. It
does not establish a three-dimensional continuum, remove the doublers, or
identify native ports with observed spacetime. Numerical execution is not
an outward-rounded or physical hardware error certificate.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.defect_field_propagation --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/defect_field_propagation
```
