# Uncertainty of the actual compensated probes

The six family/spin addition maps of the [finite Higgs source](../interacting_extension/higgs_join/PROOF.md)
give noncommuting, charge-compensated observables on its complete 438- and
1596-dimensional sectors. Their uncertainty relation follows from the
[finite-state Robertson theorem](../../Lean/EventAlgebra/Robertson.lean).
The density matrix, observables and measurement services have the same
conditional interpretation as the interacting source. This is not an
identification of physical position and momentum.

## Exact source algebra

Write a creation map as `b:H0→H1`, including its odd compensating reference.
The combined operator is even and gauge neutral. Its exact normalized matrix
is a signed partial permutation with 384 entries. Put

```
P0=b†b,  P1=bb†,  P=diag(P0,P1),
X=[[0,b†],[b,0]],  Y=[[0,-i b†],[i b,0]],  D=diag(P0,-P1).
```

Direct multiplication gives

```
X²=Y²=D²=P,  [X,Y]=2iD,  {X,Y}=0.
```

For either `A=X` or `A=Y`, the three spectral projectors are

```
E_minus=(P-A)/2,  E_zero=I-P,  E_plus=(P+A)/2.
```

They are Hermitian orthogonal idempotents summing to the identity on all
2034 states. Their ranks are 384, 1266 and 384. The zero outcome includes
54 neutral blocked states and 1212 charged states outside the probe range.
Deleting it would give an incomplete instrument.

`uncertainty.py` reads the exact numerators and positive raw basis norms,
verifies their normalized entries exactly, and compares all six resulting
maps to the executed Higgs operators. Sparse arithmetic over the rational
complex numbers checks the displayed identities on the full carrier.
The upstream Higgs certificate supplies the Gauss-kernel and reference-CAR
identification; this module does not rederive those kernels.

## Every neutral density, including entangled inputs

Let `rho` be any density supported on `H0`, possibly a reduced density of a
state entangled with another observer. Define `w=Tr(rho P0)=1-<n_j>`.
Both quadratures have mean zero, variance `w`, and complete outcome law

```
(Pr[-1],Pr[0],Pr[+1])=(w/2,1-w,w/2).
```

Their commutator expectation is `2iw`. Consequently Robertson saturates:

```
Delta X Delta Y = w = |<[X,Y]>|/2.
```

The tracial neutral state gives `w=384/438=64/73` and probabilities
`(32/73,9/73,32/73)`. No ground-state approximation is used. For a general
state on the direct sum, the Robertson theorem applies;
the displayed saturation is asserted only for neutral-supported states.

The complete maps are `rho→E_a rho E_a`, tensored with the identity on any
untouched reference. `Probe.outcomes` returns all unnormalized poststates
for a vector or purification columns. Summing their reference Gram matrices
preserves the reference marginal. Measurement disturbs the joint state;
access to the simulator specification does not make these two instruments
one simultaneous sharp measurement.

This supplies an exact projective-instrument specification and executable
finite matrix actions. A native pulse compiler for these three-outcome
projective instruments is not asserted. The existing two-outcome Ramsey
instrument measures interference correlations and is a different instrument.

## Finite registers and canonical field limits

Finite matrices cannot obey `[Q,P]=i hbar I` for nonzero `hbar`: the trace of
a commutator is zero. More strongly, in dimension `d>0`,
`|Tr A|<=d||A||` implies `||[Q,P]-i hbar I||>=|hbar|`.
Thus a uniform whole-register canonical approximation cannot improve this
norm bound. It does not obstruct restricted-state or unbounded-operator
limits, or the state-dependent commutator above.

The [shared-state fluctuation construction](../predictions/README.md#shared-state-fluctuation-field)
has a quantitative canonical limit for independent prepared copies. Its
copy count is not the physical memory-capacity parameter. The
[ordered-spin construction](../spin_canonical_field/README.md) supplies
another regular canonical field limit from its specified ground. Neither
canonical algebra nor finite uncertainty alone identifies Lorentz particles,
derives a physical Planck scale, or proves a scattering interpretation.

```sh
python -m model.particle_identification.uncertainty
python -W error -m pytest -q model/particle_identification/test_uncertainty.py
```
