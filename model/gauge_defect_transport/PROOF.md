# Complete local transport of a charged defect

The colored native repair construction gives a genuine fermionic exchange
sign and a global internal representation. Independent local gauge changes
also act on the transported flux. This calculation supplies a complete
local transport instrument on the existing fifteen-component matter
representation and sixty-state quantum link. Those extra registers and
their representation are declared inputs. The result neither selects the
extension from the unextended repair law nor identifies a physical particle.

## Literal link and endpoint actions

For each of the five matter blocks `r`, let its dimension be `d`. The link
is the complete space

```math
\mathcal L=\mathbb C|0\rangle\oplus\bigoplus_r V_{r,x}^*\otimes V_{r,y},
\qquad \dim\mathcal L=60.
```

Use either the literal raising coordinates or the coordinate completion in
[`QUANTUM_LINK_PROOF.md`](../particle_gauge_band/QUANTUM_LINK_PROOF.md).
The latter additionally assumes its compact-group rotor attachment. For
both choices the combined coordinate on `V_r tensor L` obeys

```math
A=\sum_{ab}|a\rangle\langle b|\otimes U^r_{ab},\qquad
(A^\dagger A)^2=d A^\dagger A.
```

Consequently `B=A/sqrt(d)` is a partial isometry. Set `P=B†B` and `R=BB†`.
The code independently reconstructs all link entries from the exterior
currents and checks both projector identities. It keeps all sixty link
states, including the flux states on which transport is unavailable.

For an anti-Hermitian current `g`, write its matter action as `M=g tensor I`
and the two link actions as `L` and `R_link`. Before a move from `y` to `x`,
the endpoint generators are `(L,M+R_link)`. After the move they are
`(M+L,R_link)`. The exact link covariance equations give

```math
(M+L)B=BL,\qquad R_{\rm link}B=B(M+R_{\rm link}).
```

Thus `B` intertwines the independent endpoint actions. Multiplication by its
adjoint proves that `P` commutes with the source action and `R` with the
target action. The executable checks these three equations at both
endpoints for all twelve literal currents and all five matter blocks,
with both coordinate choices. Local covariance is stronger than commuting
with the same internal transformation at both endpoints.

## Complete outcomes and reversible transport

On the two-position space put

```math
T=|x\rangle\langle y|\otimes B,\qquad
P_y=T^\dagger T,\qquad P_x=TT^\dagger.
```

The supports are orthogonal, and `T²=0`. Every endpoint generator is block
diagonal in position, so the preceding intertwining equations give
`[G_x,T]=[G_y,T]=0`. The two Kraus operators

```math
K_{\rm transport}=T,\qquad K_{\rm hold}=I-P_y
```

satisfy `sum K†K=I` exactly. Each commutes with the total endpoint action.
The reverse instrument uses `T†` and `I-P_x`. All hold outcomes remain in
the output. Omitting them fails trace preservation on every state outside
the directional source support. Finite sequences, including record-selected
choices of direction, are therefore complete and gauge covariant on this
same carrier, also in the presence of arbitrary untouched references.

The Hermitian operator

```math
W=T+T^\dagger+I-P_x-P_y
```

has `W²=I`. It exchanges the source and target supports and acts as identity
on their complement. A continuous family of gauge-invariant unitaries is

```math
U(c,s)=I+(c-1)(P_x+P_y)-is(T+T^\dagger),\qquad c^2+s^2=1.
```

This is an implementation conditional on the chosen parameter. Covariance,
normalization and the supplied representation do not determine that
parameter, a clock, a scalar potential or a Hamiltonian coupling. The
existing [Gauss construction](../native_quantum_gauss/PROOF.md) supplies
explicit neutralized charged preparations and complete gauge-invariant
detectors on this link representation.

This parameter freedom has an exact gauge-invariant experimental witness.
Adjoin a declared reference in the dual representation at `y` and prepare
`d^(-1/2) sum_b |y,b,0_link> tensor |bar b>`. Every endpoint Gauss generator
annihilates this state. Its image under `T tensor I_reference` is another
normalized endpoint singlet, orthogonal to the initial state. Measure the
target-support projector `P_x` and its complement after `U(c,s)`. The two
unconditional probabilities are `s²` and `c²`. The exact choices
`(c,s)=(3/5,4/5)` and `(5/13,12/13)` therefore give `16/25` and `144/169`
on the same preparation, reference and detector. Both choices satisfy all
of the proved finite covariance and completeness identities. The source
must specify or select the operation parameter to predict this experiment.
No measured mass or charge magnitude enters this control.

## Fermionic exchange and the domain

In a hard-core configuration, suppose the destination endpoint is empty
and the origin contains one particle. Put the occupation modes in site order
with internal labels inside each site. Applying `c†_(x,a)c_(y,b)` gives the
sign `(-1)^m`, where `m` is the number of occupied spectator modes between
the endpoints. The endpoint emptiness conditions make this sign independent
of `a,b`. Every entry of `B` therefore receives the same sign on a fixed
spectator configuration. This preserves the partial-isometry identities
and the endpoint covariance. An occupied spectator crossed once contributes
minus one, consistently with the actual `-SWAP` law in the
[colored repair construction](../defect_field_bridge/REPRESENTATION.md).
The tests apply creation and annihilation directly to full occupation masks
and compare the unchanged-internal-label case with the actual full-domain
repair rule.

This statement does not supply a full many-particle collision law. In
particular, the adjoint of a finite raising coordinate can annihilate an
existing link excitation but cannot create its reverse excitation from the
link vacuum. Hence these complete transport instruments are not asserted
to reproduce the original C/E instrument after links are added. They are
an explicit locally gauge-covariant operation on the declared extended
carrier. Its integration with protected defects needs a common preparation,
source-produced operation schedule and error-controlled propagation. No
chirality, anomaly, continuum limit, mass or laboratory species follows
from this edge calculation.

Run from the repository root:

```sh
python -m model.gauge_defect_transport.transport --check
python -m pytest -q -W error model/gauge_defect_transport/test_transport.py
```
