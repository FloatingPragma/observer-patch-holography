# The shared native Higgs excitation has a measured current response

The two preparations in [the filtered-state construction](../native_higgs_filtered_state/PROOF.md)
are `A=c_ec^dagger` and `B/2=-epsilon phi^dagger L/2`. They act on the
same complete native Hamiltonian `h=H0+2^-34 T`, with its actual ground,
all occupations and PW119 links. This result calculates the original
quantum dressed weak density W on that shared excitation. It supplies
an actual current measurement rather than transferring a gauge label.
The scalar remains the five-state quantum operator.

## Which charge and which current

The [native charge/reflection symmetry](../coherent_defect_wave/chiral_joint/interacting_charge_conjugation/PROOF.md)
is antiunitary. At the origin it sends `c -> i c^dagger` and
`phi_i -> phi_i^dagger`; every gauge current changes sign. On the full
source it also reflects positions, preserves h and reverses the time
parameter. In particular the dressed scalar density is transformed:

```
W = sum_a (sum_ij sigma_a[j,i] phi_i phi_j^dagger) j_a,
Wbar = Theta W Theta^-1
     = -sum_a (sum_ij conjugate(sigma_a[j,i]) phi_i^dagger phi_j) j_a.
```

The compiler checks `Wbar != W`, its Hermiticity, all twelve Gauss
commutators and all three spin commutators on the full local CAR domain.
It also checks `Theta Wbar Theta^-1=W`. Thus the earlier removal-channel
mean `0.0915665...` of W cannot be assigned to the addition channel by
renaming its charge. W on an added excitation equals Wbar on the
conjugated removed excitation, with the time reversed. The error bound
below is uniform for either time sign. This is not a claim about measured
CP violation: the transformed measurement is a different operator.

The complete old/dual/port reference is retained. Its common vacuum
and equal-cost conjugate ports allow this pairing. Particle-hole applied
only to the old `4N` reference would change its clock and is not used.
The selected even compensator ports have clock cost zero in this
experiment. Subtracting four from the old source-plus-reference bin
therefore identifies the native h bin; it does not change h or subtract
an unknown ground energy from the apparatus.

## Exact atomic contraction

The calculation uses all three complete charged target blocks of
dimensions 7, 17 and 7. The full CAR polynomial images close on each
block and are Hermitian in its actual Gram metric. The complete
nine-root annihilator q(H)=0 supplies the Riesz projector

```
P_e = sum_i H^i sum_(j>i) q_j e^(j-i-1) / q'(e).
```

Here e is the least root of
`e^4-2e^3-212e^2+213e+11114=0`. Exact algebraic arithmetic establishes
rank one in each magnetic block. All three contractions are equal, so
the result also holds for an arbitrary retained spin environment of
the atomic ground. No local numerical eigenvector is inserted.

For Wbar on the removal target, or W on the addition target, the mean is

```
mu(e) = (991825 e^3 +147657925 e^2-285583600 e)/100725009473
        -3967300/18125789.
```

It is strictly negative. The same calculation gives
`<W^2>=-mu/2`, probability `-2mu` at W=-1/2 and probability `1+2mu`
at W=0. The other three probabilities vanish on this atomic selected
state. The full interacting instrument still retains all five outcomes;
the atomic zeros are not assumed to survive hopping.

## Transfer to the unchanged interacting source

Use the actual spectral bin I of radius `128v`, where
`v=360*2^-34`, from [native charge resolution](../native_charge_resolution/PROOF.md).
Let psi0 be the selected atomic packet, of squared norm w0. Let eta
bound the actual neutral ground's distance from its atomic trial, d0
be the selected atomic energy's distance to the bin complement, and
rho_other the sum of all other visible atomic leakage bounds. For
`a_t=exp(-ith) P_I A g`, with the harmless atomic phase absorbed into
psi0, the complete-source estimate is

```
||a_t-psi0|| <= sqrt(w0_upper) v/d0 + rho_other + eta
               + |t| v sqrt(w0_upper) = r.
```

The first term follows from the selected packet's full residual second
moment. Every other visible root is included in rho_other. The final
term is the full-space Duhamel bound, not an evolution in a selected
charged subspace. The two actual source-clock intervals before the
current record give `|t|=2`; r is less than 0.005233. Source evolution
on all partners and failure branches continues throughout.

Let `[w_lo,w_hi]` be the actual unnormalized bin-weight enclosure. Since
`||W||=2`, the unnormalized insertion differs from `w0 mu` by at most

```
2 r (sqrt(w_hi)+sqrt(w0_upper)).
```

Interval division by the strictly positive bin weight gives a bound on
the conditional mean. For each complete W projector the same formula
uses operator norm one. This controls the entire five-outcome detector,
including the outcomes that were zero in the atomic comparison.

The executed preparation specification is the depth-three native
agreement-filter construction: eight native trials, seven filters,
all partners and all filter records. Its retained purification error
epsilon is less than `6.451e-13`. Conditional bin states differ from
their exact-ground comparison by at most `epsilon/w_lo` in trace
distance. Consequently the W-mean error is at most `3epsilon/w_lo`,
using W's spectral diameter three. Every projector probability has
error at most `epsilon/w_lo`; unnormalized effects have error epsilon.
This argument never places an exact ground projector in the apparatus.

The resulting original-W mean is enclosed by

```
A preparation:   [-0.07004020, -0.00539719],
B/2 preparation: [-0.07005609, -0.00538131].
```

The displayed decimals are outward summaries of the exact rational
receipt. Both upper endpoints are strictly negative. For the second
line, the normalized filtered A and B/2 states differ by less than
`3.971483e-6` in trace distance. Applying the same current instrument
contracts that distance. The B/2 mean differs from the A mean by at
most four times that bound. This transfers a calculated current
response to the actual quantum-Higgs preparation on the same h.
It does not equate their different preparation success rates.

## Complete apparatus and its scope

Choose A or B/2 with probability one-half. Retain the native preparation
flags, pump failure, rejected energy bin, and all five current labels
on the admitted branch. The full setting/pump/energy/current alphabet
has 80 labels, including zero-probability unused labels. If Q is the
actual preparation success and C is the selected pump, the current
records are

```
Q/2 tr(P_w U(2) P_I C rho C^dagger P_I U(2)^dagger P_w).
```

Together with both pump failures, both bin rejections and every
preparation failure, these sum to one. The receipt gives positive
original-attempt probabilities for the W=-1/2 outcome of both pumps.
They include their unequal success weights. The finite cubic reference
and charged apparatus are those of the filtered-state construction;
W is a scalar under its finite rotations and commutes with the native
Gauss constraints. Its measurement can change the energy admitted by
the earlier detector. No simultaneous exact energy/current value is
claimed.

The finite readout is the complete controlled pointer shift with a
Hermitian logarithm L. The declared source-clock service
`V(s)=exp(-iHext s) L exp(iHext s)` gives
`exp(-iHext) U_read` over one source unit. The record refers to its
interval start; h, partners and references are unchanged. This remains
a programmable full-matrix service, with the resource qualifications
of its parent construction, not an autonomous local repair compiler.
The five-state zero-energy current pointer adds `ln5` nats to the
filtered parent construction's retained-register capacity. Coefficient
and synthesis workspace remain separate resources.
Alternatively every pulse can use the parent's constant full-matrix
drive: for the final pointer unitary take
`L=i Log(exp(-iHext) U_read)` and `V=L-Hext`. The complete four-pulse
schedule ends at source time four. With
`B_H=8(600+360*2^-34)`, each constant control has norm at most
`B_H+pi`; net bare-source work has magnitude at most `2B_H`, and gross
switching work for the four constant pulses is at most `8(B_H+pi)`.
These bounds include the retained reference sectors and use the same
Hext on every failure branch.

The proved object is a finite native excitation with a shared
elementary/Higgs energy and a nonzero measured quantum current.
Relativistic chiral propagation, weak scattering, a W boson and a
Standard Model lepton mass are separate requirements.
