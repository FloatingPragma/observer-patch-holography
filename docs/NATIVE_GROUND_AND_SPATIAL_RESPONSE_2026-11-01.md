# Native ground readout and spatially local interacting response

The [continuum-field and spectral-readout extension](CONTINUUM_FIELD_AND_SPECTRAL_READOUT_2026-11-01.md)
continues this construction with a common spacing limit, complete finite
ground-family spin selection and positive finite-resolution energy POVMs.

The original cube ground can now be prepared through the declared native
observer service, and its charged response can be read by a complete
interferometer. The full eight-site response has an exact matrix spectral
resolution. On the infinite lattice, the original ground has provably
non-Gaussian field correlations, nonzero interaction memory and a spatially
local inverse response at
complex frequencies separated from the real axis.

These results extend the [ground and response construction](GROUNDED_FIELD_AND_GEOMETRY_RESPONSE_2026-11-01.md).
They close specific preparation, matrix-spectrum and spatial-response
inputs. The physical particle dictionary, interacting chiral 3+1 continuum
and joint Einstein premises remain open.

## Preparing and reading the finite ground

The [native ground experiment](../model/ground_spectral_source/README.md)
derives its 62 nonzero preparation amplitudes from the certified ground of
the original `c=1/8,U=1` action. Their squared weights are exact integers
divided by 1296. It compiles and executes 5379 preparation pulses and 184
inverse-frame pulses, retaining all export failures. Independent rational
and trigonometric bounds certify the preparation error; no unknown-state
oracle or fitted spectral amplitude is supplied.

An origin-supported Majorana operation creates the charged probe. Because
the origin's down mode is empty in the ground, this unitary produces exactly
the same vector as the down creation operator. All initial exterior
observables are unchanged, without postselection. The mean added energy is
`1/2` in model units, with variance `67/256`; the resulting state is a
superposition of energies.

The complete experiment evolves both interferometer arms under the same
original Hamiltonian for model time one. It executes the 17-bit
field-and-pointer state, all 130 scheduled repair outcome slots and all six
final pointer outcomes. Each complete setting costs 681459 native pulses,
including its own preparation. The executed probabilities are:

| Analyzer phase | Probability of pointer outcome zero |
| --- | ---: |
| 0 | 0.8849880740005236 |
| pi/2 | 0.7053341583573907 |

The certified complete source error is below `0.0006` for either setting;
the independent original-Hamiltonian comparison has output error about
`2.014e-6`. These are complete instrument probabilities, not sampled
laboratory counts or a 78-line finite-shot reconstruction.

The repair guarantee has a precise interval. The ideal state occupies
`N=4` and `N=5` during the Hamiltonian waits, where the specified repair
service acts as identity. Preparation and analyzer pulses are explicit
observer operations; no repair is inserted after the final analyzer, which
can populate other number sectors. All termination and failure branches
remain in the output law. The experiment uses the native coherent qubit
reference and initialized blank input. It does not establish autonomous
ground selection or preparation of the infinite-lattice ground.

## The complete finite matrix response

The [matrix spectral construction](../model/matrix_spectral_response/README.md)
proves that all eight hidden-frame down-addition probes have the same 78-frequency
support. Every full matrix residue has rank four, giving an observed cyclic
space of dimension 312 inside the complete 560-state charged sector.
There are 40 lines in one signed-adjacency band and 38 in the other.
The remaining 248 dimensions are orthogonal to these probes and their time
translates; their full spectrum is not claimed here.

For the signed cube adjacency `A`, with `A^2=3I`, the exact response is
`G(z)=F(z)I+J(z)A`. Its inverse and memory have the same matrix form.
This determines their static coefficients exactly over `Q(sqrt(3))`,
replacing a numerical solve with exact rational functions. On the observed
reducing space, the certified frequencies lie between `0.0299` and `1.85`.
The resulting memory expansion is controlled for `|z|<0.0299`, improving
the earlier conservative full-sector gap disk.

This exact locality uses the special two-eigenvalue cube adjacency. It is
not a spatial continuum theorem or a derivation of relativistic dispersion.

## A definitive test of the finite-gap particle dictionary

The complete visible charged spectrum has largest-to-smallest energy ratio
at most approximately `61.3027359769`. The retained PDG comparison gives
the muon/electron mass ratio above `206.7682784259`, including the declared
reported-uncertainty intervals. Consequently no assignment of these lines,
with one positive energy-to-mass conversion factor, can reproduce that
pair. At least one relative error is at least approximately
`0.54263808705`. The corresponding electron–tau lower error bound is about
`0.96535092456`.

This rejects the explicit proposal to treat the finite visible gaps as
those particle rest energies. It does not assume that a finite gap is
already a relativistic rest mass, and it does not exclude other sectors,
an appropriate continuum limit or OPH as a whole. No physical species has
been assigned by fitting these values.

## Interactions and spatial response on the infinite ground

The [spatial response theorem](../model/retarded_spatial_locality/README.md)
uses the same original infinite-lattice ground constructed by averaging
full finite-box ground projectors. Its inherited hidden SU(2) symmetry and
density bound give the onsite occupation matrix `rho=rI`, with
`1/8<=r<=1/2`.

Exact CAR identities determine the coupling to the unobserved response
space:

\[
W^\dagger W=U^2r(1-r)I,\qquad
\frac{7U^2}{64}I\le W^\dagger W\le\frac{U^2}{4}I.
\]

Thus the infinite ground has nonzero retarded interaction memory for
`U>0`. Its normal response cannot be replaced by a static,
number-conserving quadratic Hamiltonian on the same modes. This is stronger
than merely retaining a nonzero quartic operator in a faithful
representation. Two-point memory alone does not rule out descriptions
using additional modes or anomalous pairing blocks.

A stronger four-point result now does distinguish the full field. Put
`a(t)=tau_t(c_xsigma)` and center its anticommutator with the original
creation operator:

\[
D(t)=\{a(t),c_{x\sigma}^\dagger\}
 -\omega(\{a(t),c_{x\sigma}^\dagger\})I.
\]

The original equation of motion gives
`D'(0)=-iU(n_x,opposite-r)`. With
`M=9c^2+12cU+U^2`, the exact local derivative bounds imply

\[
\omega(D(t)^\dagger D(t))\ge U^2t^2/64>0,
\qquad 0<|t|\le U/(16M).
\]

This variance equals a sum of four connected fermionic four-point
correlations, with all normal and anomalous pairings retained. At least
one is therefore nonzero: the original infinite-ground field is
non-Gaussian. Linear fields under quadratic dynamics, including pairing
and auxiliary modes, have a scalar anticommutator and zero such variance.
The witness excludes that representation of these observables. It
establishes interaction at fixed spacing; it supplies no scattering or
continuum limit by itself.

The same positive-energy generator now constructs every ordered local
imaginary-time correlation. All orders use the same ground and the same
finite-volume subsequence. Positive-time words define Hilbert-space
vectors; reflecting and adjoining a second word gives their Gram matrix,
which proves time-reflection positivity. The construction preserves CAR
order, including odd operators. A Poisson integral of real-time evolution
proves convergence from the finite grounds, with an explicit integration
tail bound. This supplies the full ordered Euclidean hierarchy at fixed
spacing, without assuming Euclidean rotation symmetry or a spatial
continuum limit.

At fixed `Im z>=eta>0` and `|z|<=Z`, the original full-CAR locality bound
gives exponential spatial decay of the Green operator. A uniform inverse
bound and a bounded-weight operator argument then give exponential decay
of its inverse and memory. Their spatial Fourier symbols are analytic,
with explicit spatial and frequency Taylor remainder bounds. This
construction requires no assumed ground-state clustering.

## The remaining geometry and field join

The existing `wave_density_metric` and `source_tetrad` theorems remain the
geometry interface. The new result supplies a controlled spatial inverse
response on the actual interacting infinite ground. Its frequency domain
stays above the real axis, and the estimates deteriorate as `eta` tends to
zero. Controlled real-frequency spectral measures and the appropriate
thresholds, resonances or poles, a common hyperbolic symbol across the
matter sectors and action normalization are still required before that
response can define the physical metric.

The same-family modular, stress, area-response, coarse-graining and
null-balance premises are then needed for the Einstein conclusion. The
one-dimensional massive continuation remains a separate construction;
its chosen coupling flow is not a derived mass spectrum of this
three-dimensional ground.

The scaling limitation is explicit. Along `c=v/a` at fixed positive
`U,eta,Z`, the present bound on the non-Gaussian time interval is of order
`a^2`. The guaranteed physical momentum strip for the inverse response is
of order `a^8`. Both shrink in the continuum limit. These estimates prove
the fixed-spacing field results, but do not prove that the interaction or
derivative expansion survives at a fixed physical scale. This identifies
the needed uniform estimate; it is not a no-go theorem for the action.

## Validation and claim status

The three component packages pass 98 strict tests: 31 for native ground
preparation and readout, 21 for the matrix spectrum, and 46 for the
infinite-ground field and spatial response. Each package also passes a
fresh-process write and independent fresh-process check. Independent
source and proof reviews cover all three constructions.

The integrated tree also passes 495 integration tests and 22 subtests,
plus 136 supporting regression tests. All 1,414 model artifact hashes,
generated claim/premise/completion views, postdiction-ledger parent hashes
and the eight-paper release manifest check pass.

The integrated map records 67 realized scoped attachments among 72
reviewed domains. The custody audit preserves all 84 canonical premise
dispositions, 346 claim classifications and existing completion-gate
statuses. It adds no physical species assignment or accepted physical
completion gate. The eight updated papers build successfully, pass the
publication warning checks and pass visual review of the changed pages.

These are mathematical proofs with exact finite witnesses and executable
error checks. This extension adds no new Lean theorem and does not claim
that testing the finite witnesses by itself proves the infinite limits.
