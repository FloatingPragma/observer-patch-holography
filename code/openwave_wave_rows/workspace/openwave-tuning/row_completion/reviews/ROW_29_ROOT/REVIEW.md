# Full radiative Maxwell criterion: independent root review

**full_criterion_met = true** for the radiative observation domain of
`native-local12-dirac-noncompact-quadratic-v4-proposed`, bound by row
manifest `993744d14c2da8d66a93dc76422f81799074febad9914ae33bb60bf455cf00b6`
and size-three manifest
`90178885fcaaaba405b950a6f132fd47d519d21a130a4f982ec09423db129ec2`.
The literal criterion is Maxwell recovery with transverse waves at the
common speed. This review accepts the controlled vacuum radiative limit
and its terminal native first-moment measurements. It does not accept an
interacting continuum theory or an electromagnetic strength prediction.

## From the actual faces to Maxwell

The microscopic inputs are the translated oriented carrier faces,
noncompact quadratic electric and face energy, coefficients `a=1/2,b=3`,
spectral encoding and completed-source clock. These are declared
architecture choices. The resulting continuum curl matrix and its
coefficient are calculated from those faces. A continuum Maxwell stencil
or a list of target photon frequencies is not used to construct the face
matrix. Matching the field speed to the matter speed influenced the
coefficient choice and is retained as input ancestry.

The independent executable reconstructs the face symbol and its first
derivatives from integer translated terms. Multiplication by the incidence
symbol gives zero coefficient by coefficient as a Laurent polynomial.
This proves the complex identity for every momentum, rather than only
at sampled momenta. Nullspace bases and their exact Gram inverses give
the left, right and harmonic projectors, independently of the producer's
pseudoinverse construction. They verify

    W^T W = (2/3) I,  V^T V = I,
    P D_i Q = W cross(e_i) V^T.

This is an oriented factorization. Squared eigenvalues alone would not
fix the electric/magnetic sign. The adjoint of `i W cross(k) V^T` is
`i V cross(k) W^T`, because the real cross matrix is antisymmetric.
With the stated A, E and B dictionary this gives Faraday's equation and
Ampère's equation with coefficient `ab(2/3)=1`. The cross matrix kills
the longitudinal vector and has exactly two nonzero singular values.
Rank of the incidence symbol excludes additional positive low branches.

The executable also checks paired left/right singular frames on the
actual finite symbol for all 26 nonzero directions with components in
`{-1,0,1}`. It compares the physical vector-potential, electric and magnetic
amplitudes, including their orientation and Fourier normalization, with
the stated bounds. These numerical checks supplement the uniform proof;
they are not substituted for that proof or described as enormous lattice
executions.

## Controlled limit and common units

The entrywise translated-term bounds independently give the constants 11
and 17 for the first difference and second-order remainder. The gap
estimates `||A0^+||<4` and `||A0||<10` follow from the exact rational norm
certificates. Left and right skew rotations each have norm at most
`4 eta`. Differentiating their unitary conjugation twice bounds the
quadratic remainder by `32 eta^2 (10+eta)`; their action on the first
perturbation contributes `8 eta^2`. This establishes the stated D bound.

The Hermitian dilation separates the two positive low eigenvalues from
the zero and high clusters. Projection of their eigenvector equations
onto the comparison complement, followed by polar alignment, gives the
stated paired-frame estimate. Using a squared operator with a gap of
order `|k|^2` instead would not supply this estimate. Paired left/right
control is also what avoids an uncontrolled division by `|k|` in B.

The factors `h^(3/2)` from Fourier normalization, the microscopic
oscillator amplitudes and the declared field rescalings give equal
electric and magnetic amplitudes `sqrt(|p|/6)`. The potential amplitude
is `sqrt(1/(6|p|))`. The dictionary and frequency errors tend to zero for
each bounded nonzero momentum set. Finite Fourier sums therefore obey
the vacuum Maxwell equations with two polarizations. Constant harmonic
fields obey the zero-momentum equations separately. The source clock
uses the same h per completed macro as the matter source.

I independently recomputed the exact rational frequency, polarization
and complete native/source bounds at the three stated refinements. Their
complete record total-variation bounds are approximately `5.70437e-5`,
`1.42841e-5` and `3.59419e-6`. Physical first-moment bounds additionally
include the stated field-dictionary error. Fixed native control precision
leaves its reported floor; only the source limit is asserted to converge
at fixed compiled precision.

## Preparation, observation and spectators

The cap contains vacuum and the two one-photon polarizations. It commutes
with the exact free Hamiltonian. Thus compression preserves the free
linear-field Heisenberg equations even though it changes the canonical
commutators. Squaring compressed fields would give a different Hamiltonian
and is correctly rejected. The source vacuum scalar is included.

The selected-mode code and reset/Hadamard preparation are explicit. The
complete observed instrument is equivalent to its restriction from the
free full source by partial trace, for every normalized spectator state
and even with correlations to an external reference. This statement does
not prepare a full interacting vacuum, nor does it make spectator field
fluctuations vanish. The measured observable is the selected radiative
Fourier component. It equals the full mean when spectators have zero
first moments. This is a sufficient field observation for the stated
vacuum-wave criterion.

The native bit analysis is a bounded two-outcome field first-moment POVM
on the three-dimensional physical code. Its dilation can reach padding,
which is retained. It is not a sharp measurement of the full continuous
electric-field distribution. The position phase selects a physical
quadrature; the source's oriented-symbol proof supplies its propagation
law. The reduced word uses the derived limiting frequency and pays the
finite-source error. Its observed speed is verification of this controlled
implementation, not an independent prediction of a number supplied to
the reduced word.

The stored eight-outcome arrays are prepared-input word/pointer component
laws. Complete experiments compose them with the pinned native reset
instrument and retain its fine successes and timeouts, together with all
first export failures. Those arrays alone are not described here as a
complete unconditional reset-history distribution. The uniform reset
bound is included in the complete channel estimate; success conditioning
is not used to improve that estimate.

All three pointers are priced: `phase_clock.exact_word_bound` includes
the clock pointer, and two explicit pointer allowances cover both field
code bits. The apparent count of two in the caller is therefore correct.
The three blank and three pointer-ready reset costs are also retained.
Terminal coalescing uses exact free group composition and its own native
word bound; no per-macro history or intermediate observation is claimed.

The complete size-three implementation is additional finite evidence on
all 646 physical coexact modes with its retained padding and clock. The
huge periodic boxes in the refinement theorem are not described as fully
allocated or executed native experiments. These distinctions are explicit
parts of the accepted source/observation domain.

The packet meets the literal radiative Maxwell criterion in this domain.
An empirical coupling, charged scattering, interacting spatial limit,
or an equivalence with the literal classical record source does not follow
from this decision. Such claims require their own source transfer.
