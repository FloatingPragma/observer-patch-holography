# Klein–Gordon dispersion in the routed Dirac field construction

The construction is the routed four-spinor, noncompact quadratic field source
with its native phase clock and the three register-derived flavor masses.
An explicit coupling parameter g multiplies every route current and physical
longitudinal charge; its Coulomb energy scales as g². The size-two g=1 word is
the frozen `native_field_connection` word, including all 190 coexact modes,
three harmonic canonical pairs, source scalar phases and clock. The mass slot
has the exact substitution supplied by `species_reuse_connection`. The g=0
intervention gives an exact free factor with controlled Klein–Gordon dispersion
and native complete measurements. It is a neutral/free intervention on a
specified source family. It is not a proof of free charged-particle dispersion
at g=1 or of interacting mass poles.

The numerical and exact receipts are constructive producer evidence. They do
not issue an official row validation. The accepted Maxwell field construction
and this free matter factor have the same unit propagation speed in the
same h-time convention. The scalar walk with speed 1/sqrt(48), and the separate
90-coordinate Higgs-Wilson Hamiltonian, are not identified with this source.

## Actual routes and the free intervention

Every root-port route is a path through actual LOCAL12 seams. Positive routes
join the old root to its outgoing port, traverse its physical intercell seam,
and join the incoming port to the destination root. Negative routes reverse
the corresponding positive path. The exact incidence identity is

    B J = delta_start − delta_end.

The code checks every route for L=2,3,4, its reverse and every seam. A deleted
intracell current fails the identity. No continuum translation stencil is
inserted in place of these paths. The higher-size family uses the same local
port rule and root code, with its periodic coordinate size as a parameter.

Let P_L and P_T be the longitudinal and transverse projections of the actual
edge space. The charge-dependent particular electric field is g P_L J and
the oscillator displacement is g P_T J. Their sum is gJ, so g n+B E retains
its declared background. Its longitudinal energy is g² times the original
Coulomb energy. The source has an explicit neutral g=0 intervention. Merely
deleting the transverse kicks while retaining the Coulomb phase would leave
a position-dependent potential; the rejecting control measures the resulting
plane-wave leakage.

At g=1 `parameter_macro` delegates the complete frozen iterator. For other g
its current coefficients and Coulomb phase have exactly the displayed scales;
the spin basis, mass, destination permutation, free field word and clock are
unchanged. At g=0 all current generators are zero and the Coulomb phase is zero.
The complete source is therefore

    U_0(h) = U_D(h) tensor exp(−i H_field) tensor U_clock(h).

This identity keeps the entire field space and its phases. For every normalized
or subnormalized joint state, including an arbitrary external reference,
partial trace over the autonomous field commutes with this product evolution
and each local matter/clock instrument. The mathematical identity supplies the
reduction; an independent correlated-state test checks it. The normalized
coexact vacuum and freely spreading harmonic Gaussians are valid references.
No translation-invariant nonnormalizable field vector is used to turn off g=1.

The g=1 native field implementation is a size-two construction. Higher-L
coupled canonical sources have the same finite geometric definition; an
interacting spatial native or field-limit theorem is not claimed here. The
free factor requires only the exact autonomous-field reduction. Its field
sector is the free Maxwell source on the same actual cell complexes.

## Clock and field speed

A completed lattice macro advances physical source time by h and native clock
phase by h/8. The field free duration per lattice macro is one, as in the
frozen field word. Thus its physical generator is H_field/h. For k=hp its
low frequencies approach sqrt(ab rho)|p|, where a=1/2, b=3 and rho=2/3.
Their speed is one. There is no conflict between the lattice field-time
counter n and the physical time n h.

The free matter symbol derived from the signed routes is

    U_D(k,h) = exp(−i k_z alpha_z) exp(−i k_y alpha_y)
               exp(−i k_x alpha_x) exp(−i h beta M),

where alpha=(XX,XY,XZ), beta=ZI and M is the diagonal flavor mass operator.
The same source time h gives its continuum speed one. This matching is an
architectural consequence of the stated field coefficients and source clock,
not an independent empirical determination of c.

## Fourier transfer, spin and flavor

The native signed increment uses the diagonalized spin rail: rail zero adds
one and rail one subtracts one, modulo L=2^b. Exact reversible checks cover
all coordinate columns, clean scratch return and the defined dirty-scratch
bijection through L=16. Each dyadic increment is an explicit finite Toffoli
ladder. Its scratch is restored before the spin basis is restored.

For the normalized actual coordinate vector
phi_nu(x)=L^(−3/2) exp(2 pi i nu.x/L), translating by sign s along an axis
gives exp(−2 pi i s nu_axis/L) phi_nu. This exact intertwiner, the tested
spin diagonalizer and the unchanged mass slot give U_D. Both signed rails
and all four spin columns are retained.

At L=4 a complete 350-pulse spatial source word is executed from the full
six-level native responses on all sixteen spin/flavor columns. Its cell
wave is prepared by eighteen native pulses, including its scalar phase.
The full 1024-dimensional output agrees with the Fourier intertwiner to
8.60e−15 in the numerical diagnostic. The projection is not used to discard
an output component; the comparison uses the complete spatial output vector
and accumulates every first export complement.

The three positive masses are the exact rational midpoints of the frozen
finite register-effect intervals. Their register, analyzer and unit-kick
ancestry remain explicit. M has one additional zero entry for the retained
flavor padding label. That block still has its defined kinetic evolution.
It is not removed or declared a fourth physical family. The mass slot is the
actual native beta tensor M word, preserving coherent flavor superpositions
and arbitrary untouched references. No empirical GeV calibration is used in
this connection. A comparison to the interval-valued register coefficients
adds T times their maximum midpoint radius.

## Controlled Klein–Gordon limit

Exact Clifford identities give, for every flavor block,

    H(p)=alpha.p+beta m,
    H(p)^2=(|p|²+m²) I_4.

Both positive and negative energy branches and both spin multiplicities are
part of the source. There is no projection onto a preferred sign. The full
flavor identity is H²=|p|² I+I_spin tensor M². Consequently every component
of its continuum matter wave obeys the Klein–Gordon equation.

For a=(m,p_x,p_y,p_z), unitary product comparison yields

    ||U_D(hp,h) − exp(−i h H)|| <= h² sum_(i<j) |a_i a_j| = delta,
    ||U_D(hp,h)^n − exp(−i n h H)|| <= n delta.

The bound follows by integrating each commutator; the Clifford commutator
has norm two, which cancels the product formula's factor one-half. It is
uniform on every bounded momentum and mass set and tensors with a reference.
Any complete terminal instrument has TV error at most min(1,n delta).

When delta<1/2 and the phases are inside the stated unaliased band, unitary
spectral perturbation bounds the phase displacement by
2 arcsin(delta/2)<=2 delta. Hence the frequency error is at most
epsilon=2 delta/h and the squared dispersion error is at most
(2 E+epsilon)epsilon, using an outward upper bound for E. All constants in
the receipts use rational upper bounds, including pi<22/7.

The executable checks 120 complete four-spinor spectra, six momentum
directions, all four flavor code blocks and five spatial refinements. The
actual finite Floquet splitting of spin branches is retained. The eight
Brillouin corners also remain: their zero/π structure is not removed by a
species claim or an unannounced Wilson regulator.

For L=2^18, physical T=1/4 and nu=(1,2,3), the complete observed native
record bound is about 0.000422021, including source product error, native
controls, reset timeout and coefficient enclosure. The specified controlled
spatial program has a finite pulse upper bound of 178,026,512,664 and uses
128-bit rational controls. That word is **unexecuted**. Its finite generator,
native primitive certificates and exact count define the construction; the
literal spatial executions are the finite members identified above.

## Native source tomography and complete records

The controlled apparatus controls every I, X, Y, Z, H and J source primitive,
including its scalar phase. It probes U_D, the exact reduced source factor.
It does not substitute the controlled full interacting field evolution or
reuse that evolution's failure tape. The explicit apparatus is a fresh
complete implementation of the proved autonomous-factor observable.

Twenty-eight native spin preparations comprise four basis states and four
relative phases for each of six basis pairs. Full probe X/Y measurements of
their quadratic forms reconstruct every entry of the four-spinor propagator.
All three physical flavor labels and the padding label are measured. The
native clock is measured jointly in both quadratures. The raw twelve control
coefficients are applied to the full six-level response matrices; every first
export complement is accumulated on all 32 controlled source input columns.

Four finite source cases, one/two/four macros and both field-record laws give
5,376 component records and 908,544 explicitly evaluated fine probabilities.
The maximum reconstructed propagator discrepancy is 2.38e−13. Native first
failures, pointer outputs, padding and no-operation outcomes are retained.
The exact interval certificate supplies the control bound; floating matrix
residuals are diagnostics and do not replace it.

The native pointer is separately executed on all eight columns, including
its READY preparation from blank work bits. Its complete C/E effects define
every fine joint output by tensor contraction on unnormalized success states.
The stored component arrays precede actual reset-success weights. The full
instrument composes rank-one resets, multiplies later component laws by the
actual success weights and retains each first timeout. It never appends reset
failure masses to an otherwise unchanged unit-mass array.

Changing the actual native mass slot or clock tick is detected. Additional
controls reject an always-positive rail, retained Coulomb energy at g=0,
deletion of a controlled scalar phase, omission of a physical route leg,
discarded flavor padding and replacement of the Dirac speed by the scalar
walk speed. The scalar acoustic velocity squares have rank two and fail
the Clifford identity, so no relabeling establishes the connection.

## Full field failure cost and finite coupling

The size-two terminal full-field factor has a separate finite coalesced native
holding contract. At four lattice time units it keeps all 190 half-frequency
phases and three harmonic p² generators, its blank resets and first export
outcomes. Its additional reduced-record bound is about 6.25e−14. Its enormous
field word is specified but unexecuted. Harmonic cutoff and splitting alter
only an autonomous spectator unitary for a matter-only marginal; native
export and reset failures are still counted. Full-state or electric-field
measurements require the harmonic and spectral transfer bounds. No size-two
field cost is transferred to a larger-L field register.

A state-controlled finite-size connection also bounds g tending to zero.
On the true coexact vacuum and normalized harmonic Gaussian reference, each
free-prefix route generator has norm at most sqrt(5 v_n), where
v_n=a/(2 omega_min)+1/2+n²/8 bounds the coordinate covariance and five is the
maximum route length. The cell/rail predicates are orthogonal, so their
maximum, rather than a branch-count factor, bounds the controlled kick.
For n macros, product telescoping on free prefixes gives

    ||U_g^n psi − U_0^n psi||
      <= min(2, 3n |g| sqrt(5 v_n) + n g² max E_Coulomb).

This is a normalized Gaussian-state bound with arbitrary matter/flavor
reference. It is not operator-norm smallness for an unbounded oscillator
generator. The receipt supplies finite values for g=10^−5, 10^−3 and 1.
The g=1 bound is generally vacuous; the proof does not convert the free
Klein–Gordon result into an interacting charged-particle claim.
