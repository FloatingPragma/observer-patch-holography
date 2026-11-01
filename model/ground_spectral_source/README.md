# Native preparation and charged response of the original cube ground

The existing native qubit service now prepares the proved ground `G` of
the original three-dimensional cube and reads a charged response above it.
The experiment executes the complete 17-bit field/pointer state, with
all source export failures, all 130 repair outcome slots and all six final
pointer outcomes retained. Its amplitudes come from the original action's
exact ground vector, not from an unknown-state or fitted-spectrum oracle.

This is a bounded source execution at `c=1/8,U=1`. The native controls,
blank initialization, coherent pointer and Hamiltonian model time are
declared observer services. A physical parity reference, laboratory unit
calibration, autonomous preparation by repair, thermodynamic preparation
and relativistic particle identification are not derived here.

## Preparing the source-derived ground

The [full-Fock ground certificate](../interacting_vacuum/README.md) supplies
the highest-weight hidden-spin ground with four particles,
`E_G=-sqrt(3)/4`, and 62 nonzero occupation amplitudes. Each amplitude is
`(a+b sqrt(3))/36`, with integer `a,b` and `ab=0`. Thus its squared weight
is an exact integer divided by 1296; every hidden down mode is initially
blank. Compact the eight hidden-up modes to an eight-bit register.

`preparation.signed_weights` obtains all 256 weights and signs from those
exact determinants. The existing real-state compiler proposes a tape of
5,379 literal native pulses. Host square roots and inverse trigonometric
functions propose angles only. The independent signed-tree verifier
consumes every literal event against its Pauli template, inverts each
Walsh layer, and encloses sine and cosine with the existing rational
interval library. For every nonempty conditional node it checks the
signed ratios `sqrt(w_left/w_total)` and `sqrt(w_right/w_total)`.
The zero-weight nodes are checked as identity branches.

The rotations for a fixed tree level have disjoint control blocks, so
their operator error is bounded by the largest node error. Telescope
over the eight levels, the fixed-angle macros and the complete source
load/export maps. This bounds the initialized-state error, including
every first-failure output. A missing leaf sign fails this certificate.

Next embed all 256 rows into the even-numbered modes of the full
65,536-row Fock register, retaining eight initialized blank spectators.
Apply the inverse of the already proved onsite frame. The 184 native
frame pulses implement the exterior-power inverse of

```
R_x=(i sigma_x)^x0 (i sigma_y)^x1 (i sigma_z)^x2.
```

The exact frame check uses
`exp(i s pi P/4)=(I+i s P)/sqrt(2)`. Each local word has an even number
of quarter turns, so its complete four-column matrix is checked using
Gaussian integers and powers-of-two denominators. Vacuum and doublon
phases are included. The literal pulse errors are separately bounded by
rational/transcendental enclosures. The complete preparation state error
is below `10^-10`; its retained exact bound is about `8.67e-13`.
The numerical full-state comparison is a diagnostic, not this proof.

## Charged interferometer and energy differences

At the origin the hidden and original frames coincide, and its down mode
is empty in `G`. Therefore

```
gamma_1=Z_0 X_1=c_1+c_1^dagger,
gamma_1 G=c_1^dagger G,          ||gamma_1 G||=1.
```

The Jordan-Wigner `Z_0` is essential: omitting it produces an orthogonal
probe, since the exact origin-up occupation is `1/2`.
The creation is unitary and supported entirely at the origin. Disjoint
even observables commute with it; disjoint odd expectations vanish in both
definite-parity states. Thus the entire initial exterior CAR restriction
is unchanged, without postselection. The local down occupation changes
from zero to one and total number from four to five. Total number remains
the protecting flag during the stated Hamiltonian and repair waits. This
is an initially localized defect, not a proof that its later wave packet
remains spatially localized. Its mean energy increase is `1/2` and energy
variance is `67/256`; it is not an energy eigenstate or a rest mass.
A 27-pulse pointer word prepares
`(|0>G+|1>gamma_1 G)/sqrt(2)`. Its controlled Majorana is the product of
the four commuting Pauli rotations in

```
exp[i pi (I-Z_pointer)(I-gamma_1)/4].
```

This control explicitly uses the native qubit service's coherent phase
relation between the `N=4` and `N=5` arms. It does not derive such a
reference from a physical fermion-parity superselection rule.

Both arms evolve under the **same original full Hamiltonian**, with no
ground-energy scalar inserted into one arm. The phase difference is
therefore the original energy difference `E_N5-E_G`. A second controlled
Majorana and pointer rotations give, for analyzer phase `phi`,

```
C(t)=<gamma_1 G, exp[-it(H-E_G)] gamma_1 G>,
P(pointer=0)=(1+Re[exp(i phi) C(t)])/2.
```

The two settings use `phi=0,pi/2` and total model time `t=1`. The source
executes two half-time segments, each with 32 symmetric splitting steps
and 337,920 native pulses. The analyzer has 29 pulses per setting. A
complete setting therefore uses 681,459 pulses, counting its independent
ground preparation. The two settings are separate experiments; prefix
caching in the replay only reuses identical deterministic computations.

## Complete repairs and the precise protection interval

After the first half-time segment, execute the five C repair seams; after
the second, execute the five E seams. The complete prior repair theorem
acts as identity on `N4 directsum N5`, including cross-number coherences:
every native branch is zero and the outside branch is identity on that
ideal subspace. This supplies passive protection during the two waits.

Preparation and control pulses need not remain in this subspace at every
intermediate event. In particular, the final inverse Majorana can create
an `N=6` component after the charged particle has moved. **No repair is
applied after that analyzer.** The claim is not invariance of the complete
controlled experiment under arbitrarily inserted repairs, nor active
correction of an arbitrary disturbance.

The floating execution retains all 130 repair records, including their
small off-subspace outputs from control error. It continues on the outside
branch and terminates with each native record on the other branches.
Every export failure also terminates in its own retained record. These
are unconditioned output branches; no surviving branch is renormalized.

## Source and spectral error bounds

`bounds.py` consumes every distinct field and pointer tape against the
source Pauli templates. It counts the complete operator error of every
occurrence, and adds the existing full-Fock Strang bound for both time
segments. The bound applies to the full field with arbitrary reference,
including all intermediate number sectors and every output branch.
Contractivity through the complete repair/readout channels permits the
usual telescoping comparison to the ideal protected experiment.

The resulting complete output-vector error, and hence the error in any
unconditioned pointer probability, is below `3/5000=0.0006` for either
setting. The exact rational value is retained, separately from host
floating-point diagnostics and the host normalization admission.

The [resolved local spectrum](../retarded_matter_response/README.md)
provides all 78 positive gaps and normalized weights of this exact probe.
Rational sine/cosine enclosures evaluate
`C(t)=sum_j w_j exp(-i gap_j t)` using every retained gap/weight interval.
Adding the complete source error gives certified pointer-probability
intervals. The literal complete-state execution is compared with those
intervals and with an independent exponential of the original Hamiltonian.
The `pi/2` outcome is separated from `1/2`, whereas the `U=0` comparator
has the purely real response `cos(sqrt(3)t/16)`. A zero-time native control
checks the two analyzers as well.

These two pointer settings measure the charged characteristic function;
they do **not** independently resolve 78 energies and weights. The full
78-line resolution remains an exact source calculation. Sampling precision,
long-time spectral reconstruction, physical apparatus identification and
a relativistic mass shell require additional work.

## Replay and retained data

Run the package's strict tests and `python -W error -m
model.ground_spectral_source --check`. The JSON receipt and lossless NPZ
archive retain every first-failure weight table, the prepared full ground
state, both complete final states, all repair records and every pointer
probability. Source first-failure and repair state vectors remain
individually addressable from their actual prefixes. The source proof
and scalar bounds do not rely on a floating eigensolver or tolerance.
