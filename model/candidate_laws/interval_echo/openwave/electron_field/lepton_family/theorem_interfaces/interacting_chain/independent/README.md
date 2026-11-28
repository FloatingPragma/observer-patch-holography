# Eight-site charged spin reception

After a coherent low-energy pass on the archived eight-site Weyl/contact
chain, free evolution carries the full admitted source record halfway
around the ring. A receiver at sites 4 and 3 or 4 and 5 can read a signed
unit of excess material and a neighboring spin on the same trial.
The spin retains entanglement with the original source reference for
either charge sign. This two-site read uses a terminal occupation gate
and the local spin analyzer. It does not use a whole-chain ground-return
report. [Protocol](PROTOCOL.md) · [Replay receipt](records.json)

## Independent spectral count

The replay builds the repulsive contact Hamiltonian with all up-spin
orbitals followed by all down-spin orbitals. Its hopping is $t=1/16$,
its on-site repulsion $U=1$, and the boundary is periodic. An explicit
signed permutation equates its full integer matrices to the
archived interleaved chart. The supplied Weyl-to-hopping onsite spin
frame closes on the eight sites and agrees at sites separated by four
addresses. Both translations therefore refer to the same physical
half-circuit.

The half-filled singlet is unique at calibrated momentum index 4.
Its energy lies near −0.08697882485 in original Hamiltonian units,
and the half-filled gap lies within [0.0080272286, 0.0080272461].
The marked Néel preparation reaches this ground on **14.545735%** of
original starts. The spin-raised sector has exactly zero states below
the ground comparison cut. Spectra were counted with a Hermitian
factorization and orthogonality enclosure in each momentum block.
The calculation verifies the integer translation and spin-raising
intertwiners on the complete Fock sectors. A finite source-frequency
row inherits the exact rank statement from these counts.

The exact dyadic charged cut is 0.830078125 for addition and
−0.169921875 for removal in the original units. Both cuts isolate the
same sets as their respective midpoint cuts. The complete spin
multiplicities between the floor and the cut are shown here. Energies
use the independently specified half-filled reservoir correction.

| Corrected transfer energy | Spin | Addition momenta | Removal momenta | Local source |
| ---: | ---: | --- | --- | --- |
| 0.395474289 | ½ | 2, 6 | 2, 6 | Bright |
| 0.408396129 | 3/2 | 3, 5 | 1, 7 | Symmetry forbidden |
| 0.415816787 | 3/2 | 4 | 0 | Symmetry forbidden |
| 0.416560645 | ½ | 1, 7 | 3, 5 | Bright |

Each listed momentum occurs once in its indicated spin multiplet.
Within the fixed $S_z=1/2$ chart there are seven eigenvectors below
the cut, four bright and three dark. The dark rows are exactly
spin-3/2 by subtracting the counts in successive raised sectors.
A local spinor acting on the singlet reaches spin-½. This argument
forbids the dark overlap independently of numerical truncation.
The next distinct charged level is separated from the upper bright
line by more than 0.00102727.

On the reachable part of the window, translation by four sites has
eigenvalues $+1$ and $-1$ at the lower and upper bright energies.
Thus half of their beat evolves every admitted local spinor as the
same four-site translation up to a common phase. The proof keeps
both momentum characters of each energy and both spin components.
The phase comparison time is enclosed in
[148.98698797, 148.98698945] in native units. An independent sparse
native evolution at the fixed real wait 149 checks the same
comparison without using the sign-flip identity as the propagator.
The coherent ideal window accepts about **0.67987%** of original
Néel starts after the local pump.

## Simultaneous local charge and spin

The terminal charge observer reads a double occupation at site 4 for
addition or an empty site 4 for removal. A separately calibrated
neighbor measures its spin on single occupation. The two-site excess
is therefore exactly $+1$ or $-1$ relative to one fermion per site on each
accepted trial. The chosen spin site is either 3 or 5 in separate
banks. The observer uses the same pair of sites at every wait.
The gate is a spin scalar, and its weight is independent of the
source spin. The source and local observer retain their original
reference labels and compare their Pauli reports.

Let $p$ denote the charge-and-single gate probability within the
charged window, and let $m$ denote the mean of the spin report
$2S_z$ times that gate when the source reports spin-up. The accepted
spin map has Bloch contraction $\eta=m/p$. The exact spin
intertwiners imply covariance. For this one-spin input and output
covariance gives a rotation-invariant qubit channel with the same
contraction on all three axes. The joint reference and received-spin
Bell overlap is $(1+3\eta)/4$. It exceeds $1/2$ exactly when the
accepted two-qubit comparison is nonseparable. This inference uses
only the two neighboring material reads and the retained reference
after the source window has passed.

| Two-site test | Joint pass per original Néel start | Bloch contraction |
| --- | ---: | ---: |
| Parked | about 0.085463% | −0.092868 |
| Half beat | about 0.099872% | +0.646565 |
| Fixed wait 149 | about 0.099872% | +0.646565 |
| Fine energy resolved and forgotten | about 0.092667% | +0.3056 |

The half-beat and fixed-149 local reference records are nonseparable.
Their Bell overlap is about **0.73492**. Their unfiltered CHSH maximum
is about **1.829**, below 2. The parked and fine-energy controls are
entanglement-breaking on the same two-site receiver. Each comparison
uses one chosen neighbor and retains all original-start refusals. The
three Pauli banks in the receipt are simulated draws from the verified
instrument laws. Their simultaneous 99 percent intervals include the
spectral enclosures.

A sharp spin measurement at site 4 alone has a different outcome.
With single occupation there its Bloch contraction at the half beat
is about −0.003905 and its channel is entanglement-breaking.
The two-site read places the signed charge at 4 and the recoverable
spin on a neighboring site.

## Exact spin memory with the original ground-return detector

The supplied local pump and drain are contractions with factors
$1/\sqrt2$ and keep their refusal records. In the hole chart the two
spinor columns are $(c_\downarrow,-c_\uparrow)$. The replay checks the
raising intertwiner on the full CAR space and verifies the two drain
columns independently. After a return to the original singlet, the
accepted matrix is a single complex amplitude times the identity on
the marked spin. It gives a perfect reference Bell report regardless
of spin direction. At site 4 this exact but rarer detector succeeds
on about **0.0317774%** of original starts at the half beat, versus
**0.0037970%** parked and **0.0177872%** with the two energies
resolved and forgotten. The perfect record requires its global
ground-return flag. The simultaneous local charge-and-spin result has
the weaker, nonseparable fidelity and requires only the local gate.

The finite experiment shows a spin-dependent observable that changes
with the native free wait on a declared charged window. Its initial
coherent spectral operation is global and allows some distant
records even when parked. The finite contact ring supplies no
electromagnetic calibration, three-dimensional dispersion, or
electron rest-mass assignment.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.theorem_interfaces.interacting_chain.independent.measure
~~~
