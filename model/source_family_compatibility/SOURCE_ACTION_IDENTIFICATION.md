# Hubbard action readout with complete native records

The configured two-flavour field has c=1/8 and U=1. These are literal
`FieldGeometry` inputs. This result constructs an operational inverse from
that action to complete native repair records, including finite precision.
It also identifies an exact limitation of repair-only parameter selection.
FIFO capacity, finite record capacity, graph size and field coupling are
distinct quantities in the source implementation.

## Finite selection premises and the implemented archive

`RepairWordSchedule.selected_factorizes` applies when the complete
quotient-visible move alphabet, every finite word simplex and counting
reference are supplied. With moves (seam,law,coin), it selects independent
uniform attempts and hence equal C/E weights and fair coins. This theorem
does not contain a coherent Hamiltonian, an occupation trajectory, a
Hamiltonian time per attempt or a state-preparation rule. The fixed-point
uniqueness input of `CapacityClosurePrinciple.closure_value_eq_of_unique`
concerns a supplied source return map. Neither theorem identifies the
interpreter's FIFO slot count with a physical capacity.

For a fixed policy that uses the same complete history, changing FIFO
capacity with `export_oldest` changes only the partition of that history
into memory and archive, by the archive resegmentation theorem in
`particle_identification/NATIVE_CAPACITY.md`. At each step their concatenation is the same
ordered record list. Induction on the actual `_commit` operation proves
equality of every source/reference output block, including zero blocks.
The policy must not access a different truncated history or stop at a
capacity-dependent timeout. This is an equality before discarding any
record. The executable witness runs two actual H66 programs with capacities
one and two, a common two-attempt fixed-seam policy, an entangled input and
all 625 allocated history blocks. It is separate from the two-flavour field
extension; no undocumented embedding of that extension into H66 is used.

## An all-domain invariant of the literal repair instruments

Let D_x=n_(x,up)n_(x,down), N=sum_x n_x. Every full-domain branch K_h
commutes with every D_x and with N. A branch holds at a doubly occupied
endpoint, and a permitted hop transfers a single occupation to an empty
endpoint. The signed CAR hop therefore changes no doublon projector.
Collision holds and all zero raw slots are included.

Consequently a repair-only history with record-dependent seams and laws
commutes with any phase exp(i f(D_1,...,D_V,N)). Inserting such phases
between attempts leaves the complete classical-record/reference marginal
unchanged: on each fine history the source output is conjugated by a unitary,
which disappears only when taking the source partial trace. This statement
retains all classical records and an arbitrary quantum reference. It does
not identify the joint quantum source outputs. A coherent source observer
can detect phases, and hopping does not commute with the D_x. No
indistinguishability of two different nonzero-c Hubbard histories follows.

This proves a precise boundary for uniqueness on the repair-record menu.
The menu cannot select an onsite interaction phase from those observations
alone. A complete source return principle using a larger observable menu
has additional premises. The following construction uses the declared
hopping clock and makes its interaction coefficient observable.

## The actual full-CAR dimer and complete collision detector

Use the signed two-site source from `defect_budget.exact_source`, with
positive declared c,U. Its hopping magnitude per flavour is c/2. In the
total-number-two sector the centered action, up to a scalar, is

    H = T + U D,                 D=D_0+D_1.

The known input is |up,down>, occupation mask9. It is an equal norm
superposition of a dark vector and a coupled spin vector. The source itself
gives the exact coupled block

    H_coupled = [[0,c],[c,U]],

whose second basis vector is a doublon-hole pair. No state or sector is
postselected. This known calibration input is not a thermodynamic ground
preparation. With W=U²+4c², the original collision probability is

    p(t) = c²/W [1-cos(t sqrt(W))].                       (1)

The factor one half from the dark input component is retained. The full
16-state action is used to verify this block and its invariant complement.

For each law C/E and each fair coin, every number-two input has raw0 on
the singly occupied subspace P and raw8 on the pair subspace Q. With the
uniform law choice, the four corresponding Kraus operators are P/2 and
the other four are Q/2. The remaining 56 maps are zero on this input sector.
Thus the 64 fine records have probabilities (1-p)/4, p/4 and zero, and the
complete state/reference instrument is retained. There is no accepted-only
frequency. The detector is the literal `full_domain.branch_action`.

The original Weyl frame is also checked. For each of the three single-link
directions the exact onsite frame intertwines the full 16-state centered
Weyl action with this signed action. Independent unit hopping and contact
pieces are checked, so linearity gives the identity for every c,U. The
needed original occupation inputs
are masks5,5,9, respectively, up to a scalar phase. Every bare repair map
commutes with this frame on the number-two sector. Collision and protected
projectors depend only on local number. This narrow covariance does not
extend the detector to an arbitrary full-domain change of frame.

## Finite-time identification and error budget

Assume specified bounds c<=c_max and U<=U_max, and complete probability
intervals at t and 2t. Multiplying (1) by two gives exactly the source form
of `source_strong_coupling.readout.infer_coefficients`, with d=1, J=c and
epsilon=U. The implementation reuses that rational enclosure theorem; it
does not substitute a numerical derivative.

For clarity, write p1=p(t), p2=p(2t), W_max=U_max²+4c_max² and

    C_hat = (8p1-p2/2)/(3t²),
    Q_hat = (8p1-2p2)/t⁴.

Taylor's theorem gives the deterministic errors

    |C_hat-c²| <= c_max² W_max² t⁴/54,
    |Q_hat-c²W| <= 17 c_max² W_max² t²/90.                (2)

Independent errors at most epsilon on each complete probability add
17epsilon/(6t²) and 10epsilon/t⁴, respectively. The reused interval routine
propagates actual interval endpoints before its Taylor bounds. Once
C_->0 and Q_->=0, one obtains

    U² in [Q_-/C_+ - 4C_+, Q_+/C_- - 4C_-] intersect [0,infinity),
    (U/c)² in [Q_-/C_+² - 4, Q_+/C_-² - 4] intersect [0,infinity).

The executable reuse can give wider outward intervals, which contain
these quantities without any equality assumption about the errors.
The sign U>0 is a source-domain input; these records alone determine U².
Finite-time trigonometric inversion without a frequency domain would have
aliases. Equations (2) use a priori bounds and have no such inversion step.

For t->0 and epsilon=o(t⁴), all coefficient intervals contract to the
positive source values. This supplies a convergent finite experiment at
every stage, rather than an observable derivative at zero. Deterministic
preparation, evolution and detector errors must be added to epsilon using
an actual complete-channel bound. No large native coherent pulse tape or
hardware error certificate is executed here.

For statistical probabilities, independent identical trials at each time
and the collision indicator on every original outcome give simultaneous
frequency error epsilon with failure probability at most delta when

    samples_per_time >= log(4/delta)/(2epsilon²).          (3)

The code supplies a rational sufficient count by replacing the logarithm
with ceil(log2(4/delta)). All noncollision outcomes, zeros and implementation
failure flags remain in the observation protocol. A failure is not removed
from the denominator; a known complete implementation error is charged
separately. The exact simulated probabilities in the certificate are not
reported as sampled observations.

## Source reading versus parameter selection

At c=1/8,U=1 and deliberately larger bounds c_max=1/4,U_max=2, the ratio
squared intervals contract around64. The action is definite and its ratio
is observable with this common source, state, clock and detector.
Replacing (c,U,t) by (sc,sU,t/s) preserves every record probability, so an
absolute energy scale requires the source-time calibration. U/c survives
this clock change.

The finite schedule theorem selects the attempt distribution under its
stated premises. FIFO export preserves complete observations under the
common policy. Neither operation computes c(N), U(N), the normalized ground
prescription, lattice spacing or the duration per global attempt. A source
principle selecting those quantities must be applied on its actual enlarged
observable and return-map domain. The operational coefficient readout is
available to test that map without particle-mass inputs; it supplies no
physical species, mass, topology or relativistic dispersion assignment.
