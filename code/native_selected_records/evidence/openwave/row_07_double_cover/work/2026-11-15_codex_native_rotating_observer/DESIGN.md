# Native incidence rotation: prospective exploratory design

Candidate **A-native-incidence-rotation-v1**, descendant of
A-native-triangle-optical-holonomy-v1. Written before this candidate's first
execution. The earlier polar96 result is known; this is a different observer,
preparation and coupling, not a retrospective gate for that result.

## Object and observer

The target has one integer unit on three ports, embedded as any specified
triangle (default labels 0,1,2; proposed common-source labels 0,12,13). Its
quantum state is an arbitrary density operator on these three classical basis
configurations. An invalid fourth basis state has its own retained rejection
flag. A second one-unit register is the observer marker. A two-port vacuum
register supplies a native zero-load reference event. A path qubit, coherent
choice/coin/current memories, real orthogonal record erasers, and path
interference are declared architecture. They do not install a spin matrix,
angle-dependent sign, or polar cosine coupling.

The native endpoint puts floor((a+b)/2) and ceil((a+b)/2) on an addressed seam,
in fair coin order. The fine recorded current is old left minus new left.
The alternative swaps every addressed seam, including balanced seams, using
the same coin and current recorder. Both laws preserve the triangle symmetry.

The marker starts at vertex 2. Its controller addresses the next directed
edge of 2->0->1->2, repeating a native event until a nonzero current confirms
the hop, or a frozen cap terminates the run. Its unwrapped turn is the signed
sum of the three separately recorded 120-degree moves. The equilateral
embedding, orientation and step-to-angle dictionary are declared geometry;
neither the target phase nor a half-angle is part of that dictionary.

After each marker hop, the rotating arm addresses only the two target seams
incident to the marker, fairly, for h events. The reference arm addresses the
two seams incident to the initial marker. Both arms use identical fine flag
alphabets (choice,coin,current); the coherently computed address is uncomputed
using path, marker and choice, rather than left as unobserved which-path memory.
The complete fine record is transformed to a Helmert orthonormal basis. Flag0
is the uniform eraser; flags1..11 are all retained as termination outcomes.
The final path +/- readout is the phase measurement. No hidden target state
enters the decoder; state vectors in the receipt are explicitly audit data.

## Process preparation, including previously fatal common-mode contamination

One committed-current event addresses all three target seams uniformly. No
outcome is discarded. A fresh path compares one native event on edge(0,1)
with a native event on the auxiliary vacuum. Erasing its six fine flags and
reading path-minus gives K=(I-P01)/(2sqrt(3)). The native endpoint calculation,
not an installed contrast projector, implies I-P01=Q01. All 18 warmup flags
and all 12 subsequent (eraser,path) flags are retained. The accepted output
is the exact edge contrast for every admitted input density operator; its
probability is at least 1/72. Under the swap law it remains possible and
selects the same contrast. Complement states are explicitly rejected.

## Analytic prediction and known assumptions before pilot

On the optical plane, the native incident-star mean has eigenvalues 3/4,1/4;
the slow eigenline is the contrast between the two unaddressed leaves. No
eigenprojector is applied during rotation. Put r=3^-h. The scaled full-turn
product has determinant r^3 and trace (-1+9r+9r^2-r^3)/8. This predicts a
negative dominant eigenvalue at h=4,6, but does not predict a negative
eigenvalue at h=2. The swap star has optical eigenvalues +/-1/2, hence every
even-h frame block is a positive scalar on the whole optical plane. It should
give positive full-turn interference with the identical reader/controller.
This is a target-sign distinction, not only a magnitude distinction.

The initial contrast need not be the exact finite-h full-turn eigenvector.
Finite ray return error, unequal arm amplitudes, all rare-success costs and
finite preparation error must therefore be measured. The first pilot has no
success-dependent selection: h=1,2,4,6, one/two turns, both directions and both
laws are all retained. Gates for a later confirmatory run will be separately
frozen after this explicitly exploratory receipt. Initial pilot cap is8.

## Falsifiers, controls and interpretation

Representable failures include nonnegative native full-turn visibility,
nonpositive two-turn visibility, poor ray return, a swap control with the same
minus, preparation output retaining a constant component, incomplete Kraus
normalization, and loss of the effect after label/direction changes. Timeouts
are inconclusive outcomes and remain in unconditional probability accounting.
No-time-marker hold and deliberately planted path-minus test respectively
missing motion and reader sensitivity, not emergence. Committing fine records
instead of erasing them tests the role of coherent memory. Uniform-input and
common-mode-contaminated inputs challenge the former preparation defect.

This tests a recorded actual observer circuit with local incidence coupling,
not a freely named parameter loop. A successful result would be a finite
operational full-turn sign/two-turn return mechanism on this candidate. It
does not by itself prove arbitrary-axis SO(3) covariance or exchange statistics.
It borrows no energy or charge award, and leaves the previous separated-pair
energy excess to row2. A common-source join requires separate implementation.
No G6/G9 or upstream pin change is claimed.
