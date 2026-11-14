# Record preparation and supplied energy operators

This package independently implements finite cases of the J02 register and
J03 spectral-energy equations. It imports no worker result receipt and does
not execute a reserved acceptance configuration. The controls use the oriented
triangle `((0,1),(1,2),(2,0))`, mass 0.23, and hop time 0.37.

`record_preparation.uniform_fibers` maps each bare position to the normalized
uniform superposition of record states at that position. For the integer
cutoff-1 `JointRecords` sector this is an isometry, but `H_joint V = V H_bare`
fails. A normalized preparation alone does not license transfer of the bare
spectrum, clock, or compact-state persistence.

`modular_control` builds the complete connected sector with half-seam registers
in Z_3. The sector has 18 states, with three record states per position. The
uniform-fiber isometry obeys `H_mod V = V H_bare` and its time evolution agrees
with the bare model after that map. This is a positive control on a distinct
register algebra: cyclic shifts are unitary and uniform fibers normalizable.
Uniform unbounded integer fibers are not normalizable; a bounded integer shift
has a wall. The modular control is not a replacement for the candidate integer
state space. Definite modular records give different position probabilities
from the flux-free preparation, despite the same generator and Gauss sector.

For a port incidence matrix B, `RecordEnergy` defines the supplied functional

```text
U_phi(E) = (1/2) E^T phi(B^T B) E;
rho = B E;
P_cycle = I - B^T L^+ B,  L=B B^T;
U_phi(E) = (1/2) rho^T phi(L) L^+ rho
           + (phi(0)/2) ||P_cycle E||^2.
```

A singular value decomposition gives this identity directly: positive singular
modes contribute through rho and zero singular modes through the cycle
projection. The implementation assembles the energy kernel in seam space and
the charge kernel independently in port space. It reports both contributions.

The supplied choices are `first_order`, phi(lambda)=sqrt(lambda);
`quadratic`, phi(lambda)=1; and `massive`, phi(lambda)=sqrt(m_energy²+lambda).
For the first choice phi(0)=0 and every record with a given rho has the same
energy. This is an identity, not a relaxation or minimization prescription.
For the massive choice phi(0)=m_energy, so loops acquire a positive energy.
Even if `m_energy` is chosen equal to the propagation mass, the spectral record
functional is an additional declared operator.

The triangle control fixes rho=(1,-1,0) and compares all admissible integer
records `E=(1,0,0)+k(1,1,1)` for k=-2 through 2. The integer quadratic minimum
is 1/2 at k=0. Grounded rational inversion gives the real minimizer energy 1/3;
the integer expression is a convex quadratic in integer k, so this is also
the global integer minimum. First-order energies agree for every k, while
quadratic and massive energies distinguish their cycle content. No infinite
screen coefficient or physical coupling is inferred from this triangle.

`joint_diagonal` evaluates this functional on actual `JointRecords` states.
For a committed port state, `E_e=F_ea=-F_eb` is integer. For an occupied transit
seam, the supplied extension is `E_e=(F_ea-F_eb)/2`. Its charge identity is

```text
B E = g_port - n_port - (1/2) unsigned_incidence * n_transit.
```

The half unit at each endpoint is an effective charge of this extension;
it does not split the fundamental integer load or prove fractional particles.
The code checks the identity on every state of the finite sector.

`augment(joint,coupling)` adds the diagonal energy to H_joint. It preserves the
sector and Hermiticity. It changes the generator and its square; the free
massive-square or bare-clock claim cannot be reused. The returned experiment
compares complete free and augmented position distributions for the same
prepared state. The energy choice, coefficient 0.29, transit extension and
preparation are supplied. No calibrated force, moving two-charge dynamics,
string breaking, or physical tension follows from these checks.

Implementations: [record_preparation.py](../oph_openwave/record_preparation.py),
[record_energy.py](../oph_openwave/record_energy.py).
Independent controls: [test_record_operators.py](../tests/test_record_operators.py).
The results bear on state-transfer prerequisites for rows 3/4/30 and
interaction-energy prerequisites for rows 20/22. They establish no full criterion.
