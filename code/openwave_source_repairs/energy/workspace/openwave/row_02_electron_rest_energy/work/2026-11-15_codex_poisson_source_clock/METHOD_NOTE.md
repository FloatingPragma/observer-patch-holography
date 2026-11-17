# Method

    P_tick = sum_n exp(−1) T^n/n! = exp(T−I),
    H_count|lambda> =1260(1−a)|lambda>,
    a_single=209/210, a_pair=104/105, H_single=6, H_pair=12.

`clock_source.py:96` executes every native branch of the four single and
four product modes, retaining full output coordinates and current records.
It measures the eigenvalues used above; no additive pair energy is assigned.
`clock_source.py:25` bounds exp(−1) by exact alternating rational arithmetic;
lines47,57 and66 implement the movement-bit count decoder. Line78 implements
the full630-edge selector. The complete16-bit native prefix tree is at
line147 and the actual clocked trajectories at line183. Line229 encloses
finite rare-event energies, preserving their nonzero excess.

All9 finite gates pass: source residuals zero; all65,536 movement words
accounted; Poisson probabilities/transfer factors enclosed; eight complete
histories; zero and multiple counts observed; HOLD-auxiliary counts zero;
fixed-count discrepancy1/44100; positive finite rare-event excess with
O(epsilon) bound; finite sampler failure below10^−30.

The HOLD-source comparison is an analytic altered-law identity T=I,
not a separately executed native-source control. The active fixed-count
comparison uses the actually enumerated native one-event columns.
The HOLD-auxiliary test executes the altered auxiliary rule and retains
its records; those altered events are not counted as native calls.

`verify.py:14` independently allocates endpoint floors/ceilings.
`verify.py:25` verifies the full record archive, native call accounting,
prefix probabilities, actual seed streams and energy intervals. It imports
no producer code. Full first execution3.352seconds; separate check1.305seconds.

No complete coherent U_L program, phase reader, physical proper-time
normalization or rest-energy observable is implemented by this packet.
The full absorbing pair matrix is not called a self-adjoint Hamiltonian;
the proved scalar generator blocks are the consumer interface.

**Scope:** original-source compact generator and an executed native
event-count clock with complete finite errors; no full physical row award.
