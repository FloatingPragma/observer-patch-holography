# Same-sector selection and spatial return

Let D=(F,G) on the full240-port carrier, with the integer modes in DESIGN.md.
Their Gram matrix is [[4,2],[2,4]], their full-graph Laplacian eigenvalue is6,
and the first-cell source projector is

    R = D [[1/3,-1/6],[-1/6,1/3]] D*.

The body permutation induced by Ry=diag(-1,1,-1) sends F to -F and G to G-F.
It therefore preserves their span and R, including all off-diagonal entries.
The density R/2 has equal probabilities1/6 on the six ports4,5,6,7,8,11.
It differs from I12/12 by position TV1/2. Rx sends the last two of those ports
to9,10 and keeps the other four within the set, so its position TV is1/3.
This gives directional information that uniform reinjection does not supply.
For any nonzero y-turn smaller thanpi, the y=±1 support points do not return;
the frozen finer-angle tests use the actual coordinates.

Because R delta0=0, the initial source has no component or coherence in D.
For the exact ideal source Q, D is reducing and Q|D is the scalar phase
exp(-4iE6). Uniform reinjection is block diagonal across D and D-perp and
adds a scalar multiple of R in D. Induction therefore gives

    rho_N = (1-w_N) R/2 + sigma_N,
    sigma_N >= 0, R sigma_N=0, Tr sigma_N=w_N.

The source energy certificate bounds w_512 by4^-8. Hence the trace distance
from rho_512 to R/2 equals w_512, without the square-root loss required for
arbitrary cross-coherent inputs. Its ideal Ry trace distance is at most w_512
because only the positive bad-sector remainder changes. The finite compiled
map need not preserve these identities exactly: if its complete error is
delta, the conservative rotated-law bound is2(w_512+delta). Route failures
and terminal-reader error are added with their original probabilities.
No finite exactly symmetric output is asserted from a small numerical residual.

This invariant target alone is insufficient evidence for the mechanism:
I12 reinjection is itself Ry-invariant. The executed dephasing control retains
the same native flow, monitor and reinjection but commits every source position
after each epoch. Its full position kernel is

    K_ij = |Q_ij|² + escape_j/12.

It includes all escape/reinjection events. This memory change preserves the
geometry and reader, while removing the dark-sector interference. The packet
requires its observed Ry asymmetry to exceed the full implementation budget.
If that gate fails, this proof supplies a conditional symmetry identity only;
it does not establish the joint dynamical-return claim.

The oriented frame is the actual native current displacement x1-x0=(2,0,0)
recorded while routing0→1; a second native movement returns the source to0.
Ry reverses that arrow, while its square returns it. The ordinary vector
transformation is supplied kinematics. The field's shorter return must be
established by the pump result and its discriminating active-escape control.

All statements concern the reduced source and its complete position readings.
The controller archive retains the route, escapes and reinjection addresses;
its correlations need not share the source marginal's Ry invariance. The same
source energy sector governs every D-supported density, including R/2, so this
join can use the energy packet's restricted phase/readout bounds. No conclusion
about physical rest energy or a different source/clock is added by this lemma.
