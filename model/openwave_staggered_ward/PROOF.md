# Compatible finite curvature and the scalar difference graph

The parent linearized Ricci polynomial uses commuting derivatives. Its
collocated centered difference has an extra even-lattice checkerboard kernel.
This package changes the declared record placement and retains the same
curvature polynomial; it does not silently change the parent construction.

Store coordinates in half-spacing units. A scalar lives on the even lattice,
a vector component b on parity e_b, and a symmetric tensor component ab on
parity e_a+e_b. Define D_a f(y)=[f(y+e_a)-f(y-e_a)]/a. The two samples are one
physical spacing apart. Every derivative changes the corresponding parity;
every term of a tensor component therefore has the same target location.
The full zero extension and all derivative halos are retained.

These differences commute and are skew adjoint under the complete finite
sum. Consequently the original four-dimensional linearized Ricci/Einstein
polynomial has exact Bianchi divergence zero and annihilates
h_ab=D_a xi_b+D_b xi_a. Summation by parts makes its quadratic action
self adjoint and annihilates the complete gauge-direction variation. The
program verifies these statements on independent nonstatic compact fields,
while the parent symbolic polynomial proves the arbitrary-field identity.

On scalars D_i D_i is the original nearest-neighbour second difference:
f(x+a e_i)-2f(x)+f(x-a e_i), divided by a². The actual native dyadic Green
profile therefore has exactly its original Laplacian and stiffness, without
the parent's extra fourth-difference correction. The static formulas remain
G00=2gamma Delta phi and Gij=(gamma-1)(DiDj phi-delta_ij Delta phi).
They are linearized metric identities, not a derivation of that metric or
of pressureless physical stress.

For a connected scalar graph, D lambda=0 means equality at the two endpoints
of each actual nearest-neighbour edge. Propagation along graph paths gives
one constant. Equivalently the periodic symbol is 2i sin(ka/2)/a, with no
additional zero inside the physical Brillouin interval. The complete finite
cycle matrices retain the constant mode and expose the old stencil's extra
checkerboard mode. No boundary mode is removed to obtain the rank.

If B=lambda eta+R on the stated staggered locations, every R component is
bounded by delta and every divergence by w, then
|lambda(x+e_j)-lambda(x)|<=a w+8delta. Four differences each contain two
remainder entries. Paths sum these edge budgets on this actual connected
graph. The nine-null decoder concerns co-located tensor entries; applying
its error constant here requires a separate co-location/transport map and
its error. This package does not supply that observer instrument, nonlinear
curvature, physical coordinates, stress or an Einstein continuum.
