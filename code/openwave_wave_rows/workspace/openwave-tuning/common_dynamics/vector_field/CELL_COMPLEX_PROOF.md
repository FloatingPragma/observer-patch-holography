# Local faces on the routed carrier graph

Each cubic cell carries the twelve vertices, thirty edges and twenty triangular
faces of the icosahedron. For every positive cubic direction the routing has
two external channels. Choose a spanning tree inside each icosahedron and a
root vertex. A routed coarse edge is the internal tree path from the first
root to its outgoing port, the actual external seam, and the internal path
from the incoming port to the second root.

The difference between the two coarse-edge chains bounds an intercell ribbon.
For each pair of cubic directions, the four oriented coarse-edge chains bound
a macro plaquette. These attachments give twenty internal faces, three ribbons
and three macro plaquettes per cubic cell. Their signed edge vectors are the
columns of A. With incidence B oriented from the lower port index to the
higher, BA=0. The current convention satisfies n+B F=G.

The attachments preserve the actual carrier graph. A ribbon uses two adjacent
cells, and a macro plaquette uses at most four. The development member of size
two has 96 vertices, 288 edges and 208 faces. No face in that member uses more
than sixteen edges. The face attachments are declared geometry; they are not
deduced as a unique consequence of the bare graph.

Contract the chosen internal trees. The internal triangular boundaries kill
all nineteen integral cycle classes inside each cell: the icosahedral sphere
has integral first homology zero. Each ribbon then equates the two parallel
coarse channels. The macro plaquettes impose the square relations of the
periodic cubic two-skeleton. Its first integral homology is Z cubed. Changing
the internal tree paths changes the coarse chains by internal boundaries and
leaves this quotient unchanged.

The executable checks the integral assertion directly for the size-two member.
Deleting a global spanning tree supplies 193 integral graph-cycle coordinates.
In those coordinates the face matrix has 190 nonzero Smith factors, all one.
The quotient is free of rank three. The three coordinate-period functionals
annihilate every face. Three explicit winding cycles have the identity period
matrix, so each primitive torus winding is present.

Keeping only the internal triangular faces gives rank 152 and leaves 41 cycle
classes. The intercell attachments remove 38 of these classes. Thus the
single-carrier charge-lattice statement cannot be extended to the federation
by ignoring its glue seams. The completed geometric quotient retains the
three noncontractible global windings required by the periodic domain.

The attachment supplies local divergence-free field shifts and a well-defined
distinction between face cycles and global periods. It supplies no autonomous
field update, field energy, wave speed or particle-force measurement. Such
dynamics requires an explicit extension of the same source and its complete
instrument. The exact finite checks and input hashes are in
`CELL_COMPLEX.json`.
