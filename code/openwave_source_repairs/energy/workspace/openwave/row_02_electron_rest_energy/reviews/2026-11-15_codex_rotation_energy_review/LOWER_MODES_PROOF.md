# Exactly eighty lower source modes

This is a retrospective proof about the already frozen OLD240 source graph.
The independent integer checker establishes that it consists of twenty
identical twelve-port body graphs and thirty exterior edges forming a matching.

Let B be one body Laplacian. Its connectedness makes the zero eigenspace
one-dimensional. The checked integer matrix identity

    B^4-16 B^3+80 B^2-120 B=0

restricts its eigenvalues to0,6,5-sqrt(5),5+sqrt(5). The conjugate irrational
roots have the same multiplicity because the characteristic polynomial has
integer coefficients. Dimension12 and trace60 therefore give multiplicities
1,5,3,3 respectively. Traces360 and2280 of B² and B³ supply further exact
checks, without relying on a floating eigensolver.

Write L0 for the direct sum of twenty copies of B. The exterior matching
Laplacian E has disjoint blocks [[1,-1],[-1,1]], so0<=E<=2I exactly. Hence

    L0 <= L=L0+E <= L0+2I.

The minmax theorem places the first eighty eigenvalues of L at most
7-sqrt(5)<6 and all remaining160 at least6. Thus exactly eighty full-source
modes lie strictly below the compact value6. The logarithmic preparation
energy is increasing in L, so its count of modes below the compact sector is
the same.

None is deleted from the source flow. The regional pump's independently
proved contraction removes its admitted input's non-dark population. This is
a dissipative selection process with declared apparatus, not an assertion
that the compact sector is the global lowest charged eigenvalue. The lower
mode count and the selection theorem answer different, compatible questions.
