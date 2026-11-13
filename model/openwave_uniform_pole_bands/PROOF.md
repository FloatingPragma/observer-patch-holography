# Complete momentum-circle pole isolation

The source is the unchanged native strong action with the positive
inverse-rest-matrix correlated propagation rule. The fixed charge, spin and
isospin source has five compact coordinates with Gram diag(6,1,1,1,1),
plus the complete countable relative-position complement. All calibrations,
GeV conversion, supplied physical speed and strong length are inherited.

## A smaller exact finite theorem

Split the compact frame into L=(0,2,4) and R=(1,3). The two R coordinates
have diagonal c and mutual native entry -t/2; the added inverse-rest term
is positive. Their source restriction is therefore at least c-t/2.
The old infinite Q has lower bound C and its coupling to R is at most
sqrt(5)t. Thus the complete new high block on R plus Q has lower bound

    D=min(c-t/2,C)-sqrt(5)t.

Let B0 be the compact L-to-R native hopping block and JLR the corresponding
inverse-rest block. The phase only multiplies entries of B0 by unit complex
numbers; its Gram-normalized Frobenius norm is constant. Since
0<=f(k)=2(1-cos k)<=4, the entire L-to-high squared norm is bounded by

    b²=(||B0||_F+4||JLR||_F)²+5t².

No high coordinate is deleted. These two exact inequalities allow reuse of
the same complete block Schur theorem with a three-dimensional finite part.

## Scalar subtraction and all-real-momentum transport

Write the light symbol as A_L(k)=K(k)+f(k)jG_L, choosing the exact scalar
j=J00/6. It has no effect on gaps or eigenvectors. In normalized coordinates
the remaining inverse-rest matrix is J_L-jI. Its exact Frobenius bound is
small because the three light levels have nearby inverse energies. Only
one native light hopping pair depends on phase; its derivative norm is
sqrt(3/2)t. Consequently

    ||K'(k)|| <= sqrt(3/2)t+2||J_L-jI|| = L.

The charts z=(1+it)/(1-it), -z and their conjugates, 0<=t<=1, cover the
entire circle. Each chart is split into 32 closed intervals. Its angle has
|dk/dt|<=2, so a midpoint symbol differs by at most L times the interval
width. The midpoint symbol and Gram matrix are exact Gaussian-rational
matrices. Realification duplicates eigenvalues; exact Sturm enclosures
retain one member of each pair. The scalar f(k)j is monotone on each chart
and is bounded by its two exact rational endpoints.

For each interval the finite three-level ceiling is below D. The Schur
error epsilon=b²/(D-upper_level_2) controls every low level. If midpoint
reduced eigenvalue intervals are [l_i,u_i] and the transport radius is r,
then the complete source ground gap is at least

    l_1-u_0-2r-epsilon > 0.

The common scalar cancels in this same-momentum difference. The full high
probability is at most b²/(D-upper_level_0)². The three native interpolators
therefore have a rank-one ground spectral residue with trace at least
one minus that probability. All other low and high spectral outcomes
remain. Every positive isolation denominator is checked for every interval,
for both charges. Complex conjugation completes the lower half-circle.

The source is an infinite one-dimensional chain at fixed spacing. Uniform
momentum-fibre poles are an input to normalizable direct-integral packet
records; they do not by themselves identify physical relativistic particles,
a four-dimensional continuum, scattering amplitudes or decay lifetimes.
