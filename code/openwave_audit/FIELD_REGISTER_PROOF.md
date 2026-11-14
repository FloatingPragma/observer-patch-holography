# Integer field records from finite cyclic registers

Consider a finite graph with integer incidence matrix B and integer face
columns A_f satisfying B A_f=0 and |(A_f)_e|<=1. Matter transport carries
integer edge current and changes each edge by a total absolute amount at
most J over the specified observation interval, including every intermediate
prefix. The field factors are electric diagonal unitaries and magnetic
factors exp(i beta_f (U_f+U_f*)/2), where U_f adds A_f. Every coefficient,
update order, preparation and clock is shared between the two constructions.

The integer construction uses the field Hilbert space ell²(Z^E). The cyclic
construction uses (Z/mZ)^E for even m. Its centered representatives are
-m/2,...,m/2-1. Its electric phase agrees with the integer electric phase
inside that window; for example both use alpha sum_e F_e²/2. The centered
basis map is an isometric embedding into the integer space.

Suppose the initial field has support |F_e|<=E0 and set d=m/2-E0-J>0.
Let Lambda be the sum of |beta_f| over every face and every magnetic update
in the observation interval. The norm distance between the embedded cyclic
output and the integer output is at most

delta = min(2, 2 sum_{k=d}^infinity Lambda^k/k!).

The bound holds uniformly on the admitted initial subspace with any quantum
reference. Every common complete terminal instrument has total-variation
distance at most min(1,delta). This includes the complete joint field and
matter record, provided its output labels use the same centered embedding.

Expand each bounded magnetic exponential in its shift generator and retain
all interspersed matter and electric unitaries exactly. A term of total
magnetic order k<d has |F_e|<=E0+J+k<m/2 at every prefix. The cyclic and
integer executions therefore agree term by term, including electric phases.
At order k the total operator-norm bound is Lambda^k/k!, by the multinomial
identity for the ordered product of exponential series. The two omitted
tails give the stated factor two. Each electric exponential is a diagonal
unitary; this argument does not require a bounded integer electric energy
or a domain theorem for a sum of unbounded generators.

For initial leakage probability epsilon outside the admitted support,
decompose a purification into its two orthogonal support components. The
additional vector-norm error is at most 2 sqrt(epsilon). All errors are
capped at two. Native pulse synthesis, reset and readout errors must be
added using their own admitted domains and complete instruments.

The bound is state dependent. A cyclic shift maps |m/2-1> to |-m/2>, while
the integer shift maps it to |m/2>. The embedded outputs are orthogonal.
Thus increasing m does not produce a uniform operator approximation over
all initial states. This boundary example also prevents replacing the
integer electric/shift commutator with a finite-dimensional canonical
commutator without a state-dependent estimate.

The executable evaluates the positive tail with exact rationals. For an
integer k0>=d with k0+1>Lambda, sum orders d,...,k0-1 exactly and use the
upper remainder (Lambda^k0/k0!)/(1-Lambda/(k0+1)). The source coefficient
can be passed as the exact rational represented by its actual floating
control. No downward decimal approximation enters the certificate.

For three updates with beta equal to binary64 .41 on one face, E0=0 and
J=3, modulus32 gives delta<5.2e-9. Activating208faces at every update
changes Lambda by208; the same modulus32 bound is then vacuous. The face
count, observation interval and initial support are part of the theorem.

This is a finite-time compact-rotor transfer. An integer field register
does not by itself supply a Maxwell vacuum, a canonical real-field limit,
physical units or particle species. A spatial continuum claim must specify
its state, clock, coefficients, field normalization, retained observable
weight and error along a common refinement family.
