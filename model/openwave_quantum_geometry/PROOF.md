# Autonomous native quantum geometry

## Declared source

Each of two cells carries the complete native eight-level matrix
`openwave_energy_lapse.source.native`. Write its positive excitation operator
as K=M−E0 I. Both cells retain all eight levels, including their unique zero
state. No fitted masses enter this construction.

The spatial Dirichlet matrix is L=[[2,−1],[−1,2]]. Declare equal inertia and
stiffness κ=10000. Its normal frequencies are 1 and √3. At product occupation
cap C, both normal modes have occupations 0,...,C. With truncated lowering
operators a,b, put A=(a+a*)/(2√κ), B=(b+b*)/(2√κ 3^(1/4)), Q0=A+B,
Q1=A−B and Hg=N_a+√3 N_b. The complete source is

    H = Hg + exp(Q0) K0 + exp(Q1) K1.

The independent product cap gives [A,B]=[Q0,Q1]=0 exactly. A cap on total
occupation would not preserve this identity at its boundary. Hg includes
the finite oscillator boundary completion: for a normal mode of frequency
ω, the projected quadratic oscillator, less ω/2, equals
ωN−ω(C+1)P_C/2. Adding ω(C+1)P_C/2 gives the declared Hg. No exact finite
canonical commutator is asserted. The spatial matrix, κ, occupation caps,
exponential coupling and time unit are architecture inputs.

## Positivity, energy and autonomy

All three terms in H are positive. Its unique zero state is native00 times
the geometry vacuum. Every nonzero native energy makes its conditional
geometry block strictly positive. Both native K operators commute with H;
their populations are conserved. H generates a linear unitary group on
all 64(C+1)^2 states, including arbitrary external references. Hg does not
commute with the coupled source, so geometry energy can change while total
energy is conserved. No minimization service, refreshed geometry, bath or
irreversible relaxation appears in the dynamics. Finite unitary recurrence
precludes a general convergence claim to a static minimum. Matter occupies
two fixed cells; particle translation is not implemented by this source.

## Complete spectral certificate

The native gaps have the rational interval certificates imported from
`openwave_energy_lapse.source`. At C≤4, each of A and B has norm below 1/50:
the maximal absolute row sum of a+a* is at most √3+2<4, and its coefficient
is at most 1/200. The degree12 exponential series has operator remainder
at most 2(1/50)^13/13!. Applying this bound to every matrix entry is
conservative. Because A and B act on independent factors,
exp(Q0)=exp(A) tensor exp(B), and exp(Q1)=exp(A) tensor exp(−B).
Interval arithmetic at 80 decimal digits encloses all algebraic and
exponential source entries independently of the floating eigensolver.

For each of all 64 native sectors, the eigensolver proposes O,D, rounded
to dyadics with denominator 2^46. Exact integer products bound
η=||O*O−I|| and R=||HO−OD||, including the interval source-to-dyadic error.
The certificate uses Frobenius upper bounds from an integer square-root
ceiling. η<1/2; the polar factor W is unitary and ||O−W||≤η, ||O||≤2.
Every exact conditional H and proposed D has norm below30: Hg≤4+4√3<11,
the native gaps are below4, and ||Qx||<1/25. Consequently

    ||H−W D W*|| ≤ 2 R + 4·30 η = ε,
    ||exp(−itH)−O exp(−itD) O*|| ≤ |t|ε + 3η.

The first inequality follows by adding and subtracting ODO* and using
||OO*−I||=η. Weyl's bound encloses every source eigenvalue by D_j±ε.
The propagator bound follows from the Duhamel identity and the polar
replacement. The maximum block error covers the full native source and
any reference extension. The producer encloses all 49 mutual ground
energies Eij−Ei0−E0j+E00, i,j>0, strictly below zero for C=1,2,4.
These are conditional finite interaction energies, not calibrated GeV masses.

## Mediated clock entanglement and complete observations

Prepare both native clocks in (|0>+|1>)/√2 and the geometry in its vacuum.
This uses four native sectors while retaining the other 60 sectors in H
and in every instrument. The four geometry vectors evolve under their
own complete conditional blocks. Their Gram matrix divided by four is
the reduced two-clock density matrix. At source time10000 its partial
transpose has a negative eigenvalue.

A rational complex witness is proposed from that eigenvector. The
producer evaluates its expectation using interval phases of the exact
dyadic eigenframe. If δ bounds the global state-vector error, then
the density trace-norm error is at most 2δ+δ². The partial transpose of
a normalized rank-one projector has norm at most1, so the same bound
encloses its expectation. Every reported upper endpoint is negative;
this certifies entanglement of the two native clocks in the exact finite
source, not only in a floating eigensolver approximation.

Nine separate product X/Y/Z clock measurements, followed by all geometry
occupation projectors, form complete orthogonal instruments. Each retains
4(C+1)^2 active outcomes plus the 60(C+1)^2 complementary native outcomes,
whose prepared-state probabilities are exactly zero. The same δ bound
encloses every active outcome. Orthogonal record-copy isometries preserve
arbitrary quantum references and copy only the measured classical record.
The producer records measurement work for complete XX/number dephasing;
no conservation of H through measurement is assumed.

Deleting the spatial edge replaces L by diag(2,2). The exact source then
factorizes into independent native-clock/geometry systems. A product
preparation cannot entangle the clocks at any time, for any cap or native
levels. The finite calculation checks this control. Frozen geometry also
gives a sum of independent clock Hamiltonians and cannot entangle them.

## Scope

This is an autonomous quantum extension of the exponential native-energy
coupling, with its own explicitly declared spatial graph, geometry state,
clock and apparatus. It supplies a finite mediator and backreaction
example. It is not the static31-cube source of the planar wave experiment.
The construction does not derive physical G, metric universality,
macroscopic gravity, a continuum Einstein equation, an observed gravity
entanglement experiment or cap-independent dynamics. Cap comparisons are
finite source comparisons; their agreement is not a limit theorem.
