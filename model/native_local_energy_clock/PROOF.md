# Local exchange energy on the fair source

Consider the declared connected binary federation, with V ports, M seams,
maximum port degree Delta at most six, and fixed occupation q. Its exact
population transfer and counted transfer clock are

    G = sum_e h_e,     h_e = I-S_e,
    P = I-G/(2M),      H = -2M log P,      t = N/(2M).

S_e exchanges the endpoint occupations. The source's triangle cover proves
alpha I <= P <= I with alpha=10/21 in every occupation sector. On the original
twelve-port q=2 source alpha=2/3 suffices. These are the operators from
`native_exchange_clock`, not fitted Hamiltonian coefficients. The reference
used throughout the controlled-clock comparison is

    tau=I/D,   c=Tr(tau P),   rho=P/(Dc),   D=binom(V,q).

The original twelve-port preparation of rho=P/56 belongs to
`native_phase_experiment`. This formula specifies one faithful finite
reference at each federation level; it does not prepare that entire growing
reference. Native C/E repair, its record clocks and controlled unitary
continuation are distinct channels on this same algebra.

The full binary transfer floor can be checked without a spectral fit. On
three binary ports the triangle operator T=sum_(three edges)(I-S_e)
satisfies T(3I-T)=0. Since T is positive, T<=3I. Each twelve-port cell has
twenty triangular faces and each internal seam belongs to two faces, so
its thirty internal seam terms sum to at most 30I. Each cell has three
external half-seams. Their total number is 3C/2 for C cells, and each has
norm at most two. Thus G<=33C I while 2M=63C, proving P>=10I/21.
The binary local dimension is essential: three distinct local labels admit
the fully antisymmetric representation, where T=6. The triangle argument
does not extend to an enlarged local matter carrier by notation alone.

## Exact finite moments without a global matrix

Write A=sum_i deg(i)(deg(i)-1), the number of ordered distinct adjacent seam
pairs, and B=M(M-1)-A. Set

    u = q(V-q)/(V(V-1)),
    v = q(q-1)(V-q)(V-q-1)/(V(V-1)(V-2)(V-3)).

For V>=4, direct trace counting gives

    Tr(tau G)=2Mu,
    Tr(tau G^2)=4Mu+A u+4B v,
    v_tau=Var_tau(G)=4Mu+A u+4B v-4M^2 u^2,
    c=1-u.                                                       (1)

Indeed h_e^2=2h_e and Tr(tau h_e)=2u. Two distinct adjacent swaps multiply
to a three-cycle, giving Tr(tau h_e h_f)=u. Two disjoint swaps give
Tr(tau h_e h_f)=4v. This includes every occupation sector, with the
vanishing factors at its extremes. Only four-site counting is needed.

Since rho=tau P/c, the same reference has

    mu=Tr(rho G)=2Mu-v_tau/(2Mc),
    Var_rho(G)<=v_tau/c.                                         (2)

For the second inequality, center first at Tr(tau G), use P<=I, and then
minimize the centered variance over real centers. No independence of seam
energies is assumed. On the actual federation

    V=240*4^ell,   M=630*4^ell,   A=5400*4^ell.

The last count follows from three degree-six ports and nine degree-five
ports in each of the 20*4^ell cells.

A simple uniform estimate establishes the rate. For V>=6,
v<=u^2 V(V-1)/((V-2)(V-3)), u<=1/3, and
(4V-6)/((V-2)(V-3))<=16/V. Thus v-u^2<=16/(9V).
Using A<=2M(Delta-1), M<=V Delta/2, and dropping the nonpositive term
-4A v in (1) gives

    v_tau <= (4/3)M+(2/3)M(Delta-1)+(64/9)M^2/V
           <= ((38 Delta+6)/9)M.                                (3)

The exact expression (1), rather than this loose majorant, enters the
executable bound.

## One local controlled clock at finite density

Use ||X||_(2,rho)^2=Tr(rho X*X). Define

    a = 1/(1-mu/(2M)),       alpha_H^t(X)=exp(itH)X exp(-itH).

For every matrix X on the complete sector,

    ||[H,X]-a[G,X]||_(2,rho)
      <= ( ||[G,X]|| sqrt(v_tau/c)
                      + (1/2)||[G,[G,X]]|| )/(2M alpha^2).       (4)

To prove (4), diagonalize G with eigenvalues g_i. For
h(x)=-2M log(1-x/(2M)), h'(mu)=a and
0<=h''(x)<=L=1/(2M alpha^2) on the entire actual spectrum. Its divided
difference obeys

    |(h(x)-h(y))/(x-y)-a| <= (L/2)(|x-mu|+|y-mu|),

including x=y by continuity. Put Y=[G,X]. Entrywise comparison followed
by the triangle inequality for the weighted Euclidean matrix norm gives

    ||[H,X]-a[G,X]||_(2,rho)
      <= (L/2)(||(G-mu)Y||_(2,rho)+||Y(G-mu)||_(2,rho)).

The second term is at most ||Y||sqrt(Var_rho G). Commuting G past Y in
the first term adds ||[G,Y]||_(2,rho)<=||[G,Y]||. Equation (2) proves (4).
The transfer floor covers the complete spectrum, so no discarded spectral
tail or bounded-correlation-degree restriction is needed.

Because H and G commute and rho commutes with both, their adjoint groups
are isometries in this norm and commute with each other. Integrating their
difference therefore gives

    ||alpha_H^t(X)-alpha_(aG)^t(X)||_(2,rho)
      <= |t| times the right side of (4).                        (5)

All commutators on the right act on the original X. Evolving X does not
introduce an unaccounted support-growth factor in this estimate.

For ||X||<=1 supported on s fixed ports, at most s Delta seams meet X.
Each h_e has norm at most two, so

    ||[G,X]|| <= 4s Delta,
    ||[G,[G,X]]|| <= 16s(s+1)Delta^2.                            (6)

The second inequality expands over the first seam and then seams meeting
its union with the support of X. Compression to a global number sector
preserves these bounds for number-preserving local X. Equations (3)-(6)
give a dimension-independent O(M^-1/2) bound at fixed s, Delta and t.
The independently computed exact rational square-root upper bound in the
code avoids rounded numerical moments as theorem premises.

If q/V tends to p, (1)-(3) imply

    a -> a_p=1/(1-p(1-p)).                                       (7)

The code retains the exact finite correction |a-a_(q/V)| and adds
|t| |a-a_(q/V)| ||[G,X]|| to (5). This compares the actual logarithm with
one density-dependent local exchange clock. The coefficient is common
to all local observables at that density. Different densities supply
different states and clocks; no density choice is selected by the formula.

The result is stronger than a bounded-correlation-degree operator estimate
but weaker than operator-norm or complete-instrument convergence. It
controls local expectation differences for any state eta<=C rho by
sqrt(C) times (5), since Cauchy-Schwarz gives
|Tr(eta Y)|^2<=Tr(eta Y*Y)<=C Tr(rho Y*Y). In particular every unitarily
rotated rho has C<=1/alpha. It supplies no uniform comparison for arbitrary
adversarial states or every complete detector history.

## Complete fixed finite controlled detector laws

There is a finite complete-observation extension with explicit outcome cost.
Fix n local Kraus instruments, each with its entire declared finite outcome
menu, acting after controlled intervals t_1,...,t_n. Every branch K_(j,o)
is supported on s_j ports, preserves the fixed source number, and has norm
at most one by instrument completeness. Zero-probability outcomes, failures
and erasures remain in the menu. Set

    b_j=4s_j Delta,       c_j=16s_j(s_j+1)Delta^2.

For any fixed future record starting at instrument j, its Heisenberg effect
E_j has norm at most one. Leibniz's rule, and commutation of G with either
controlled evolution, give

    ||[G,E_j]|| <= 2 sum_(k>=j)b_k = 2B_j,
    ||[G,[G,E_j]]|| <= 2 sum_(k>=j)c_k+4B_j^2 = 2C_j+4B_j^2.     (12)

To verify (12), write the history effect as a product of the Kraus factors
and their adjoints, with intervening unitary conjugations. Every Kraus
factor appears twice, each factor has norm at most one, and each unitary
commutes with G. Differentiating twice gives the sum of each factor's
second derivative plus every ordered pair of first derivatives; the
displayed square safely includes all those terms.

Hybridize the n intervals between H and aG. In each hybrid the earlier
unnormalized source branch sigma is a product of Kraus contractions and
unitaries applied to rho. Since rho<=I/(Dc),

    0<=sigma<=I/(Dc)<=rho/alpha,       Tr sigma<=1.

Thus a difference Y in the future effect has
|Tr(sigma Y)|<=alpha^(-1/2)||Y||_(2,rho). Applying (4)-(5) and (12)
to the interval being changed bounds every complete history probability by

    epsilon_word <= alpha^(-1/2) sum_j |t_j|
       [2B_j sqrt(v_tau/c)+C_j+2B_j^2]/(2M alpha^2).             (13)

If L is the product of the complete outcome-menu sizes, the total variation
distance between the two complete classical laws is at most
min(1,L epsilon_word/2). All outcome probabilities are used in this bound;
there is no postselection or discarded failure. At fixed n, support, times
and finite menus it tends to zero at O(M^-1/2). The bound can be expensive
when the outcome menu grows, and says nothing uniform about a growing
archive. It compares classical detector records from initial rho; it is
not a diamond bound on arbitrary inputs or a joint terminal-reference-state
comparison. The source floor and the declared complete local instruments
are premises, both testable independently of this scalar bound.

## Conserved local energy for this controlled exchange

Define the Hermitian seam energy current

    J_(e->f)=i[h_e,h_f].

It is antisymmetric and vanishes for disjoint source seams. Its support is
at most three ports. Since [h_e,h_f]=[S_e,S_f], its norm is at most two.
For controlled exchange evolution, the exact identities are

    d h_e/dt=i[a_p G,h_e]=-a_p sum_f J_(e->f),
    d G/dt=0.                                                   (8)

Summing over any seam set leaves exactly its complete outward current.
Site energies epsilon_x=(1/2)sum_(e incident x)h_e sum to G. Their currents
i[epsilon_x,epsilon_y] vanish beyond port-graph distance two. These are
actual energy currents for the source-derived exchange interaction. No
component is discarded by imposing a nearest-neighbor current ansatz.
Equation (5) relates their bounded local observables to the transfer
logarithm in the same reference-weighted sense. It does not assert exact
local continuity for finite-M H with the h_e densities.

This is an interacting hard-core exchange energy on the source graph.
Its spatial refinement, dispersion and tensor transformation law are
additional physical questions. A graph energy continuity identity does
not supply symmetric Lorentz stress or a gravitational field equation.

## Actual repair has a positive energy input on the same reference

Let Phi be one unread fair uniform native C/E repair, retaining the same
raw Kraus maps K_b before summing the record. It is trace preserving and
unital. The faithful tracial commutator identity used by
`native_fixed_algebra` gives

    Tr[G(G-Phi(G))]=(1/2)sum_b ||G K_b-K_b G||_HS^2.

Consequently, on the same rho=P/TrP,

    Tr[G(Phi(rho)-rho)]
      = Tr[G(G-Phi(G))]/(2M TrP)
      = sum_b ||[G,K_b]||_HS^2/(4M TrP) >= 0.                   (9)

The equality follows by substituting rho and using Phi(I)=I. Expanding
the squared commutators proves the identity using trace preservation and
unitality. Self-adjointness of the complete channel is not required.

For a nontrivial connected fixed-number sector G is nonscalar and the
owning fixed-algebra theorem makes (9) strictly positive. It cannot be
written as a source-only boundary flux on a closed graph. It is the
explicit energy input of the dissipative repair operation when energy is
read by the derived exchange observable. Source-plus-apparatus energy
conservation requires a separately specified apparatus energy operator;
retaining records alone does not select that operator.

Thus (8) and (9) are two different dynamical statements on one carrier and
one reference. The positive local comparison (5) provides a controlled
transfer-clock route to (8). It does not replace passive repair by a
unitary or infer an autonomous physical clock from service pulses.

The same repair raises the exact logarithmic energy. Diagonalize G and
write W_ij=sum_b |(K_b)_ij|^2. Completeness and unitality give unit row
and column sums. Write h(x)=-2M log(1-x/(2M)). Both g and h(g)-g are
increasing lists. A doubly stochastic mixing cannot increase their sorted
pairing, hence

    Delta<H>-Delta<G>
      = sum_ij (h(g_i)-g_i)(delta_ij-W_ij)g_j/(2M TrP) >= 0.     (10)

For completeness, write the increasing list h(g_i)-g_i as its initial
value plus a sum of nonnegative successive increments. Each increment
multiplies a difference between the sum of the largest k entries of g
and a weighted sum with weights in [0,1] adding to k. That difference
is nonnegative. This proves the mixing inequality without symmetry of W
or a selected Kraus representation. The independently enumerated H66
complete instrument gives

    TrP=56,       Var_tau(G)=1220/121,       Delta<G>=37/840.

This is visible in principle to the compiled complete energy instrument
of `native_energy_measurement`. Assign every processor failure and final
residual the score zero. Its 66 energy labels lie in [0,30], its complete
half-diamond implementation error is less than 1/12000, and its nearby
energy operator differs from H by at most 27/10^8. With rational labels
accurate to eta, comparing independent pre-repair and post-repair settings
therefore has a systematic expectation error at most

    4*30/12000 + 2*(27/10^8+eta).                              (11)

The factor four retains the bounded-gain stored-map convention; it does
not use a cancellation requiring exact normalization. For eta<=10^-30,
the expected complete score increase exceeds 0.034. This statement has
the ideal rho preparation premise. A separately certified half-trace
preparation error epsilon adds 120 epsilon to (11), and statistical error
requires independent repetitions. The owning numerical phase-preparation
admission is not promoted into that certified preparation error.
