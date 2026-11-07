# Weighted repair across connected cells

The endpoint translation and complete C/E raw codes are the literal maps in
[`native.py`](../observer/native.py),
[`local_bias.py`](../local_repair_family/local_bias.py) and
[`interpreter.py`](../native_program/interpreter.py). Positive endpoint weights
are inputs to the local law. The connected extension assigns such weights
to every vertex of a specified graph. Its source graph can be the retained
`port_pair` federation, without changing its seams or routing.

The native program class admits weighted local cells and fair connected
federations. Its admission function rejects this weighted connected extension.
No contract or interpreter law is changed here. The construction establishes
consistency and exact comparisons for this declared extension; it does not
select the extension, a physical capacity trajectory, or an interacting
relativistic field theory.

## Complete endpoint instrument

Let G=(V,E) be a finite connected simple graph with an ordered vertex set,
positive integer vertex weights w_i and positive seam clock weights c_e.
The carrier is the complete binary tensor space; each fixed-number sector
H_m is invariant. Select e=(a,b), a<b, with probability
q_e=c_e/sum_f c_f. A placement coin c has probabilities

\[
 p_e(1)=\frac{w_a}{w_a+w_b},\qquad
 p_e(0)=\frac{w_b}{w_a+w_b}.                                      \tag{1}
\]

Call the existing native endpoint function on (n_a,n_b). It conserves their
sum, places a singleton at a for coin 1 or at b for coin 0, and fixes 00 and 11.
C records the signed current n_a-n'_a. E records n_a-n_b and decodes the same
current from that contrast and the retained coin. The supplied C/E mixture
weight affects record backaction, while the population kernel is identical.

For each coin and raw code use the existing decoding table. An admitted code
defines the partial injection from compatible binary input configurations to
their native outputs; other codes give the zero operator. Multiplication by
the square root of the clock, mixture and coin probabilities gives its Kraus
operator. All four C codes and eight E codes per coin are retained, together
with the identity hold outcome, including zero weights. A fixed-seam mixture
thus has 25 blocks; the twelve-port thirty-seam mixture has 721 blocks. These
codes include the impossible binary contrasts and the reserved abort codes.

Each input belongs to exactly one code for each selected law and coin. Thus
the sum of K†K is I. Every branch preserves total number and untouched factors.
The same identities hold after tensoring every branch with I_R for any finite
reference R. The executable comparison checks every native H66 branch domain,
image and positive coefficient against the unmodified interpreter. Since the
branches are scalar partial injections with no extra phases, this determines
their action on every complex matrix unit, including entangled inputs.

## Reference selected by the specified law

Write e_m(w) for the degree-m elementary symmetric polynomial. For a binary
configuration S with |S|=m, define

\[
 \pi_m(S)=\frac{\prod_{i\in S}w_i}{e_m(w)},\qquad
 \rho_m=\sum_{|S|=m}\pi_m(S)|S\rangle\langle S|.                    \tag{2}
\]

All weights are positive. For two configurations differing by a singleton
transfer from a to b, their probability ratio is w_b/w_a; their transition
probability ratio is the same. Therefore each seam obeys detailed balance
and fixes (2). Every nonadaptive fixed mixture of seams also fixes it. This
claim concerns the unread source channel: retaining the record produces
source-record correlations.

For a connected graph with positive seam probabilities, adjacent occupied/
empty swaps connect all m-subsets. For 0<m<|V| this follows by moving holes
along graph paths; the empty and full sectors have one configuration. All
configurations have positive self-loop probability. Consequently (2) is the
unique stationary population law on H_m and the population chain converges
to it. No uniform mixing time or unique stationary coherent state is asserted.
An outcome-adaptive selector is a different joint source-controller chain;
stationarity of (2) alone does not follow for it.

This is the existing finite conditional-repair theorem with the fiber map
given by the occupations outside e and the endpoint total. Each fiber has
one configuration or the two possible singleton placements. Their conditional
weights are precisely (1). The results `heatBath_detailedBalance`,
`heatBath_stationary`, `mixture_stationary` and `heatBath_secondLaw` in
[`FiniteConditionalRepair.lean`](../../Lean/Thermodynamics/FiniteConditionalRepair.lean)
therefore apply once this finite identification is made. The source attachment
and general weighted-graph proof here are prose and executable checks; no
additional Lean theorem is claimed.

## Exact local reference and restriction

For a retained subset A and complement B, tracing B gives the entire number
mixture

\[
 \operatorname{Tr}_B\rho_m
   =\bigoplus_{k=0}^{|A|}h_A(k)\rho_{A,k},\qquad
 h_A(k)=\frac{e_k(w_A)e_{m-k}(w_B)}{e_m(w)}.                        \tag{3}
\]

The convention e_j=0 outside the possible degree range includes every zero
sector. The identity follows by grouping configurations according to S∩A.
For a nested subset C⊂A, the convolution
e_j(w_A)=sum_k e_k(w_C)e_{j-k}(w_{A\C}) proves exact composition of restrictions.
The executable evaluator uses polynomial coefficient recursion rather than
enumerating the global binomial carrier. Equal weights reduce to the exact
hypergeometric restriction in
[`native_capacity_restriction`](../native_capacity_restriction/PROOF.md).

If h_A(q)>0, conditioning on local number q gives exactly rho_{A,q}. The
half trace distance to that normalized conditional state is 1-h_A(q); for an
arbitrary normalized state supported in number q this is only a lower bound.
The conditioning operation and its probability are explicit. Ignoring other
number sectors is not the original restriction or a deterministic preparation.

On the literal federation, choose the repeated local weight word
(25,25,8,8,8,8,8,8,8,8,8,8). Each cell retains its original 30 internal seams.
The conditional q=2 state in every cell is exactly the existing weighted H66
reference. All cross-cell seams use (1), so they exchange loads with adjacent
cells and the tensor product of local fixed-two sectors is not invariant.
The repeated weight word is a declared assignment to existing cell addresses;
the source graph alone does not select it.

There is also an exactly restriction-compatible static family on full binary
tensor spaces. For one declared fugacity z>0,

\[
 \sigma_z=\bigotimes_i
   \frac{|0\rangle\langle0|+zw_i|1\rangle\langle1|}{1+zw_i}
 =\sum_m\frac{z^m e_m(w)}{\prod_i(1+zw_i)}\rho_m.                  \tag{4}
\]

Every seam fixes sigma_z, and tracing vertices simply deletes their product
factors. Consistent vertex weights and the same z are required across a
restriction. The source does not select z, nor does (4) identify cell ancestry
with a physical spatial refinement or construct a dynamical intertwiner.

### Modular current gaps and restriction

The diagonal modular generator on H_m is determined by (2):

\[
 K_m=-\log\rho_m=\log e_m(w)-\sum_i n_i\log w_i.
\]

For a recorded current j from a to b, every nonzero matrix element of that
branch has the same reference ratio and modular gap,

\[
 \frac{\pi_m(S')}{\pi_m(S)}=(w_b/w_a)^j,\qquad
 K_m(S')-K_m(S)=j\log(w_a/w_b).                                  \tag{4a}
\]

The normalization e_m(w), graph size and conserved occupation cancel. For
the repeated native weight word, every possible transfer gap is zero or
plus/minus log(25/8), when the occupation sector permits the transfer. Thus
these dimensionless modular gaps have an exact capacity-independent
continuation within the declared family. They are not identified with
physical energies or particle masses.

For the faithful full-tensor reference (4), the modular generator is
sum_i log(1+zw_i)-sum_i n_i log(zw_i). The fugacity contribution also cancels
in number-preserving repairs. With inherited weights and the same z,
sigma_{z,V}=sigma_{z,A} tensor sigma_{z,B}. Consequently the algebraic
unitary U_V(t)=sigma_{z,V}^{-it} factors exactly as U_A(t) tensor U_B(t), and
for every source-reference operator X,

\[
 \operatorname{Tr}_B[(U_V(t)\otimes I_R)X(U_V(t)^\dagger\otimes I_R)]
 =(U_A(t)\otimes I_R)\operatorname{Tr}_B(X)
                 (U_A(t)^\dagger\otimes I_R).                    \tag{4b}
\]

No factorization of X is assumed. The same statement on embedded local
observables follows by conjugating A tensor I_B. Exact rational probability
ratios on every matrix unit verify the phase identities for arbitrary real t;
no floating exponential is used. Equation (4b) is algebraic modular
covariance of the declared product reference. It is not an execution of a
coherent clock service, a repair-attempt equivalence, or a physical clock
calibration.

## Complete-record comparison at a boundary

Fix a window A and the independent seam selector q_e. Let beta be the total
selection probability of seams crossing from A to B. Use full repair attempts
with repair probability one. Compare the actual process with a control on
the same full carrier, selector and outcome alphabet: on a crossing seam,
the control keeps the source fixed and writes raw zero with the same law and
coin probabilities; the other raw slots remain zero. The other seams are
identical. Neither process suppresses a boundary record.

Attach an auxiliary classical bit indicating whether a crossing seam has
ever been selected. On its zero block the entire products of Kraus operators
and complete raw histories agree. Under either process this block has total
probability (1-beta)^t, independent of the input, even for arbitrary initial
source-exterior-reference correlations. On the complementary block each
process is completely positive and has trace 1-(1-beta)^t. The trace-norm
triangle inequality gives, for every joint input state,

\[
 D(\mathcal I_t(\rho),\mathcal I_t^{\rm cut}(\rho))
 \le 1-(1-\beta)^t.                                             \tag{5}
\]

Here D is half trace distance and each output includes every classical raw
record and the terminal source/reference. The bound also holds in half
diamond norm. Forgetting the auxiliary bit or applying a common readout
cannot enlarge it. The formula is an upper bound, not a claim of equality:
a selected boundary seam can act as a hold on a particular state. The control
is a comparison instrument and is not an admitted substitute source program.

In the cut control, after discarding exterior record values, operations on
the exterior leave the local/reference marginal unchanged. Internal records
therefore implement the ordinary local weighted instrument at its own repair
count. Keeping exterior records can expose correlations; no factorization of
the complete archive is asserted. This local partial-trace step uses the
entry algebra of `liftLeft_conj_apply` in
[`NoSignalling.lean`](../../Lean/ObserverPatchHolography/Locality/NoSignalling.lean).
The theorem `ptraceFst_local_kraus` supplies the complementary statement that
an unread complete local operation preserves the untouched marginal; its
partial trace is over the active factor.

Let a be the probability of an internal seam and b=beta that of a crossing
seam. Discarding independent exterior-only waiting choices gives internal
probability a/(a+b) among relevant choices. Starting before the first choice
and stopping immediately at the kth internal choice, the probability of no
crossing choice through that stopping event is

\[
 \left(\frac{a}{a+b}\right)^k.                                  \tag{6}
\]

Thus the corresponding complete-history coupling error is at most one minus
(6). There is no waiting interval after the kth internal record in this
convention. The stopping time is almost surely finite, although its record
space is countable; the finite-time inequality extends by summing the
nonnegative history probabilities. Matching finite timeout records gives a
finite evaluator and the simpler bound (5). Arbitrary record-adaptive
selectors or other adaptive stopping conditions do not satisfy (5) or (6)
without a separate bound on their conditional boundary hazards.

For the specified federation every cell has 30 internal and 3 crossing seams.
Its cells are triangular faces of the retained subdivision; the dual graph
has one neighbour per side and the declared routing assigns one port pair
to each dual edge. The executable census checks every cell at each tested level.
At level ell the total seam count is 630*4^ell. Hence a=30/|E|,
b=3/|E|, and (6) is (10/11)^k at every level. Increasing the total cell count
reduces boundary probability per global attempt and reduces the internal
repair rate by the same factor. It does not establish vanishing interaction
error at a fixed local repair count.

## Local positive generator

Let T_e be the population transition matrix in column convention, D the
diagonal matrix of pi_m, and S_e=D^(-1/2)T_e D^(1/2). Detailed balance makes
S_e symmetric. On the endpoint singleton basis |10>,|01>,

\[
 h_e=I-S_e=\frac1{w_a+w_b}
 \begin{pmatrix}w_b&-\sqrt{w_aw_b}\\-\sqrt{w_aw_b}&w_a\end{pmatrix}; \tag{7}
\]

h_e vanishes on 00 and 11. It is an orthogonal rank-one projector. The specified
per-seam clock therefore gives the finite positive local operator
L=sum_e q_e h_e=I-D^(-1/2)T D^(1/2). Its zero vector has components
Omega_m(S)=sqrt(pi_m(S)), and connectedness makes the zero eigenspace in each
number sector one-dimensional. On the full tensor space there is one such
vector per conserved number. Assigning rate c_e to each seam instead of one
global attempt multiplies L by sum_e c_e; this is a clock convention. Equal
unit rates per seam correspond to c_e=1 for every seam.

Equation (7) is the symmetric representation of classical population
relaxation generated by the literal endpoint law. Omega_m is its square-root
representation, not a source-prepared coherent state. Neither its spectrum
nor exp(-itL) is identified with a physical Hamiltonian, masses, modular
evolution, or a selected Higgs action.

## Verification and domain

The tests compare the complete H66 instrument against the independent native
interpreter, enumerate every fixed-number sector of small weighted graphs,
check actual partial traces for every retained window, and verify nested
restriction and fugacity factorization. They check the exact generator by
independent conjugation of the population matrix and check the complete
boundary comparison on a source entangled with a full reference. Source
levels 0, 1, 2 supply 240, 960, 3840 ports; their full geometry and seam inventories
are checked without constructing their binomial Hilbert spaces. The receipt
binds all imported source files and these proof/test fixtures.

The result supplies no selected cosmological N, no law for total occupation
as N varies, no universal population mixing bound, no complete coherent
service embedding, no interacting continuum limit and no measured particle
identification. All finite parameters and preparation conditions appear in
the statement rather than being inferred from consistency alone.
