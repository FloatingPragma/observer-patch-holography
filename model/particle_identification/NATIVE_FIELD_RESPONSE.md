# One native source for preparation, spatial response and exchange time

The [preparation](../native_correlation_preparation/PROOF.md),
[spatial pair](../native_pair_response/PROOF.md) and
[exchange clock](../native_exchange_clock/PROOF.md) theorems compose on one
specified source. Their conclusion is a controlled finite correlation
experiment. The physical particle-identification requirements in the
[research route](ROADMAP.md) are separate from this conclusion.

## Shared hypotheses and complete experiment

Use a connected binary federation graph with V ports, M seams and the
unchanged fair native endpoint law. Every unread attempt selects a seam
independently and uniformly. Repair has probability one and the endpoint
coin is fair. Fix a conserved occupation 2<=q<=V-2. The state to which the
observed expectation is compared is tau_q=I/binomial(V,q). Neither q nor
this reference is inferred from an external physical vacuum. The default
biased twelve-port law, adaptive seam selection and a weighted connected
extension have different contracts.

Choose two disjoint native seams at each of the past and future settings.
Their contrast fields are phi_i=(n_a-n_b)(n_c-n_d), with primitive
two-subset vectors h_i of squared norm four. Distinct settings may overlap.
Starting from any source/reference state rho_QR supported on the fixed q
sector, perform n warmup attempts, two past E reads, T unread attempts,
and two future E reads. Retain every native record and FIFO export.

For the comparison observation, ignore warmup and body raw values without
conditioning. Multiply the two past outgoing contrasts and the two future
incoming contrasts to give a score S in [-1,1]. Let Z_R be the signed sum
of terminal reference blocks weighted by this score. The preparation
theorem includes source coherences and arbitrary source/reference
correlations. It does not assert independence of the complete warmup tape
and the future experiment.

The same source has unread transfer and exchange operator

    P_q=I-G_q/(2M),       G_q=sum_e(I-S_e),
    v_q=4 binomial(V-4,q-2)/binomial(V,q)>0.

Here S_e swaps endpoint occupations. Under tau_q the exact complete-record
expectation is

    C_ij(T)/v_q = h_i^*P_2^T h_j/4.

The two-subset expression is an exact coefficient-space representation of
the experiment on q occupied ports. It does not change the prepared sector.
The two past probes' backaction is included in the effect
H_(past1)H_(past2)[phi_i P_q^T phi_j], of degree at most four.

## A combined preparation and action bound

Let gamma_0>0 be a certified lower bound on the unread exclusion gap,
N=binomial(V,min(4,q,V-q))-1, and

    epsilon_prep=sqrt(N/v_q)*(1-gamma_0)^n.

The preparation theorem proves

    ||Z_R/v_q - [C_ij(T)/v_q] rho_R||_1 <= epsilon_prep.

This is a bound on one signed reference-valued observation. It is not a
trace-distance bound on the full source state or on the complete output
instrument. The rational graph certificate gamma_0=2/(M V D), for a proved
diameter bound D, and integer half-contraction blocks give sufficient n
without uncertified eigenvalues or logarithms.

For T>=1 define the counted exchange time t=T/(2M) and

    B_ij^ex(t)=h_i^*exp(-tG_2)h_j/4,
    epsilon_clock=min(1,3/(5T)).

The native positivity bound P_q>=(10/21)I and spectral calculus establish
the clock error in every occupation sector. Since ||h_i/2||=1, triangle
inequality gives the composed result

    ||Z_R/v_q - B_ij^ex(t) rho_R||_1
      <= epsilon_prep + epsilon_clock.                    (1)

At T=0 the clock error is zero. For a requested t with
T=floor(2Mt)>=1 and kappa=2Mt-T, replace the clock allowance by
min(1,(3/5+kappa/2)/T). These constants compare native counted attempts
with the exchange operator of the same graph. They supply no time in
seconds, physical Planck constant or relativistic rest energy.

For a finite linear combination of entries with weights a_i,b_j, the
preparation allowance in (1) is multiplied by ||a||_1||b||_1. The sharper
exchange allowance is epsilon_clock times
||sum_i a_i h_i/2||*||sum_j b_j h_j/2||. The probe Gram matrix is retained;
the settings need not be orthogonal or jointly measurable.

If an expanded execution stops on prefixes with total probability p_stop,
assign score zero to those returned stop outcomes. Complete positive branch
extension bounds the extra normalized signed-score error by p_stop/v_q.
Conditioning on successful returns would define a different experiment.
For R independent completed repetitions, Chebyshev's inequality gives
probability at most 1/(R v_q^2 epsilon_stat^2) of sampling error exceeding
epsilon_stat in the normalized scalar response. This sampling allowance
does not bound a full reference-valued output estimate.

## Spatial observation and continuum scope

The spatial pair theorem supplies exact entries for several probe supports.
For the observation consisting of those probes and a final local state,
let S include every observed seam endpoint and every retained final-state
vertex. A radius-r ball B around S has escape allowance

    p_escape <= min(1, |S| binomial(T,r+1)(Delta/M)^(r+1)).

The matching local control uses the actual restriction of tau_q, including
all local occupation sectors. Replacing the stationary response by that
control adds at most 2 p_escape/v_q. This estimate traces remote raw values;
an observer who also reads those values or uses them in feedback requires
their dependencies in S. It cannot be applied to a global final detector.

A different complete observation accesses arbitrary quadratic occupation
snapshots. One pure-E uniform body followed by a final spanning-tree E
readout determines its terminal configuration from the public total q.
Reversing all retained body currents determines each earlier configuration.
It uses one final tree, with its disturbance after the body, and does not
assume that local seam contrasts span every quadratic field. For quadratic
snapshots F_i,F_j its initial effect is F_i P_q^T F_j, of degree at most
four. The generic polynomial preparation theorem applies using an explicit
centered L2 bound and normalization for this effect. Its error cannot be
replaced by the unit-contrast bound without those norm checks.

For the bounded smooth profiles, occupation q_l=V_l/3 and spatial
subsequence of the owning density limit, exact removal of constant and
one-body components gives a primitive quadratic field. Its one-body
correction has a controlled L2 error of order V_l^-1/2. Uniform canonical
moments justify convergence of quadratic covariances to the Wick products
of the inherited Gaussian diffusive field. The quadratic variables
themselves need not have Gaussian distributions. Their response generator
is the symmetric restriction of A tensor I+I tensor A. This limit has no
additional interaction supplied by the finite exclusion term.

Finite exclusion responses can differ from the collision-retaining
independent-pair control. The spatial proof bounds that difference without
normalizing away collisions. That finite difference alone does not prove
an interacting relativistic continuum sector, and the smooth limit does
not identify unsmeared microscopic seam contrasts with continuum particles.

## Resource and physical-identification boundaries

The four-E experiment costs n+T+4 attempts. Its finite description, attempt
ceiling 2^32, FIFO rule and expanded evaluator limits are unchanged.
Sufficient preparation budgets beyond single-program admission are
mathematical bounds, not returned source runs. The body-plus-tree experiment
has its own description and attempt budgets. Sparse marginal propagation,
authenticated raw-word queries and expanded history execution remain
distinct computations.

Equation (1) gives a common finite source, reference target, observable and
counted clock. It identifies neither the tracial working reference with a
transfer ground vector nor either with an interacting Higgs vacuum. A
physical capacity selection, admissible preparation family, growing quantum
measurement services, physical spatial and time dictionary, interacting
chiral field limit, observable particle sectors and independent calibrated
predictions require their own proofs. No measured mass ratio follows from
these response bounds.
