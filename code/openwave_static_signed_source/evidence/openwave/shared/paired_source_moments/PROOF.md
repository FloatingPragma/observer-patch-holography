# A primitive paired source with apolar and polar current records

Let v_p be the twelve vectors in the finite current quotient, with
phi^2=phi+1 and |v_p|^2=phi+2. The explicit table is reindexed to the sealed
relative-source graph. Its adjacency is v_a dot v_b=phi. The Lean theorem
`SeamCurrentEdge30Moment.carrierSeamDifference_eq_endpointDifference` proves
that the finite quotient of a unit current is v_b-v_a; no physical translation
or clock is inferred by that theorem.

## Measured source and complete instrument

Use the primitive source from the relative-source construction:
x=delta_p+delta_s, y=delta_r+delta_s, with r antipodal to p. A current-only
restoring census independently reads both inputs. For each tree edge its
first current and first opposite-coin current give the initial edge contrast
2j+j*, after correction by the accumulated incidence. Integrating the tree
and using the binary total-two promise fixes the configuration. Native routing
restores the decoded input, with explicit timeout outcomes. This is the
existing source recovery construction, not a discovery claimed here.

For a body record let k=jX-jY and delta=k(v_b-v_a). Continuity gives
q_next=q-k B_e for q=x-y. The body enumeration computes all sixty elementary
branches from each restored input independently. It is an exact instrument
calculation conditional on a read preparation. It neither copies an unknown
state nor represents sixty branches as one continuing trajectory.

The two endpoint cases establish independence of the common background.
An edge containing neither discrepancy gives k=0. An edge incident to exactly
one discrepancy gives one moving coin and one holding coin. The antipodes are
not adjacent. Each of five edges at p gives displacement w-v_p for one of its
coins. At r, write its neighbour as -w; the negative discrepancy gives
-((-w)-(-v_p))=w-v_p. Thus each of the five vectors w-v_p occurs with probability
2/60=1/30, and the zero vector has probability50/60. The moving coin may depend
on background. Distributional equality is claimed; same-coin pathwise equality
across backgrounds is not claimed.

## Exact mean and covariance

Let m1,C1 denote the one-unit moments derived by finite endpoint enumeration in
the bound current-covariance proof. They obey

    m1=((phi-3)/30)v,
    C1_parallel=11(3-phi)/180,
    C1_transverse=(phi+2)/30.

The paired law doubles each nonzero probability, so its raw second moment is
twice the one-unit second moment, while m=2m1. Therefore

    C=2C1-2m1 m1^T,
    |m1|^2=(3-phi)/180,
    m=((phi-3)/15)v,
    lambda_parallel=(3-phi)/9,
    lambda_transverse=(phi+2)/15.

The covariance gap is (9-8phi)/45, which is nonzero. Exact arithmetic in Q(phi)
checks these formulas independently of the formula values used by the gates.
The measured mean is nonzero and points along the same undirected axis that
is selected by the measured covariance. The covariance is computed from its
raw second moment; it is never defined as the square of the mean.

## Proper rotations and reader identities

R=diag(-1,-1,1) is a proper half-turn and a permutation of the source graph.
It sends v_0=(-1,phi,0) to v_3=-v_0, and R^2=I. Background-independent native
statistics and the exact uniaxial form yield

    m_3=R m_0=-m_0,
    C_3=R C_0 R^T=C_0.

The nonzero polar mean returns under R^2. Under a continuous rotation of the
fitted vector/tensor representation about this same axis, the tensor period
is pi and the vector period is2pi. Intermediate angles are a mathematical
interpolation of the fitted representation. Native inputs here implement only
the finite half-turn and its square.

Reversing the ordered source pair negates k for every identical body law.
Its odd mean and even covariance alone are kinematic. Native law dependence
enters the independently computed nonzero uniaxial shape, its attachment to
the independently read source axis, and its covariance under the fixed R.
The deterministic fixed-seam (0,5) exchange has a single nonzero jump direction
v_5-v_0=(1,1-phi,phi). It has both transverse-plane and z components, so its
nonzero rank-one covariance does not equal R C R^T. Pair reversal continues to
leave that control covariance unchanged, illustrating the distinction. Frozen
body erases both observables. The planted off-axis two-vector distribution
passes through the identical moment accumulator and rotation detector.

Pooling one complete law from every source orientation gives zero mean and
covariance2I/9. The polar axis is a preparation-conditioned observable. This
is a static witness on the same primitive source class as the row1 response;
it proves neither free persistence nor a physical spin angular momentum.
No central phase or exchange theorem is required by the literal row7 test.

## Finite stack and scope

The consumed stack includes the exact native binary endpoint cases in
`Lean/EventAlgebra/NativeBinarySourceEndpoints.lean`, the regional incidence
identity in `Lean/Screen/RegionalContinuity.lean`, the signed-load quotient in
`Lean/Screen/SeamCurrentCarrierQuotient.lean`, and the port/seam identities in
`Lean/Screen/SeamCurrentEdge30Moment.lean`. The source recovery theorem and
`model/native_clock/configuration.py` supply the active census context.
The sealed row1 source proof supplies the primitive relative class and its
native static response. This experiment supplies the moment join on that
source; it does not rerun every row1 stopping response or topology certificate.
No new Lean build or formalization of this joined theorem is claimed.

The declared apparatus contains the graph and its finite geometric quotient,
paired sources, binary preparation domain, common scheduler, current archive,
census/reset controller and moment reader. It has no calibrated nature value.
The finite proper port action needs an accepted physical rotation/source
dictionary for a full row interpretation. Physical angular momentum, magnetic
coupling, rest mass and clock are not consequences of these moment identities.
