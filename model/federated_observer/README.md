# Elementary-observer federations and their observable limits

The elementary observer is twelve-port. This construction reuses the
retained screen federation: 20, 80 and 320 identical elementary carriers,
with their original thirty internal seams and the existing `port_pair`
connections between neighbouring screen cells. It executes the native
binary endpoint repair across both kinds of seam. The routing convention
and uniform seam schedule are declared inputs; neither is uniquely selected
by this calculation.

This closes a finite composition gap and predicts accessible correlations
on the coupled screen. It does not recover relativistic spacetime. In
particular, coarse screen blocks do not automatically form a closed Markov
process, and the density modes of the reconstructed native transfer are
not automatically relativistic modes.

## What actually executes

The [source adapter](source.py) checks the retained source files and builds
the existing federations. It verifies local instrument identities, supplies
an explicit cross-carrier continuing experiment with complete outcomes,
and checks current conservation. The joint occupation is global: particles
can move between elementary carriers. Keeping exactly two particles in
*each* carrier would discard valid inter-carrier events. The binary local
endpoint operation itself is unchanged.

The continuing experiment uses four specified cross-carrier and internal
seam attempts. It retains all 81 signed-current words, with the 16 hidden
coin trajectories summed inside their CP outcomes. Its six-dimensional
support is an exactly invariant subspace of the full 28,680-dimensional
global total-two sector. The initially empty second carrier ends with
0, 1 or 2 particles with probabilities `5/16`, `5/8` and `1/16`; the mean
is `3/4`. Removing the connecting seams leaves that carrier empty. The
observer obtains the same law from boundary-current sums, and its serialized
checkpoint continues the same experiment.

The [correlation computation](correlations.py) predicts the records of
detectors that sum their incident signed currents. It retains the supplied
uniform fixed-number preparation and the declared attempt windows. These
increment reads require no access to an initial hidden occupation. Absolute
occupation reads require an initial local preparation record.

The [refinement calculation](REFINEMENT.md) checks actual parent/child
screen maps, static restrictions and dynamical defects. The derived
[coarse-memory law](MEMORY.md) reproduces the fine process's coarse density
means exactly, retaining both memory and the initial unresolved profile. The
[scaling calculation](SCALING.md) connects native density modes to the
transfer logarithm and bounds the influence of distant data over a finite
auxiliary evolution interval. These calculations have distinct contracts;
their conjunction is not a proof of a continuum QFT.

```sh
python -m model.federated_observer --check model/federated_observer/results.json
python -m pytest -q -W error model/federated_observer
```

## A closed detector law from the native interactions

Let `V` count ports, `M` count seams and `L` be the graph Laplacian of the
actual federation. For a native occupation configuration `n`, one fair
seam attempt has the exact conditional mean

```text
E[n(next) | n] = P n,                    P = I-L/(2M).
```

This holds for every binary particle number, including interacting
exclusion sectors. It closes the first moment without replacing the full
many-body history or quantum instrument by a one-particle process.

Each endpoint update preserves the uniform distribution on configurations
with a specified global number `k`. Write `rho=k/V`. Its static covariance
is `C(0)=rho(1-rho)V/(V-1) (I-11^T/V)`. Conditional expectation then gives

```text
C_ij(n) = Cov(n_i(0), n_j(n))
        = rho(1-rho)/(V-1) [V(P^n)_ij-1].
```

The readable signed incident-current sum over a window is
`J_i[0,n]=n_i(n)-n_i(0)`. Stationarity and reversibility therefore give

```text
Cov(J_i[0,n], J_j[0,n]) = 2[C_ij(0)-C_ij(n)],
Cov(J_i[0,a], J_j[a,a+b]) = C_ij(a)+C_ij(b)-C_ij(a+b)-C_ij(0).
```

The receipt computes these quantities as exact rational numbers for all
three screen levels, two disclosed preparations and windows of 1, 2, 4 and
8 attempts. Integer sparse propagation replaces an exponentially large
configuration matrix. Independent tests compare the formulas with complete
native histories on small graphs, including the actual current records.
This reduction concerns these correlations only, not all higher moments.

For example, at the two endpoints of an inter-carrier seam, the one-attempt
cross-current covariance divided by either local variance is **`-1/6`**.
Both endpoints have degree six. This ratio holds at every executed screen
level and both supplied particle densities. It is a prediction of the
specified detector protocol, not a comparison with laboratory data.

## Existing results used, with their hypotheses

| Existing result | Attachment here | Boundary |
| --- | --- | --- |
| [Native repair four-cycle](../../Lean/Thermodynamics/NativeRepairFourCycle.lean), `coordinate_drift` and `static_covariance` | The same local fair-swap identity and fixed-number covariance give the graph formulas above. | The general graph argument is supplied here and tested; the existing Lean file formalizes its four-cycle case. |
| [Regional continuity](../../Lean/Screen/RegionalContinuity.lean) | Interior currents cancel; oriented boundary currents determine a regional occupation change. | Complete boundary records are needed. A region's count need not itself be a Markov state. |
| [Native transfer construction](../primitive_observer/DERIVATION_TRANSFER.md) and [finite fields](../native_effective_fields/README.md) | The federation exchange generator has the same seam terms, with its own normalization and positive-transfer bound. | The history transfer is not the native channel on density matrices. |
| [Slab continuum connection](../../Lean/SourceLaw/ContinuumConnection.lean) | Identifies the additional covariance and response hypotheses of the Gaussian slab branch. | Its square-root frequency follows from that branch's action. It cannot be substituted for the logarithm of this native transfer. |

The refinement and scaling notes identify their further theorem mappings.
Exact integer identities, analytic bounds, numerical diagnostics and the
existing formally checked special cases are labelled separately. This packet
adds no Lean theorem and does not rewrite the parent receipts.

## What a spacetime description must reproduce

A literal infinite screen need not be a physical constituent of OPH. A
sequence of federations can be a mathematical approximation tool. On one
finite screen, an effective spacetime claim can instead specify a range of
detectors and interventions with a controlled error. Either route must
match observable joint laws, propagation and clock comparisons on its
declared domain.

The present calculation gives two-point predictions and a quantitative
locality bound for a named sector. It does not identify the screen's
two-dimensional mesh with physical space, establish a bulk dimension, or
derive a Lorentz metric. Native density relaxation can also occur within
relativistic theories; its diffusion-like spectrum is an obstruction only
to identifying these particular modes directly with relativistic vacuum
excitations. It is not a no-go theorem for OPH.

The next constructive test is to find a jointly compatible observable
sector and state whose *own* correlation and intervention laws satisfy the
existing relativistic reconstruction hypotheses. Changing a mode's
frequency by hand, or renaming the attempt counter as modular time, would
not supply that connection. Observer-relative modular time is the
intended type of description; no preferred universal clock is required.

These are model predictions, not empirical matches. The retained
anharmonic comparison is unchanged.
