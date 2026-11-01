# Charged motion under the uniform interacting action

The existing charged dimer candidates move under the **unchanged Weyl/contact
Hamiltonian** when two dimers are joined. The complete four-site calculation
connects their local residues to hopping, physical number current, spectral
continuation and a source-executed preparation/evolution/read history. It also
quantifies why an isolated band is insufficient to identify a localized
electron or positron.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/charged_dynamics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.charged_dynamics.packet --check model/charged_dynamics/results.json
```

The [contract](contract.json) fixes the carrier, action, references, source
services and remaining joins. The compact [receipt](results.json) is rebuilt
from the proofs and complete executions, with recursive source bindings.
No physical comparison or completion gate is accepted by this package.

## One action and two distinct references

Use the open four-site chain, two CAR modes per site, `c=1/8`, contact `U=1`
and all `2^8=256` occupation states. Let

\[
K(\lambda)=H(\lambda)-N/2=K_A+K_B+\lambda V,\qquad 0\leq\lambda\leq1.
\]

`V` is the actual middle Weyl bond; the two outer bonds retain their original
strength. The endpoint `lambda=1` is exactly the uniform action in
[`full_fock_field`](../full_fock_field/README.md). Intermediate values certify
continuation from the known dimers; they are not fitted physical couplings.
The site-local monomial gauge `d_j=(i sigma_axis)^j c_j` turns all hopping
amplitudes into the real Hubbard value `t=1/16`, without dropping states.
Integer CAR, spin and particle-hole identities are checked on the full carrier
and the gauge is bound to all three coordinate-axis presentations.

For the preparation experiment the input is a lower charged doublet on A
and the isolated dimer ground on B. This is not a stationary state of the
joined chain. Stationary addition/removal spectroscopy instead uses the
unique coupled ground, whose energy is approximately `-2.0364147301` in K.
The source experiment does not claim to prepare that coupled ground. Within
each fixed-N history, H and K differ only by a scalar phase, so their charge
and candidate-location probabilities agree.

## Exact continuation and positive transfer

For every `lambda` in `[0,1]`, the lowest four states in **each complete
N=3 and N=5 sector** form a band separated from their sector complement by
at least `1/256`. The coupled half-filled ground is unique and its excitation
gap is at least `1/128`. Other number sectors remain in the model.

The proof in [`certification.py`](certification.py) uses the fixed-spin
blocks of the full SU(2) representation. Every odd-N multiplet meets
`Sz=1/2`; checking `Sz=3/2` excludes an uncounted lower quartet. Similarly,
the `Sz=0` and `Sz=1` blocks distinguish the half-filled singlet ground from
every higher-spin multiplet. Symmetric integer congruence counts eigenvalues
on a 64-cell grid, at dyadic thresholds with denominator 65536. Floating
eigenvalues propose brackets; exact inertia must validate every bracket.
The middle bond has norm at most `1/8`, because its two commuting spin hops
each have norm at most `1/16`. Weyl's inequality therefore costs at most
`1/512` in a gap between neighboring grid points. The surviving uniform
margins are the bounds above. No sampled gap is promoted to an interval proof.

Finally, positive contact energy and the one-particle hopping row bound give
`K >= -5N/8 >= -15/8` for `N<=3`; particle-hole symmetry covers `N>=5`.
The half-filled product-ground trial has energy at most `-2`. Thus the
half-filled ground is also the global ground, separated from other number
sectors by at least `1/8`.

Let `W_A,W_B` be the phase-fixed bare charged doublets, with the other dimer
in its ground. The [dimer residue](../particle_response/README.md) fixes

\[
W_B^*VW_A=-jI_2,\qquad
j=t\frac{1+1/\sqrt{17}}4=\frac{1+1/\sqrt{17}}{64}.
\]

This is an action-derived hopping coefficient, with no independently chosen
translation operator. `K(0)` preserves A's number parity and V reverses it.
Consequently only odd Dyson terms connect A to B. Both bare doublets have
the same `K(0)` energy, so the first term has magnitude `jT`, and

\[
\sigma_{\min}(W_B^*e^{-iK(1)T}W_A)
\ge jT-[\sinh(T/8)-T/8].
\]

At `T=4` an exact rational radical enclosure and Taylor-tail bound give a
probability lower bound greater than `0.00319`, for **every internal spin
state and untouched reference**. This is a full-Hamiltonian transition
bound, without a band projection, conditional normalization or fitted gap.

## Actual charge and the localization cost

The physical observables are `Q_A=N_A-2`, `Q_B=N_B-2` and outward current
`J_A=-i[K,Q_A]`. They obey `dQ_A/dT=-J_A`, `dQ_B/dT=J_A` on the full carrier.
The particle-hole unitary preserves K and reverses both charges and current.

For the isolated band projector P, use the canonical polar chart
`W=P W_0 (W_0* P W_0)^(-1/2)`. This is a coordinate choice inside an invariant
subspace, not a new physical definition of position. At the uniform endpoint,
with `q=-1` for holes and `q=+1` for particles, the numerical matrices give

\[
W^*KW=\bar E I-J\tau_x,\quad
W^*Q_AW=\frac q2(I+\eta\tau_z),\quad
W^*J_AW=qJ\eta\tau_y,
\]

where `J=0.0064624483` and `eta=0.0791310376`; each Pauli matrix acts on
location and is tensored with the two-dimensional internal identity.
The current follows from the same K and charge, rather than a guessed band
current. Charge readout also has a nonzero off-band part, which is retained.
Polar effective-Hamiltonian coordinates are a standard finite-subspace
construction; see [Bravyi, DiVincenzo and Loss](https://arxiv.org/abs/1105.0675).

| Experiment | Result in model units | Meaning |
| --- | --- | --- |
| Bare charged input, uniform action, T=4 | Opposite bare-candidate probability about 0.00585153 | Genuine motion of the original local excitation, with all other outcomes retained |
| Same input, T=pi/(2J), about 243.0652 | Opposite bare-candidate probability about 0.773790 | Ideal full-Hamiltonian calculation; the source certificate below covers T=4 |
| Dressed band location swap at pi/(2J) | Logical swap probability 1; physical charge contrast only 0.0791310 | A logical location swap does not move a sharply localized unit charge |
| Bare input's isolated-band weight | About 0.510988 | Almost half of the input is outside the retained spectral quartet |

The eigenvalues of `W_0* P W_0` are approximately 0.270479 and 0.751497,
each twice. The dressing is therefore substantial. Complete occupation
histories, block-charge distributions, current, energy, exterior reduced
states and band-complement weights are available from
[`observables.py`](observables.py). The coupled-ground response retains all
256 final states for each of eight local addition and removal probes.
Every per-mode addition/removal weight totals 1/2; band residues and their
complements are reported separately, with CAR sums and moments through
order two checked against the original operators.

## Complete source history and budget

[`source.py`](source.py) prepares an unknown internal qubit, optionally
entangled with an untouched reference (up to 16 explicit reference columns
for complete pointer execution; the analytic bound is reference-uniform).
Counted X/CNOT apparatus operations
supply the seed occupation. A number-preserving complete dimer basis encoder
prepares the actual candidate and ground columns. Symmetric split evolution
then applies the uniform action for `T=4` with 256 steps.

There are two separate terminal refinements: all 256 occupation outcomes,
or the inverse dimer basis followed by all 256 outcomes. The latter resolves
candidate A, candidate B and every other outcome; neither read discards the
remaining addresses. Both also retain the 25 block-number bins. Their full
raw records are different experiments, not interchangeable histories.

Each local source word is executed on its entire input algebra, including
every six-level export complement. An immutable factorized archive owns the
primitive matrices and input, and reconstructs every first-failure quantum
output by its exact prefix and literal target remapping. Failure mass is
summed without clipping or renormalizing success. Archive consumers bind the
executed schedule, input, outputs and resource counts; dropping a step or
replacing a retained field cannot produce a smaller accepted error budget.
The eight pointer interactions are counted. Their eight fresh blank bits and
the seven field
blanks accompanying the unknown input qubit are supplied initialization
resources; their physical preparation cost is uncalibrated.

The analytic split error in operator norm, hence half diamond distance,
is at most `841/1572864`, about `0.0005346934`. Local floating source replay
uses fixed complete-instrument admission gates: `1e-10` per small macro and
`1e-7` per dimer encoder/decoder. If these admitted errors are epsilon_j,
the complete composition error is bounded by
`product(1+epsilon_j)-1 <= s/(1-s)`, where `s=sum epsilon_j<1`.
This includes possible norm growth; it is not an outward interval proof of
floating arithmetic or a hardware-error guarantee. The receipt gives the
exact rational propagated budget and the resulting positive transfer
margin after subtracting errors from the **amplitude** lower bound.

Encoder control synthesis uses a fixed 12-decimal convention and `1e-12`
empty-pivot threshold, retaining every pulse and service request. This avoids
arbitrary rotations generated from roundoff at empty elimination pivots.
The resulting complete source instrument is still checked against the
original analytic encoder; quantum states and outcome weights are never
rounded into agreement. Fresh receipt replay reloads the bound source-control
data and recomputes cached primitive and encoder matrices.

The occupation histories use 283,877 hole or 283,881 particle native pulses;
the candidate-read histories use 307,625 or 307,629 respectively. All load,
export, scalar-phase and decoder operations are included. Fresh source
execution gives opposite-candidate probability about `0.00585147`; the
analytic split and admitted numerical source budgets leave a probability
lower bound above `0.0031`. Neither a native repair count nor this supplied
Hamiltonian duration has been identified with SI time.

## Particle and topology status

This completes a finite same-action connection from the charged local
candidate to motion, current, coupled-ground response and a reproducible
source history. The auxiliary interpolation proves persistence of spectral
subspaces; it does **not** classify a nontrivial topological defect. CAR
number is conserved by this Hamiltonian, but is not thereby conserved under
every native C/E repair branch. The [native conservation obstruction](
../completion/defect_algebra.md) and the separate actively recovered memory
constructions keep their own domains.

The existing [periodic dispersion benchmark](../particle_comparison/README.md)
still fails under its fixed particle dictionary. The [isolated-dimer
family exclusion](../dimer_family/README.md) also remains in force: the
four-site construction does not supply a three-generation lepton spectrum
or repair that dimer dictionary. This open-chain result
neither changes that benchmark nor establishes electromagnetic charge,
Lorentz spin, a chiral continuum, native passive protection or a physical
electron/positron assignment. Those are explicit remaining joins in the
[completion program](../../docs/COMPLETION_PROGRAM.md).
