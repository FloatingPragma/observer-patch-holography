# Larger-volume controls for the native colour loop coupling

The declared source is the same four-dimensional periodic pure-colour Wilson
ensemble as `openwave_scale_coupling`. The five-level seam trace restricts
exactly to three fifths of the fundamental SU(3) Wilson cost. The supplied
cost multiplier ten gives beta6 and bare g²1 on every ensemble. No observed
strong coupling or imported perturbative beta function enters this test.
The change is volume only: sides6,8,10 at the same lattice spacing and
bare action. It is not a continuum limit or a physical energy scan.

The loop readout and its normalization are unchanged. All five rectangle
means enter the nonlinear Creutz ratios. The finite Gaussian Green covariance
gives k_R(N), and `g_C²(R,N)=chi_R/k_R(N)` at lengths0.5,1.5,2.5. In the
Gaussian control this ratio is independent of scale; finite geometric loop
shape alone therefore cannot produce the reported coupling trend. The
normalization is a declared finite observable scheme, not a claim that all
interacting finite-volume zero modes admit that perturbative expansion.

Two additional productions, independent cold/hot starts on side10,
retain400 records each. Each record follows two compound heatbath sweeps
with two overrelaxations; warm-up lengths400,500 are declared and
retained. Final links, independent plaquettes and every rectangle record
are saved. The producer uses the pinned runtime and source kernel from the
own-ensemble package. A separate production command regenerates the raw
histories; analysis replay does not silently replace Monte Carlo production.

## Correlated nonlinear estimates

For a common equal-block jackknife, let w be the five ensemble means and
f(w) the three normalized couplings. Delete each block in all five loops
simultaneously, evaluate f on that deletion, and use their full covariance.
Common deletions preserve correlations between distances. Block lengths
5,10,20,40,80 divide each400-record chain; no block crosses an independent
chain boundary. Equal chain weights are retained. Hot/cold differences and
first/second-half differences are displayed, not selected away.

The infinitesimal influence matrix is the actual derivative of the logarithmic
loop ratio. For example the middle row is
`(-1/w11,2/w21,-1/w22,0,0)/k2` and the final row is
`(0,0,-1/w22,2/w32,-1/w33)/k3`. Multiplying each centered loop record by
this matrix gives the local influence histories. Their autocovariances use
one common denominator n. The initial-positive monotone adjacent-pair sum
supplies an integrated-autocorrelation diagnostic and effective record count.
Independent numerical differentiation verifies the influence matrix, and a
correlated autoregressive control rejects treating every record as independent.
These finite-history estimates neither prove Markov-chain mixing nor supply
rigorous confidence coverage.

Across independent boxes the covariance of a difference is the sum of the
two covariance estimates. The receipt keeps every adjacent-volume difference,
its error diagnostic and all three scales. A small difference relative to a
sampling estimate is not a rigorous bound on finite-volume bias; a visible
hot/cold or time-split difference is retained even when inconvenient.

The experiment tests whether the finite-scale trend survives larger boxes.
It supplies no source-selected core radius, physical GeV conversion,
renormalized continuum beta function or coupling to the interacting particle
Hamiltonian. Matching a reference coupling across independently varied bare
actions would be a different step-scaling experiment.
