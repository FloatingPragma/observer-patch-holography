# Native clocks and the phase resource

The source is the complete native four-occupation reduction of the charged
branch, its Gauss reference and neutral fermionic memory. The native addition
operator A obeys [H,A]=mA on the full retained carrier. The collision uses the
actual one-copy Kraus operators K0, K+ and K- of the autonomous-clock package,
with p=delta gamma=1/1048576. The inactive occupations are retained, and K-=0
is an explicit recorded output. The source mass and all96 native probe rows
are inherited separately for electron, muon and tau. Their GeV units retain
the electron-anchor and analyzer ancestry.

## Conditional phase obstruction

Let a finite source S and battery B have Hamiltonian Hs+Hb. Assume their
initial joint state commutes with Hs+Hb, U commutes with Hs+Hb, and every
recorded battery effect Eb commutes with Hb. The unnormalized conditional
state is sigma_b=Tr_B[(I tensor Eb) U rho U*]. Conjugating this expression by
exp(-itHs), inserting exp(-itHb), and using invariance of the partial trace
gives exp(-itHs) sigma_b exp(itHs)=sigma_b. This proof permits correlations,
degenerate energy blocks, zero-probability branches and arbitrary finite
ancillas. It iterates for adaptive branches satisfying the same premises.
Every fixed source observable therefore has a constant free expectation in
each such branch. This is stronger than cancellation of an ensemble mean.

For the actual collision [H,K_b]=omega_b K_b, with omega_b/m equal to0,1,-2.
The three battery output energies are m,0,3m, from an input of energy m.
Thus (Hs+Hb)V=V(Hs+m), and the Gram sum is exactly the identity. Exact
arithmetic executes all3^4 battery words for each of the four native input
occupations. Every allowed branch is energy invariant; every forbidden
word has explicit zero probability. Source energy gain and battery loss
cancel on every nonzero word. These finite records instantiate the general
theorem without a collective-spin or success-only projection.

Measuring the battery in the superpositions of its m and0 states changes
the source instrument to (K0 plus or minus K+)/sqrt(2). Each of its two
outcomes has probability1/2 on the source vacuum, and its unnormalized
native X expectation is plus or minus sqrt[p(1-p)]. The third battery
outcome is retained with probability zero. This readout fails to commute
with battery energy. Its mean work is m(p-1/2), including the battery;
source energy is mp before and after the readout. The resulting conditional
phase uses an energy-incompatible measurement resource. No zero-resource
onset follows from these records.

## Relational control on the same native pairs

Two full pairs have dimension16. Their energy is H1+H2 and the bounded
exchange is lambda(A1 A2*+A1* A2). The exchange commutes with total energy.
Starting in the first occupied pair, it creates a relative coherence between
|11,00> and |00,11>, retaining all other occupations. The complete terminal
occupation records are evaluated at interaction phases0,pi/4,pi/2 and pi.
Each state commutes with H1+H2, so this relative coherence has free total-energy
frequency zero. With the exchange left on, populations oscillate at2lambda.
The interaction coefficient is a declared architectural parameter; the
source supplies no identity that sets it to the particle rest gap.

The exchange starts with zero expectation and conserves that expectation;
its on/off mean work is zero. Total native energy is m throughout. Each
single-pair energy coherence is zero even when the relative coherence is
nonzero. These complete controls do not establish an autonomous single
particle rest-mass clock. They exclude the finite energy-invariant collision
route under the stated premises, while permitting other declared sources,
phaseful initial states and separately proved infinite-system limits.

The algebraic energy-shift branch and finite coarse-graining identities are
proved in `ClockPhaseCovariance.lean`, with only the standard propositional
and quotient axioms. The partial-trace argument and physical native
attachment have the explicit hypotheses stated above.
