# What the first pair reading changes

One observer on the 20-cell sphere can find an adjacent pair after the two loads start separated and evolve freely. Its probability of finding a pair anywhere in its own cell at duration five is 0.499425% under Interval Echo’s free κ = 8 interaction. Inserting a sharp reading of one nominated seam at duration one raises the late probability to 0.512724%, counting both outcomes of the inserted reading. The choice to read changes the observed population. [Frozen protocol](PROTOCOL.md) · [Complete patch probabilities](result.json)

| First service at duration one | Late pair inside the observer cell, κ = 8 | Exchange only, κ = 0 |
| --- | ---: | ---: |
| Bypass the reader | 0.00499425 | 0.23581465 |
| Choose read or bypass with a fair bit | 0.00506074 | 0.23591358 |
| Sharp read, keeping every trial | 0.00512724 | 0.23601251 |
| Condition on a sharp positive | 0.99322582 | 0.33867101 |

The first three rows use the same denominator of all enrolled starts. The fourth row conditions on the loading of one particular seam at the first reading. That seam gave a positive with probability 0.00018914 under κ = 8 and 0.01106690 under exchange. The freely evolving exchange-only control gives the larger terminal adjacency yield at this wait. After a sharp positive, the contrast interaction gives the larger probability of a repeated pair observation. The unconditional rows and the postselected row answer different experimental questions. [Prior two-read observation](../../evaluation/observer_encounter.json) · [Controls and quiet outcomes](result.json)

A separate experiment retains the native coherent repair on the same seam and tests the terminal yield away from its endpoints. [Retained-write experiment](MEMORY.md)

## What the patch can read

A patch of one echosahedron occupies cell 0. The 36 preparations put one load on each of two distinct, nonadjacent ports in that cell with equal weight. All 28,680 basis states of the 20-cell two-load sector enter the free propagation. The first binary reading tests seam (0, 1). The final local read resolves the 30 internal seams and a complementary quiet outcome. Quiet includes two separated loads and a pair beyond this observer’s aperture. Preparations, measurement windows, service decisions and waiting durations are supplied. The clock is idle. Each quoted probability is calculated from the complete finite record distribution. [Native sector](../../candidate_laws/contrast_memory_spherical_v1/sectors.py) · [Local records and compiler](../../candidate_laws/contrast_memory_spherical_v1/presence.py)

Let E denote the projector onto occupation of the nominated seam and Q = I − E. A sharp reading with both outcomes retained sends the input to EρE + QρQ. Deleting the record also loses the relative coherence between E and Q. Applying the local relative sign Q − E to half an ensemble and identity to the other half gives the same measured channel. With a sign on a quarter of the ensemble the channel matches the fair read-or-bypass service. Every final outcome in both interactions agrees to within 10⁻¹⁰. The sign comparison checks the channel identity mathematically. The random read-or-bypass service requires an ordinary sharp read and a fair service choice. [Readout comparison](__main__.py)

## What the free binding theorem explains

The dressed bound band is the isolated low-energy subspace of the contrast Hamiltonian. Its projector commutes with the free Hamiltonian. Its mean population therefore stays constant between interventions, even as the probability of adjacency on a particular seam changes. The simulator can evaluate the band weight. The patch detector directly measures local adjacency. [Binding theorem](../../candidate_laws/contrast_memory_spherical_v1/PROOF.md) · [Band diagnostic](result.json)

| κ = 8 ensemble at duration one | Dressed bound-band weight |
| --- | ---: |
| Separated input, then bypass | 0.00291210 |
| Sharp positive | 0.99556927 |
| Sharp quiet | 0.00290052 |
| Both sharp outcomes, weighted by their probabilities | 0.00308827 |

The all-outcome increase is 0.00017618 per enrolled start. A sharp positive prepares a state with 99.5569% weight in the bound band. The initial separated preparation also has a small bound component, which propagates freely. An unconditional disturbance and conditioning on a rare reading both enter the two-read result. An encounter-created durable record that predates the nominated read is an additional condition in OBS-04. [Observer criterion](../../evaluation/EVIDENCE.md#obs-04)

For the shifted κ = 8 Hamiltonian the binding theorem places the lower band at energies at most −10 and the upper band at energies at least 0. Positivity and the maximum hopping row sum put the full spectrum in [−16, 12]. With X = (H + 5I)/17, the spectrum obeys 5/17 ≤ |X| ≤ 1 and the projector is (I − sign X)/2. For a = 5/17, q = (1 − a)/(1 + a) and T = (1 + a² − 2X²)/(1 − a²), the Legendre generating function gives

$$
(X^2)^{-1/2}=\frac{2}{1+a}\sum_{k=0}^{\infty} q^k P_k(T).
$$

The spectrum of T is in [−1, 1], so truncation at k = 44 bounds the absolute operator error of the projector by q⁴⁵/(2a), below 2.43 × 10⁻¹² in exact arithmetic. A dense one-cell control independently diagonalizes the Hamiltonian and directly evolves the split density matrices. This checks both the numerical projection and the branch weights. [Dense check](test_encounter.py)

## Reproduce

From the RER repository root with the pinned model dependencies, run

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m model.experiments.interval_echo_encounter
python -m pytest -q model/experiments/interval_echo_encounter
```

The replay reconstructs the state amplitudes and compares all final seam probabilities, quiet outcomes and band diagnostics with the archived values. The protocol digest is checked before the run. [Implementation](__main__.py) · [Receipt](result.json)
