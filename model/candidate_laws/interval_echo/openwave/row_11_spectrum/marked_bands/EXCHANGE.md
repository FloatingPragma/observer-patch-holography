# Which exchange phase do the three transported orbitals display?

The matched two-rail device distinguishes the **positive exchange** of two native vacancies from a negative exchange phase with identical single-vacancy motion and identical contact interactions. The distinction requires a path that lets the two vacancies pass around each other. Opening two existing physical contacts closes the circuit through the second rail. An independent one-vacancy experiment measures the phase around that same circuit. [Complete port records](exchange.json) · [Native circuit and readout](exchange.py) · [Independent tensor-product and time-evolution checks](test_exchange.py)

## What do two arriving vacancies do at a beam splitter?

On the original two disjoint rails, prepare the same one of the three [well orbitals](README.md) in the left and right wells of rail A. Both charges occupy the same rail coordinate that previously carried a qubit. A preparation bank containing one vacancy selects the exposure closest to half transfer from the 257 published pilot times. A fresh bank prepares two vacancies and reads both port addresses by calibrated ternary interrogation. For levels in increasing order of launch energy, a single vacancy transfers with probabilities 0.49687, 0.50035 and 0.50204. Two vacancies then emerge one per original well in **8,383,229, 8,380,131 and 8,383,500** of 8,388,608 starts. Preparing one vacancy on each disjoint rail gives one per longitudinal well in **4,194,027, 4,189,162 and 4,193,151** starts. Each bank includes every address pair. The stored controls also change the right-hand orbital and compare the independent symmetric and antisymmetric predictions.

Two native vacancies cannot occupy the same ternary port. On adjacent ports their exact energy relative to two independent one-vacancy generators is −17 units. This follows by comparing the blank-seam energies for 00, 01, 10 and 11. Keeping this term gives identical generators in the ordered hard-core boson and fermion position bases on each *open* rail. The apparent exclusion at a balanced splitter therefore admits both signs for every possible position read. It supplies an experimental control for the loop test.

## The loop and its independent flux calibration

Both neighboring pairs of buffer ports contain unit loads in every prepared local orbital. Reversing their zero-contrast writes reconnects the two rails at each end. The closed route has 26 ports; the left and right sources and all twelve well-side marks retain their original preparations. The three pairs of input orbitals all start on rail A. The two compensating excesses are held separately in two remote isolated dimers. Their physical cuts can be written on filled contacts by parking each excess at its other dimer endpoint. The [stock-pair preparation](material.py) supplies the marked well sides; the two vacancies under test are additional source loads.

At the same exposure **134,217,728** for all three preparations, the terminal observers record both vacancies anywhere on the original rail. The table counts this event against all **8,388,608** starts per orbital. The last column gives the probability in the corresponding interacting fermion model. Both models retain the −17 contact energy, all 26 physical hopping amplitudes in their one-vacancy sector, and the measured well and reference potentials.

| Local orbital | Both vacancies on A | Observed probability, simultaneous 99% interval | Matched negative-exchange probability |
|---|---:|---:|---:|
| Lower | 8,386,748 | 0.999762–0.999794 | 0.354368 |
| Middle | 8,385,851 | 0.999651–0.999691 | 0.007166 |
| Upper | 8,387,785 | 0.999891–0.999912 | 0.470818 |

The exposure was selected by inspecting the predicted separation on a power-of-two ladder. The six subsequent one- and two-vacancy outcome banks use independent seeded draws at that common exposure. Exact binomial intervals cover all six rail-count probabilities simultaneously. The receipt also retains the 26 single-address and 325 two-address outcomes. An independent calculation squares the short-time full propagator 34 times and reproduces the middle two-vacancy separation.

In the ordered coordinate basis, crossing the cut between the last and first sites takes one vacancy past the other’s index. A fermionic exterior product contributes a negative sign on just these transitions. A symmetric hard-core product contributes a positive sign, as does the measured native seam. Independent tensor-product construction reproduces both complete matrices. The choice of cut only changes the coordinate phase convention. A one-particle π flux would reverse this two-particle comparison, so observers also launch a single vacancy from the left well. The three measured counts on A are **4,995,091, 695,014 and 5,754,201**. The predictions with an added π flux exceed 99.99% on A in each case.

A second instrument measures the phase directly. A reference path holds a localized vacancy while the other path makes 26 adjacent native π contacts, going once around the whole loop. The one-particle flux comparison returns 8,388,608 real-quadrature bright outcomes in 8,388,608 starts. Next place a vacancy in each well. A sequence of 26 contacts moves them around each other, winds once, returns the initial occupation word and exchanges their tracked paths without ever hopping through an occupied port. This comparison gives the same bright count. Dividing by the independently calibrated one-particle phase gives exchange phase **+1**. A quarter-turn of the phase reference produces close to half bright in both runs. For each of the three common source orbitals, physically compiled inverse well analyzers localize its two copies before the exchange, and their inverses restore them. Each phase comparison again gives all bright. Inserting π in the reference arm sends the same inputs to the dark output.

The coherent reference, addressed isolated contacts and supplied exposure clock belong to the phase apparatus. In the separate autonomous correlation test, every active physical seam operates throughout release. Taken together the two instruments identify positive exchange for these orbital vacancies on this calibrated material circuit. Other OPH sources with additional mobile field degrees of freedom or statistical phases require their own comparison. The three measured orbitals therefore fail the same-label fermionic exchange test for an electron, muon or tau.

From the repository root, run the experiment with an unused output path for a new receipt.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.marked_bands.exchange
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/marked_bands/test_exchange.py
```
