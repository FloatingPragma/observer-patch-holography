# Stretching the contact and measuring its work

The reciprocal first-arrival length of the contact between ports 1 and 38 can be stretched by a measured amount. A pulse on that contact changes the speed of its isolated exchange clock. A separate Ramsey measurement determines whether the dimensionless energy scale of the contact material has stayed fixed. Together these two readings select a longitudinal work operator on the enrolled contact. A further experiment measures how an impulse on the contact changes the center-of-energy velocity of each adjacent cell.

The work and impulse readings use the same two stationary heads, four-tick free flight and equiprobable forward/reverse tagged launches as the energy-current experiment. The actuator acts *after* that appointment. Ordinary and pairing variations retain the independently prescribed radar volume. Clock and material calibration use separate preparations. [Source, tags and current](../README.md) · [Complete simulated work and calibration records](strain.json) · [Impulse readings](results.json)

## An observable strain and its conjugate work

Let $D=P_1+P_{38}$ select contact occupation and let $P_-$ project onto the antisymmetric combination of the two ports. Both act as identity on the coin factors. Direct restriction of the active blank bridge gives $b=P_-+8D$, measured relative to the outside level. Write $s=e^{-\epsilon}$. Consider two actuators
$$
b_A(s)=s(P_-+8D),\qquad b_B(s)=sP_-+8D .
$$
The first scales the entire contact pulse. The second adds the contact-occupation phase $8(1-s)D$ to that pulse. The two controls commute on the isolated contact. Each port occupancy is a local load-word projector in the enrolled single excess or vacancy sector. The compensation can be implemented as separate phases of those two words. Both actuators agree with the unperturbed law at $s=1$. [Native seam meters and port phases](../../../../../model/experiments/interval_echo_null_balance/README.md)

Prepare a sharp occupied endpoint of the isolated bridge and read the opposite endpoint after time $t$. Both actuators give $\sin^2(st/2)$ at every time because their Hamiltonians differ by a scalar on the occupied two-port space. Relative to the frozen blank clock the reciprocal first-arrival length is $1/s=e^\epsilon$. On the full board, every coefficient at the minimum graph order for a given source and destination is also identical for A and B. Every factor of the Hamiltonian along such a path must change position. The two generators have identical off-diagonal position blocks. This proves equality of these coefficients even for a selected coin outcome whose coefficient vanishes. An integer replay compares 6,016 matrix entries at graph orders one through three for four launches in both sign and reversal sectors. Longer-exposure interference with outside paths can resolve the occupation phase. [Strain check and first-arrival replay](strain.py)

The work released per unit logarithmic stretch, $W=-\partial_\epsilon H|_0$, is $b$ for A and $P_-$ for B. For a sudden finite change, $H(1)-H(s)=(1-s)W$ exactly. One fresh contact projection at the appointment distinguishes the outside, symmetric and antisymmetric outcomes. Their work scores are $(0,8,9)$ and $(0,0,1)$ for A and B. The outcomes can be read by load inspection and two-port interference. Their difference is the measured contact-occupation energy $8D$. The ordinary-calibrated held ratio divides B's pairing/ordinary work ratio by A's. It would equal one for two work prescriptions related by a single normalization across both source arms. The bounds cover the certified stationary source, the volume adjustment and propagation error. [Bounded computation](certificate.py) · [Rigorous intervals](certificate.json)

| Mobile charge | A ordinary work slope | A pairing work slope | B ordinary work slope | B pairing work slope | Ordinary-calibrated held ratio |
|---|---:|---:|---:|---:|---:|
| Vacancy | −0.369872 | −0.117916 | −0.0343753 | −0.0102204 | [0.93062, 0.93456] |
| Excess | −0.113442 | −0.043910 | −0.000026661 | −0.000068003 | [5.1149, 9.3696] |

The four slope columns are nominal model expectations. A near-zero ordinary slope makes B's excess ratio sensitive to uncertainty. The signed determinant $W'_{A,o}W'_{B,p}-W'_{A,p}W'_{B,o}$ avoids this division. It lies in $[-2.764,-2.696]\,10^{-4}$ for a vacancy and $[6.383,6.701]\,10^{-6}$ for an excess. The independent work-count experiment samples a multinomial with $2^{56}$ trials per preparation. Its certificate and the deterministic enclosure answer different questions. [Counts and simultaneous sampling errors](strain.json)

### Selecting a fixed material

The symmetric and antisymmetric energies are $E_+,E_-$. A separate outside rail fixes $E_{\rm out}$. An interference measurement of the symmetric and outside states gives the dwell frequency $E_+-E_{\rm out}$. A second interference measurement between the symmetric and antisymmetric states gives the exchange frequency $E_--E_+$. Each experiment recombines a balanced superposition after a calibrated exposure. The initial source heads are spectators of the isolated blank contact. The outside reference is essential to the first measurement. Implementation requires coherent splitting and recombination with this rail.

The unstrained contact has the dimensionless index
$$
m=\frac{E_+-E_{\rm out}}{E_--E_+}=8 .
$$
At $s=3/4$ the two isolated end-to-end transfer curves coincide. A gives $m=8$ and B gives $m=32/3$. Independent simulated Ramsey counts give simultaneous 99% intervals $[7.9953,8.0092]$ and $[10.6571,10.6727]$. The declared frequency ranges $[0,16]$ for dwell and $[0,2]$ for exchange remove the cosine alias at the selected exposures. [Calibration counts and frequencies](strain.json)

This gives a uniqueness statement on a precisely specified contact. Assume a local, swap-symmetric, head-independent actuator with fixed outside rail, no hidden levels, and calibrated exchange gap $s>0$. The symmetric and antisymmetric sectors are each one-dimensional in position. Holding the measured material index at 8 forces $E_+-E_{\rm out}=8s$ and $E_--E_{\rm out}=9s$, hence $b_A(s)$ up to a common energy origin. Thus **A's work is the stress conjugate to this contact strain at fixed material index.** Its units are energy per unit logarithmic length. Identifying it as a continuum normal pressure requires a proper stressed volume and a demonstration that this actuator is the local variation of the reconstructed physical metric.

There is a second operational ambiguity even after selecting A. If local lapse $N$ and proper longitudinal stretch $e^\epsilon$ act through $H_e=N e^{-\epsilon}b$, every isolated-contact record depends only on $N e^{-\epsilon}$. Stretching by $4/3$ at fixed lapse and slowing the local clock to $3/4$ at fixed length give the same contact operator and the same material index. Distinguishing these two changes requires an independently calibrated clock or ruler exposed to a different combination of the two factors. This applies even with unlimited samples and complete tomography of this single contact.

A nearby-clock control distinguishes the imposed stretch from a *common regional* slowing on the same post-flight mixture. At the four-tick appointment freeze the flight, read the mobile position, and retain both the initial endpoint and the outside outcome for the three disjoint blank dimers $(5,9),(1,38),(40,47)$. Each inside result starts a sharp endpoint. Pulse the dimers in isolation and read whether the carrier changes endpoint. For each dimer $j$ the conditional probability is $\sin^2(\omega_j t/2)$ for every initial joint head state. A bridge stretch of $4/3$ gives native frequencies $(1,3/4,1)$. A common regional lapse of $3/4$ gives $(3/4,3/4,3/4)$. For the bridge stretch, the two side-dimer frequencies lie between 0.9988 and 1.0014 and the central one between 0.7464 and 0.7535 across the two mobile species in the simultaneous 99% count intervals. A lapse confined to the bridge still agrees with the local stretch experiment. [Initial enrollments, misses, endpoint flips and comparison clocks](strain.json)

### The even null comparison

The specified bridge action $H_e=N e^{-\epsilon}b$ has lapse derivative $b$ and negative longitudinal strain derivative $b$ at the reference settings. If these represent the integrated energy and normal pressure on one common proper element, their sum is $2b$ in the corresponding units. Summing the forward and reverse null contractions cancels the directional energy current. The independent native-time optical quantity is the sum of the two contractions $-\mathrm{tr}(\ddot J J^{-1})$ of the measured image maps. Fix the material conversion on the ordinary source, keep the affine-time and pressure normalization fixed, then test the pairing arm. The measured focusing divided by its prediction from the work meter lies in **[0.98491, 0.99490] for a vacancy** and **[0.77909, 0.78841] for an excess**. Both exclude one even when each species has its own ordinary calibration. These results reject this bridge-only local source and clock assignment for the screen reading. Extended regional energy, independent physical time or different proper-volume normalization lead to distinct tests. [Forward-plus-reverse source analysis](../../jacobi_screen/source_response.json) · [Work calibration](certificate.json)

## Independent impulse on a material cell

Let $x_p$ measure distance inward from the gateway in one of the two adjacent cells. To define the internal energy moment, allocate the energy of a native internal seam equally to its two endpoints and keep each degree-reference energy at its port. In units with the blank speed scale set to one,
$$
B=\sum_{e=(p,q)\ \mathrm{internal}}\frac{x_p+x_q}{2}h_e
+\sum_{p}x_p h^{\rm ref}_p,\qquad P=i[h_c,B].
$$
This $P$ is the rate of change of the cell's internal energy moment during isolated evolution. Adding a constant to every $x_p$ adds that constant times $h_c$ to $B$, which leaves $P$ unchanged. Constant parked energies give the same cancellation within each reversal sector. The bridge is the external actuator. Energy stored on the other incident bridges belongs to the interrupted free-flight apparatus. The meter reads internal-cell motion, with no claim of vanishing surface momentum. [Operator construction](measure.py)

At the four-tick appointment apply only the bridge pulse for $u$ and then only the internal cell action for $v$. Read the resulting energy moment $B$. Taking the mixed derivative at $(0,0)$ gives
$$
F=\partial_u\partial_v\langle B\rangle=\langle i[b,P]\rangle
=-\langle[b,[h_c,B]]\rangle .
$$
It is the rate at which the bridge pulse changes the subsequently measured internal energy-moment velocity. The finite replay also reads the four interrupted histories $(u,v),(u,0),(0,v),(0,0)$ using positive durations. Richardson extrapolation of the mixed forward difference at steps $1/1024$ and $1/2048$ recovers every nominal force or slope within 0.0020 native units. These impulse results are simulated population expectations. The work and Ramsey records above include simulated detector counts. [Interrupted pulse readings](results.json)

Two coordinate calibrations use the same unit blank link. The hop ruler assigns each traversed internal seam one unit. The front ruler uses reciprocal charge-pooled first-arrival coefficients measured on the reference heads. The left and right marked links read 0.172207 and 0.176777 blank-link units. Calibrating to the nearest $1/4096$ blank-link unit selects numerators 705 and 724 across the full source enclosure. All other internal links have weight one. Shortest-path distances start at each gateway. Both coordinate tables are frozen across source variations. [Reference-ruler proof](certificate.json) · [All site coordinates](results.json)

| Mobile charge | Receiving ruler | Ordinary impulse derivative | Pairing impulse derivative |
|---|---|---:|---:|
| Vacancy | Hop | [−84.10, −82.87] | [−25.40, −24.89] |
| Vacancy | Front | [3.905, 3.981] | [0.8249, 0.8518] |
| Excess | Hop | [71.44, 72.46] | [22.22, 22.65] |
| Excess | Front | [−13.50, −13.27] | [−3.893, −3.801] |

The sign change excludes replacing one ruler by the other with a fixed positive unit conversion. The results for both sides use coordinates directed from each gateway into its own cell. A mechanical identification requires an independently justified spatial coordinate, any material recoil and the storage on the contact. Separate impulse and strain records distinguish these requirements experimentally.

There is also an exact limit on a closed finite-system interpretation. For a bounded energy dipole $B$, $P=i[H,B]$ has matrix elements $i(E_m-E_n)B_{mn}$. Conservation $[H,P]=0$ kills every element between different energies, and $P$ vanishes within each degenerate block by construction. Thus a nonzero $P$ of this form cannot be an exactly conserved global translation charge of the finite Hamiltonian. A local momentum balance can include boundaries, winding transport or forces on additional material. Establishing that local balance on reconstructed balls is separate from the contact measurement.

From the repository root run `PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.directional_stress.traction.measure --write`. The module names `strain` and `certificate` write the work, calibration and rigorous-enclosure receipts. Run `PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/traction` for the focused regressions.
