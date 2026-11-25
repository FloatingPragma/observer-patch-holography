# What a material-written strain register conveys

A mobile load and two correlated seam memories can write a detectable change into a qubit on the bridge between them. On the 80-cell federation, calibrating that change on ordinary rail energy predicts an advance of a long-wave density signal when the memories are decorrelated. The signal is delayed. An independent clock gives the same sign disagreement. A finite rank calculation identifies information that the terminal bridge qubit loses. A second experiment recovers the missing modular information from repeated preparations with a directional launch and multiple appointments. [Source and pump](source/README.md) · [Wave and clock responses](results.json) · [Recovery from simulated counts](history.json)

When the next probe keeps both seam memories and the strain, the registered change becomes visible. On the same live source, two marginal-preserving rotations keep five-station area and volume fixed while the retained entropy increases. [Live geometry, held stations and fixed-volume proof](LIVE_GEOMETRY.md)

## One change to the contact action

The pump uses the native Interval Echo generator with a single extra two-state register on one blank bridge. The experiment sets its frequency parameter to $\Omega=1/2$, its coupling to $\epsilon=\log(5/4)$, and its initial state to $|+X\rangle$. With $D$ the projection onto the two bridge endpoints and $X_{ab}$ their exchange, define

$$
 b=8D+(D-X_{ab})/2,\qquad
 H=I\otimes H_{\rm IE}-\Omega X\otimes I+(e^{-\epsilon Z}-I)\otimes b .
$$

The native outside, symmetric and antisymmetric bridge energies are 0, 8 and 9. The coupling multiplies this whole action, with the same outside energy zero for every source. Ordinary energy changes the strain, and the strain changes the subsequent phase and hopping of a material probe. The register, its parameters, the selected bridge and the preparation and read appointments are design choices of this prototype. [Native contact energy and its outside reference](../../source_implication/directional_stress/traction/README.md) · [Active law](../../../../model/candidate_laws/interval_echo/SPEC.md)

One mobile excess or vacancy forms a joint two-head state $\sigma$ in four native ticks. The common pump reference uses these heads and a 2:1 lottery of symmetric and antisymmetric mobile bridge rails. The ordinary arm reverses the lottery. The collar arm independently pairs the two heads from the same formation ensemble, preserving both single-head marginals and the rail lottery. The ordinary energy increases by $1/3$ in the ordinary arm and is constant in the collar arm. A separate modular read uses $K=-\log(\tau_{\rm rail}\otimes\sigma)$ on the prepared rail and head algebra. The ordinary change is $(\log 2)/3$ nats. For vacancy formation at 80 cells the collar change is $+0.200825$ nats. All three preparations propagate freely under the same pump Hamiltonian. [Source energy, modular measurements and certified transfer errors](source/README.md)

The read stage discards the pump mobile and head degrees of freedom, retains their reduced strain state, and supplies an independent blank probe on the same port board. This operation isolates the information that passed through that single register. A probe that continues to interact with the original joint heads can access additional information. The finite records below use supplied release and analysis operations.

## Independently reading the phase and the spatial signal

A single load on a blank, empty host has Hamiltonian $48I+L/2$, with $L$ the port graph Laplacian. The constant 48 comes from the native degree reference. Let $g$ be the normalized constant mode. Project the oriented endpoint difference $|b\rangle-|a\rangle$ into the first three nonconstant modes of $L/2$ and normalize to $u$. The board fixes this choice independently of the pumped states. The three eigenvalues may differ. Set $E=\langle u|L/2|u\rangle$ and $F_p=u_p/\max|u|$.

* The transport probe launches $(g+iu)/\sqrt2$ and $(g-iu)/\sqrt2$ separately. Every terminal port reports its score $F_p$. Half the difference between these two means measures the spatial movement of the density disturbance.
* The band clock launches $(g+u)/\sqrt2$. Two phase analyzers on $g,u$ give a signed Ramsey quadrature. Amplitude outside their subspace scores zero and every preparation enters the average.
* A local clock compares the symmetric bridge level to the outside rail. A separate directed bridge transfer compares opposite launch phases at the endpoints. Both read the local multiplier at release after normalization by their respective reference frequencies.

Write $c=\operatorname{Tr}(e^{-\epsilon Z}\rho_{\rm strain})$. The initial normalized response of each band read is $1+(c-1)R$. The transport coefficient follows from the commutator with the port score. Diagonal contact energy commutes with that score. The band clock reads the difference of diagonal energies on $u$ and $g$, which includes the contact occupation term. The respective gains are

$$
 R_{\rm travel}=\frac{(u_b-u_a)^2}{2E},\qquad
 R_{\rm clock}=R_{\rm travel}
             +\frac{8(u_a^2+u_b^2-2/n)}{E},
$$

where $n$ is the number of ports. The two isolated-contact gains are each one. These are fractional changes of finite instrument rates per unit change of the bridge multiplier.

| Cells | Mean slow gap $E$ | Travel gain | Band clock gain |
| ---: | ---: | ---: | ---: |
| 20 | 0.022142 | +0.069635 | −2.383895 |
| 80 | 0.0051765 | +0.019903 | −3.030101 |
| 320 | 0.0012809 | +0.004095 | −3.167303 |

The discrepancy has an exact polarization test. Every real three-dimensional slow subspace contains a unit vector $u$ vanishing at both marked endpoints. On this vector the contact has zero contribution to the spatial current, while its contribution to the normalized gap against $g$ is $-16/(nE)$. A universal ratio between those two global sensitivities cannot hold over all polarizations. A proper clock and a local ray in an independently reconstructed metric require their own identification. [Independent population and control tests](test_readout.py)

## A held prediction at finite travel time

For a signal or clock record $S(t)$, use the reference level at each of two appointments $t_*=\pi/(6E),\pi/(3E)$. On the local crossing branch the change in log inverse duration per unit admixture of a source arm is $\delta S(t_*)/[t_* S'_0(t_*)]$. This follows by differentiating the fixed-level condition and requires a nonzero reference slope. It allows a comparison at a finite fraction of the unperturbed beat without assigning a speed to a graph edge. The simulator evolves the full qubit and port state throughout the read. It retains scattering into all port modes. The receipt includes every reference slope, including small slopes that amplify the crossing response.

The following example uses 80 cells, a vacancy for formation and for the pump, and a four-tick pump. The reference ratio of collar to ordinary modular changes is 0.869187. Multiplying the measured ordinary response by this ratio fixes the collar prediction for each instrument independently. The first line is an absolute change in the dimensionless contact multiplier. The last two are changes in log inverse duration at the first band appointment, 101.149 native time units.

| Read | Ordinary | Collar from ordinary calibration | Collar from independent evolution |
| --- | ---: | ---: | ---: |
| Bridge multiplier at release | +0.00224697 | +0.00195303 | +0.00047401 |
| Long-wave travel | +0.00731885 | +0.00636145 | −0.00035393 |
| Band clock | +0.00183847 | +0.00159798 | −0.00008937 |

The experiment covers 32 combinations of board size (20 or 80 cells), two formation charges, two pump charges and four pump times. At the two band appointments the travel prediction has the opposite sign to the computed response in 55 of 64 cases, and the clock prediction in 54 of 64. These are simulated population means. The 320-cell entries are release gains. Independent dense evolution checks the 80-cell example against the sparse evolution, and the uncoupled control reproduces the sum of the three free modal beats. [All source comparisons and read effects](results.json) · [Native collective control on stationary sources](../../../../model/experiments/seam_program/collective_rate/long_wave/README.md)

## What one final strain state loses

Fix the same 80-cell, four-tick, vacancy source. A change $D$ of the two-head density that preserves both marginals lies in the nine-dimensional real span of $P\otimes Q/4$, with $P,Q\in\{X,Y,Z\}$. Holding the mobile lottery fixed, let $M D$ be the three changes in the terminal strain Bloch vector and $hD=\operatorname{Tr}[(-\log\sigma)D]$ the independent modular slope. If all such slopes could be read from the terminal qubit, $h$ would vanish on the kernel of $M$.

The three qubit rows have rank 3. Adding the modular row gives rank 4. The rational certificate encloses a four-column minor with $\|A_0^{-1}(A-A_0)\|_\infty<0.000058$, so the true minor is invertible. There exists an exact $D$ with $MD=0$ and $hD\ne0$. Faithfulness of $\sigma$ allows both $\sigma+\eta D$ and $\sigma-\eta D$ at sufficiently small positive $\eta$. They give exactly the same final qubit. Every subsequent measurement receiving only that qubit, including an arbitrarily precise clock or flight measurement, has identical outcome laws on this pair. The kernel statement concerns the full algebra of fixed-marginal perturbations. The table above separately tests decorrelation prepared directly by re-pairing actual histories. [Rank certificate](hidden.json) · [Enclosure calculation](hidden.py)

## Recovering the information with a directional measurement

On independently prepared copies of the same formed head source, use two mobile launches, the reference rail lottery and a sharp launch from the left bridge endpoint. Read all three strain components after $1/2,1,2,4$ pump ticks. The resulting 24 real measurements span all nine fixed-marginal correlation directions. A rationally enclosed nine-row minor has inverse error below 0.00565. Repeated preparations and a directional launch suffice for this finite tomography experiment. The measurement costs 24 settings and an ensemble of copies. [Settings and certified decoder](history.py)

Calibrate a linear decoder for the vacancy-formed modular observable on the nine correlation basis directions. The worst-case modular decoding error for any difference of two densities with matching marginals is below 0.000015 nats at exact populations. The decoder has Euclidean weight norm 185.28. With $2^{36}$ independent simulated binary trials per arm and setting, its simultaneous 99% counting allowance for two comparisons is 0.00347 nats. The second comparison generates its heads with an excess and uses the same vacancy-calibrated logarithm and decoder.

| Forming load | Direct charge change, nats | Decoded from simulated counts, nats |
| --- | ---: | ---: |
| Vacancy | +0.200825 | +0.199619 |
| Excess | +0.259754 | +0.260321 |

The recovered charge can enter the finite-time test without reading a simulator density matrix for the held arm. Substituting the vacancy count estimate and its two allowances into the ordinary calibration gives a predicted travel response between +0.00621 and +0.00644 and a predicted clock response between +0.00156 and +0.00162 at 101.149 ticks. Both computed collar responses in the finite-time table are negative. These intervals propagate the modular measurement uncertainty with the ordinary response taken at its calculated population value.

This procedure supplies an independent material measurement of a specified regional charge. Relating it to gravity requires a persistent source and an independent geometric response on the same proper balls. In particular the bulk modular term and the fixed-volume area term must cancel with the ordinary-matter coupling and errors smaller than the leading fourth-order area signal. The supplied read host, global phase analyzer and pump clock in this experiment provide no proper-ball volume or pressure measurement. They therefore furnish no Newton constant or galaxy prediction. [Einstein premise chain](../../PROOF.md) · [Proper-volume clock and area test](../README.md) · [Complete synthetic counts](history.json)

Run from the repository root with the model requirements installed.

~~~sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -m matched_einstein_limit.rate_geometry.seam_strain.source.measure --write
python -m matched_einstein_limit.rate_geometry.seam_strain.source.certificate --write
python -m matched_einstein_limit.rate_geometry.seam_strain.measure --write
python -m matched_einstein_limit.rate_geometry.seam_strain.hidden --write
python -m matched_einstein_limit.rate_geometry.seam_strain.history --write
python -m pytest -q -W error code/matched_einstein_limit/rate_geometry/seam_strain
~~~
