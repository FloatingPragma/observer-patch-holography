# The same twelve pulse strengths on two physical carriers

The registered one-load excitation and the two-load bound pair each have a twelve-dimensional current response. When exactly the same signed strength at each of the twelve addresses acts on both preparations, their joint connected response has **dimension 23**. A full central turn calibrated on the pair is exactly the identity there. It sends the one-load launch to ten other ports, disjoint from its two input ports. The two results follow from the existing Interval Echo current service at its common unit gain. [Current specification](../../../candidate_laws/interval_echo/CURRENTS.md) · [Executable proof and detector](measure.py) · [Complete records](results.json)

## The central turn as a physical test

In a sector with $N$ loads, let $P_N$ be its retained first triplet projector. The circulation terms sum to zero because the two neighbor rings using each directed edge have opposite orientations there. The native occupations sum to $N I$. The positive class-sum projector acts as identity on the first triplet and annihilates the second. The protected address generators therefore satisfy

$$
\sum_{p=0}^{11} H^{(N)}_p=N P_N .
$$

The same summed exposure $\pi$ acts by $I-2P_1$ on a single load and by the identity on the entire pair sector, including its 60 spectator modes. The one-load odd launch is $(|p\rangle-|\bar p\rangle)/\sqrt2$. On this subspace of the icosahedral cell, the adjacency $A$ has eigenvalues $\pm\sqrt5$ and $P_1=(I+A/\sqrt5)/2$. The changed output is $-A/\sqrt5$ times the launch. The two original sites have distinct, disjoint sets of five neighbors. The sham therefore has mass one on the original antipodal pair and the central word has mass zero there. This distinguishes the central word by measurements within one load sector. No phase reference between different total loads is required.

The physical address schedule uses complete symmetric ascending and descending sweeps. Let $L=\pi\sum_p\|H_p\|_{\rm op}$. The first two coefficients of one symmetric slice agree with those of the summed exponential. The individual exponential series and the ordered product each have remainder bounded by $e^{L/n}(L/n)^3/6$ after their quadratic terms. Telescoping the $n$ unitary slices bounds the operator error by $e^{L/n}L^3/(3n^2)$. Proper rotations give the same current norms at all ports. The measured trace squares certified below and the exact first-triplet trace $N/4$ bound every current norm by 1.4 for one load and 0.81 for the pair. With 2,048 sweeps the resulting rigorous error bounds are 0.0121 and 0.0024 respectively. The original two sites have total one-load probability below $0.0121^2$ in the physical pulse arm.

| One-load record | Sham | Central word | Occupation gain zero |
| --- | ---: | ---: | ---: |
| Detections on the input antipodal pair before free flight, 8,192 simulated starts per setting | 8,192 | 0 | 8,192 |

At the registered free exposure of 16, the numerical one-load total variation for the same 2,048 sweeps is 0.34405114. The ideal integrated turn gives 0.34405098. The changed gain deletes the real occupation term in each protected band, keeping the same circulation, geometry, preparations and addresses. Its integrated sum vanishes on both carriers. The receipt includes all twelve one-load outcomes, all thirty pair seam outcomes and quiet. [Counts, numerical pulse errors and analytic upper bounds](results.json)

## The two independent simple-factor tests

For a real triplet the current on the second band is imaginary skew, so its anti-Hermitian generator is in $\mathfrak{so}(3)$. On the first band it is in $\mathfrak u(3)$. In either carrier, the twelve generators jointly span these spaces. To see whether two such responses could attach to one Lie algebra with the same addresses, remove the trace from the port-0 generator in each simple factor and compute $\operatorname{Tr}(H_0^2)$ there. The entries below are certified enclosures.

| Simple factor of the local current | One load | Bound pair |
| --- | ---: | ---: |
| $\mathfrak{su}(3)$ on first triplet | [1.849, 1.852] | [0.556, 0.559] |
| $\mathfrak{so}(3)$ on second triplet | [0.689, 0.692] | [0.114, 0.117] |

The compressed rotations are anti-Hermitian after multiplication by $i$. On each given simple type these trace squares are the negative Killing form times the same constant in the two carriers. Any Lie automorphism preserves the Killing form, including the conjugating outer automorphism of $\mathfrak{su}(3)$. The disjoint intervals rule out identifying either pair of simple factors under the registered address action.

More explicitly, the compact joint algebra projects onto every single simple factor. Each simple ideal of its derived algebra maps either isomorphically onto a simple image factor or to zero. Surjectivity of the whole projection makes each nonzero ideal image an ideal of that factor, and two nonzero commuting images cannot both fill a simple factor. Thus the only possible identifications between image factors come from a common simple ideal mapping to both. The unequal trace squares exclude this for the two copies of each type. The different types have unequal dimensions. Consequently the derived algebra contains two independent copies each of $\mathfrak{su}(3)$ and $\mathfrak{so}(3)$.

Port transitivity and the summed-current identity give first-triplet trace $N/4$ at each port. Every commutator is traceless in each block. The shared central generator thus has weights 1 and 2 and is the entire central line. The full connected matrix image on the two registered sectors is

$$
\{(U_1,U_2)\in U(3)^2:\det U_2=(\det U_1)^2\}
 \times SO(3)^2.
$$

Its dimension is $8+8+3+3+1=23$. For covering coordinates $(A_1,B_1,A_2,B_2,z)$ in $SU(3)\times SU(2)\times SU(3)\times SU(2)\times U(1)$ the map is $(zA_1,\operatorname{Ad} B_1,z^2 A_2,\operatorname{Ad} B_2)$. Its kernel consists of $(z^{-1}I_3,\pm I_2,z^{-2}I_3,\pm I_2,z)$ with $z^3=1$ and independent signs. This rational circle and the compact simple factors give a closed connected image. Both weak centers act trivially. The finite group statement concerns the guarded address service viewed on the direct sum and makes no assumption of coherent preparation across sectors.

## Certificate for the native bands

The source sectors use the integer symmetric operators $2H$. For one load the source is annihilated by $x^2-10x+20$ evaluated on $2H$. For the pair the polynomial is $x^4+34x^3-128x^2-6904x+41936$. The verifier obtains the factors from the full characteristic polynomials and multiplies them into every integer odd launch. Exact real-root counts isolate three eigenvalues, with multiplicity, in each interval of radius $2/25$ about the enrolled energies of $2H$. The other sourced pair triplets have agreement below one half.

For each isolated triple, round the three eigenvectors to the columns of a rational matrix $B$ with denominator $10^{11}$. Direct integer products bound both $\|B^TB-I\|_F$ and $\|(2H-t)B\|_F$ below $10^{-8}$ at its recorded rational center $t$. The polar isometry $W$ of $B$ has $\|W-B\|<10^{-8}$. For the exact spectral projection $P$, the complement gap $2/25$ implies $\|(I-P)W\|\le 2(10^{-8})/(2/25)$. Taking the polar isometry of $PW$ gives an exact basis $U$ of the whole triple with $\|U-B\|<10^{-6}$. This follows by bounding the tilt by $\|(I-P)W\|$ and the polar normalization correction by a second copy of that bound.

All local operators before compression have row norm at most five and are Hermitian. Their true and rational compressions in these matched bases differ in operator norm by less than $1.1\times10^{-5}$. Traceless projection is an orthogonal projection for the Hilbert–Schmidt norm, so each squared norm in the table differs from its rational center by less than 0.001. The rank check uses an exact positive-definiteness test on the rational twelve-column Gram matrix minus $I/20$. The twelve columns together perturb its singular-value map by less than $\sqrt{12}\sqrt6(1.1\times10^{-5})$. This is smaller than $\sqrt{1/20}$ and proves full rank. The commuting integer class sum has exact positive or negative sign on each enrolled triple by the same basis estimates and rational three-by-three definiteness checks. For the source filter, the computed nonzero signed agreement scores exceed $10^{-5}$ in magnitude, larger than the projection error in the normalized quadratic form. [Rational verification](measure.py) · [Regression](test_joint.py)

This places an obstruction on the FT-GLOBAL transfer from the measured pair to the independently accessible one-load source at the same physical pulse strengths. Any common twelve-dimensional symmetry acting on a proposed larger matter family through these same representations would project to an image of dimension at most twelve. These two restrictions alone have image dimension 23. A separate transport or flux action could supply other generators, but it requires its own common-source measurements. The measured sectors do not realize a weak doublet.

From the repository root, replay the calculation and its regression with

    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m model.experiments.interval_echo_gauge.joint_currents.measure
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/joint_currents
