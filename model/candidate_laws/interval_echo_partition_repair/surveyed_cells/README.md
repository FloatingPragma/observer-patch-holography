# Area and entropy in overlapping cubical regions

A material region has one freely readable binary register at the center of each internal face. A separate survey measures the lengths of cube edges, the proper area of exposed faces and the proper volume enclosed by them. In the reversible repair described below, **one equilibrium makes bulk entropy plus an area term stationary for every union of enrolled material cells simultaneously**. An ordinary excitation calibrates the coefficient. A change in three-register correlation passes the held prediction on an overlapping region without changing the ordinary energy. These are proper areas in a supplied flat cubical geometry. [Frozen comparison](PROTOCOL.md) · [Terminal records](results.json) · [Implementation](measure.py)

## A single repair and its regional readers

Use 216 identical unit cubes in a flat periodic $6^3$ box. The center of each of the 648 shared faces carries a bit. Fix three other bits $h,y,z$ at the origin. A separate excitation meter reads $h$ in a standard with energies 0 and 1. The assembly control uses $g=h\mathbin{\mathrm{XOR}}(y\mathbin{\mathrm{AND}}z)$. In sector 1 two declared pairs of material cubes exchange their positions. The geometry can be read from the labeled cubes and their rods. Each update compares the completed assemblies. The actuator prescribes the exchanges; the model of the survey contains no trajectory of a moving cube.

A cycle offers an independent fair flip to every register. On the controller bits a change of $g$ actuates the indicated exchange, leaving the face memories at their geometric positions. The maps on the complete classical records are reversible. After one full sweep the actual $3+648$ bits are independent and fair from any initial distribution: at each register the two possible outcomes have equal conditional weight. The sector is then automatically fair. The same equilibrium supplies every observer.

A cut is any fixed set of material labels, including disconnected and overlapping sets. Its bulk observer reads its labeled assembly, every face register strictly inside it, and all three controller bits whenever the origin belongs to the cut. A register exactly on the region boundary is outside this bulk algebra. All observed cells and words are retained. A second instrument constructs the eight corners from each cube's rods, pairs faces with identical vertices across periodic faces and sums the unpaired parallelogram areas. The volume meter sums the proper rod determinants. Neither geometric meter reads a face bit or an entropy.

## Exact balance and stability under arbitrary source variations

Let $n,b,m$ denote the numbers of enrolled cubes, exposed faces and internal faces. Counting face incidences gives $6n=2m+b$. In rods of length $h_0$, write $a=h_0^2$ and $v=h_0^3$, so $A=ab$ and $V=vn$. Every regional geometry/controller outcome has the same marginal probability $1/c_R$ under the common equilibrium. If the trigger cube is included, the whole tag triple is seen and $c_R=8$. Otherwise the local geometry distinguishes either both sectors or neither, hence $c_R=2$ or $1$. This uses fixed material membership and the binary geometry selector. Retaining only some controller bits or admitting the trigger in a sector-dependent way requires a separate multiplicity check. The conditional equilibrium on the $m$ face registers is uniform. For a complete output $w$ the regional reference law is thus $P^*_R(w)=1/(c_R2^{m(w)})$.

Set $\eta=(\log2)/(2a)$ and $\lambda=3\log2/v$. For any global preparation with regional law $Q_R$, direct expansion of relative entropy gives

$$
H(Q_R)+\eta\,\mathbb E_Q A_R
=\log c_R+\lambda\,\mathbb E_Q V_R-D(Q_R\Vert P_R^*).
$$

All cuts use the same global distribution and the same $\eta$. The overlap observer can obtain its record by restriction from either larger window. Every two-sided differentiable source variation through the equilibrium at fixed expected proper volume therefore has $\delta S+\eta\,\delta A=0$. This includes arbitrary changes of the joint bit law, not just changes of the prepared tag frequencies. In each fixed-label experiment here the volume is even constant on individual trials. Replacing a bit by a uniform $q$-state face memory replaces $\log2$ by $\log q$. A boundary capacity per unit area and this choice of the strictly interior bulk instrument are physical inputs to the comparison. [Probability and restriction checks](test_survey.py)

## Independent ordinary and held measurements

The ordinary variation gives $h$ probability $1/2+t$ and leaves $y,z$ fair and independent. Its slope in sector probability is $1/2$. The held law gives the triple $(h,y,z)$ weight $[1+t(-1)^{h+y+z}]/8$. Every one- and two-bit marginal stays fixed, including the ordinary excitation, while the sector slope is $1/4$. The ordinary and held regions share three cells. The third cut is a two-cell subset of their intersection and the same regional outcome can be recovered from either set of records.

| Enrolled material region and source | $V$ | $A$ in sectors 0,1 | $\delta E$ | $\delta S$ (nats) | $\delta A$ |
|---|---:|---:|---:|---:|---:|
| Four-cell line, ordinary | 4 | 18, 16 | 1 | $(\log2)/2$ | −1 |
| Six-cell sheet, held correlation | 6 | 22, 26 | 0 | $-(\log2)/2$ | +1 |
| Two-cell overlap, held correlation | 2 | 10, 12 | 0 | $-(\log2)/4$ | +1/2 |

At $t=\pm1/8$ each bank independently samples $2^{32}$ terminal copies from the specified finite law. The entropy reader uses the complete regional histogram. The face/volume reader uses the measured sector rods, and the excitation reader uses the separate binary meter. The two-point secants equal the reference slopes: the entropy of the tag triple on the ordinary variation and the entropy of the sector on the held variation are even in $t$, while the mean number of internal faces is affine. The ordinary gain inferred from these counts lies in $[0.34006,0.35310]$ nats per square unit. Applied without refitting to the held region it gives an entropy-plus-area residual in $[-0.01319,+0.01327]$ nats. The ordinary energy slope on the held intervention is in $[-0.000278,+0.000269]$.

These simultaneous intervals have at least 99% coverage. Changing one draw changes empirical entropy by at most $(1+\log N)/N$. McDiarmid's inequality bounds the sampling error about the mean by $(1+\log N)\sqrt{5/N}$ except with probability $2e^{-10}$. With $k$ possible words the downward bias is between zero and $\log(1+(k-1)/N)$ by the expected empirical relative entropy, its chi-squared bound and Jensen's inequality. Each geometric or energy sector frequency uses Hoeffding radius $\sqrt{5/N}$. A union bound over eighteen marginal intervals covers the transported held gain. This simulation draws compressed terminal multiplicities; it does not evolve billions of separate lattice arrays.

Biasing acceptance multiplies the equilibrium odds of sector 1 by two. The same tag and geometry proposals obey detailed balance with respect to these biased weights. For a region resolving the selector the extra term in the fixed-volume differential is exactly $-(\log2)\,\delta\Pr(g=1)$. The ordinary residual becomes $-(\log2)/2$ and the held residual $-(\log2)/4$. The sign also depends on the readout convention: counting all boundary face registers as bulk words adds $b\log2$ and reverses the sign of the area contribution. [Changed kernel and tests](test_survey.py)

## Cropping a genuinely smooth ball

A Euclidean ball of radius $r$ has true area $4\pi r^2$. With a uniform shift of the lattice by one grid period $h_0$, the expected number of actual face centers strictly inside any fixed smooth region is $3V/h_0^3$ by integration over the shift. That smooth-domain extension of the point-register reader samples volume alone. A different extension accepts cells by their centers and counts faces crossing the digital cut. Each axis gives two oriented contributions from points in the ball whose translate by one spacing lies outside it. For $0<h_0\le2r$ the intersection of the two smooth balls is two spherical caps, with volume $V-\pi r^2h_0+\pi h_0^3/12$. Consequently

$$
\mathbb E[\text{exposed cubical area}]=\tfrac32(4\pi r^2)-\tfrac{\pi}{2}h_0^2.
$$

For a plane with unit normal $n$ the limiting projected area factor is instead $|n_1|+|n_2|+|n_3|$. These facts distinguish smooth cropping from measurement on actual exposed cube faces. A surface-dependent conversion or mesh-dependent gain does not establish a universal ordinary coupling. At fixed alphabet size $\eta$ itself varies as $h_0^{-2}$ if the physical spacing is reduced. A shrinking curved proper-ball realization would have to control all these conversions below the fourth-order signal at a fixed physical gain, locate the boosted regions and independently measure the action stress there. This model establishes none of those identifications. [Exact crop calculation](cropped.py)
