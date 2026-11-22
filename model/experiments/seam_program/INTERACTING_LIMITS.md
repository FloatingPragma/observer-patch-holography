# Which interactions survive spatial refinement?

A **local, scale-invariant two-load scattering sector** exists on a programmed cubic guide. Its threshold occurs at effective contact strength $\kappa_*=0.9682678127$. Time scales as distance squared. On the propagating eight-slot chart a different program maintains an interaction at ballistic clock scaling by surveying a fixed physical neighborhood per opportunity. Its range is a controller input. A local, interacting relativistic conformal vacuum inferred by Interval Echo observers remains open. [Cubic scattering and packet records](cubic_resonance.json) · [Ballistic range comparison](range_scaling.json)

## Local resonance on a three-dimensional cubic guide

Place binary loads on a homogeneous periodic cubic address graph. Give its sites spacing $h$ and use $h^2$ units of physical time per native exposure. Isolation, routing, clock conversion and the uniform distinct-pair launch are supplied. The blank-seam exchange and contrast give the exact binary pair operator $K-2\kappa_{\rm eff}P_{\rm adj}$ up to a constant. $K$ is half the pair-configuration graph Laplacian. At zero total momentum the relative coordinate $r$ evolves by the cubic graph Laplacian on $r\ne0$, with no edge into the occupied origin, and the potential $-2\kappa_{\rm eff}$ on its six neighbors. Comparing the full native pair matrix with this reduction on a complex symmetric state checks the signs, hopping rates and excluded site. [Native sector](../../candidate_laws/contrast_memory_spherical_v1/sectors.py) · [Independent controls](test_kinetic_resonance.py)

The potential can be adjusted with the *same* physical seam at $\kappa=8$. A native exposure $\delta$ followed by the occupation phase below of angle $2(\kappa_{\rm eff}-8)\delta$ yields the selected generator on one edge, up to its singleton reference phase. On a regular guide the sum of these references is constant at fixed load number. The two edge terms commute on that edge. Even periodic cubes admit six disjoint edge matchings. At fixed two-load number each layer has bounded operator norm independent of side length. Splitting the layers into slices of native duration $\delta=h^3$ gives accumulated operator error $O(h)$ over physical time of order one by the product-formula commutator bound. The computation evolves the exact logical generator. A physical 20-cell run of the counterpulse for one resonance-adjusted contact at native slice $1/64$ returned the helpers and agreed on all four binary input columns within $3.34\times10^{-12}$. The counterpulse and the threshold setting are prescribed apparatus choices.

### Finding the threshold from the native operator

Let $G_z(r)$ be the full cubic Green function of $L-z$, with symbol $\epsilon(p)=2\sum_j(1-\cos p_j)$. Write $G_0=G_{z=0}(0)$ and $G_1=G_{z=0}(e)=G_0-1/6$ for a nearest neighbor $e$. At zero energy extend a cubic-symmetric unit-incident wave to the excluded origin by its common neighbor value $B$. This keeps the same Laplacian on the six physical neighbors. The extended equation is $L\psi=2\kappa_{\rm eff}B\,1_{\{|r|_1=1\}}$. Summing the six translated Green functions at a neighbor gives $6G_1$, so $B=(1-12\kappa_{\rm eff}G_1)^{-1}$. At large radius $G_{z=0}(r)\sim(4\pi |r|)^{-1}$. In the convention $\psi(r)\sim1-a/(h|r|)$ the physical scattering length and the pole are

$$
a=-\frac{3h\kappa_{\rm eff}}{\pi(1-\kappa_{\rm eff}/\kappa_*)},
\qquad \kappa_*=\frac{1}{12(G_0-1/6)} .
$$

The onsite heat-kernel integral gives $G_0=0.252731009858663$. The pole is selected by the two-body threshold equation before any packet comparison. On the uncorrected $\kappa=8$ guide the scattering length tends to zero. The cubic-symmetric channel also has a bound level below threshold at that strength. The below-threshold denominator crosses zero as the attraction grows through $\kappa_*$, and the attractive potential is monotone. The statement about vanishing scattering concerns fixed physical momenta on the open scattering branch. [Enrolled arms](cubic_resonance_protocol.json) · [Green functions and finite-cube release](cubic_resonance.py)

At energy $z$, write $u=1-z/6$, $\Delta=1+2\kappa_{\rm eff}$, and $g_0=G_z(0)$, $g_1=G_z(e)$. Then $g_1=u g_0-1/6$. Imposing the missing origin by a Dirichlet projection puts the uniform shell eigenvalue of the projected free resolvent at $g_1/g_0$. The residual shell potential is $-\Delta$. For a cubic-symmetrized incident plane wave of value one at the origin and $u$ at each neighbor, eliminating the origin source gives the on-shell source strength

$$
Q(z)=\frac{\Delta u-1}{\Delta/6-(\Delta u-1)g_0}.
$$

Indeed, a uniform shell response $B$ has origin source $b_0=-(1+6\Delta B g_1)/g_0$. At a neighbor the response is $B=u-g_1/g_0+\Delta B g_1/g_0$, and its outgoing on-shell source is $b_0+6u\Delta B=Q(z)$. An independent sparse solution on an $8^3$ punctured finite torus agrees with the corresponding finite-volume shell resolvent to $6\times10^{-16}$ at negative energy.

For fixed outgoing physical momentum $k>0$, $g_0(h^2 k^2+i0)=G_0+i h k/(4\pi)+O(h^2)$. The coefficient follows from the quadratic minimum of $\epsilon(p)$. Equivalently the low-energy density of states is $\sqrt{E}/(4\pi^2)+O(E^{3/2})$. At negative energy the first correction is $-h\mu/(4\pi)$ for $z=-h^2\mu^2$. Thus $hQ/(4\pi)$ tends to $f(k)=(-1/a-ik)^{-1}$ whenever the inverse scattering length has a finite limit. At the pole this is $i/k$, independent of a scattering-length scale. This is the three-dimensional two-body unitary limit with dynamical exponent two. On the finest mesh doubling the physical relative momentum from one to two multiplies the kinetic energy by 3.999985. The exponent follows from the one-load dispersion and is independent of contact tuning. Additional particle sectors and a relativistic vacuum have separate requirements. [Universal two-body amplitude and few-body scope](https://arxiv.org/abs/cond-mat/0410417) · [Nonrelativistic conformal symmetry](https://arxiv.org/abs/0706.3746)

| Effective contact family | Physical scattering length | Limiting outcome |
| --- | --- | --- |
| Unmodified 8 | $O(h)$ | Vanishing fixed-momentum scattering |
| Exchange only | 0 | Vanishing fixed-momentum scattering |
| $\kappa_{\rm eff}^{-1}=\kappa_*^{-1}+3h^{1/2}/\pi$ | $-h^{1/2}$ | Vanishing fixed-momentum scattering |
| $\kappa_{\rm eff}^{-1}=\kappa_*^{-1}+3h/\pi$ | $-1$ | Finite momentum-dependent scattering |
| $\kappa_{\rm eff}^{-1}=\kappa_*^{-1}+3h^2/\pi$ | $-1/h$ | Unitary limit |
| Exact pole | Infinite | Unitary limit |

At $k=1$ the outgoing pole response was $0.040552+0.998027i$ at $h=1/16$ and $0.002538+0.999992i$ at $h=1/256$. Linear detuning gave $0.499998+0.498729i$ at the finest mesh against the limit $(1+i)/2$. The calculation uses the outgoing lattice Green function at real energy and the decaying Green function at negative energy independently. The packet check uses neither function. On a unit periodic $32^3$ guide at physical time $1/8$, the probability of separation below one quarter in the maximum-coordinate distance was **0.55714 at the pole**, **0.44645 with linear detuning**, **0.08710 at unmodified contrast eight** and **0.10297 with exchange only**. These are conditional finite population probabilities. All ten separation bins are retained, including the boundary bin at half the circumference. [Scattering amplitudes, packet bins and controls](cubic_resonance.json)

## Local one-dimensional guide

For the degree-two guide at the same kinetic clock scaling, put $z_j=e^{ihp_j}$ and $\Delta=1+2\kappa_{\rm eff}$. The adjacent-state equation is $2\Delta\psi(x,x+1)=\psi(x,x)+\psi(x+1,x+1)$ for the continued two-wave solution on $x<y$. The coefficient of the exchanged plane wave is

$$
S=-\frac{1+z_1z_2-2\Delta z_2}{1+z_1z_2-2\Delta z_1}.
$$

Unmodified contrast eight tends to $S=-1$ on the scattering branch. Exchange only tends to $+1$ at unequal physical momenta. A programmed repulsive detuning $\kappa_{\rm eff}=-h/4$ gives $(q-i)/(q+i)$ with $q=p_1-p_2$. This is the exchange phase of the unit-strength one-dimensional delta interaction in these kinetic units. The $-h^{1/2}/4$ and $-h^2/4$ controls tend to $-1$ and $+1$ respectively. At $q=1$ the linear-detuning phase on 256 sites was $-0.002686-0.999996i$ against $-i$. The independent ring release retains the separation distribution and the original attractive law's below-band weight. [Guide data](guide_scaling.json) · [Ring calculation](guide_scaling.py) · [Spin-chain scattering and continuum Bose gas](https://arxiv.org/abs/0705.3569)

## Fixed-range interaction on the ballistic scalar chart

Use one mobile excitation in each of two independently addressable eight-slot scalar registers, with one private site for the order-zero implementation. At a mesh $h=1/N$ let an opportunity include all pairs of slots with relative station coordinate in a specified cube, followed by two original walk periods in each register. Its declared chart duration is $2h$. The centered cube has width $R$ lattice sites per axis. Each selected contact receives native exposure $\tau=2h/8$. One independently serviced fair clock after the whole opportunity supplies the stopped records. The larger charts use an exact translation-sector reduction. At $N=32$, the half-width survey requires 8,589,934,592 logical pair contacts per opportunity before routing. This clock is a programmed cadence, with no bounded microscopic time per field period. [Protocol and resource definition](range_scaling_protocol.json) · [Opposite-momentum readout](range_scaling.py)

**Native contact.** On binary rails $a,b$ the blank-seam generator is $(8+1/2)(n_a+n_b)-17n_an_b-(b_a^\dagger b_b+b_b^\dagger b_a)/2$. A sign reversal on $a$ between two half exposures reverses the hopping. The diagonal and hopping commute on this seam. The sign word uses $a$ and 31 vacant helpers. For each nonzero five-bit XOR mask, apply a matching of native $\pi/8$ pulses. In the single-occupation bank the matching generators commute, with sum $31\kappa I+16(I-P)$ for the equal-slot projector $P$. At $\kappa=8$ their product is $-I$. The vacuum is fixed and all work rails return vacant. Full swaps conjugate each addressed pulse onto its native seam and restore intervening occupations. Inverse refocused singleton services against two extra vacant rails remove the two one-body phases. The resulting word is exactly $\exp(i17\tau n_an_b)$. The sham matches the singleton services on all inputs. A routed 20-cell check used 30,906 physical seam pulses for one coupled block. Its maximum four-column error was $3.34\times10^{-12}$. [Compiled seam word](range_contact.py) · [Microscopic and Fourier tests](test_range_scaling.py)

Let $W_h(r)$ indicate membership in the cube. The pair gate per opportunity is $\exp(i\,2h g W_h)$ with $g=17/8$. Holding transport gives joint zero-momentum return probability $1-4p(1-p)\sin^2(gt/2)$ with $p=(R/N)^3$. The source is the normalized uniform packet in each register. With transport active, translation invariance gives opposite outgoing momenta and the detector sums over every coin polarization. At deterministic chart time one the probability of a nonzero transfer was as follows.

| Width or control | $N=4$ | $N=8$ | $N=16$ | $N=32$ |
| --- | ---: | ---: | ---: | ---: |
| One site | 0.04105 | 0.00463 | 0.00054 | 0.000064 |
| $\lfloor\sqrt N\rfloor$ | 0.11229 | 0.03044 | 0.02866 | 0.00700 |
| Half the chart width | 0.23985 | 0.20605 | 0.20004 | 0.19878 |
| Full width or sham | 0 | 0 | 0 | 0 |
| Half width, exposure divided by $N$ | 0.01861 | 0.00404 | 0.00098 | 0.000245 |

At $N=32$ the half-width experiment gave probability 0.19925 of nonzero momentum before the clock alarm, with timeout probability 0.00105 on the original launch denominator. A particular opposite pair of physical wave numbers, $k_B=(2\pi,0,0)=-k_A$, retained probability 0.01130 at deterministic time one. A full-width survey supplies only a common phase at fixed occupations. The shrinking widths and the vanishing dose have free continuum evolution for smooth input. [All bins and held-transport controls](range_scaling.json)

**Continuum operator.** Write $U(k)$ for the original one-period reflection and shift on the full coin at private order zero. $U(0)^2=I$. The derivatives $v_j=(i/2)\partial_{k_j}U(k)^2|_0$ are Hermitian. On each fixed Fourier mode the doubled free step is $I-2ih\sum_j v_j k_j+O(h^2|k|^2)$. On the two-register space the limiting generator is
$$
H=\sum_j(v_j^A p_{A,j}+v_j^B p_{B,j})-gW(x_B-x_A).
$$
The kinetic operator is the self-adjoint Fourier multiplier on the torus. Trigonometric polynomials form a core, preserved as a core under bounded potential perturbation. Identify the sampled Fourier modes isometrically with their continuum modes and extend the step by identity on the omitted modes. On each core vector, sampled multiplication by the window converges in norm to $W\psi$ because its Fourier coefficients and squared norm converge by Riemann sums and its faces have zero measure. The first difference quotient therefore converges to $-iH$ on the core. Stability of unitary steps and the contraction product formula give strong convergence at bounded times. The stopped clock has chart-time variance tending to zero and its cap has vanishing timeout probability. In the half-width limit the initial derivative into a fixed opposite nonzero Fourier pair is $ig\widehat W(k)$ times the input coin. It is nonzero on each first axis mode. Translation symmetry confines this transfer to correlated momenta, so the evolution can entangle initially independent packets. The full coin is retained through collisions. This derivation permits mixing out of the one-particle equal-coin scalar sector.

The fixed width supplies a physical length and an instantaneous interaction in chart coordinates. For example a relative displacement with first component $3/16$ is inside whereas its dilation by two is outside. This precludes a scale-free local CFT interpretation of this program. On the original resonant word, homogeneous fixed-number symmetric binary states are invariant under every partial swap. At filling $p$ and $S$ slots, every two distinct occupation probes have covariance $-p(1-p)/(S-1)$ independently of distance. The three tested fillings yield no spatial density power law in that stationary family. Other states and operators require their own measurements.

## Reproduction

Run from the RER root with the model dependencies. Each producer writes a separate result file.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m model.experiments.seam_program.cubic_resonance --output /tmp/echo-cubic-resonance.json
python -m model.experiments.seam_program.guide_scaling --output /tmp/echo-guide-scaling.json
python -m model.experiments.seam_program.range_scaling --output /tmp/echo-range-scaling.json
python -m pytest -q -W error model/experiments/seam_program/test_kinetic_resonance.py model/experiments/seam_program/test_range_scaling.py
```
