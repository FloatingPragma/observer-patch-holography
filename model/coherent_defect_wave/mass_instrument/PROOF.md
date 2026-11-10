# Complete particle and hole mass records on the graded flux clock

The source is the [graded flux construction](../many_body_route/PROOF.md),
with its declared Wilson Hamiltonian, Gauss encoding, fermionic grading and
filled negative-band reference. No spectrum is transferred from the ordered
routing tick or its native interaction. The rest parameter is the literal
paired-flux face count, with a chosen unit coefficient in this Hamiltonian.

## The same state, clock and momentum

On a finite periodic cube with spacing a, physical momenta are
`p=2*pi*n/(L*a)` in the chosen Brillouin zone. The four-component symbol,
identical for each of the five face channels, is

\[
 h_a(p)=\sum_i\alpha_i\frac{\sin(ap_i)}a+
 \beta\left[2+\frac1a\sum_i(1-\cos(ap_i))\right],\qquad
 E_a(p)^2=\sum_i\frac{\sin^2(ap_i)}{a^2}+
 \left[2+\frac1a\sum_i(1-\cos(ap_i))\right]^2.
\]

The Clifford identities give `h_a^2=E_a^2 I`. Both energy signs have
multiplicity two. Write `P_+/-=(I+/-h_a/E_a)/2` and fill the negative
band. For one momentum and face block, the full Fock space has sixteen
states and the reference has two occupied modes. Its pure density is the
second exterior power of P_minus. The code constructs its six-by-six
two-occupation block from the actual projector minors, embeds it in all
sixteen states and verifies its trace, purity, number and zero energy.

The ground-relative Hamiltonian is `H=dGamma(h_a)+2E_a I`. Its entire
spectrum is `0,E_a,2E_a,3E_a,4E_a`. Its spectral projectors are computed as
polynomials in that same H and checked exactly. The other momentum/face
blocks remain in their filled reference. The global ground contains `10V`
source pairs. Relative charge means change in that count, not an assumed
electromagnetic charge. A hole removes a pair from an occupied reference
mode and has positive excitation energy. It is not a negative-energy
particle in the empty routing reference.

The source encoding intertwines this declared Fock Hamiltonian with its
encoded link program. Position and momentum refer to its occupied source
slots and translations. A particle created at p carries relative momentum
p; removing an occupied p mode gives a hole of relative momentum -p.
`E_a(-p)=E_a(p)`. These are free particle/hole records on this specified
source. The physical charge and Standard Model identification need their
own current and interaction attachment.

## An apparatus with every outcome retained

Let v be normalized and `c=c(v)` its CAR annihilator. Set
`w=<v,P_minus v>`. An extra fermionic reservoir mode is prepared occupied
for addition and empty for removal. A genuine fermionic mode SWAP exchanges
it with c. It conserves their total number, including its exchange signs.
The implementation checks all thirty-two basis states and all five CAR
mode transforms. It retains the terminal reservoir flag.

For the specified even-parity reference, the addition outcomes have reduced
source states `c* rho c` and `n rho n`; removal has `c rho c*` and
`(I-n)rho(I-n)`, where `n=c*c`. An energy measurement follows the reservoir
readout. A failed addition or removal need not leave the reference unchanged:

\[
 n\Omega=w\Omega+c^*(P_+v)c(P_-v)\Omega.
\]

The second term is orthogonal to the ground, has energy `2E_a`, relative
charge zero, and squared norm `w(1-w)`. Consequently the complete nonzero
conditional records are

| Setting | Flag | Relative charge | Energy | Probability |
|---|---|---:|---:|---:|
| Add | Success | +1 | E_a | 1-w |
| Add | Failure | 0 | 0 | w squared |
| Add | Failure | 0 | 2E_a | w(1-w) |
| Remove | Success | -1 | E_a | w |
| Remove | Failure | 0 | 0 | (1-w) squared |
| Remove | Failure | 0 | 2E_a | w(1-w) |

The setting is selected with recorded probability one half. The receipt
retains all thirty-six combinations of two settings, two reservoir flags
and nine joint energy/charge effects. Those nine effects span the entire
sixteen-state block: they are the labels `energy=(N_particle+N_hole)E_a`,
`charge=N_particle-N_hole`, with both occupation counts in `{0,1,2}`.
The other sixteen nominal energy/charge products are identically zero
operators. All thirty zero-probability records in the prepared experiment
are retained. This eliminates redundant labels without postselection. The
probabilities are computed from the full thirty-two-state SWAP output and
the source energy/number projectors. They sum to one without conditioning
on a particle-creation success. At p=0, v equal to the first literal spin
basis vector has w=1/2. The setting-weighted probabilities are
`1/4,1/8,1/8,1/4,1/8,1/8`; the detected energies are `2` and `4`, in source
units, as well as the zero-energy failure outcome.

The nonzero-momentum execution uses `a=1`, `p=(0,0,pi/2)` on a periodic
side-four source, giving `E_a=sqrt(10)`. Every populated mode is measured
with the same reference and Hamiltonian convention. An ideal spectral
projector is an explicitly supplied finite detector, not an assertion that
a finite Fourier bin implements an exact projector.

## A finite energy-conserving apparatus completion

The driven SWAP can transfer energy to the source. Its work is included by
a finite nine-level battery with energies `b E_a`, b=0,...,8. The reservoir
has zero energy in this apparatus convention. Let V be the thirty-two-state
SWAP and P_j the source's energy-jE_a projector, extended over the reservoir.
For each n=4,...,8 define an isometry into the 288-dimensional source,
reservoir and battery space by

\[
 F_n\psi=\sum_{j=0}^4 P_j\psi\otimes|n-j\rangle.
\]

The verified projector identities imply `F_n*F_m=delta_nm I`. Each range
is the complete total-energy nE_a shell. Thus

\[
 U=\sum_{n=4}^8 F_n V F_n^*+
       I-\sum_{n=4}^8 F_nF_n^*
\]

is a finite unitary commuting with total energy. Partial edge shells are
retained with the identity action. Preparing the ground, reservoir setting
and battery level four uses the complete n=4 shell. An energy-jE_a result
leaves the battery at level `4-j`. The reservoir/energy/charge probabilities
above are exactly unchanged. The battery records energy transfer; tracing
it out can remove coherence between different energies, so equivalence is
claimed for this energy-resolved observation grammar, not arbitrary later
interference. The code instantiates this isometry theorem using the actual
projectors and SWAP; it does not materialize a 288-by-288 matrix.

## Mass readout, cutoff and units

The exact zero-momentum particle and hole energies are two at every volume
and spacing. The lowest neutral pair costs four. Away from zero, the
observable dispersion is E_a, including its declared regulator term. Its
algebraic inversion is

\[
 \sqrt{E_a(p)^2-\sum_i\sin^2(ap_i)/a^2}
       -\sum_i(1-\cos(ap_i))/a=2.
\]

This identity uses the specified Wilson action. It is not a correction
selected to improve a measured mass. On `|p_i|<=P`, the owning source proof
gives

\[
 |E_a(p)-\sqrt{|p|^2+4}|\le b(a,P)
 =a^2P^3/2+3aP^2/2.
\]

The same bound controls a finite-time one-particle evolution by `|t|b`.
Spectral lines on a finite torus are exact; replacing its discrete momentum
grid by continuum wave packets needs the stated volume/refinement limit.
A stable free pole has zero width in this architecture. No nonzero
interaction correction or physical decay width is inferred from these free
records. Its five equal face channels are not five observed species.

If the energy of one source unit is `u_E` GeV, the particle rest energy is
`2 u_E` GeV, or mass `2 u_E` GeV/c^2. The source calculation does not select
u_E. Assigning a measured particle's mass to fix u_E is calibration of that
comparison. Dimensionless mass ratios and independent current/scattering
properties remain necessary tests. No measured mass or empirical scale is
used in this packet.
