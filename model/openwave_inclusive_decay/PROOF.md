# Inclusive charged-current emission and finite survival records

## Declared source and exact invariant sector

The internal coefficient is the actual fifteen-component weak raising current
from `model/common_source/calibrated_matter.py`. Its three generators are
`(J+J†)/2`, `-i(J-J†)/2`, and `(JJ†-J†J)/2`. The charged masses are the complete
native analyzer's three positive effects, in the middle-branch unit. The light
mass is `m=e_1/e_2`, enclosed by the rational intervals in the producer; the
coarser enclosure `3/625 < m < 49/10000` suffices for the proof. The phase `2/9`
and family labels are architecture choices. Two distinct massless neutral
species accompany the charged transition. The spectator native fields have
their free dynamics and their declared filled-sea reference.

This member specifies continuum relativistic quasiparticle dispersion, the
positive-frequency middle-to-light charged-current vertex and its Hermitian
reverse. Its smooth form factor is a supplied ultraviolet choice. It is not
identified with the unfiltered three-momentum contact Hamiltonian of
`model/openwave_weak_decay`, whose scattering and vacuum-response words have
different action. No rotating-wave error bound for that unfiltered source is
asserted. No theorem for the separate massive-neutral family architecture is
transferred to these massless neutral modes.

In a finite momentum box the interaction has the complete form

`V = G sum_channel (g_channel b_e† c_nu† d_antinu† a_parent + h.c.)`.

All allowed momenta, daughter spins and neutral helicities occur. The parent
has two spins. Each channel satisfies total momentum conservation. The initial
state is one rest parent and daughter vacuum relative to the free sea, with an
idle opposite-charge reference. Spectators are retained. Three exact conserved
number combinations are

`N_parent+N_e = N_parent+N_nu = N_parent+N_antinu = 1`.

Consequently the admitted full Fock sector consists of one parent or three
daughters. The reverse vertex returns three daughters to the parent. This is
an invariant sector proved from the full quartic source, not a selected
two-level approximation. States orthogonal to the current's angular spectral
vectors are retained dark states. An arbitrary external quantum reference can
be tensored to the evolution and the complete record channel.

## Native coefficient, spin and helicity normalization

Writing `T_a` for the three generators, the two ordered cross-family current
products give

`2 sum_a (T_a)_(e,nu) (T_a)_(nu,e) = 1`.

The native raising sign cancels. Thus the left-Weyl amplitude coefficient is
exactly `G`, with no inserted Fermi constant:

`M = G (u_e† bar_sigma^mu v_antinu) (u_nu† bar_sigma_mu u_parent)`.

Use `bar_sigma=(I,-sigma)`, metric `(+---)`, and spinor normalization
`u†u=2E`. The left-component spin sum is
`D(p)=p^0 I-p dot sigma`. For each massless neutral the nonzero left component
has one helicity; its rank-one `D(p)` sums that helicity once. The electron spin
sum is `D(p_e)`. The unpolarized rest-parent density is `D(P)/2`. Therefore the
complete spin sum and parent average are

`sum_(mu,nu) eta_mu eta_nu tr[D_e bar_sigma_mu D_anti bar_sigma_nu]
 tr[D_nu bar_sigma_mu (D_P/2) bar_sigma_nu]`

`= 8 (P dot p_anti) (p_e dot p_nu)`.

This identity follows by expanding the four Pauli matrices, equivalently from
`tr(sigma_i sigma_j)=2 delta_ij` and
`sigma_i sigma_j=delta_ij I+i epsilon_ijk sigma_k`. The producer independently
evaluates all sixteen trace contractions at nine distinct three-body
kinematic points and verifies the invariant expression. Omitting the parent
spin average doubles every trace. The two neutrals have distinct family and
particle/antiparticle labels; there is no identical-particle `1/2!`. The
electron is distinct from both neutrals.

Let `s=(p_e+p_nu)^2`, parent mass `M`, and light mass `m`. Then the spin-averaged
amplitude squared is

`|M|^2 = 2 G^2 (M^2-s)(s-m^2)`.

## Box, Fourier and phase-space factors

For box volume `V`, normalized creation operators obey the discrete canonical
anticommutation relation. A free field contributes
`V^(-1/2) u(p)/sqrt(2E_p)` per mode. The spatial integral of the four-field
vertex gives a Kronecker momentum delta and a factor `V`. Thus its finite-box
matrix element is

`G A_spin delta_(P,sum p) / [V sqrt(2M 2E_e 2E_nu 2E_anti)]`.

The two independent momentum sums have continuum density
`V^2/(2pi)^6`. Equivalently, with continuum operators normalized by
`{a(p),a†(p')}=(2pi)^3 delta^3(p-p')`, each final particle contributes
`d^3p/[(2pi)^3 2E]`. Energy and momentum conservation give

`dPhi_3=(2pi)^4 delta^4(P-sum p) product_j d^3p_j/[(2pi)^3 2E_j]`.

Integrating the two-body radial delta in its rest frame gives
`Phi_2(S;a,b)=sqrt(lambda(S,a,b))/(8pi S)`: the radial Jacobian is
`p E_1 E_2/(sqrt(S) E_1 E_2)=p/sqrt(S)`, with the angular integral `4pi` and
the displayed Fourier factors. Factoring the three-body measure at invariant
mass `s` gives `ds/(2pi) Phi_2(M^2;0,s) Phi_2(s;m^2,0)`. Hence

`dPhi_3/ds=(M^2-s)(s-m^2)/(128pi^3 M^2 s)`.

The rate's parent flux factor is `1/(2M)`, so

`Gamma=G^2/(128pi^3 M^3) int_(m^2)^(M^2)
 (M^2-s)^2(s-m^2)^2/s ds`

`=G^2 M^5 f(m/M)/(1536pi^3)`,

where

`f(r)=1-8r^2+8r^6-r^8-12r^4 log(r^2)`.

The coefficient `1536=128*12` is reconstructed by the native current factor,
spin trace, two two-body integrals and the exact final integral. The producer
symbolically verifies `f(r)=12 int_(r^2)^1 (1-z)^2(z-r^2)^2/z dz`. In particular,
the massless integral is `1/12`. These factors are not imported from a measured
lifetime or a tabulated decay formula. `G` has dimension energy to power minus
two; the rate has dimension energy. The box Riemann sums have this continuum
normalization. A quantitative approximation by that physical momentum box is
not part of the finite energy-bin theorem below.

## Complete reservoir and finite observation window

Set `M=1`. The chosen ultraviolet amplitude is `u(E)=(2-E)^2` on `m<=E<=2`,
zero elsewhere. It equals one at the parent energy. For a daughter total energy
`E`, the parent spin density contributes `M/E` times that of a rest spinor of
mass `E`. The finite-box external-parent normalization `1/(2M)` cancels this
factor to `1/(2E)`. The spectral delta has one fewer factor `2pi` than the
on-shell rate. The resulting spectral density without `G^2` is therefore

`rho(E)=E^5 f(m/E) (2-E)^4/(3072pi^4)` on `[m,2]`.

Rotational covariance makes the angular spectral Gram matrix proportional to
the identity on parent spin; its scalar is this density. Spectral coarea gives
two identical parent-reservoir blocks. It is an isometry for the admitted
states, with the entire angular dark complement retained. Every daughter
configuration and its correlated spins occur in the native vertex; inclusive
emission sums their complete probabilities.

In either block the Hamiltonian is a bounded selfadjoint Friedrichs operator
on `C plus L2([m,2],dE)`:

`H(a,b)=(a+G int sqrt(rho) b, E b+G sqrt(rho) a)`.

The omitted notation denotes its orthogonal angular dark factors and free
spectators, not discarded outcomes. Initial parent state `(1,0)` has energy
expectation one. In the rotating frame its survival amplitude satisfies the
exact Volterra equation

`a'(t)=-G^2 int_0^t C(s)a(t-s) ds`,
`C(s)=int rho(E) exp[-i(E-1)s] dE`.

The model is unitary, so `|a|<=1`. Suppose
`int_0^infty |C| <= B0` and `int_0^infty s|C| <= B1`. Then
`|a'|<=G^2 B0`. Put
`alpha=int_0^infty C=pi rho(1)-i PV int rho(E)/(E-1) dE`.
The latter identity follows by inserting an exponential regulator and taking
its limit; smooth compact support and the integrable correlation justify it.
Subtracting `-G^2 alpha a(t)` from the exact equation gives an error bounded by
`G^4 B0 int_0^t s|C(s)| ds + G^2 int_t^infty |C(s)| ds`.
Integration with the contractive scalar semigroup yields

`|a(t)-exp(-G^2 alpha t)| <= G^2 B1 (1+G^2 B0 t)`.

Thus the parent survival probability differs from `exp(-Gamma t)` by at most
twice this bound, where `Gamma=2pi G^2 rho(1)`. This is a finite-time comparison,
not exact exponential decay for all times.

## Rational correlation certificate

Write `rho=S q/(3072pi^4)`, `q=(2-E)^4`, and

`S=E^5-8m^2 E^3+8m^6/E-m^8/E^3-12m^4 E log(m^2/E^2)`.

The threshold integral shows `S` vanishes to fifth order at `E=m`; the cutoff
vanishes to fourth order at `E=2`. The first two integration-by-parts boundary
terms are zero, and the third derivative is integrable. Since `0<=f<=1`,
`int rho <= [int_0^2 E^5(2-E)^4 dE]/(3072*3^4) = (256/315)/(3072*3^4)`.
For `m_hi=49/10000`, `m_lo=3/625` and `log(2/m)<7`, the following bounds on
`int|S^(j)|` follow by termwise differentiation and integration:

`S0=32/3`,
`S1=32+64m_hi^2+9m_hi^5+384m_hi^4`,
`S2=80+96m_hi^2+179m_hi^4`,
`S3=160+96m_hi^2+52m_hi^3`.

The logarithmic bound is certified by a rational positive Taylor partial sum
for `exp(7)>2/m_lo`. For example, the logarithm and constant in `S'` contribute
at most `384m_hi^4`; the reciprocal terms contribute `9m_hi^5`. The cutoff
derivative suprema are `16,32,48,48`. With `D=3072*3^4`, define

`A0=(256/315)/D`,
`A2=(16S2+64S1+48S0)/D`,
`A3=(16S3+96S2+144S1+48S0)/D`.

Fourier integration by parts gives
`|C(t)|<=min(A0,A2/t^2,A3/t^3)`. Splitting the first integral at
`sqrt(A2/A0)` and the weighted integral at `(A3/A0)^(1/3)` gives

`int|C| <= 2sqrt(A0 A2) < 1/2000`,
`int t|C| <= (3/2) A0^(1/3) A3^(2/3) < 1/250`.

The producer checks the equivalent squared and cubed inequalities using exact
rational arithmetic. It also checks
`Gamma/G^2 > (1-8m_hi^2-m_hi^8)/(1536*(22/7)^3) > 1/50000`.
At `G=1/64`, time `t=k/Gamma` has amplitude error at most
`(1+25k)/(4096*250)`. The one-lifetime survival is `exp(-1)` within
`0.00005078125`; the four-lifetime survival is `exp(-4)` within
`0.000197265625`. The complementary inclusive emission probabilities have the
same error bounds.

## Finite source, complete records and work

Round each reservoir energy to its containing bin's midpoint, leaving the
coupling function and all angular/spin vectors unchanged. With maximum bin
width `h`, the full source operator difference is at most `h/2`. Duhamel's
identity gives a unitary error at most `t h/2`, including arbitrary external
references. Each bin's coupled vector has weight `int_bin rho`. Those bright
vectors span a finite invariant parent-plus-all-bins sector of the rounded
source. Every orthogonal bin/angular mode is retained and has zero support
for the admitted parent preparation. Inclusive and energy-bin records on
this source have survival probability error at most `t h` relative to the
continuum comparison.

The same finite bright-sector Hamiltonian has an exact full fermionic
realization: one parent mode and three daughter modes per bin, with positive
individual daughter energies summing to the bin energy and quartic vertices
of strength `G sqrt(int_bin rho)`. Its full Fock complement is retained.
These are orthogonal encoded three-daughter channels. This realization proves
the finite inclusive record law, not a finite-box reconstruction of every
individual continuum daughter momentum. Angularly resolved physical-box
observables require their own discretization. No two-state surrogate is used.

For `h=2^-60`, the rational time bound gives a finite-bin probability error
below `7.11e-10` through four lifetimes. The required number of bins is about
`2.30e18`. This finite-carrier result is analytical. The producer executes
32-, 64- and 128-bin source matrices, records every emission-bin probability
and the parent complement, and verifies full energy and normalization. These
small matrices are diagnostics; their large-time recurrence is not an
irreversible decay law.

The complete measurement is parent survival versus all daughter energy-bin
outcomes. It destroys parent/daughter interaction coherence. Its mean energy
work is

`W=E_free_after-1=-2G Re[a* int sqrt(rho)b]`,

with `|W|<=G sqrt(int rho)`, since `2|a| ||b||<=1`. Full autonomous energy is
conserved before readout. The executable retains this measurement work and
every no-emission outcome. A clock or detector implementing this service is
part of the declared observer architecture.

## Physical interpretation

The native ratio, current and phase-space calculation determine this member's
dimensionless rate once its contact coefficient and cutoff are declared.
Electron anchoring converts the compiled middle mass to about
`0.1056594144959 GeV`; multiplying the dimensionless rate by that number gives
a conditional width in GeV. Under this conversion the contact coefficient is
`(1/64)/M_GeV^2`, not an empirically fitted weak coefficient. No observed muon
lifetime, weak-scale match, massive-neutrino identification, or interacting
vacuum is claimed. The exact finite spectral theorem, invariant-sector proof
and finite unitary comparison replace any assumption of irreversible decay
on a finite source.
