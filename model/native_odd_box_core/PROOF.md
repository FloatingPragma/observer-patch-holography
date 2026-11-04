# The unmodified odd-box charged core

The action is the centered source clock on the open cube
`{0,...,L-1}^3`, with two fermion modes per site:

\[
 H=d\Gamma(h)+W,\qquad h_{x,x+e_j}=ic\sigma_j/2,
 \qquad W=U\sum_x(n_{x\uparrow}-1/2)(n_{x\downarrow}-1/2).
 \tag{1}
\]

Here `L` is odd, `L>=3`, and `c,U>0`. The clock is the declared
`H_uncentered-U N/2+UV/4`; the scalar does not change dynamics. Thus the
interaction, hopping and centered number clock are the source laws. No mass
term, order parameter or auxiliary defect background occurs in (1). The
`FieldGeometry.single_h` entries and their separation below are checked
exactly at the retained 27-site parameters. The centered contact convention
is bound to `centered_source_ground`. Every number sector is retained.

## Exact finite quadratic spectrum

On one coordinate let
`K=(ic/2)(upper_shift-lower_shift)`. A diagonal phase converts it to the
real path hopping matrix. Its eigenvalues, with multiplicity, are
`c cos(pi k/(L+1))`, `k=1,...,L`; reversing their order removes a possible
overall sign. Its kernel is one dimensional, with normalized vector

\[
 z_j=\sqrt{2/(L+1)}\,1_{j\ {even}}.
 \tag{2}
\]

The equal nonzero entries have the same sign. In three dimensions
`h=sum_j K_j tensor sigma_j`, where the three `K_j` commute. Pauli
anticommutation gives `h^2=sum_j K_j^2 tensor I`. For each simultaneous
spatial eigenvector the spin block is `k dot sigma`, with eigenvalues
`+/-|k|`. Consequently the only zero space is
`phi tensor C^2`, `phi=z tensor z tensor z`, of dimension two. All other
one-body energies have absolute value at least

\[
 \Delta_L=c\sin(\pi/(L+1)),\qquad
 I_L=\sum_x|\phi_x|^4=(2/(L+1))^3.
 \tag{3}
\]

There are `V-1` negative modes and `V-1` positive modes, `V=L^3`.
After subtracting the filled negative-sea energy, the complete quadratic
Fock Hamiltonian has a four-dimensional kernel `P`, with particle numbers
`V-1,V,V,V+1`, and satisfies `H0 >= Delta_L (1-P)`. Excitations can empty
any negative mode and fill any positive mode, so this inequality holds on
the entire `2^(2V)` carrier. The corresponding one-dimensional chain has
the same conclusions with `V=L` and `I_L=2/(L+1)`.

The executable bounds use `Delta_lower=2c/(L+1)`, following
`sin t>=2t/pi` on `[0,pi/2]`. For `L=3` they use the sharper rational
`Delta_lower=2c/3 < c/sqrt(2)`. The proof applies to every finite odd
length; the executable geometry API limits explicit inputs to the owning
source's 4096-site domain.

## Filled sea and the original interaction

Spatial inversion `S:x -> (L-1)-x` anticommutes with `h` and fixes `phi`.
If `R_-` and `R_+` are the negative and positive spectral projections,
then `S R_- S=R_+` and `R_-+R_+=I-|phi><phi| tensor I_2`. The weight
`rho_x=|phi_x|^2` is inversion invariant. Taking diagonal spin blocks and
summing with this weight gives the exact matrix identity

\[
 \sum_x\rho_x R_-(x,x)=(1-I_L)I_2/2.
 \tag{4}
\]

Write `f_1,f_2` for the two zero-mode annihilators and `n_a=f_a^* f_a`.
Wick expansion in the filled sea, including its exchange term, gives

\[
 PWP=C+U I_L\left(n_1n_2-\frac{n_1+n_2}{2}\right),\qquad
 C=U\sum_x\left[\det R_-(x,x)-\frac{\operatorname{tr}R_-(x,x)}2+\frac14\right].
 \tag{5}
\]

For example, the linear zero-mode coefficient from the local contact is
`rho_x f^*[(tr R_-(x,x)-1/2)I-R_-(x,x)]f`. Equation (4) sums it to
`-I_L(n_1+n_2)/2`. The quartic term sums to `I_L n_1 n_2`. Keeping only
the quartic zero-mode term would give the wrong charged-level ordering.

The four eigenvalues in (5) are `C-g/2,C-g/2,C,C`, where `g=U I_L`.
Since each centered local term has eigenvalues `+/-U/4`, and the empty
configuration achieves the simultaneous upper endpoint,

\[
 w=\|W\|=UV/4.
 \tag{6}
\]

## Complete-Fock interacting gap and actual reference

Assume `w<Delta_lower/4` and put

\[
 \epsilon=\frac{w^2}{\Delta_{\rm lower}-2w},\qquad
 \delta=\frac{w}{\Delta_{\rm lower}-2w}.
 \tag{7}
\]

Min-max places exactly four low levels in `[-w,w]`, and every other
level above `Delta_lower-w`. For a low eigenvalue `lambda`, the full
complement resolvent exists and the Feshbach matrix is

\[
 PWP-PWQ\bigl(Q(H_0+W)Q-\lambda\bigr)^{-1}QWP.
 \tag{8}
\]

The second term is positive with norm at most `epsilon`. Inertia, or the
monotonic eigenvalue branches of (8), puts each ordered low energy within
`epsilon` of the corresponding ordered eigenvalue of `PWP` (in fact on
its lower side). The same argument applies separately in each conserved
number sector, so sector labels are retained. No virtual positive or
negative mode is removed from the resolvent.

The physical spin antiunitary `Theta=i sigma_2 conjugation` commutes with
the original hopping: complex conjugation reverses `i`, while spin time
reversal reverses each Pauli matrix. It also preserves `W`. Its Fock lift
squares to `(-1)^N`, hence the two low `N=V` levels are an exact Kramers
pair because `V` is odd. For particle-hole symmetry take the canonical
CAR automorphism

\[
 c_x\longmapsto(-1)^{x_1+x_2+x_3}i\sigma_2 c_x^*.
 \tag{9}
\]

The staggered sign reverses each bond; together with spin reversal this
gives `J h^T J^*=-h` in the particle-hole quadratic identity. The
zero-trace hopping and centered contact are invariant. The automorphism
maps `N` to `2V-N`, so the `N=V-1` and `N=V+1` low levels agree exactly.

If `g/2>2 epsilon`, the actual full-source ground is the `N=V` Kramers
doublet, and its charge-conjugate charged excitation gap `m_ch` obeys

\[
 0<g/2-2\epsilon\le m_{\rm ch}\le g/2+2\epsilon.
 \tag{10}
\]

The associated rank-two ground projection `P_G` differs from the filled
sea with one zero fermion projection `P_1` in operator norm by at most
`delta`. Indeed the `Q` component of every ground vector follows from
the resolvent in (8), and the equal-rank principal-angle identity converts
`||Q P_G||<=delta` to `||P_G-P_1||<=delta`. The analogous bound holds
for the rank-one low charged projection in each of its number sectors.
Thus the actual source reference is the normalized exact ground mixture
`P_G/2`; a preparation of `P_1/2` has unhalved trace distance at most
`2 delta` from it. This last estimate does not assert an executed sea
preparation or an exact measurement of interacting spectral projections.

At `L=3,c=1,U=1/65536`, `I_L=1/8`, `g/2=U/16`, and (10) is strictly
positive by rational arithmetic. The complete 54-mode problem is bounded
analytically; the complete 64-dimensional Fock chain is an independent
executed comparison, including its filled sea and all number sectors.

## Odd-length scaling and boundary meaning

An explicit controlled diagonal family is

\[
 U_L=\frac{c}{2^{16}(L+1)^{10}}.
 \tag{11}
\]

Using `Delta_lower>=2c/(L+1)` gives `2w<=Delta_lower/2` and
`epsilon<=U_L^2 L^6(L+1)/(16c)`. Therefore
`2 epsilon/(g/2)<=U_L L^6(L+1)^4/(32c)<=2^-21`. The family has a
positive gap in each finite box and `g/2=4c/[2^16(L+1)^13]`, tending to
zero. The small-coupling condition is volume dependent.

The normalized zero orbital has uniform weight on one parity sublattice,
with `IPR=8/(L+1)^3`. At `L=3` it has probability `1/8` at every corner
and zero at the central site. Its existence uses the odd open boundary.
It is a finite standing wave; this argument provides no topological
localization, nonzero thermodynamic charged mass, or relativistic
interacting continuum spectrum.

## Complete native occupation readout

For `L=3` apply on each line the real one-body rotations
`G_01(-pi/2)` followed by `G_12(pi/4)`, with
`G(theta)=[[cos theta,sin theta],[-sin theta,cos theta]]`. Their product
sends `(1,0,1)/sqrt(2)` to the central basis vector. Applying this to both
spin components along all three axes gives 108 fermionic Givens rotations.
On the entire CAR carrier, for `p<q`,

\[
 a_p^*a_q-a_q^*a_p=
 \frac i2(X_pZ_{p+1}\cdots Z_{q-1}Y_q-
          Y_pZ_{p+1}\cdots Z_{q-1}X_q).
 \tag{12}
\]

The two Pauli strings commute. Their exact native compilations therefore
implement the full Gaussian unitary, including double occupations and
all spectator parity signs. Each distance-`d` Givens costs `10+12d`
primitive pulses. The three axial mode distances are `18,6,2`, with 36
Givens each, so one decode costs 12,312 pulses. Copy the two central
occupation bits to separate blank memory bits by two native CNOTs, then
undo the Gaussian decoder. This costs 24,630 pulses and implements the
four-outcome instrument with projectors
`D^* Pi_(n_center_up,n_center_down) D`. All four occupation outcomes,
every native complement failure and the independent reference are retained.

The original six-level primitive identities and scalar phases are used.
The exact rational-port interface bounds a `p`-bit command's complete
isometry error by `29*2^-p`. The complete unhalved trace error for `N`
pulses is at most `58N*2^-p`, for arbitrary reference entanglement.
The producer emits every distinct dyadic command and its certified source
interval at a precision giving cube error below `10^-9`. These exact
commands are a separate interface from legacy binary64 controls.

The chain uses 278 pulses on all 64 field rows with 64 reference columns
and two blank pointers. Its full six-level native action is executed,
including every first-failure weight; the small binary64 error has an
independent exact source/angle certificate. The cube compilation uses
exact local identities and a 54-dimensional orbital check, without
executing a `2^54` native matrix. Scheduled addressed pulses, blank memory,
and the declared integrated control precision are apparatus inputs. The
theorem does not protect this detector against simultaneous Hamiltonian
drift or interleaved C/E repair; no quiet interval or successful outcome
is silently conditioned upon. Native observation of the Hamiltonian pole
itself requires the separate response/readout interface.
