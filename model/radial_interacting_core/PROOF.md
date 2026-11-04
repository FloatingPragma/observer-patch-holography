# Interacting charged states of the bounded radial source

On the symmetric open cube `{-R,...,R}^3`, with `V=(2R+1)^3` sites,
retain both spin modes and every Fock state. The one-body action is

\[
 h_{x,x+e_j}=\frac{iv}{2a}\sigma_j,\qquad
 B_x=m\sum_j(-1)^{x_j}\frac{x_j}{|x|}\sigma_j,\quad B_0=0.
 \tag{1}
\]

This is the same bounded radial profile as
[the exact radial gap theorem](../radial_defect_gap/PROOF.md).
The interacting clock is

\[
 H=d\Gamma(h)+W,\qquad
 W=U\sum_x(n_{x\uparrow}-1/2)(n_{x\downarrow}-1/2),\quad U>0.
 \tag{2}
\]

It is the original centered contact on the same carrier, with the supplied
bounded quadratic profile (1). Its autonomous preparation or selection is
not a consequence of this theorem. No finite linear mass texture is used.

## The exact zero frame and its source symmetries

Choose a rational `k>=sqrt(3)v/a` with `m>k`, and put `d=m-k`.
The parent theorem gives exactly two inversion-even zero modes and no
other one-body spectrum in `(-d,d)`, for this finite cube and the infinite
lattice. Its finite normalized origin frame is

\[
 F=\begin{pmatrix}I_2\\-A^{-1}V_0\end{pmatrix},\qquad
 A=QhQ,\quad V_0=QhP_{\rm origin},\qquad
 W_0^{\rm frame}=F(F^*F)^{-1/2}.
 \tag{3}
\]

The origin complement inverse obeys `||A^-1||<=1/d`, and
`V_0^*V_0=3v^2 I_2/(2a^2)`. Write `W_x` for the two by two block of
the normalized frame at site x. The symbol `W_x` in this section denotes
an orbital frame, rather than the contact operator W in (2).

Let `S_x=(-1)^(x_1+x_2+x_3)`, `J=i sigma_2`, and let `R_inv` be spatial
inversion. Complex conjugation is K. Plain spin time reversal `J K`
preserves `h0` and reverses B. The staggered S reverses `h0` and preserves
B; inversion reverses both. Hence

\[
 (SJK)h(SJK)^{-1}=-h,\qquad
 \Theta=R_{\rm inv}SJK,\quad \Theta h\Theta^{-1}=h,
 \quad\Theta^2=-I.
 \tag{4}
\]

The real spatial operators commute with each other and J, which proves
the square in (4). Neither plain `J K` nor `R_inv J K` is a commuting
time reversal of this radial source. Its Fock lift fixes the vacuum and
squares to `(-1)^N`, so every odd-number sector has Kramers degeneracy.
The centered onsite contact is invariant under these spin rotations,
site inversion and staggered phases.

The canonical finite-CAR particle-hole map

\[
 c_x\longmapsto S_xJ c_x^*
 \tag{5}
\]

also preserves (2). The matrix condition is `S J h^T J^* S=-h`; the
hopping has zero trace, the onsite mass is traceless, and the centered
contact is invariant under particle-hole conjugation and spin rotation.
This maps N to `2V-N`. Equation (5) is derived for the actual radial
profile; it is not imported from a different mass source.

Inversion acts as identity on the zero space. Its frame in (3) has origin
block equal to a positive scalar times `I_2`: the positive two by two
normalization matrix commutes with the origin antiunitary `J K`, hence
is scalar. Since `Theta` preserves the zero space and its origin basis,

\[
 W_{-x}=W_x,\qquad S_xJ\overline{W_x}=W_xJ.
 \tag{6}
\]

For either sign of S, the second identity gives
`W_x=[[z,-S conjugate(w)],[w,S conjugate(z)]]`.
Consequently

\[
 W_x^*W_x=W_xW_x^*=\rho_xI_2,\quad
 \sum_x\rho_x=1,\quad W_0=\sqrt\alpha I_2,\quad
 \alpha\ge\alpha_*=\frac1{1+3v^2/(2a^2d^2)}.
 \tag{7}
\]

In particular the effective contact coefficient satisfies
`I=sum_x rho_x^2`, `alpha_*^2<=I<=1`. The spin frame need not be a scalar
spatial orbital. Equation (6), checked independently over exact algebra,
is the needed replacement for that assumption.

## The filled sea and the original contact

The one-body nonzero spectrum occurs in plus/minus pairs by inversion.
Fill its `V-1` negative modes. Let `P` be the resulting four-dimensional
Fock space with zero-mode occupations `00,10,01,11`. Its particle numbers
are `V-1,V,V,V+1`. Subtracting the fixed filled-sea energy gives
`H0>=d(1-P)` on the complete `2^(2V)` Fock carrier.

If `R_-` is the negative one-body projection, inversion gives

\[
 R_-(-x,-x)=I_2-W_xW_x^*-R_-(x,x).
 \tag{8}
\]

Pairing x and -x, including the half-weight identity at the fixed origin,
and using (7), yields

\[
 \sum_xW_x^*R_-(x,x)W_x=(1-I)I_2/2.
 \tag{9}
\]

For zero annihilators `f_1,f_2`, Wick contraction of a local centered
contact has linear term
`f^*[(tr R_-(x,x)-1/2)W_x^*W_x-W_x^*R_-(x,x)W_x]f` and quartic
term `|det W_x|^2 n_1 n_2`. Since `|det W_x|^2=rho_x^2`, (9) proves

\[
 PWP=C+UI\left(n_1n_2-\frac{n_1+n_2}{2}\right).
 \tag{10}
\]

Here `C=U sum_x[det R_-(x,x)-tr R_-(x,x)/2+1/4]`. The filled-sea
exchange and centering terms are retained. The two singly occupied core
energies are `C-UI/2`, and the empty/full core energies are C.

## Controlled full-Fock interaction and spectral poles

Exactly `||W||=w=UV/4`. Assume `w<d/4` and define

\[
 \epsilon=\frac{w^2}{d-2w},\qquad t=\frac{w}{d-2w},\qquad
 g_*=\frac{U\alpha_*^2}{2}-2\epsilon>0.
 \tag{11}
\]

The full complementary Feshbach inverse at every low energy
`lambda in [-w,w]` has norm at most `1/(d-2w)`. Its correction to (10)
is positive and has norm at most epsilon. Each of the four ordered low
energies therefore differs from its compressed energy by at most epsilon.
All other states have energy at least `d-w` above the subtracted free
sea reference. The same argument in each conserved number sector retains
the labels of the four low states.

Kramers symmetry forces the actual N=V pair to be exactly degenerate,
and (11) makes it the full-source ground. Particle-hole symmetry equates
the lowest N=V-1 and N=V+1 excitation energies. Their common charged gap g
satisfies

\[
 0<g_*\le g\le U/2+2\epsilon.
 \tag{12}
\]

Every other charged excitation lies above `d-2w` relative to the actual
ground. For each of these three number sectors, the resolvent equation
gives `||(1-P)psi||<=t||psi||`. Equal-rank principal angles give
`||G-G0||<=t`, `||C_+/- - C_+/-0||<=t`, where G has rank two and each
charged projector C has rank one. Thus the normalized ground mixtures
and the respective charged pure states have unhalved trace distance at
most `2t` from their free-sea core references.

The positive addition/removal spectral measures of `c(f)` on `G/2`
have atoms at g. On the entire one-body space their residue matrices obey

\[
 \|M_\pm-Z/2\|\le3t,
 \tag{13}
\]

where Z is the exact two-dimensional zero projection of the same bounded
radial h. Changing the ground density costs `2t` and the charged
projection costs t in a normalized CAR quadratic form. This applies
the complete finite pole argument in
[defect_core_response](../defect_core_response/PROOF.md) with the actual
symmetries, frame and interacting reference proved above.

Each cube residue has rank exactly two. Indeed rank is at most two from the
two ground-to-one charged transitions, while (11) gives `t<I/V<=1/V`:
`UI/2=2wI/V>2wt`. Since each cube has V>=27, this implies `3t<1/2`.
The independent chain controls use that same bound explicitly, and the
API does not certify rank two when it fails. Kramers symmetry makes the zero-frame compression
scalar, and particle-hole symmetry equates its two signs. For each
normalized zero-frame smearing, the low weight is in `[1/2-3t,1/2]`,
and the total weight above `d-2w` is at most `3t` for each sign. The
anticommutator correlation is `2z cos(gs)+R(s)`, `|R(s)|<=6t`.
At the actual origin and for each spin, the pole weight is at least
`alpha_*/2-3t`. These are positive spectral measures on (2); no native
long-time estimator is asserted to have been executed.

## An actual local matter observable

At the inversion-fixed origin, (8) gives
`R_-(0,0)=(1-alpha)I_2/2`. The free references have local covariance
eigenvalues `(s,s)`, `(s+alpha,s)` and `(s+alpha,s+alpha)`, respectively,
where `s=(1-alpha)/2`. Wick's rule for these Slater states gives the
expectation of the actual occupation-parity projector

\[
 D_0=(n_{0\uparrow}+n_{0\downarrow}-1)^2
 \tag{14}
\]

as `(1+alpha^2)/2` in the empty/full core and `(1-alpha^2)/2` in each
singly occupied core state. Their equal mixture has the same value.
This does not identify D0 with a zero-orbital occupation. It is the native
local occupation operator, and its local centered interaction is exactly
`W_center=(U/2)D0-U/4`.

For an effect `0<=D0<=I`, the unhalved trace bound `2t` changes its
expectation by at most t. Consequently the actual ground and charged
states of (2) satisfy

\[
 \langle D_0\rangle_G\le(1-\alpha_*^2)/2+t,\qquad
 \langle D_0\rangle_{C_\pm}\ge(1+\alpha_*^2)/2-t,
 \tag{15}
\]

with charged-minus-ground contrast at least `alpha_*^2-2t`. A source
instrument receiving this local matter observable therefore has a typed
state/algebra input; a geometry or Einstein conclusion needs its own
operation and refinement theorem.

## A faithful finite reference and its modular clock

The same full Hamiltonian has faithful stationary Gibbs state
`rho_beta=exp(-beta H)/Tr exp(-beta H)` for every finite beta>0.
Its complete ground gap is at least g_*: the only four low states are the
ground doublet and charged pair, and the rest lie above `d-2w>g_*`.
For an integer n>=1 choose `beta=(2V+n)/g_*`. The ground degeneracy is
two, so its excited Gibbs probability is bounded by

\[
 p_{\rm exc}\le\frac{2^{2V}-2}{2}e^{-\beta g_*}
 <2^{-n-1}.
 \tag{16}
\]

Here `e>2` supplies an exact elementary bound. Since the Gibbs density
is uniform on the ground doublet, `||rho_beta-G/2||_1=2p_exc<=2^-n`.
Any effect expectation incurs at most `2^(-n-1)` additional error. On
the full matrix/CAR algebra its modular automorphism is exactly
`rho_beta^(it) A rho_beta^(-it)=tau_(-beta t)(A)` for Hamiltonian flow
`tau_s(A)=exp(isH)A exp(-isH)`. The finite faithful input and clock are
explicit; preparation, equilibration under repair and geometric modular
identification are separate requirements.

## Retained parameters and boundaries

The retained cube uses `R=1,a=v=1,m=64,k=2,U=1/65536`. It has
`d=62`, `alpha_*=7688/7691`, and `t=27/16252874`. Rational bounds give
`7.62310e-6<g<7.62974e-6`, origin-spin pole weight greater than0.49979,
ground parity below0.000392, and charged parity above0.999608. Exact
fractions, rather than these display decimals, determine admission.

The independent three-site chain uses the same source restricted to one
coordinate axis and all 64 Fock states. Its kinetic norm and origin-coupling
norm are smaller than the three-dimensional bounds used by the executable
API. Those conservative bounds therefore apply to the chain controls;
no three-dimensional many-body matrix is replaced by a chain calculation.

The one-body strong-mass gap is uniform in cube size. The contact estimate
here is global, with `w=UV/4`; it proves no volume-uniform interacting
gap at fixed U. The strong-mass premise `m>sqrt(3)v/a` also does not
provide a fixed-mass small-spacing continuum. The supplied bounded mass
profile, finite reference preparations and source units are explicit.
No Lorentzian interacting continuum, observed particle species, repair
protection, regional entropy/area law or Einstein dynamics is inferred.
