# Scalar source response at the magnetic node

The proper-cubic averaged centered full-CAR ground of `MOMENTUM_RESPONSE.md`
has two scalar signed response channels at its magnetic kinetic node. Source
particle-hole symmetry relates them by energy reflection, leaving one
independent scalar measure. The eight transported original Weyl-valley probes
fall into these two channels according to their kinetic winding. Their nodal
Lorentzian scores are equal on the interacting ground. A scalar measure is
not assumed to be a particle pole, and this reduction supplies no mass value.

The every-momentum extension in `SPECTRAL_RESIDUE.md` gives the nodal measures
their precise meaning. Statements at one momentum do not by themselves give
uniform estimates for growing packets with shrinking energy resolution.

## Actual little-group matrices

Use the physical-cell Bloch convention, parity components `s in {0,1}^3`,
and the same two-step magnetic translations as the source construction.
Write

    k_*=(pi/2,pi/2,pi/2),
    h(k)=c sum_j cos(k_j) Gamma_j,
    J=Gamma_1 Gamma_2 Gamma_3,   C=iJ.

The actual signed Gamma matrices satisfy `J*=-J`, `J^2=-I` and `Tr J=0`.
Consequently `P_+=(I+C)/2` and `P_-=(I-C)/2` are orthogonal rank-four
projectors in the eight-dimensional parity space. Hidden spin contributes
an independent identity factor of dimension two.

At the node a one-step magnetic translation is `-i A_j`, where

    (A_j v)_s=(-1)^(sum_(l<j) s_l) v_(s xor e_j).

For a proper signed permutation `R e_j=epsilon_j e_(pi(j))`, set
`m_(pi(j))=(epsilon_j-1)/2`. The actual lifted rotation, expressed back in
the same Bloch chart, is the signed permutation

    (W_R)_(R s mod2,s)
      =(-1)^(m.(R s mod2))
         (-1)^(sum_(i<j,pi(i)>pi(j)) s_i s_j).          (1)

The first sign transports `R k_*=k_*+pi m` across the reciprocal-cell
boundary. Omitting it gives a different representation. Direct source
identities are

    [A_j,Gamma_l]=0,
    W_R Gamma_j W_R*=epsilon_j Gamma_(pi(j)),
    Comm({A_j,W_R})=span_C{I,J}.                       (2)

The last equality is finite exact linear algebra. Conjugation by each
signed permutation identifies entries of an arbitrary eight-by-eight
matrix with signed entries. Signed-orbit elimination leaves exactly two
independent matrices. Independent rational elimination of all 1,728
commutator equations has nullity two. Magnetic translations alone have
an eight-dimensional commutant; the cubic lifts provide the reduction.

An independent periodic four-by-four-by-four source calculation verifies
the translation and rotation matrices directly on literal hopping links
and Bloch columns. Their second quantizations also preserve the onsite
contact action: they permute sites and multiply both local flavours by
the same sign. No replacement of bare repair under a changed chart is
inferred.

## The same interacting spectral measure

The averaged ground, action and retarded clock are invariant under the
stated symmetries. Hidden-spin invariance removes the hidden-spin matrix
factor. Each original matrix spectral measure at the node commutes with
the operators in (2). Hence, for every Borel set F,

    nu_(k_*)(F)=[mu_+(F)P_+ + mu_-(F)P_-] tensor I_2. (3)

Here each mu is a positive scalar probability measure: positivity follows
by testing a unit vector in its projector, and the total matrix measure
is the identity. This holds for the actual interacting reference, without
factorization of its correlations.

There is an additional source identity. In the real signed occupation
chart, K and its normalized ground projector are real, and this property
survives the translation and proper-cubic averages. Define
`delta a=[a,K]` and the original moments

    H_n(x,y)=omega({delta^n a_x,a_y*}).

They are real. Particle-hole conjugation sends
`a_x -> (-1)^|x| a_x*`, preserves K, and obeys
`delta^n(a_x*)=(-1)^n(delta^n a_x)*`. Invariance of omega gives

    H_n(x,y)=(-1)^(n+|x|+|y|) H_n(x,y).

With `D_ss=(-1)^|s|`, the Bloch matrices therefore obey
`D H_n(k)D=(-1)^n H_n(k)`. The uniform exponential-moment determinacy of
`UV_TAIL.md` extends this equality from moments to spectral measures:

    nu_k(-F)=D nu_k(F)D.                              (4)

Since `D C D=-C`, equations (3)-(4) imply

    mu_-(F)=mu_+(-F).                                 (5)

This concerns the combined signed retarded measure. Addition and removal
are retained as their separate positive physical-energy records; the
negative sign in retarded time is not a negative physical excitation
energy. Equation (5) does not assert equality of the two measures at
each signed energy. Every even spectral score, including
`eta^2/(eta^2+E^2)`, has equal expectation in both channels.

## Transport of the original eight valleys

Let U_s be the existing exact onsite frame from the original Weyl chart
to the signed chart. The two columns associated with an original node
`pi b`, `b in {0,1}^3`, have parity blocks

    V_b(s)=8^(-1/2) exp(-i k_*.s) (-1)^(b.s) U_s.      (6)

The sixteen columns of all V_b form a unitary matrix. Source arithmetic
gives, with `chi_b=-(-1)^|b|`,

    (C tensor I_2)V_b=chi_b V_b,
    -(Gamma_j tensor I_2)V_b=V_b[-(-1)^b_j sigma_j].    (7)

Thus the exact transported probes have nodal response

    V_b* nu_(k_*)(F) V_(b')=delta_(b,b') mu_(chi_b)(F) I_2.

There are four original valleys in each channel. Their kinetic degrees
are the windings in `KINETIC_NODE_MAP.md`; the mapping in (7) holds as a
probe identity for the interacting source too. It does not establish that
these probes create isolated relativistic particles. In particular the
eight valley labels do not supply eight independently selected mass laws.

## Consequence for the unresolved memory

Both original and composite probes transform by the same source symmetry
matrices. The full two-block response and its exact Schur memory inherit
the two-channel decomposition. At `z=i eta`, `eta>0`, selfadjointness and
(4), with the composite probe's opposite particle-hole sign, give

    Sigma(i eta)=a_eta C-i b_eta I,
    a_eta real, b_eta>=0.

Define `m_eta=eta b_eta`, `j_eta=eta a_eta`, and
`H_eta=alpha^2+eta^2+m_eta`, with `alpha=U/2`. The original Lorentzian
score is exactly scalar:

    -eta Im G_00(i eta,k_*)
        =[1-alpha^2 H_eta/(H_eta^2+j_eta^2)] I.        (8)

The nodal kernel theorem gives the common zero-energy residue
`m_0/(alpha^2+m_0)`, where
`m_0 I=R*1_{0}(L_perp)R`. Its value is undetermined by the positive
total residual Gram. The finite-frequency reduction in (8) helps compute
that missing source quantity; it does not prove a positive m_0 or a
controlled joint limit `eta~a`, `k-k_*~a`. Physical spin, relativistic
covariance, scattering and Standard Model assignments require the same
state's further dynamical response.
