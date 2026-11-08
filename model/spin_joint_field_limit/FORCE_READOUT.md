# A filtered current of the actual joint-family clock

Use the state, physical spin frame and H_B clock in `PROOF.md`, with
p=n^-2, B=Jn^-5, n>=2. The symbol P_p below denotes the canonical field
built from the uniform internal spin. It is distinct from the generator
of spatial translations.

For a real compact smooth G and positively oriented cubic bonds put

\[
 P_p(G)=-p\sum_xG(px)S_x^2,\quad
 Q_p(G)=p^2\sum_x\epsilon_xG(px)S_x^3,\quad
 C_{x,i}^2=(S_x\mathbin\times S_{x+e_i})^2,
\]

\[
 D_i^pG(X)={G(X+pe_i)-G(X)\over p},\qquad
 Z_p(G)=p^3\sum_{x,i}D_i^pG(px)C_{x,i}^2.
\]

The current has norm 1/2 on a bond. Exact source continuity and the
selection torque from `CLOCK_COMPARISON.md` give

\[
 i[L_{B,p},P_p(G)]=-p^{-2}Z_p(G)-\lambda_p Q_p(G),
 \quad L_{B,p}={H_B-E_B\over Jp},\quad
 \lambda_p={B\over Jp^2}={1\over n}.                    \tag{1}
\]

This is an operator identity on every finite box containing the support
and its incident bonds. The same local identity holds in the positive-field
ground representation. It does not replace the source current by a linear
field force.

## An energy filter supplies the natural current scale

Let rho be a complex spectral multiplier with |rho|<=1 supported in
[0,M], and apply it to this same L_B,p. The ground equation and (1) imply

\[
 \rho(L_{B,p})Z_p(G)\Omega_B
 =-ip^2L_{B,p}\rho(L_{B,p})P_p(G)\Omega_B
  -p^2\lambda_p\rho(L_{B,p})Q_p(G)\Omega_B.              \tag{2}
\]

Consequently, writing r_p=||P_p(G)Omega_B|| and q_p=||Q_p(G)Omega_B||,

\[
 \|\rho(L_{B,p})Z_p(G)\Omega_B\|
       \le p^2(Mr_p+\lambda_p q_p).                     \tag{3}
\]

The joint-family transverse bounds make r_p and q_p uniformly bounded.
Thus (3) gives the quantitative natural O(p^2) scale of the actual
energy-filtered current. No estimate on an H_0-evolved state is used.
For the Fourier normalization of `PROOF.md`,

\[
 \limsup_{n\to\infty}p^{-2}
 \|\rho(L_{B,p})Z_p(G)\Omega_B\|
 \le {M\over2}\left(\int |k||\widehat G(k)|^2\right)^{1/2}.
                                                               \tag{4}
\]

If a limiting vector is needed, use a continuous compactly supported
rho and the existing enlargement by the uniformly bounded vectors
P_p(G)Omega_B and Q_p(G)Omega_B. The projection onto the original
bounded-history Hilbert space identifies their projected limits with
the actual Stone-field vectors P(G)Omega and Q(G)Omega, by the weak
derivative argument in `spin_canonical_energy/README.md`. That projection
commutes with time translations. Both rho(E) and E rho(E) are bounded
continuous functions. Passing (2) to this common Gram limit and projecting
therefore identifies the current vector in the original representation as

\[
             j_\rho(G)=-iL\rho(L)P(G)\Omega.            \tag{5}
\]

It is on the domain of L because the multiplier has compact energy
support. Any orthogonal norm in the enlarged source-vector limit is
excluded from the identification (5). No equality of all source-vector
norms, common unbounded-product domain or full-history force closure
is inferred. At fixed n, (2)-(3) also allow sharp finite spectral filters;
passing their discontinuities through a continuum limit requires a
separate boundary-weight estimate.

## The stronger cancellation required for a linear force

For a proposed real coefficient a and the ordinary continuum Laplacian
of G define

\[
 Y_p(G,a)=Z_p(G)+p^2aQ_p(\Delta G).
\]

The filtered defect from the proposed equation
i[L_B,p,P_p(G)]=aQ_p(Delta G) is exactly

\[
 \rho(L_{B,p})\{i[L_{B,p},P_p(G)]-aQ_p(\Delta G)\}\Omega_B
 =-p^{-2}\rho(L_{B,p})Y_p(G,a)\Omega_B
  -\lambda_p\rho(L_{B,p})Q_p(G)\Omega_B.                \tag{6}
\]

The second term vanishes. The norm of the finite-source filtered force
defect on the left of (6) tends to zero if and only if

\[
            \|\rho(L_{B,p})Y_p(G,a)\Omega_B\|=o(p^2).  \tag{7}
\]

This equivalence is stronger than a statement about its projected
limiting field vector: the projection in (5) can discard orthogonal
source-vector norm. Equation (7) is sufficient for the projected linear
equation; no necessity for that weaker conclusion is asserted.

The natural O(p^2) bound in (3) supplies no cancellation in (7).
The generic gap estimate for an unamplified derivative-weighted current
is at most (792/n+24/n^6)a_0^2 V_0, with the coefficient and support bounds
defined in `PROOF.md`. Its mean is zero because an internal pi rotation
about axis one reverses the transverse current and preserves omega_B.
Multiplication by p^-4 changes that bound into
(792n^7+24n^2)a_0^2 V_0. Applying the energy filter first, as in (3),
avoids this divergent bound and gives bounded force vectors. It gives
neither a vanishing residual nor a spectral atom. The required estimate
concerns current-current and current-field cross correlations in the
specified ground and spectral window. Scalar displacement-energy
coefficients and two canonical pairs do not determine those correlations.

## Complete finite records of the current norm and residual

For a finite source the specification determines H_B, the actual source
currents, and every chosen finite spectral multiplier. Set

\[
 A=\rho(L_{B,p})Y_p,\quad c\ge\|A\|,\quad c>0,
 \quad F=A^*A/c^2.                                    \tag{8}
\]

The two effects F and I-F are complete. The first outcome has original
probability ||rho Y_p Omega_B||^2/c^2. The other outcome is retained.
For the current alone take Y_p=Z_p. A finite exact or certified matrix
calculation can use the source norm sum as c; it does not condition the
prepared state on a spectral window or divide by its acceptance weight.

There is also a direct complete Pauli readout that avoids assuming an
additional spectral measurement primitive. On N spins expand
F=sum_u f_u sigma_u, with real f_u=2^-N Tr(F sigma_u). Choose one of
the 4^N Pauli settings uniformly and measure its nonidentity axes using
the existing native local probes. Identity axes have a deterministic
plus outcome and a retained zero-probability minus label. For record
(u,s), output the score

\[
           4^N f_u\prod_{x:u_x\ne I}s_x.               \tag{9}
\]

Its expectation is Tr(rho F); every setting and every pointer record is
included. If S=4^N max|f_u|, the score lies in [-S,S]. An error delta in
the complete output half trace distance changes the expectation by at
most 2S delta. Native failure records have score zero and are included
in that output distance. Multiplying the score by c^2 gives the unnormalized
filtered-current or residual squared norm, with error at most
2c^2 S delta. Dividing by p^4 to test (7) requires that scaled error
to tend to zero. The positive-field source preparation and native probe
planners permit any specified finite error; no practical or uniform
continuum cost is inferred. A certified operator approximation to F
adds its operator-norm error to the expectation error.
For R independent fresh preparations, Hoeffding's inequality bounds the
probability of an additional score-mean error u by
2 exp(-R u^2/(2S^2)). Thus a finite requested confidence and accuracy
specify a finite repetition budget, with preparation and compilation
errors charged separately. This is a finite observable protocol, not
a nonvanishing continuum probability for an artificially narrowed band.

The finite control uses the actual two-spin source at n=2, p=1/4 and
B=1/32. Its full H_B spectrum is

\[
 E_-=-(16+\sqrt{257})/32,\quad0,\quad0,\quad
 E_+=(-16+\sqrt{257})/32.
\]

The rational projector Pi_0=I-H_B(H_B+I)/B^2 is the complete zero-H_B-energy
doublet. Relative to the ground its scaled excitation is
(16+sqrt(257))/8<5, so rho=Pi_0 is an admitted [0,5] filter. Exact algebra
checks (1)-(3), all four source rows, the full spectrum, the complete
two-effect current record and the 64 labelled Pauli records in (9).
The record scores reproduce the original success probability exactly.
This doublet is a finite source control. Its single finite energy is
not a sharp continuum branch or a measured particle mass.
