# Weighted reversal of complete repair records

The biased ordered-cover source has a complete finite quantum instrument
whose Petz reverse consists of its own raw branches. This statement uses
the same 132-dimensional carrier, seam selector, tagged C/E laws, endpoint
coins and persistent ordering memory as [the native reference](PROOF.md).
The reverse does not require a replacement repair law or a fitted energy.

## The reference and the full instrument

Write the 66 occupation configurations as (x), with positive weights

\[
 \pi_x=Z^{-1}\prod_{i:n_i(x)=1}\theta_i,
 \qquad \theta=(25,25,8,\ldots,8),\qquad Z=7505.
\]

The memory flip is (X), its projectors are (P_\pm=(I\pm X)/2), and

\[
 \rho_\alpha=\operatorname{diag}(\pi)\otimes\sigma_\alpha,
 \qquad \sigma_\alpha=\alpha P_++(1-\alpha)P_-.
\]

For (0<\alpha<1), this reference is faithful. For \(\alpha=0\) or
\(1\), every operation preserves the supported charge corner; all inverse
powers below are restricted to that corner. No inverse is assigned to an
absent charge sector.

A complete label is (h=(\ell,e,c,r)), where \(\ell\) is C or E,
\(e=(a,b)\) is one of the 30 seams, \(c\) is the endpoint coin, and
\(r\) is its full raw address. There are 720 labels, including impossible
and abort addresses. The source produces

\[
 K_h=\sqrt{w_h}\,V_h,\qquad
 V_h=\sum_{x\in D_h}|f_h(x)\rangle\langle x|
           \otimes X^{p_h(x)}.
\]

Here (f_h) is an injective partial map, (p_h(x)) is the actual ordering
change, and (w_h) is the squared Kraus amplitude supplied by the native
selector. It includes the seam and law probabilities and the biased coin
probability. The (w_h) are not outcome probabilities: the probability of
an outcome also contains the weight of its input domain.

## Exact weighted branch reversal

At the same seam and for the same tagged law, define

\[
 (C,0,1)\longleftrightarrow(C,1,2),\qquad
 (E,0,5)\longleftrightarrow(E,1,6).
\]

The displayed triples give law, coin and raw address. Every other label is
fixed. In particular every holding branch is self-dual, and each zero or
abort branch retains its address. Fixing a zero operator is a declared
completion of the record involution; the zero operator alone cannot select
a unique label.

For every nonzero column (x\mapsto y), source enumeration gives

\[
 V_{\bar h}=V_h^\dagger,\qquad
 w_h\pi_x=w_{\bar h}\pi_y.
 \tag{1}
\]

The reversal of a move uses the opposite endpoint coin. Its ordering flip
equals the forward flip. The two biased amplitudes in (1) are therefore
different when the endpoint weights differ. Holding columns have
\(x=y\), flip zero and the same weight on both sides. Zero and abort
slots have empty support. The executable certificate checks (1), exact
completeness and the full involution on all 720 labels and both memory
orientations.

The memory factor \(\sigma_\alpha\) commutes with (X). For the column
of \(V_h^\dagger\) from \(y\) to \(x\), equation (1) consequently gives

\[
 \rho_\alpha^{1/2}K_h^\dagger\rho_\alpha^{-1/2}
 =\sqrt{w_h}\sqrt{\frac{\pi_x}{\pi_y}}V_h^\dagger
 =K_{\bar h}.
 \tag{2}
\]

This is an operator identity on the full supported quantum carrier.
It applies to coherences and remains an identity after tensoring every
operator with the identity of an arbitrary untouched reference system.
No diagonal-input restriction enters (2).

## Detailed balance and the clock

Let \(\Phi(A)=\sum_hK_hAK_h^\dagger\), and let \(\Phi^\dagger\) be
its trace adjoint. The reference is stationary. The Petz reverse with that
reference is

\[
 \mathcal R_{\rho_\alpha,\Phi}(A)
 =\rho_\alpha^{1/2}
   \Phi^\dagger(\rho_\alpha^{-1/2}A\rho_\alpha^{-1/2})
   \rho_\alpha^{1/2}
 =\sum_hK_{\bar h}AK_{\bar h}^\dagger
 =\Phi(A).
 \tag{3}
\]

Thus the actual unread native channel is its weighted Petz reverse. The
Heisenberg channel \(T=\Phi^\dagger\) is self-adjoint for the finite
Kubo-Martin-Schwinger inner product

\[
 \langle A,B\rangle_{\rho,1/2}
 =\operatorname{Tr}(\rho^{1/2}A^\dagger\rho^{1/2}B).
\]

Indeed, if \(\Gamma(A)=\rho^{1/2}A\rho^{1/2}\), (3) states
\(\Gamma T=T^\dagger\Gamma\), which is precisely the weighted
self-adjointness identity. The declared reference and the actual instrument
therefore supply quantum detailed balance without a population-only
replacement.

For \(K_\alpha=-\log\rho_\alpha\), an actual recorded current (j)
on seam \((a,b)\) obeys

\[
 [K_\alpha,K_h]=\Delta_hK_h,
 \qquad \Delta_h=j\log(\theta_a/\theta_b).
\]

Equation (1) gives \(w_{\bar h}/w_h=e^{\Delta_h}\), and reversal
changes the modular transfer to \(-\Delta_h\). The clock and the
weighted branch reversal use the same reference. The modular parameter,
repair-attempt count and physical time retain their separate definitions.

## Complete histories and the control law

For a finite complete record word \(\boldsymbol h=(h_1,\ldots,h_n)\), set
\(A_{\boldsymbol h}=K_{h_n}\cdots K_{h_1}\). Reverse both label and order:
\(\overline{\boldsymbol h}=(\bar h_n,\ldots,\bar h_1)\). Consecutive
reference factors in (2) cancel, proving for every finite word

\[
 A_{\overline{\boldsymbol h}}
 =\rho_\alpha^{1/2}A_{\boldsymbol h}^\dagger\rho_\alpha^{-1/2}.
 \tag{4}
\]

The identity includes holding and zero histories. With stationary initial
state \(\rho_\alpha\), cyclicity of the finite trace gives

\[
 \Pr_{\rho_\alpha}(\boldsymbol h)
 =\operatorname{Tr}(A_{\boldsymbol h}\rho_\alpha A_{\boldsymbol h}^\dagger)
 =\Pr_{\rho_\alpha}(\overline{\boldsymbol h}).
 \tag{5}
\]

This equality uses the fixed native state-independent selector and its
biased coin probabilities. It compares complete records under that
stationary experiment. For an arbitrary input, (4) continues to hold as
an operator identity; the two probabilities in (5) need not agree.
Outcome-dependent feedback requires its own reversed controller and
record probabilities. Equation (5) supplies no reversal theorem for an
arbitrary adaptive policy.

The certificate checks all two-event raw histories on two declared seams,
including empty supports, and both three-event odd exchange loops. The
all-length statement follows from (2), rather than extrapolation from
these finite examples. The odd loops return the occupation configuration,
have zero total modular transfer, and retain the memory action (X).

## Reversal and irreversible repair

Petz self-duality is a relation between channels with a specified reference.
It does not undo an observed sample path. Two orthogonal cover basis
states have strictly overlapping unread outputs; the certificate computes
their exact trace-distance contraction. A completely positive trace-
preserving left inverse on all states would contradict that contraction.
Tracing the complete repair tape is part of this unread comparison.

Equations (2)-(5) establish a finite quantum reversal and detailed-balance
interface for the declared weighted source. They select no coherent
interacting Hamiltonian, antiunitary physical time reversal, laboratory
temperature, stress tensor or relativistic field limit. The ordered-cover
memory is an explicit part of the source construction. Its consistency
does not identify that extension with a physical particle sector.

Run the exact certificate through
`python -m model.protected_reference.reference --check` after its source
binding, and the direct tests with
`python -m pytest -q model/protected_reference/test_reversal.py`.
