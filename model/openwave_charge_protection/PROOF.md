# Protected integer winding and the missing current relation

## The actual compact texture law

Use the unchanged source and complete four-proposal law of
[`compact_defect_transport`](../compact_defect_transport/PROOF.md).
Its two zero cores are fixed at \((0,0,\pm1)\), and the noncore direction
is the normalization of \(f(x,y,z)=(x,y,z^2-1)\) inside its compact
support. The exterior is north. Each outward oriented unit cube shell
has 48 triangles. A signed count of regular-target preimages gives
degrees \(-1\) and \(+1\). The executor independently solves each
triangle's rational barycentric linear system, and agrees with the
owning determinant implementation on all three continuing texture states.

Every accepted vertex replacement obeys the native relative budget and
the nonnegative scalar-product guard. The normalized affine interpolation
between the old and proposed vertex is nonsingular on each incident
triangle: the squared norm of a positive linear combination is the sum
of positive diagonal terms and nonnegative pair terms. The interpolations
agree on shared edges. They are therefore a homotopy of the complete
oriented shell. Rejected records are identity maps. Induction preserves
both degrees on every finite history, including repeated states and
rejections. This proves a complete-law result without enumerating an
exponentially growing list of histories.

The executor retains all twelve source-state/proposal records and their
exact probabilities. An additional finite guard sweep tests all shell
vertices against signed coordinate perturbations, zero directions and
sign reversals. It uses the same acceptance law and does not replace the
declared stochastic proposal distribution.

## Current quantization on the complete fermion carrier

There is a distinct exact quantization route through the actual native
current module. Write \(Q=T_3+Y\), and let \(C\) be the native colour
Cartan generator induced by
\(\operatorname{diag}(1,1,-2)/3\) on the three colour slots. The executor
checks all fifteen internal weights, obtaining \(3Q_i\in\mathbb Z\)
and \(Q_i+C_i\in\mathbb Z\). Second quantization gives, on every
occupation sector and every number of spatial or family copies,

\[
\operatorname{spec}(d\Gamma(Q))\subset\tfrac13\mathbb Z,
\qquad e^{2\pi i d\Gamma(Q)}\Gamma(Z_{\rm colour})=I.
\]

On the colour-singlet subspace \(\Gamma(Z_{\rm colour})=I\), so the
charge spectrum is integral. In particular a half-unit charge eigenvalue
is absent. This extension does not require enumeration of the full Fock
carrier: exterior products multiply the verified one-particle phase
identities. Every actual texture-record field map is
\(\sqrt{p_r}\,\Gamma(U_r\otimes I_{45})\). It commutes with this current,
so each complete accepted or rejected record preserves its charge sectors.

This is a source-current theorem with the declared generator and its
normalization. It does not need an equality between current and winding.
An expectation value can be a half: the equal mixture of charge zero and
charge one has mean \(1/2\), while its retained measurement outcomes are
zero and one. Thus the absence of a half eigenvalue is not a restriction
on mean charge or a proof that electric charge has been physically
identified with the native current.

## What unit winding does and does not forbid

The degree of a closed continuous map from an oriented sphere to the
direction sphere is an integer. If two equally charged daughter shells
had degrees \(q,q\), additivity would require \(2q=1\) or \(2q=-1\).
Neither equation has an integer solution. This is an integer-degree
obstruction, independent of the modular-centre obstruction that fails
for the electron-like singlet. The supplied law additionally fixes both
cores and forbids a noncore zero, so creation of an extra core is outside
its transition alphabet.

Integer degree alone does not prohibit the unequal signed decomposition
\(1=2+(-1)\), nor does it supply an energy cost for such a deformation in
a larger law. Removing the guards permits a texture with all noncore
directions north and degrees \((0,0)\), while keeping both zero cores.
The four zero onsite fermion coordinates survive this erasure. Thus a
low spectral-band count is not a proof of winding protection.

The arithmetic lemmas in `Protection.lean` check the no-equal-halves
claim and the signed counterexample. They do not formalize the sphere
degree theorem or the native shell homotopy; those hypotheses and their
finite source checks are stated above.

## The integer register is a different observable

The simulator's law-C register attaches an integer \(n_s\) to each seam
and reads its oriented face curl \(m_f\). An accepted repair of current
\(j\) gives \((m_+,m_-)\mapsto(m_++j,m_--j)\). The exact oriented sum is
conserved, and \(m_f\bmod6\) equals the stored class. The update rule
consumes that class and the load current; it does not consume the
unreduced lift.

All \(6^2\cdot2\cdot6^2=2592\) local load/coin/class cases are checked,
with nine independent integer-sheet choices in each case. A fresh
80-carrier source run retains 200 successive cuts and checks Stokes and
the modulus identity after each sweep. The lift is not a conserved
individual particle label. For example, loads \((4,0)\), coin zero and
face values \((4,0)\) give an accepted \(j=2\) update:

\[
(4,0)\longmapsto(6,-2),\qquad (4,0)\bmod6\longmapsto(0,4).
\]

The nonzero class has moved to the neighboring face, while its integer
value changes from \(+4\) to \(-2\); the original face retains the
hidden value \(+6\). Consequently identifying this lift with the
unchanged particle hypercharge is not licensed by class conservation.
It is also not the degree of the compact texture's direction map.

## Exact obstruction to an electric-current identification

The compact source retains its spatial-spin field and the actual
fifteen-state native internal current module, with three declared family
copies. Its uniform candidate electric current is the native
\(Q=T_3+Y\). The texture is a classical control register; the fermionic
occupation algebra is a separate factor. On a fixed plus-one texture,
take the same spatial band state and either a neutral internal mode or
the native charge-one internal mode. In their two-dimensional span,

\[
W=\begin{pmatrix}1&0\\0&1\end{pmatrix},\qquad
Q=\begin{pmatrix}0&0\\0&1\end{pmatrix}.
\]

No function \(f\) can satisfy \(Q=f(W)\) there. This rules out not only
a guessed normalization \(Q=cW\) but every universal dictionary from
winding alone to the existing current on the unchanged full carrier.
The vacuum supplies another state with the same texture and zero current.
The native spatial transport is \(U\otimes I_{45}\), so it preserves
the uniform current while transporting the same winding; current
conservation does not enforce their equality.

The compact construction supplies no local Gauss constraint tying its
classical degree to those fermion occupations. A restriction or coupling
could define a different source on which such an equality holds, but it
would require its own invariance, complete-record and energy tests. The
present calculation establishes protected unit topology on its declared
law and an exact attachment obstruction. It does not prove electromagnetic
electron indivisibility or derive the electric-charge normalization.
