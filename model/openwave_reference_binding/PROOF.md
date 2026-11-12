# Constrained binding in the complete finite charged source

## Carrier, energy and measurement

The source is the eight-cell, 240-fermion-mode construction with eight
five-state scalars and 24 complete 119-state Peter-Weyl links. Its local
Yukawa coefficients and scalar/link penalties equal one. The Hamiltonian
is the same

\[
H_{\rm tot}=H_0+2^{-34}T+4N_{\rm ref},\qquad \|T\|\leq360,
\]

used by [charged response](../coherent_defect_wave/chiral_joint/interacting_charged_response/PROOF.md).
The complete 30-mode reference is retained at cell zero. The selected
sector contains one reference particle of native type \(e^c\) or \(L\).
Its charge is compensated by source matter and links; the label is a
native representation, without a physical particle identification.

Let \(\Pi_y\) impose local matter singlets at every cell other than \(y\)
on the joint physical source/reference sector. It retains every link and
the full local matter space at \(y\). The eight projectors are orthogonal;
their complete complement contains screened clouds. Define the constrained
minimum

\[
\mathcal E_r(y)=\inf_{\psi\in\operatorname{ran}\Pi_y,\ \|\psi\|=1}
\langle\psi,H_{\rm tot}\psi\rangle.
\]

These are expectation minima of compressed operators. The subspaces are
not invariant under the full Hamiltonian. No stationary charged pole or
energy-resolved transport follows from these minima.

## Exact physical strings and their local endpoints

The link basis and all twelve endpoint currents are the actual
[Peter-Weyl source](../coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md).
On a path of \(\ell\) edges, the normalized state of a reference in
representation \(r\), the path links, and a conjugate endpoint is

\[
|\mathcal W_{r,p}\rangle=d_r^{-(\ell+1)/2}
\sum_{a_0,\ldots,a_\ell}
|a_0\rangle_{\rm ref}\otimes
\bigotimes_{j=1}^{\ell}|r,a_{j-1},a_j\rangle
\otimes|\bar a_\ell\rangle_y .
\]

For a backwards traversal the stored link uses the dual orientation and
reversed indices. At the origin, \(T\) on the reference cancels \(-T^T\)
on the outgoing link. At an interior vertex the two endpoint actions
cancel; at the destination the matter action is \(-T^T\). These identities
are checked for all twelve actual currents on every shortest stored path
to all eight cells, including parallel oriented edges. There are 6,648
exact local Gauss checks. The zero-length state is the normalized local
reference/matter singlet. Reversing the matter conjugation fails Gauss.

The endpoint is attached to the actual local CAR/scalar module. Its
highest-weight basis is the existing complete charged response block.
For \(L\), applying the actual weak lowering operator constructs its
second component. The executable verifies the four nontrivial current
intertwiners, the positive Gram matrix, and

\[
hB=B(I_{d_r}\otimes h_r),\qquad
B^\dagger GB=I_{d_r}\otimes G_r.
\]

The inherited quark module is a colour singlet. Thus the other eight
currents vanish on the entire endpoint embedding. A real isolated root
of the complete local characteristic polynomial selects an admissible
normalized eigenspace. This gives a variational state, with no assertion
that the dynamical ground-preparation apparatus produces it automatically.
The seven spectator cells occupy actual local ground doublets. Every
path link is nonblank and has energy one; all other links are blank.

## Matching lower and upper bounds

Write \(E_0\) for the local physical ground energy, \(e_r\) for the
selected local charged trial energy, and \(\Delta_r^-\) for the complete
compensated local gap lower bound. The inherited complete-sector proof
gives \(h\geq l_0=-11.71\), and
\(\delta=E_0-l_0\), with \(\delta+\Delta_r^-<1\).

For \(y\ne0\), Gauss implies that every nonblank support connects the
reference to \(y\), so it has at least \(\ell=d(0,y)\) edges. A support
of exactly \(\ell\) edges is a simple path. Its endpoint carries the
conjugate reference representation. Supports with additional edges obey
\(H_0\geq8E_0+\ell+1-\delta\), which exceeds the bound
\(8E_0+\ell+\Delta_r^-\). This is the complete-support argument of
[charge localization](../coherent_defect_wave/chiral_joint/interacting_charge_localization/PROOF.md).

The origin requires its own comparison. With no excited links, the local
matter compensates the reference and has gap at least \(\Delta_r^-\).
With one or more excited links, all other cells remain matter singlets,
so \(H_0\geq8E_0+1-\delta\). The same lower bound therefore covers
\(\ell=0\). No flux or local matter state is discarded.

Every hopping word has two distinct matter endpoints. At least one is a
local singlet under \(\Pi_y\), and a single nontrivial fermion operator
has zero singlet compression. Hence \(\Pi_yT\Pi_y=0\) exactly. Together
with the physical string trial this proves

\[
8E_0+\ell+\Delta_r^-+\mu
\leq\mathcal E_r(y)
\leq7E_0+e_r+\ell+\mu,\qquad\mu=4.                 \tag{1}
\]

For \(L\), the complete local comparison identifies its actual charged
minimum \(E_L\), below every excluded quark module. An exact polynomial
greatest-common-divisor and root-isolation check identifies it with the
trial root. The stronger inequality \(\delta+E_L-E_0<1\) excludes all
additional-flux sectors. Therefore

\[
\mathcal E_L(y)=7E_0+E_L+\ell+\mu,\qquad
\mathcal E_L(y)-\mathcal E_L(0)=\ell.                \tag{2}
\]

For \(e^c\), the visible lowest pole is an upper bound on the complete
charged minimum. Its gap is at most 0.691656680939; the complete local
gap is at least 0.100585294247. Equation (1) gives the conservative binding
intervals below. The receipt retains exact rational endpoints.

| Graph distance | \(L\) constrained binding cost | \(e^c\) constrained binding cost |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 1 | [0.408928613309, 1.591071386691] |
| 2 | 2 | [1.408928613309, 2.591071386691] |
| 3 | 3 | [2.408928613309, 3.591071386691] |

These costs use the declared source energy unit. There is no empirical
mass or length input. The actual interacting neutral ground energy lies
in \([8E_0-v,8E_0+v]\), where \(v=360\,2^{-34}\). The receipt includes
absolute total minima and excitation intervals relative to that same
ground. The constant reference cost is explicitly subtracted in the
source-only excitation table. Both this cost and the unknown common
ground shift cancel exactly from the binding differences (2). Subtracting
an energy zero cannot remove a spatial binding cost.

## Retained screened configurations

The complete complement \(\Pi_{\rm cloud}\) is separate. The executable
retains the existing physical origin pump followed by the native remote
hop on edge nine. It has nonzero norm, one excited link and a charged
cell at distance two without a flux path from the reference. Every
\(\Pi_y\) annihilates it. Thus (1) and (2) cannot be extended to arbitrary
remote charged cells or used as an unrestricted confinement theorem.

No screened-cloud minimum, string-breaking threshold, freely propagating
charged band or Coulomb limit is calculated. The result identifies
reference binding on the specified single-unscreened-charge preparation
domain of the complete finite source.
