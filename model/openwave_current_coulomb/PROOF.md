# Native current, Gauss projection and a finite Coulomb energy choice

## Actual source current and complete flux readout

Use the periodic PW119 carrier, its actual quantum Yukawa operators and
all native fermion, scalar, link and reference coordinates. At each site
let \(\rho_x^Y\) be matter hypercharge plus the scalar hypercharge. Add
the complete reference current at cell zero. The reference is retained
as a quantum system; restricting a displayed experiment to one \(e^c\)
reference does not replace its full thirty-mode Fock carrier.

Orient every stored nearest-neighbor edge from its left to its right
endpoint and define the counting-metric incidence matrix by
\((Df)_e=f_{\rm right}-f_{\rm left}\). Write
\(b_e=L_{e,0}^Y\) for the actual native origin endpoint current.
The full 119-state identity is \(L_{e,1}^Y=-b_e\). On the joint physical
Gauss carrier,

\[
D^*b=\rho,
\qquad \rho_x=\rho_x^Y+\delta_{x0}Y_{\rm ref},
\qquad \sum_x\rho_x=0.
\tag{1}
\]

Every coordinate in (1) is an actual commuting central-current operator.
The current includes the scalar and reference contributions. Its action
is compared with the native Gauss compiler on arbitrary complete
occupation/link/scalar inputs. The existing charged-transport proof
preserves the physical domain on its preparation, failure and clock
branches.

Let \(L=D^*D\), and let \(G=L^+\) annihilate the constant mode. Then

\[
E_{\rm long}=DG\rho=(DGD^*)b,
\qquad D^*E_{\rm long}=\rho.
\tag{2}
\]

The matrix \(P_{\rm long}=DGD^*\) is an orthogonal projector. Consequently

\[
\|b\|^2=\|E_{\rm long}\|^2+
\|(I-P_{\rm long})b\|^2,
\quad
\frac12\|E_{\rm long}\|^2=\frac12\rho^*G\rho.
\tag{3}
\]

For every other flux \(E\) with the same divergence,
\(E=E_{\rm long}+h\), \(D^*h=0\), so
\(\|E\|^2=\|E_{\rm long}\|^2+\|h\|^2\). This proves the exact
minimum-norm Gauss field. It uses the actual source current and flux,
without identifying two unrelated labels.

Equations (2)–(3) define simultaneous spectral readouts and classical
postprocessing of the native current outcomes. The transverse native
flux outcome is retained. No quantum state is replaced by its
longitudinal projection. On the original Hamiltonian, the quadratic
quantity in (3) is an observable, not the original unit-link energy.

On each of the six actual \(e^c\) arrival effects,
\(\rho=\delta_0-\delta_y\), with the reference fixed at zero. The
native current proof fixes the coefficient to one. In the side-two
counting geometry,

\[
\frac12\rho^*G\rho=G_{00}-G_{0y}=\frac7{48}
\]

for every neighboring arrival. Thus all six positive joint
arrival/charge/flux probabilities also support this exact quadratic
readout. Every preparation failure and the complete nonarrival output
retains its own spectral value of \(\rho^*G\rho/2\). The ten coarse
probabilities are unchanged and sum to one. Nonneutral kinematic inputs
fail the physical Gauss readout: the pseudoinverse returns
\(\rho-\bar\rho\), not \(\rho\).

## What the unchanged link penalty proves

The actual unit-penalty source has

\[
H=H_{\rm local}+N_{\rm nonblank}+\lambda T+4N_{\rm ref}.
\]

Let \(\Pi_y\) impose local matter-plus-scalar singlets at every site
other than \(y\), on the joint physical sector with one fixed reference.
The entire link carrier and target local carrier are retained. Its
complement contains screened clouds.

The local constants and physical endpoint intertwiners of
[`openwave_reference_binding`](../openwave_reference_binding/PROOF.md)
extend over every even periodic volume \(V=N^3\). A nonblank support
in \(\Pi_y\) must connect the reference to \(y\): summing the central
Gauss law over its reference component otherwise leaves an uncompensated
nonzero reference charge. It therefore has at least
\(\ell=d(0,y)\) edges. At least one endpoint of every hopping word is a
local singlet, and a single nontrivial fermion operator has zero singlet
compression. This gives \(\Pi_yT\Pi_y=0\) for all link occupations and
all axial path lengths, independently of \(\lambda\).

The complete local bounds give, for every such volume,

\[
\mathcal E_L(y)-\mathcal E_L(0)=\ell,
\qquad
\ell-w\le\mathcal E_{e^c}(y)-\mathcal E_{e^c}(0)\le\ell+w,
\quad w\simeq0.591071386691.
\tag{4}
\]

Only local spectral bounds, Gauss support, and the exact singlet
compression enter this extension. It does not assume a larger-volume
interacting ground or transfer a side-two ground gap. The source-word
checks cover sides two, four and six. The same finite word argument
holds for every allowed even side.

Suppose the positive term \((\kappa/2)\rho^*G\rho\), \(\kappa>0\),
is added while the original nonblank-link penalty is retained. On
\(\Pi_y\) it is exactly the scalar

\[
\kappa q_{\rm ref}^2(G_{00}-G_{0y}),
\tag{5}
\]

and on \(\Pi_0\) it is zero. Equation (5) therefore shifts each
compressed minimum by a nonnegative amount. It cannot cancel (4).
A change of energy zero or reference subtraction cannot do so either.
This is a uniform obstruction on the complete single-unscreened-charge
family. It is not an unrestricted confinement theorem or an exclusion
of every screened Coulomb-like sector of the full source.

## Explicit replacement of the energy law

A separate finite construction replaces, rather than supplements, the
unit nonblank-link penalty:

\[
H'=H_{\rm local}+\lambda T+4N_{\rm ref}
+\frac\kappa2\rho^*G\rho,
\qquad \lambda=2^{-34},\quad\kappa>0.
\tag{6}
\]

The executable gives the actual full-carrier column action of (6).
All original matter/scalar modes, all reference occupations and all
119 link states are retained. It uses the original quantum link
coordinates and Yukawa words, including the occupied-link transitions.
The full kinematic operator uses the mean-zero pseudoinverse; its
physical Gauss sector has the neutrality condition (1).

Every local \(Y\) density commutes with every local gauge generator,
and all such densities commute with each other. Hence the added term
preserves every Gauss constraint and commutes with every local density.
The exact continuity equation has the same native hopping commutator,

\[
\dot\rho_x=i\lambda[T,\rho_x].
\tag{7}
\]

The projected field obeys \(D^*\dot E_{\rm long}=\dot\rho\) on the
physical domain. Independent full-source column tests check all 96
cell/current constraints on a physical matter/scalar/link/reference
input. No hopping coefficient, scalar charge or endpoint orientation
is replaced by a Coulomb ansatz.

Equation (6) changes the Hamiltonian and its clock. The original unique
interacting ground, its preparation fidelity, energy bins and transport
probabilities are not transferred to \(e^{-itH'}\). The original
readout result in (2)–(3) and the alternative energy law in (6) are
separate statements.

## Exact static trial energies and their limits

Take the actual native local charged endpoint eigenstate with energy
\(e_{e^c}\), an actual local ground of energy \(E_0\) at every other
site, and any physical \(e^c\) flux path connecting its endpoint to the
reference. The complete local-current intertwiner identifies this
endpoint inside the original CAR/scalar space. All twelve Gauss
constraints hold on every path vertex. On this normalized physical
trial state,

\[
\langle H'\rangle=(V-1)E_0+e_{e^c}+4
+\kappa(G_{00}-G_{0y}).
\tag{8}
\]

The expression is independent of the path length. At \(\lambda=0\)
these are exact static eigenstates. At nonzero \(\lambda\), (8) is an
expectation value because \(\Pi_yT\Pi_y=0\); the trial need not be a
stationary eigenstate. One-edge and three-edge physical paths to the
same destination have different native flux and different original
link costs, but the same (8). Their divergence-free difference is
explicitly retained in the receipt.

The unrestricted local bound \(h\ge-11.71\) also gives

\[
(V-1)E_0-11.71+4+\kappa C_y
\le\inf\Pi_yH'\Pi_y
\le(V-1)E_0+e_{e^c}+4+\kappa C_y,
\quad C_y=G_{00}-G_{0y}.
\tag{9}
\]

The width is independent of distance and volume. This removes the
linear penalty from this constrained static family. It does not
identify the exact alternative constrained minimum, a deconfined
vacuum phase, a scattering particle, or an autonomous spectrum.
In particular closed native flux loops receive no added energy when
their density is unchanged. No transverse/link-loop Hamiltonian or
Lorentz completion is included in (6). The original neutral product
state is not asserted to be its actual ground.

## Green profile and independent coupling

The same native unit charge is placed on the periodic counting cochain
for the static trial family. The finite Green solve gives a potential
and displacement-energy differences; no \(1/r\) potential is inserted.
The computed large-box regime has
\(V_{\rm interaction}\simeq-\kappa/(4\pi r_{\rm sites})\), with
near-inverse-square differences and explicit lattice/image errors.
At grid side 128 and radius eight, the potential-difference amplitude
ratio is approximately 1.00407 and the force amplitude ratio 1.02842.
The finite directional controls and distance slopes are retained.
These compare the static trial/field energy, not the dynamics of an
identified free particle above the alternative interacting ground.

The geometry uses counting incidence and counting link norm. The
native matter kinetic convention has spacing \(\pi\); this calculation
assigns no physical length, GeV conversion or common measured light
speed. The coupling \(\kappa\) is an independently declared energy
coefficient. Changing it from one to two leaves every native current,
Gauss and Ward identity intact while doubling the static pair-energy
and force coefficients. Those identities therefore do not fix a
physical fine-structure constant. Any empirical electromagnetic
normalization is an additional input, not a prediction of this result.
