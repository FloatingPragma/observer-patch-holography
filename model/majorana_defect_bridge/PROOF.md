# Same-current Majorana masses and additive defect charges

The declared native twelve-current architecture, its five-mode spin lift,
the nonvacuum even exterior sector, the family factor C3, the canonical
Weyl kinetic form and the complex Higgs doublet are inputs. This construction
extends their conditional charged action by a dimension-five Weinberg
coefficient. The six complex entries of that coefficient are inputs. The
result is an exact field-action and finite matrix theorem, with a separate
conditional constraint on additive defect labels.

## The native-current LH line and its full family lift

Let g_a be the fifteen-row matter currents and h_a their two-row Higgs
currents. The charge -3 subspace of the primitive q=6Y operator has dimension
two. A bilinear `chi^T T H` uses a fifteen-by-two matrix T supported there.
It is invariant precisely when

```
g_a^T T + T h_a = 0, for every one of the twelve currents.       (1)
```

`native.dictionary` constructs the native spin/exterior currents, solves
all 360 coefficient equations in the four entries, and obtains rank three.
Thus the LH invariant line is one-dimensional. Its neutral column has unit
norm. The occupation-subset permutation P obeys `g_old P=P g_native` for
every port. It sends this LH tensor to the retained Weinberg tensor up to
a sign. The sign cancels in its square. This comparison uses the actual
Higgs/current dictionary, without supplying a species table to the solve.

For arbitrary H put `ell=T H`. Given any symmetric complex three-by-three
matrix kappa, define the all-left kernel

```
K_nu(H)=(ell ell^T) tensor kappa.                             (2)
```

For `delta H=h_a H`, equation (1) gives `delta ell=-g_a^T ell`.
Consequently

```
(g_a tensor I3)^T K_nu + K_nu (g_a tensor I3)
  + (delta ell ell^T+ell delta ell^T) tensor kappa = 0.        (3)
```

The code checks these identities for both symbolic Higgs components and
six independent complex family coefficients, then checks the complete
45-row tensor transport with `P tensor I3`. Fermi antisymmetry of the
two-spin epsilon contraction requires the internal kernel to be symmetric.
Family currents are the identity, so the remaining coefficient space is
exactly `Sym3(C)`, of complex dimension six. At `H=(0,w)`, the neutral
column occupies the existing left-neutrino slot and has electric charge
zero. No independent right-neutrino field is added.

The convention is `L_mass=1/2 chi^T(epsilon tensor K)chi+h.c.`. In (2),
kappa has inverse-energy units and includes any Wilson coefficient divided
by its effective-theory scale. This is a supplied dimension-five interaction,
outside the minimal dimension-four charged action.

## Actual Takagi, CAR and momentum maps

Choose a positive reference norm w0 and supplied exact unitary factors whose
product is U. For nonnegative masses m_i define

```
M0=conjugate(U) diag(m1,m2,m3) U^dagger,
kappa=M0/w0^2,
M_nu(w)=(w/w0)^2 M0.                                       (4)
```

Exact unitarity implies `U^T M0 U=diag(m_i)`. Every symmetric complex
matrix admits a Takagi factorization of this form; the executable interface
accepts the factors as exact supplied data and verifies their unitarity and
the full resulting congruence. Zero and repeated masses are retained.
A repeated eigenspace does not select a unique mixing basis.

Add (2) to the actual charged 45-row Higgs action, with its three arbitrary
complex Yukawa matrices and the prescribed conjugate-Higgs channels. In
the neutral background the charged rows and the neutrino rows are disjoint.
The embedded W45 which applies U only to the neutrino family slot gives

```
W45^T K45 W45 = K_charged + diag_neutrino((w/w0)^2*m_i).      (5)
```

`takagi_block_certificate` checks this full matrix equality, including the
unchanged charged block. In particular,

```
rank K45 = 6 rank Yu + 6 rank Yd + 2 rank Ye + rank M0.       (6)
```

The neutrino contribution ranges from zero to three. With the calibrated
charged rank 42, these four cases give ranks 42,43,44,45. A zero lightest
neutrino mass is not silently replaced by a positive one.

The physical ninety-row creation coefficient is
`Delta=epsilon tensor conjugate(K45)`, which is antisymmetric. The
annihilator basis map is `W90=I2 tensor W45`. Its creation-pairing map is

```
W90^dagger Delta conjugate(W90)
   = epsilon tensor conjugate(W45^T K45 W45).                (7)
```

The code checks (7) on the complete ninety rows and checks commutation
with all three kinetic principal matrices. It uses the actual source
matrix, rather than substituting a list of scalar masses for the action.

For `h(p)=-v sigma.p tensor I45`, retain the full opposite-momentum Nambu
symbol on `(chi(p),chi^dagger(-p))`:

```
B(p)=[[h(p), Delta], [Delta^dagger, -h(-p)^T]].              (8)
```

The identities `sigma_j epsilon=-epsilon sigma_j^T` cancel the off-diagonal
blocks of B squared. Since K is symmetric, the squared diagonal blocks are

```
v^2|p|^2 I90 + I2 tensor (conjugate(K) K),
v^2|p|^2 I90 + I2 tensor (K conjugate(K)).                   (9)
```

The symbolic control retains all six complex entries of a generic symmetric
family matrix and arbitrary real momentum. The identity is independent of
the internal dimension and applies to the full K45 in (8). It proves
`E_i(p)^2=v^2|p|^2+((w/w0)^2 m_i)^2` in the neutrino block. Each nonzero
signed frequency has Nambu multiplicity two per Takagi column. A massive
Majorana field has two physical spin states. A zero column is a massless
left Weyl field with its conjugate helicity state. At p=0 its four Nambu
zero rows encode two complex CAR zero modes, not four particles.

Independent controls diagonalize the entire 180-row symbol with exact
complex input matrices and execute the complete 64-state CAR Hamiltonian
for three mixed neutrino families at zero momentum. The phenomenologically
small neutrino block is also checked separately after a common rescaling,
so the large quark masses do not conceal its frequencies in a floating
eigensolver. These controls do not execute a microscopic native interaction
word or construct an interacting continuum quantum field theory.

## The calibrated instance and the actual weak current

The retained [calibrated neutrino packet](../common_source/calibrated_neutrino_parameters.json)
supplies normal-ordering low-energy oscillation splittings and mixing angles.
It chooses the lightest mass and both Majorana phases to be zero. The
charged and Higgs parameters use the separate calibrated common-scale
action. Setting kappa by (4) with its Higgs norm is a phenomenological tree
matching prescription. The transport of the neutrino effective coefficient
to that common scale, its loop error and its ultraviolet origin are not
certified by this prescription.

The bridge reconstructs exactly the retained full K45 with the optional
Weinberg coefficient. Its masses are

```
m1=0,
m2=sqrt(0.00007537) eV,
m3=sqrt(0.002511) eV.                                     (10)
```

The receipt supplies exact rational enclosures obtained by integer square
root and checked by squaring their endpoints. These are enclosures of the
supplied action parameters, not observational confidence intervals or
predictions of pole masses. Scaling H by t multiplies charged action masses
by t and the Majorana coefficient by t squared with fixed kappa.

The retained raising current has identity family block from charged lepton
to neutrino in the weak basis. Applying the actual neutrino mass-basis map
therefore gives `U_PMNS^dagger` in that block. Every calibrated PMNS entry
is nonzero. To certify the support without a floating threshold, set
`a=sin^2(theta12), b=sin^2(theta23), c=sin^2(theta13)`. The unmixed entries
have their positive product squares. Each mixed entry is `A+B exp(i delta)`
with |A|,|B| at most one; if their rational squares differ,

```
|A+B exp(i delta)|^2 >= (|A|-|B|)^2
  = (A^2-B^2)^2/(|A|+|B|)^2 >= (A^2-B^2)^2/4 > 0.         (11)
```

The nine exact lower bounds appear in the receipt and hold for every Dirac
and Majorana phase with those angles. No arbitrary unitary matrix is assumed
to have this support. The topology code accepts a supplied support graph
and includes a disconnected countercontrol.

## Additive charges for Majorana families

Assume an abelian additive conserved label, inverse antiparticle labels,
and a definite label for each listed mass eigenfield. Each admitted nonzero
interaction preserves that label. These are conditional defect hypotheses,
not consequences of the finite current matrices. A non-diagonal action of
a conserved charge in a degenerate family space requires a separate analysis.

For each admitted charged-current entry,

```
q(W-)-q(e_alpha)+q(nu_i)=0.                               (12)
```

On a connected charged-lepton/neutrino bipartite graph, all neutrinos have
one charge N, all charged leptons one charge E, and `W-=E-N`. Full nonzero
quark mixing similarly gives `down=D` and `up=D-E+N`. A massive Majorana
field is its own antiparticle, hence

```
2N=0.                                                    (13)
```

Neutral conjugate-pair vertices set the Z and Higgs labels to zero. Thus,
when all three neutrinos are massive, the universal perturbative relation
group is `Z^2 direct_sum Z2`. In an integer winding target, (13) forces
N=0. An order-two label can survive in a group with two-torsion.

Under the requested rule excluding massless fields from defect assignments,
any zero-mass neutrino has charge zero. Connected mixing then forces N=0
for the other two as well, even if they have positive Majorana masses.
The calibrated branch therefore has universal group `Z^2`. It would be
incorrect to assign different defect winding to its two massive columns
merely because the lightest column was excluded from a particle list.

The optional three-generation electroweak process `(udd nu)^3` adds the
full field relation

```
3*up+6*down+3*N=9D-3E+6N=0.                              (14)
```

With (13), this reduces to `3(3D-E)=0`. The all-positive branch has
`Z direct_sum Z3 direct_sum Z2`, equivalently Smith invariants
`Z direct_sum Z6`. The branch with an excluded massless neutrino has
`Z direct_sum Z3`. The unreduced 6N coefficient in (14) preserves the
fermion-parity assignment in the all-positive branch. Constructing this
electroweak process in the microscopic source is an additional premise.

`topology.py` checks the integer presentations and Smith normal forms with
all fifteen field labels retained. Independent finite cyclic-group controls
check assignment counts. A diagonal PMNS graph admits three independent
order-two neutrino labels, so the connected-graph hypothesis cannot be
discarded. Neither these relations nor a Majorana mass select a topological
texture, its energy, or a physical defect-to-field intertwiner.

## Parameter and physical boundary

The concrete [fixed-repair source](../common_source/native_repair_mass.py)
does have a coefficient readout independent of the empirical mass packet.
Its [same-current attachment](../common_source/native_repair_mass_attachment.py)
uses the identical real symmetric Y on all three charged Higgs channels.
The [retained native readout](../common_source/native_repair_mass_results.json)
keeps all native syndromes and pointer outcomes in nine fixed probe settings.
This package verifies that receipt's full source manifest, recomputes all
128 native C/E record weights, and obtains advance probability p=1/12.
For its specified order-five rotation U and seed projector P, it reconstructs

```
pi_k=(1-p)*p^k/(1-p^5),
Y=sum_(k=0..4) pi_k U^k P (U^T)^k.                       (15)
```

The exact stationary transition, rotation and projector identities are
checked. The complete native-frame K45 reconstructed from this Y agrees
with the retained source hash. Its positive spectrum has trace one and
characteristic polynomial

```
x^3-x^2+(1392/22621)*x
 -100224/511709641+418176*sqrt(5)/12792741025.               (16)
```

The retained rational sign brackets are rechecked. They give ordered ratios
approximately 30.8406 and 452.6194, disjoint from the separately supplied
comparison windows [206,207] and [3470,3485]. These broad windows are not
experimental standard uncertainties. The common Y gives identical charged
channel masses and CKM=I in a common mass basis. Its neutrino rows are zero:
this source specifies kappa=0 and three massless neutrinos. No expensive
detector execution is repeated by this algebraic replay.

Thus reading the complete specified repair/feedback law determines that
particular coefficient and its excluded comparison spectrum. The separate
empirical neutrino packet is not selected by this law. The exclusion does
not cover other interacting source laws, loop corrections, or an unspecified
physical identification beyond the retained comparison's hypotheses.

The exact current equations select the LH tensor line, not its six complex
family coefficients. For example, with the same currents, Higgs norm and
canonical kinetic form, diagonal kappa with entries `(1,2,3)` and `(1,2,4)`
both satisfy (3) and produce different neutral spectra. The analogous
charged-Yukawa freedom is established by the same-current matter attachment.
Reading a fully specified parameter packet fixes those values as source
inputs. It does not turn an empirical calibration into a prediction or
identify that declared source with the physical source.

This package proves the finite action joins, quadratic mass shells and the
stated conditional additive-charge constraints. It supplies no repair-law
selection of kappa, no measured absolute neutrino mass or Majorana phases,
no error-controlled interacting pole masses, and no stable topological
particle assignment.
