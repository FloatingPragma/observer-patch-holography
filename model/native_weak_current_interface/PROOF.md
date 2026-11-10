# Actual weak currents and the scalar-background transfer boundary

## A reusable same-source canonical current map

The exact map in
`model/coherent_defect_wave/chiral_particle_map/model.py` uses the actual
fifteen native internal states and all thirty spin/internal CAR modes.
Its charged-lepton variables are the L annihilator and the epsilon-conjugated
e^c creator, with the signs fixed by the native Yukawa tensor. Together
with the corresponding colored pairs and the unpaired neutral field,
it is a bijective canonical map on the complete occupation algebra:

```
Psi=A c+B c^dagger,
AA^dagger+BB^dagger=I,       AB^T+BA^T=0.
```

The inverse identities also hold. All scalar and link factors are identity
spectators for this algebraic change of variables. Therefore it can be
applied literally to the unchanged quantum native Hamiltonian and its
actual ground. The state and Hamiltonian are not replaced by a free
vacuum or a neutral scalar background.

For every actual native current G, normal ordering gives the number
current matrix `A G A^dagger-B G^T B^dagger`. All anomalous current terms
vanish, and all twelve normal-ordering constants vanish on the complete
family. The three weak-current matrices have support only on the
designated left components. The chosen generator Q=T3+Y has value -1
on both charged-lepton halves. This is an exact current interface,
independent of scalar state, v or Yukawa coefficients. The physical
electromagnetic interpretation of that chosen generator and Lorentz
chirality of those components require their dynamical and reference
joins. Gauge-covariant fields require a compensating reference or dressing
before they define physical preparations or measurements.
The gauge-singlet character of a dressed fermion is compatible with
nonzero physical weak-current three-point responses. A vanishing gauge
generator on physical singlets is a constraint and does not exclude weak
scattering. The dressed current and its complete source
correlators must be evaluated.

The same package also proves a different, explicitly conditional identity:
with a prescribed constant neutral scalar background, its Nambu Hamiltonian
reduces to Dirac blocks of mass v*y. That identity uses the scalar
background in the Hamiltonian. It is not a reducing identity for the
unchanged quantum scalar/source Hamiltonian merely because the canonical
current map is valid.

## Exact obstruction at the five-state scalar cap

The actual scalar fields are

```
phi_k=|0><1+k|+|3+k><0|,    k=0,1.
```

Exact multiplication gives phi_k^3=0, phi_k^2!=0 and ||phi_k||=1.
These are the literal matrices returned by `chiral_higgs.higgs_memory`.
For a normalized vector eta in this scalar carrier, with any additional
entangled reference retained, write
`epsilon=||(phi_k-v)eta||`. The identity

```
phi_k^3-v^3=(phi_k^2+v phi_k+v^2)(phi_k-v)
```

and the operator norm bound imply

```
epsilon >= |v|^3/(1+|v|+|v|^2).
```

The argument applies to phi_k tensor I on any retained reference and to
its adjoint. A nonzero exact scalar eigenbackground is therefore excluded
at this fixed cap. Stacking phi0, phi0^dagger, phi1 and phi1^dagger gives
a rank-five matrix, so there is no nonzero common zero eigenvector either.
Thus there is no one-dimensional reducing scalar state that replaces
all quantum field matrices by constants on the complete native CAR
domain. Independence of the fermion pair monomials makes this a direct
obstruction to that product-state background intertwiner.

If one tries to match the source-only E_H interval from the actual
Higgs spectral join by setting |v|=E_H at the actual unit Yukawa
coefficient, the residual floor is greater than 0.15199448. This is a
conditional error test of that proposed replacement, not a choice of v
or an empirical calibration.

The obstruction does not exclude nonzero expectations, controlled
collective limits, entangled dressed low-energy codes or Higgs physics.
For an explicit control take
`eta=(sqrt(2)|0>+|1>+|3>)/2`. Then

```
<phi0>=<phi0^dagger>=1/sqrt(2),
<phi1>=<phi1^dagger>=0,
||(phi0-<phi0>)eta||²=1/4,
||phi1 eta||²=1/2,
```

with the same variances for the adjoints. For n independent copies in
eta^tensor n, the averaged fields have these squared residuals divided
by n, since centered cross terms factor to zero. Every member is finite.
This is an explicit controlled expectation/average example, not the
interacting native ground or a stationary scalar-background theorem.
A collective or relative-reference construction must specify its actual
state, clock, coupling and capacity. It cannot silently substitute those
averaged fields for each local source phi_k.

The same-H E_H result supplies the justified dressed alternative here:
the literal quantum H†L operator and the elementary e^c probe have
nonzero cross residue in a common actual eigenspace, with their scalar
operators and all complementary sectors retained. That statement does
not require a scalar eigenbackground or make the composite canonical
on all occupations.

## An actual dressed weak-current operator

Write `J_beta=sum_i phi_i^dagger c_(beta,L_i)` and define

```
K_ij=phi_i phi_j^dagger,       S=tr K,
K_a=tr(sigma_a K),
j_a=(1/2)sum_beta c_(beta,L)^dagger sigma_a c_(beta,L),
W=sum_a K_a j_a.
```

The executed full normal-ordered CAR/scalar polynomial identity is

```
sum_beta J_beta^dagger J_beta=(S/2)N_L+W.
```

All twelve Gauss-current commutators of W vanish. It is a spin scalar,
even under fermion parity and neutral under conserved raw L. The Higgs
factors provide the dressing of the lepton weak-current density.
Consequently W is not the bare Gauss generator, whose action vanishes
on singlets. The Gauss constraint does not imply W=0.

The full norm is calculated on the complete five-scalar by sixteen-state
Fock space of the four L modes. The other twenty-six native fermionic
modes are identity spectators. The actual Hermitian eighty-state matrix
is annihilated by

```
(W+I)(W+I/2)W(W-3I/2)(W-2I)=0.
```

The exact projector at eigenvalue2 has rank two, so ||W||=2 on the full
native carrier. This is an operator bound, with no ground substitution.
The complete effects `E_plus/minus=(I plus/minus W/2)/2` are positive
and sum to identity. Their square-root instrument commutes with every
Gauss generator and parity; a retained neutral pointer records both
outcomes. At a declared apparatus location, the unchanged source h
supplies all intervening evolution. This is a physical dressed-current
readout interface. Its actual eight-cell ground expectation, inelastic
spectral weight or charged-state form factor is not calculated by this
audit. No W-boson identity, Lorentz-vector interpretation or weak
scattering rate follows from its nonzero operator norm.

## Copy symmetry and generation scope

The native current map contains one R15 family. The charged free code
built from e^c states and their conjugates does not add an independent
L representation. Its macrocell copies carry spatial support; their
number is not the number of observed generations.

For three identical independent family copies, the one-excitation
Hamiltonian is the source operator tensor I3. All three copy energies
are equal. Even allowing a family matrix invariant under every copy
permutation gives only
`M=(d-b)I3+b ones(3,3)`, with real d,b when M is Hermitian. The exact
commutator kernel has dimension two, and its energies are
`d+2b,d-b,d-b`. It cannot provide three distinct simple masses while
retaining this permutation symmetry. A declared family-breaking operator
can remove the degeneracy; the copy count supplies neither its entries
nor an electron/muon/tau assignment.

The existing three-family gauge/Higgs calculation in
`model/interacting_extension/higgs_join` explicitly supplies its family
matrix, scalar cap, graph and coefficients. Its finite response is an
admissible separate architecture, not the same eight-cell h and g used
for E_H. These distinctions impose no uniqueness requirement on a
declared family construction.
