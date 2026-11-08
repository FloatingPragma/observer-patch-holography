# Free-branch census and physical mass admission

The generated [census](table.md) is a control on one declared source family.
It is **not a physical particle spectrum**. Its physical particle mass table
is empty. The eight entries are the kinetic valleys of the centered
bounded-contact free limit, with both particle-addition and hole branches
retained. The ordered-spin construction is a separate branch whose
relativistic mass readout is undetermined.

The producer replays the retained free-limit and canonical-spin-energy
receipts against their complete producers. It binds their actual proof and
code inputs, reproduces all eight Weyl-symbol identities and checks the
charged-clock bound below. This verifies the declared model results; it
does not supply experimental particle identities or native detector
attachment.

## What the zero invariant parameter means

The [free-limit proof](../continuum_free_limit/README.md) uses the centered
ground of `K=H-(U/2)N`, `c(a)=v/a`, fixed `v>0`, and uniformly bounded
`0<=U(a)<=Umax`. One common limiting state carries all eight valley fields
and their common time and spatial translation groups. Their one-particle
and hole generators satisfy

```
h_b(p)=-v sum_j (-1)^b_j sigma_j p_j,
h_b(p)^2=v^2 |p|^2 I,       E(p)=v|p|,
E^2-v^2|p|^2=0.
```

The last equality gives zero invariant rest-energy parameter for these
free one-particle and hole subspaces, in the model convention `hbar=1`.
A massless excitation has no rest frame. This parameter does not assert a
rest-frame measurement, a calibrated SI mass, or a physical species. It
does not apply to the entire Fock space: multiparticle states have a
continuum of total invariant energies. Ratios of two zero parameters are
undefined.

The free-limit theorem fixes its complete field hierarchy, including the
absence of surviving interactions. Its state differs from the original
uncentered ground and from the ordered-spin reference. These branches
cannot be combined into one spectrum. The source-to-detector preparation
and error bounds needed for an operational particle prediction are separate
premises. A finite `H_B` selection gap is an excitation gap in that finite
problem and supplies no relativistic mass by itself.

## A common clock cannot preserve both scales on the protected trajectory

The [protected-spin theorem](../protected_spin_limit/README.md) declares
`c_E=J E^4`, `U_E=J E^8`, fixed `J>0`, with graph edge count `E>=2`.
Its charged excitation lower bound is

```
Delta_charge >= J (E^8/2-E^5) >= J E^8/4.
```

The second inequality follows from `E^3>=4`. Under any common positive
clock rescaling `K -> q_E K`, both excitation energies and the exchange
coefficient acquire the same factor. If some charged excitation has
rescaled energy at most a fixed `C>0`, necessarily

```
q_E J <= C/(E^8/2-E^5) <= 4C/E^8 -> 0.
```

Thus keeping a charged excitation at bounded energy on this trajectory
forces the fixed-spacing local spin dynamics to freeze on that clock.
Conversely, `liminf q_E J>0` makes every charged excitation escape to
infinite energy. This restricts this specified trajectory; it excludes no
different source family, phase or simultaneous scaling construction.

There is also a state-weight bound on the actual positive, ground-subtracted
source generator `L=K-E_ground`. Since charge commutes with K,
`L>=Delta_charge Pi_(N!=|V|)`. For every normalized density operator with
`Tr(rho L)<=C`, positivity gives

```
Tr(rho Pi_(N!=|V|)) <= C/[J(E^8/2-E^5)] <= 4C/(J E^8).
```

This bounds all charged-sector weight, including coherent or mixed inputs;
it assumes neither a thermal state nor a state-selection rate. The relevant
energy is excitation energy above the ground of this same K. An uncentered
Hamiltonian, a selected finite `H_B` reference or a different clock cannot
replace it. The charge is occupation excess, with no electromagnetic
identification supplied.

Unbounded contact alone is insufficient to avoid a free limit: the
[growing-contact theorem](../strong_coupling_refinement/README.md) also
covers `U(c)=o((log c)^(3/2))` with `c=v/a`. A massive interacting result
therefore requires a justified capacity-to-coupling and length/clock map,
a controlled selected state and a finite nonzero mass-shell component on
the same clock. The cited source families declare their coupling laws;
they supply no theorem selecting a physical massive trajectory from the
observer's capacity law. Choosing a law to fit observed masses would not
provide that missing implication.

## Conserved charge and its creating operator on the same finite source

The two-flavour full-domain C/E extension specifies every raw branch on
the CAR carrier, including collisions. Its basis is the declared per-site
quasiparticle basis. Every branch either holds its input or moves one
occupied mode to an empty endpoint. Thus every Kraus operator commutes
with the finite occupation excess `Q=N-|V|`. The centered Weyl/contact
Hamiltonian also commutes with Q. For its local CAR fields,

```
[Q,c_j^*]=c_j^*,       [Q,c_j]=-c_j,
P_(q+1)c_j^*=c_j^* P_q,       P_(q-1)c_j=c_j P_q.
```

These identities hold on every input charge sector, before selecting a
state. They supply a conserved charge label and an explicit field operator
that creates it. Q is not a real-space topological invariant, a winding
number, or a momentum-space Berry charge. A spin operator `c_up^* c_down`
commutes with Q; excitations inside the protected single-occupancy carrier
therefore retain occupation excess zero. No topological defect label for
those spin excitations is supplied by this argument.

The centered dimer provides a complete same-source witness on all 16
Fock states. Its exact ground has sharp `N=2`, proved by its vector and
Hamiltonian, and every onsite mode has addition and removal weight `1/2`.
The proof does not infer sharp total charge from a one-point mean.
For an infinite gauge-invariant state, the appropriate charge grading is
instead the GNS implementation of U(1), which fixes the cyclic vector.
Applying a creator or remover gives its character `exp(+/-i theta)`;
an infinite total occupation operator is not assumed.

There is a complete donor readout for the creating operator. Start an
extra occupation mode in state 1 and use the isometry

```
V psi=c_j^* psi tensor |0> + n_j psi tensor |1>.
V^* V=(1-n_j)+n_j=I,
(Q+n_donor)V=V(Q+I).
```

The donor's 0 record marks creation; its 1 record is the retained failure
branch. Both have probability `1/2` in the actual dimer reference. The code
completes V to a signed permutation on all 32 system/donor addresses,
separately within each total-charge sector. This is an exact finite
charge-conserving Stinespring construction, including arbitrary spectator
references. No native pulse execution of this donor word is claimed.

Conservation alone does not preserve a prepared reference. The complete
unread C repair changes this dimer's centered ground: its exact ground-return
probability is `33/34` at the retained `t=1/16,U=1`. Consequently this charge
intertwiner does not transfer the Hamiltonian vacuum to an uncontrolled
repair schedule or prove a sharp relativistic particle. Finite preparation,
complete native detector synthesis and the same-history limit remain
separate requirements.

Run `python -m model.particle_mass_readout_gate --check model/particle_mass_readout_gate/table.json`.
