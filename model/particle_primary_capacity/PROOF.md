# Capacity of the selected source family

Use the unchanged [primary construction](../particle_common_contract/PROOF.md).
Member n has side `L=2^(2n+1)`, `V=L^3` cells, sixty CAR modes per cell,
three owned link factors per cell, and per-link word degree `ell=n+1`.
Dangling links remain allocated even when the physical Gauss condition
restricts them. Capacity here means `N=ln D` nats for the declared carrier
dimension D. It is not a number of cells or bytes.

## Ambient capacity and the physical constraint

The ambient Fock factor has dimension `2^(60V)`. The native degree-one
link space contains the constant coefficient and the matrix coefficients
of the five mutually inequivalent native irreducible modules and their
conjugates. Their dimensions are 6,3,3,2,1. The exact existing
[degree-one construction](../particle_physical_boundary_state/PROOF.md)
therefore gives

```
d_1 = 1+2*(6^2+3^2+3^2+2^2+1^2) = 119,
N_ambient,0 = 8*(60 ln2+3 ln119) = approximately447.4096105 nats.
```

At arbitrary ell, pad each word with identities to length ell. The 225
native matrix entries, their conjugates and the identity supply at most
451 letters. Their products span the link space, so `d_ell<=451^ell`.
Linear dependencies can make this bound strict. Thus

```
N_ambient,n <= V*(60 ln2+3 ell ln451)
             <= V*(60+27 ell) ln2.
```

The last integer is the existing source-plan qubit bound. Its separate
`log_dimension_upper` also uses `ln2<=1` as a conservative cooling bound;
that bound must not be described as an exact capacity in nats.

Let `d_phys=rank P_G`. An independent lower bound uses the actual traceless
anti-Hermitian native M15 generators. The derived action on the top
exterior power of one complete M15 copy is multiplication by the generator
trace, hence zero. On the declared connected gauge group, both that filled
state and the empty state are singlets. The four Dirac rows therefore
supply sixteen orthogonal local singlet occupation states. Constant link
functions are invariant at both endpoints. Independently choosing these
sixteen states at every cell gives `16^V` orthogonal physical vectors,
including all ghost-endpoint constraints. Consequently

```
4 V ln2 <= N_phys,n = ln d_phys <= N_ambient,n.
```

This is a dimension witness, not a claim that its local span is invariant
under hopping or contains the selected ground. All occupations and parity
sectors remain allocated; the capacity bound does not claim coherent
operational access between every superselection sector. The selected
polynomial density is faithful on `ran P_G`, but its entropy is generally
less than `ln d_phys` and has not been evaluated here.

## Different capacities and execution resources

The contract allocates the full conjugate source carrier as reference,
giving joint ambient log dimension `2 N_ambient`. Its selected purification
is pure; entropy zero of that pure vector is not the capacity of its carrier.
Apparatus selectors, complete detector histories, program tables and
workspace are additional allocations.

The existing [public-record capacity](../observer/CAPACITY.md) is
`ln M0<=ln D` for a specified correctable reachable code. Its two-event
native benchmark has eleven terminal records, hence `ln11` nats, in a
different3354-dimensional allocated carrier. Neither that record count nor
the Wilson ambient bound is a same-carrier cosmological capacity theorem.

A capacity of `10^12` nats corresponds to `10^12/ln2` binary degrees of
freedom. A generic quantum state on that carrier has `exp(10^12)` complex
amplitudes, rather than `10^12` stored amplitudes. The smallest ambient
Wilson member alone has about645.48 binary degrees of freedom. None of
these statements claims dense classical execution of the specified state.
The exact small matrix controls in other packages have their own capacities.

## Capacity dependence of masses

The selected sequence changes volume, spacing, link cutoff and cooling
together. Its physical source extent is `L a=2^(n+1)`, while `a=2^-n`.
An n-dependent response therefore cannot be attributed to capacity alone
without separating these effects or proving a joint scaling law. The
source energy unit and any declared conversion to GeV must also be held
explicit. Changing that conversion with n can manufacture apparent mass
convergence and is not a source result.

Every term of the selected Hamiltonian carries the common factor `1/a_n`.
Writing `H_n=Hhat_n/a_n` gives the exact source-energy relation
`E_n=2^n Ehat_n`. The dimensionless eigenvalues of Hhat_n also change with
the box and link cutoff; they are not known constants. With a declared
time unit tau_n seconds, `E_GeV=hbar_GeV_s*2^n Ehat_n/tau_n`.
For a fixed time unit, a finite nonzero rest-energy limit would require
the relevant dimensionless excitation to scale like a_n. This is a
necessary scaling condition, not a derived mass-capacity law. Keeping one
lattice spacing at a fixed number of metres at every n changes the unit
policy and is different from shrinking spacing in a fixed coordinate unit.

No target capacity is selected here. The level-five ambient upper bound
is about `1.302193221e12` nats, with more than eight billion cells; it is
neither an exact physical capacity nor a computationally practical dense
simulation. The mass task needs controlled response estimates and error
bounds across the family. A cosmological `N` additionally needs its
same-carrier horizon and public-record closure identifications.
