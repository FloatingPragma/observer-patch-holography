# Electric sectors of the complete native link channel

The source is the configured finite quantum link in
`model/native_quantum_gauss/algebra.py`, with its fifteen matter labels,
static dual reference and all five native current blocks. Its link is
`C |vac> direct-sum_r V_r* tensor V_r`, of dimension 60. The complete
one-matter/reference carrier has dimension 27000. The admitted input space
is the ten-dimensional invariant Gauss-singlet span of the literal states
A_r and B_r returned by `reference_states`. The executable applies all
12 current generators at both endpoints, and the ambient Hamiltonian,
to these states before reducing the calculation. This is a restriction to
an invariant subspace of the configured source, not an approximation to
an infinite SU(2) rotor. No omitted rotor tail is claimed or needed.

## Actual electric projectors and complete records

Write `a_r=A_r/sqrt(d_r)` and `b_r=B_r/d_r`. Their norms and ambient action
are `||A_r||²=d_r`, `||B_r||²=d_r²`, `H A_r=B_r`,
`H B_r=d_r A_r+B_r`. Thus in each normalized pair

```
H_r = [[0,sqrt(d_r)],[sqrt(d_r),1]].
U_r = (I+i h H_r/2)^(-1) (I-i h H_r/2),   h=1/8.       (1)
```

This is the exact Cayley service and default duration of the parent
complete-history calculation. It is not identified with `exp(-ihH)`.
The full channel is `E(rho)=sum_z P_z U rho U* P_z`, where z denotes
endpoint a or b and `P_z=sum_r |z_r><z_r|`. The two branches include
return/no-motion records. `sum_z (P_z U)* (P_z U)=I` proves trace
preservation on all ten-dimensional inputs with every idle reference.
The program retains all eight three-round records for each a_r input
and replays the literal parent instrument in its independent Gram basis.

The link's left weak Casimir is inherited from the actual native current
representation. It is zero on its vacuum and on the d=1,3,3 blocks, and
3/4 on the d=2,6 blocks. Consequently the genuine weak electric
half-spin projector on the admitted space is

```
P_half = |b_(d=2)><b_(d=2)| + |b_(d=6)><b_(d=6)|.
P_zero = I-P_half.                                      (2)
```

These are SU(2) Casimir sectors, with no finite-rotation-group or classical
distance identification. The link centre shift is `A=I-2P_half`.

## Failure of a two-sector quantum transition law

Because P_half commutes with both position projectors,
`E*(P_half)=U* P_half U`. In a weak-doublet block its off-diagonal
entry is `conjugate(U_10) U_11`, which is nonzero. Therefore it cannot
be a linear combination of P_zero and P_half. The stored exact witness
uses `(a_2+b_2)/sqrt(2)` and `(a_2-b_2)/sqrt(2)`: their initial electric
probabilities coincide, but their next half-spin probabilities differ.
Even diagonal inputs do not close on these two labels: a weak-block a
state and a weak-singlet-block a state both have j=0, while only the
former can enter j=1/2. Discarding coherences alone cannot supply the
missing family information.

## Closed family and every stationary density

For the declared dephased family `rho=sum_r x_r |a_r><a_r|+y_r |b_r><b_r|`,
the exact transition law in each block is

```
K_r = [[1-p_r,p_r],[p_r,1-p_r]],
p_r = 16 d_r h² / ((d_r h²+4)²+4h²)
    = 1024 d_r / ((d_r+256)²+256).                       (3)
```

Every branch multiplicity is included in the ambient `sqrt(d_r)` matrix
element; no additional dimension weight is inserted. All p_r lie strictly
between zero and one. Every stationary diagonal density therefore has
`x_r=y_r=w_r/2`, with arbitrary nonnegative `w_r` summing to one.
Its endpoint imbalance decays by `(1-2p_r)^n` after n attempts.

The full quantum stationary family can also be solved exactly. Every
channel output is endpoint-diagonal. For the r,s block, its two stationary
entries obey the homogeneous two-by-two system

```
([[U_r00 conjugate(U_s00), U_r01 conjugate(U_s01)],
  [U_r10 conjugate(U_s10), U_r11 conjugate(U_s11)]]-I) [x,y]^T=0. (4)
```

All 25 systems are checked over exact algebraic numbers. The rank is one
if d_r=d_s, with kernel span(1,1), and two otherwise. There are three
single-copy dimensions d=1,2,6 and two distinct d=3 blocks. Hence every
stationary density is `direct-sum_d M_d tensor I_endpoint/2`, where
M_1,M_2,M_6 are nonnegative scalars, M_3 is an arbitrary positive
2-by-2 matrix on the two d=3 copies, and `sum_d Tr M_d=1`. The fixed
operator space has complex dimension seven; its densities have real
affine dimension six. Coherences between the two distinct d=3 blocks
are retained. The channel does not select their weights or coherences.

## Coupling-chain consequence

For every stationary density in (4),

```
<A> = 1-w_(d=2)-w_(d=6),       0 <= <A> <= 1.           (5)
```

A declared maximally mixed preparation on the admitted ten-dimensional
space is stationary and gives `p_half=1/5`, `<A>=3/5`. The parent's d=6
preparation a_6 converges instead to `p_half=1/2`, `<A>=0`. These different
values follow from actual source preparations, with no measured target.

The seam calculation supplies the explicitly imported comparison formula
`v=E_cell b <A>`. Its dimensionless input is therefore 3/5 for the tracial
preparation and zero for the limiting native d=6 preparation. No physical
E_cell, b, electroweak scale, thermal Casimir law or running coupling is
selected. A claim of a unique seam measure for this channel is refuted;
a complete quantum characterization and conditional chain input are proved.
