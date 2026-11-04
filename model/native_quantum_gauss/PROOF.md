# Native matter with finite local Gauss constraints

## Declared register and exact local algebra

Let V be the actual fifteen-state native matter representation, computed
by `native_matter_attachment` from port-admitted generators and its stated
spin/exterior lift. Its central current separates invariant blocks V_r of
dimensions 6, 3, 3, 2 and 1. Write T for any anti-Hermitian generator. For
an oriented edge x to y choose

```
L_xy = C|0> direct-sum sum_r (V*_{r,x} tensor V_{r,y}),
U_ab = |bar a,b><0|             for a,b in the same block r,
L_x(T) = 0 direct-sum sum_r (-T_r^T tensor I),
L_y(T) = 0 direct-sum sum_r (I tensor T_r).
```

Its dimension is 1+36+9+9+4+1=60. Cross-block U_ab are zero, so transport
preserves each native multiplet. These matrices are part of the simulator design.
They obey the same Lie brackets as T at either endpoint and commute
between endpoints. Their literal matrix commutators are

```
[L_x(T),U_ab] = -sum_i T_ai U_ib,
[L_y(T),U_ab] =  sum_j T_jb U_aj.
```

On any finite graph, attach the physical canonical anticommutation
relations at vertices and define the local generator G_x(T) as the
fermion charge `dGamma_x(T)` plus its incident link generators. Then

```
H_hop = J sum_ab c_xa^dag U_ab c_yb + adjoint,
H_link = epsilon (I-|0><0|),       J real, epsilon >= 0,
[G_v(T),H_hop+H_link] = 0         for every vertex v and native T.
```

To prove the hopping identity, use
`[dGamma(T),c_a^dag]=sum_i T_ia c_i^dag` and
`[dGamma(T),c_b]=-sum_j T_bj c_j`. The two endpoint link commutators
cancel these terms after exchanging dummy indices. This is an operator
identity on the entire finite Fock space, including every occupation
sector. The link projector commutes with its generators. Terms on other
edges obey the same identity, so arbitrary finite sums retain local Gauss.
Contracting link indices around an oriented cycle also gives a
gauge-invariant cycle term by the same cancellation, if that additional
term is included in the declared action.

Every spectral function of H preserves the joint kernel of the local
Gauss generators. Complete measurements whose Kraus operators commute
with those generators preserve it on every retained outcome, even when
the measurements exchange energy with the field. This supplies a finite
local Ward identity for this explicit extension. The `U_ab` are raising
operators and are not unitary Wilson links. No Yang--Mills continuum,
physical chiral measure or regulator anomaly cancellation follows.

Finite link registers with continuous gauge covariance have established
precedents in [Chandrasekharan and Wiese](https://arxiv.org/abs/hep-lat/9609042).
The argument here proves the displayed register's identities directly;
no continuum result from that separate construction is imported.

## A nontrivial same-Hamiltonian transport and readout

Supply one static dual reference V* at endpoint y. It belongs to the
total Gauss charge. In each block r of the one-fermion sector define
unnormalized states, with both sums restricted to that block:

```
A = sum_b |y,b> |bar b>_reference |0>,
B = sum_ab |x,a> |bar b>_reference |bar a,b>_link.
```

They have squared norms d and d^2, d=dim(V_r). Both are annihilated by every
local G_x(T), G_y(T): each fundamental factor is contracted with its dual
at the same vertex. The source checks this for all twelve actual native
generators, and detects failures if the reference, link, or orientation
is omitted. Direct operation gives

```
H A = J B,
H B = d J A + epsilon B.
```

On the orthonormal basis A/sqrt(d), B/d, the Hamiltonian is

```
H_reduced = [[0, J sqrt(d)], [J sqrt(d), epsilon]],
Omega^2 = epsilon^2 + 4d J^2.
```

This is an exact reducing two-state subspace with an explicit static
reference. Its level separation Omega is a finite transport frequency,
not a particle rest mass. Starting in A/sqrt(d), continuous evolution
has left-site probability
`4d J^2 sin^2(Omega t/2)/Omega^2`.

The executable certificate uses the exact Cayley step
`(I+i tau H/2)^(-1)(I-i tau H/2)` followed by the complete binary
position instrument. An initial complete block measurement on the
maximally entangled matter/reference singlet selects r with probability
dim(V_r)/15. The certificate retains that record and all sixteen four-round
histories in each of the five blocks for J=epsilon=1 and tau=1/8; the common
hopping coefficient is not changed between multiplets. Each Cayley step preserves the energy of H.
Each position projection commutes with all Gauss generators; it generally
changes the hopping energy. For every round the certificate calculates
that detector work as the energy after readout minus the energy before
it. The sum equals the total field-plus-link energy change. Every such
change caused by the readout is included in this work term. A Hamiltonian
model of the detector reservoir or its reset cost is not supplied.

The coherent step, detector, static reference, clock parameter and link
energy are declared services. They have no claimed bare C/E compilation.
The native internal currents are shared with the compact-defect
construction; the physical states and interaction law are different.
An undressed spatial band is generally changed by independent local
current actions. Relating this Gauss sector to those localized modes
requires an explicit dressing and a compatible interacting Hamiltonian.
Finite non-Abelian conservation therefore closes this register's local
current premise without closing the observed particle, chiral continuum
or responsive-gravity requirements.
