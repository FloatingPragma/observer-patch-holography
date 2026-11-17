# Local bounded energy and geometry coupling

This proof fixes the finite construction before its witness. It is not a
derivation of a physical gravitational field or a continuum limit.

## Complete local Hamiltonian

Let a finite graph have n vertices. On one finite three-label coordinate at
each vertex define Phi_v=diag(-1,0,1). These Hermitian geometry operators
commute. For the weighted graph Laplacian L plus specified nonnegative
diagonal boundary/regulator, define

    H_g=j sum_v L_phi,v+(kappa/2) sum_vw Phi_v L_vw Phi_w,
    N_v=I+g Phi_v,
    H=H_g tensor I+sum_v N_v tensor(h_v-e_v I)+(sum_v e_v)I.

Each h_v is a Hermitian operator on the **same complete source**; sum_v h_v
is the supplied full source Hamiltonian. The graph, shard allocation,
regulator, kinetic, stiffness, lapse coefficient and reference energies are
declared inputs. With j>=0, kappa>0 and positive L, H_g is nonnegative.
For 0<=g<1 the lapses have eigenvalues 1-g,1,1+g and are positive. The finite
joint H is bounded and Hermitian even when the source shards do not commute.
It therefore generates unitary evolution and conserves total energy.

If [G_a,h_v]=0 for **every** source Gauss generator and every shard, then

    [I_geometry tensor G_a,H]=sum_v N_v tensor[G_a,h_v]=0.

Commutation of only the shard sum with G_a would be insufficient. The dense
builder checks the stronger hypothesis; the full-source callback caller
must prove and bind it for its actual term allocation. No projection onto
Gauss agreement, a low-energy sector or a selected particle branch is used.
The source scalar is retained through the shards. This coupling makes the
allocation of that scalar among vertices part of the local law.

For h_v -> h_v+c_v I, e_v -> e_v+c_v, the operator changes only by
(sum_v c_v)I. At fixed references it instead changes by sum_v c_v N_v,
which can be nonzero when sum_v c_v=0. The zero-sum control therefore detects
an otherwise hidden local reference change while keeping bare H_source fixed.

## Locality and complete columns

The path-Laplacian kinetic changes exactly one Phi_v label by one adjacent
step. The geometry potential couples adjacent graph coordinates and local
pinning only. A lapse N_v is diagonal in the local geometry coordinate and
couples only its assigned source shard. The spatial support of h_v must be
declared by the source builder; this module cannot make an extended source
operator local by renaming it a vertex energy.

For a complete basis |a,s>, the matrix-free routine computes

    sum_b (H_g)_ba |b,s>
    +sum_v,b,t (N_v)_ba (h_v)_ts |b,t>
    -sum_v,b e_v(N_v)_ba |b,s>+(sum_v e_v)|a,s>.

It calls each shard oracle once, merges equal labels symbolically, and drops
only exact zero coefficients. All source labels are arbitrary hashable
objects. Thus the formula holds for nested keys as well as native occupation
and flux labels. Complete source-column correctness is the caller's
obligation, not an inference from the small diagnostic.

The complete Hamiltonian uses the local energies individually:
`H=sum_v[A_v+M_v(k_v-F)+F]`, with `k_v=g_v+N_v(h_v-e_v)+e_v`.
The exact shard `g_v` keeps its vertex path kinetic and regulator/boundary
potential, and one half of each incident gradient energy. Summing these
shards recovers `H_Phi`. The [complete column proof](JOINT_GEOMETRY_PROOF.md)
states the local tensor action and disconnected-component factorization.
