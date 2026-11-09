# Occupation-sensitive electric bounds

The actual primary Wilson endpoint admits a full physical ground floor and
a bound on matter defects in its selected polynomial state. At member
`n=0`, the ground enclosure is approximately
\([-727.498195,-723.7678199]\) in the declared source clock. The selected
state has fewer than five beta defects across eight cells. Consequently a
uniform choice among sixteen matched local neutral pair probes has original
success weight greater than \(11/16\). No postselection normalization enters
these statements. A particle rest mass, an isolated endpoint level and a
refinement-uniform weight bound do not follow.

## 1. Source and inherited electric reduction

The source, product link cutoff, full 60-mode CAR carrier, physical Gauss
restriction, polynomial state, idle conjugate reference and clock are those
of [the primary contract](../particle_common_contract/PROOF.md). Write
\(H=\widehat H/a\), \(a=2^{-n}\), volume \(V=L^3\), internal bond count
\(B\) and plaquette count \(P_{\square}\). The onsite term is
\(3d\Gamma(\beta_D)\) per cell. The defect operator is

\[
 D_x=d\Gamma(\beta_D)_x+30I,\qquad D=\sum_xD_x\ge0.
\]

It counts occupied positive-beta modes and empty negative-beta modes. Its
ambient range is \([0,60V]\). The filled-negative-beta product vector with
constant links has energy \(E_{00}=-90V\), as established in
[the source spectral calculation](../particle_primary_spectral_scale/PROOF.md).
This vector is an operator-analysis reference; the preparation is the
selected \(\rho=Q^2/\operatorname{Tr}Q^2\).

For one internal bond allocate onsite coefficient \(m>0\) at each end and
electric coefficient \(e>0\). The shifted cluster is

\[
 K=m(D_x+D_y)+\frac e2 C_{xy}+T_{xy}.
\]

[The endpoint electric bound](../particle_primary_endpoint_bounds/PROOF.md)
proves, on the full Fock carrier including all four Dirac rows,

\[
 m(D_x+D_y)+T_{xy}\ge-\delta(m)I,
 \qquad\delta(m)=60\left(\sqrt{m^2+\tfrac14}-m\right).
 \tag{1}
\]

Let \(P_0\) project the link to its constant word and let \(Q_0=I-P_0\).
The exact all-word Casimir bound is \(C\ge3Q_0/10\). Under
\(3e/20>\delta\), the nonconstant block of \(K\) is positive. Completion
of its block square and inverse order reduce the constant-link self-energy
to

\[
 P_0TQ_0(Q_0KQ_0)^{-1}Q_0TP_0
 \le\sum_r w_r\mathcal Q_r,\qquad
 w_r=\frac1{eC_r/2-\delta(m)}.                     \tag{2}
\]

All five native representations and their conjugates are retained. Their
dimensions and Casimirs are
\((3,4/5),(2,9/20),(6,21/20),(3,7/10),(1,3/10)\).
Each single-entry packet from a constant link lies in every admitted product
cutoff. Haar orthogonality supplies (2) before the Gauss restriction.

For a block of native dimension \(d\), two active Dirac rows participate in
the rank-two Wilson hopping. With their flavor currents
\(E^x_{ab}=\sum_{i=1}^d a_{ai}^{\dagger}a_{bi}\), \(a,b=1,2\), and the
corresponding currents at \(y\), the exact Haar quartic is

\[
 \mathcal Q=N_x+N_y-\frac2d\sum_{a,b}E^x_{ab}E^y_{ba}
 =N_x+N_y-\frac{N_xN_y}{d}-\frac4d\mathbf S_x\!\cdot\!\mathbf S_y.
 \tag{3}
\]

Its sharp full-occupation norm is \(q_d=2(d+1)\). That cap controls excited
blocks; retaining the vacuum and its crossblock gives a smaller lower shift.

## 2. Active and inactive rows

The active subspace at each endpoint is a two-dimensional eigenspace of
the relevant \(\gamma_\mu\). Since \(\{\beta_D,\gamma_\mu\}=0\), its
orthogonal inactive subspace is its image under \(\beta_D\). Choose paired
bases in these two subspaces. In this basis
\(\beta_D=\left(\begin{smallmatrix}0&I_2\\I_2&0\end{smallmatrix}\right)\).
Both active and inactive rows remain in the local Fock space. Per native
block and site let

\[
 D_x^{(d)}=d\Gamma(\beta_D)_x+2dI.
\]

Its unique zero vector \(\Omega_x\) fills the \(2d\) negative-beta modes.
The active reduced state is the normalized identity on its Fock space of
dimension \(2^{2d}\). Set \(X_{ab}=E_{ab}-\tfrac d2\delta_{ab}I\). The CAR give

\[
 \langle X_{ab}\rangle_{\Omega_x}=0,\qquad
 \langle X_{ab}X_{cd}\rangle_{\Omega_x}
       =\frac d4\delta_{ad}\delta_{bc}.
 \tag{4}
\]

Each vector \(X_{ab}\Omega_x\) is a positive-beta particle and a
negative-beta hole, hence has exactly two defects. Equation (3) becomes

\[
 \mathcal Q=dI-\frac2d\sum_{a,b}X^x_{ab}X^y_{ba}.
\]

For \(\Omega=\Omega_x\otimes\Omega_y\), (4) therefore proves

\[
 \langle\mathcal Q\rangle_\Omega=d,\qquad
 v=(\mathcal Q-d)\Omega,\qquad
 \|v\|^2=1,\qquad
 (D_x^{(d)}+D_y^{(d)})v=4v.                      \tag{5}
\]

In particular \(\mathcal Q\Omega\ne d\Omega\). Replacing (3) by its
vacuum expectation as an operator inequality would omit a nonzero coupling.

Both \(\mathcal Q\) and the onsite operator conserve the total number of
all active and inactive modes separately at each endpoint. In the sector
with total number \(2d\) at each endpoint, the complement of \(\Omega\)
has at least two defects. Every other pair of number sectors has at least
one defect. These are full Fock decompositions, without a Gauss assumption
or a restriction to a fixed active occupation.

## 3. A small local Schur problem

Consider \(L=m(D_x^{(d)}+D_y^{(d)})-w\mathcal Q\), \(w>0\). Its vacuum
diagonal is \(-wd\) and its off-diagonal vector is \(-wv\). On the
balanced-number vacuum complement,

\[
 L_\perp+cI\ge m(D_x^{(d)}+D_y^{(d)})-wq_d+c.
\]

If \(2m-wq_d+c>0\), inverse order and (5) bound the vacuum self-energy by
\(w^2/(4m-wq_d+c)\). Consequently

\[
 c\ge wd,\quad
 (c-wd)(4m-wq_d+c)\ge w^2,\quad
 2m-wq_d+c>0,\quad c\ge wq_d-m                 \tag{6}
\]

imply \(L\ge-cI\) on the full Fock space. The last condition covers every
other number sector. A convenient value is the maximum of \(0\),
\(wq_d-m\), and

\[
 c_4=\frac{wd-b_4+\sqrt{(b_4+wd)^2+4w^2}}2,
 \qquad b_4=4m-wq_d.                            \tag{7}
\]

If the complement-positivity condition in (6) fails, use the same formula
with \(b_2=2m-wq_d\). The two-defect block estimate then applies directly.
The code checks the exact rational Schur inequalities after rounding each
square root upward. Increasing \(c\) preserves them in the admitted domain.

The native blocks act on disjoint even CAR factors, so (6) sums to
\(m(D_x+D_y)-\sum_rw_r\mathcal Q_r\ge-c_{\rm occ}I\), where
\(c_{\rm occ}=\sum_rc_r\). For the full link cluster, apply its Schur
complement to \(K+c_{\rm occ}I\). The positive nonconstant block has an
inverse no greater than the one in (2), so

\[
 K\ge-c_{\rm occ}(m,e)I.                       \tag{8}
\]

The proof treats the electric link quantum mechanically. It does not
replace it by a classical link value or declare the constant-link space
invariant under hopping.

## 4. Spatial allocation and the selected state

Choose independent retained fractions \(0\le r,s<1\) for onsite and
electric energy. With maximum internal vertex degree \(z\), take

\[
 m=\frac{3(1-r)}z,\qquad e=105(1-s).
\]

Here \(z=3\) for the side-two member and \(z=6\) for the larger cubes.
Every bond receives its own electric term and \(mD_x\) at each endpoint.
The unallocated onsite terms are positive because \(m\deg(x)\le3(1-r)\).
Exterior and other unallocated electric terms are positive. Bounding the
magnetic term by \(-P_{\square}/7\) yields the ambient operator inequality

\[
 \widehat H\ge E_{00}-\mathcal L(r,s)+3rD
       +\frac{105s}2\sum_\ell C_\ell,
 \quad
 \mathcal L=Bc_{\rm occ}(m,e)+\frac{P_{\square}}7.  \tag{9}
\]

All terms preserve the physical Gauss space, so (9) holds on that actual
carrier, including boundary constraints. The source reference is idle;
tensoring (9) with its identity preserves the inequality.

Let \(U_0\) be the authenticated full-source variational ground upper from
[the endpoint trial calculation](../particle_primary_variational_readout/PROOF.md).
The selected preparation has
\(\operatorname{Tr}\rho H\le E_0+\epsilon\le U_0+\epsilon\), with the
same member's \(\epsilon=2^{-n-4}\). For \(r>0\), (9) gives

\[
 \operatorname{Tr}\rho D
 \le\min\left\{60V,
 \frac{a(U_0+\epsilon)-E_{00}+\mathcal L(r,s)}{3r}\right\}.
 \tag{10}
\]

No exact-ground preparation or spectral gap is used. In particular, retaining
onsite energy alone, \(s=0\), suffices for (10). This choice supplies no
positive retained electric penalty and cannot be read as \(rH_0\).

For a ground enclosure, (9) also yields
\(E_0\ge(E_{00}-\mathcal L)/a\). A separate ordered second-level floor
adds \(6\min(r,s)/a\), using the actual physical anchor gap. Its positivity
relative to an upper bound would require an independent numerical check.

## 5. An original neutral-pair weight

The native \(Q_6=6\) block is one-dimensional. At each cell choose paired
orthonormal positive- and negative-beta bases, with two matched rows
\(\alpha=1,2\). Define

\[
 B_{x\alpha}=c^{\dagger}_{x,+,\alpha}c_{x,-,\alpha},\qquad
 B_{x\alpha}^{\dagger}B_{x\alpha}
       =n_{x,-,\alpha}(1-n_{x,+,\alpha}).
\]

Each \(B\) is an even contraction, gauge invariant and number preserving.
The CAR imply

\[
 n_-(1-n_+)\ge1-[(1-n_-)+n_+],\qquad
 \sum_{x,\alpha}B_{x\alpha}^{\dagger}B_{x\alpha}
       \ge2V-D_{Q_6=6}\ge2V-D.                  \tag{11}
\]

Uniformly choosing one of the \(2V\) addresses independently of its outcome
therefore gives original averaged weight

\[
 \overline w=\frac1{2V}\sum_{x,\alpha}
      \operatorname{Tr}\rho B_{x\alpha}^{\dagger}B_{x\alpha}
 \ge\max\left\{0,1-\frac{\operatorname{Tr}\rho D}{2V}\right\}.
 \tag{12}
\]

The statement concerns the uniform address average; no spatial invariance
or identical weight at each site is assumed. The complete success and
failure instrument retains the original normalization. A transition band
requires a dynamical commutator or another spectral estimate in addition
to (12).

## 6. Actual endpoint values and scope

At \(n=0\), \(a=1\), \(V=8\), \(B=12\), \(P_{\square}=6\).
With \(r=s=0\), the certified bond loss is approximately
\(0.553420963619\), giving a full physical ground floor
\(-727.498194421\). The authenticated variational upper is approximately
\(-723.767819943\); the enclosure width is below \(3.731\). The inherited
occupation-independent floor is approximately \(-740.229067566\).

With \(r=2/5\), \(s=0\), \(m=3/5\), \(e=105\), the bond loss is
approximately \(0.736465810093\). Equation (10), including
\(\epsilon=1/16\), gives

\[
 \operatorname{Tr}\rho D\le4.991177197<5,
 \qquad \overline w\ge0.6880514252>\frac{11}{16}. \tag{13}
\]

The displayed decimals are outward bounds where used as inequalities; the
receipt retains exact rational values and input hashes. The values use the
actual endpoint \(\lambda=1\), all native modes and the selected polynomial
state. A diagonalization of the full spatial source is not claimed.

The retained full-anchor trials do not establish an isolated endpoint
branch. The onsite-only bound establishes occupation weight without that
premise. At \(n=1\), the same fractions \((2/5,0)\) fail the positive
nonconstant-block condition because the maximum vertex degree is six.
Thus (13) is a finite-member result, with no refinement-uniform extension.
The source-energy origin and clock are unchanged; converting declared
units is distinct from identifying a physical rest-energy branch.

`model.py` regenerates the actual member and variational bound before
issuing a receipt. The tests execute full active/inactive CAR controls,
the independent quantum-link coupling, all-occupation pair inequalities,
exact rounding checks and source/clock mutation rejection. Those finite
controls verify operator steps; their spectra are not substituted for the
primary spectrum.
