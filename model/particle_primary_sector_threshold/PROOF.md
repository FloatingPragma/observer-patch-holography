# Positive thresholds in physical occupation sectors

The unchanged eight-cell primary Wilson source has explicit nonempty
physical sectors whose entire energy support lies strictly above its full
ground. Two sector pairs have odd fermion parity and one pair has even
parity. Their lower support edges and variational upper bounds use the
actual interaction coefficient one. A complete preparation and energy
clock makes each sector measurable while retaining the original source
state, reference and a compensating opposite sector in physical memory.
No isolated line, relativistic spin, physical momentum or rest mass is
inferred from these finite thresholds.

## 1. Conserved occupations and the physical witnesses

Use member zero of the [primary contract](../particle_common_contract/PROOF.md),
including all sixty CAR modes per cell, degree-one quantum links, physical
Gauss support, selected polynomial density, idle conjugate reference and
`exp(-itH)`. The [native current calculation](../particle_primary_channel_labels/PROOF.md)
fixes the five conserved occupations in charge order

\[
 q=(-4,-3,1,2,6),\qquad d=(3,2,6,3,1),\qquad
 N^0=2Vd=(48,32,96,48,16).
\]

The actual exterior-basis central generator has eigenvalues `-q/2`;
its sorted block order is reversed before using these labels. This chart
conversion is checked against the source Casimirs and dimensions. The
beta-filled constant-link anchor has occupations \(N^0\) and energy
\(-720\). It is an analysis vector, not the selected preparation.

For \(\Delta=N-N^0\), write the number of positive-beta particles and
negative-beta holes in block r as \(n_r^+,n_r^-\). Then
\(D_r=n_r^++n_r^-\), \(\Delta_r=n_r^+-n_r^-\), so the exact operator
inequality on the fixed-occupation sector is

\[
                         D_r\ge|\Delta_r|I.              \tag{1}
\]

The Hamiltonian and physical Gauss projection preserve every \(N_r\).
The following six sectors are therefore reducing subspaces of that same
Hamiltonian:

| Occupation increment, with either sign | Local gauge contraction | CAR parity |
| --- | --- | --- |
| `(0,2,0,0,1)` | Two weak doublets contracted by the weak epsilon tensor, times the charge-six singlet. | Odd |
| `(1,0,0,2,0)` | One charge-minus-four and two charge-two color antifundamentals contracted by the color epsilon tensor. | Odd |
| `(0,0,0,3,-1)` | The charge-two color determinant and a charge-six hole. | Even |

All have zero total native gauge charge \(q\cdot\Delta\). This condition
alone does not prove a gauge singlet. The code constructs the actual
local particle/hole occupation basis and solves the simultaneous kernel
of all twelve native generators. The three positive-sign kernels have
dimensions one in bases of sizes 1, 9 and 1. Their squared integer-vector
norms are respectively 1, 3 and 1. Negative-sign kernels are independently
constructed and checked. These are exact statements on the source
representation, with both inactive Dirac rows retained.

Place the local singlet at the origin corner, the filled-beta vector at
every other cell, and constants on all links. Each local gauge generator
annihilates the resulting normalized physical vector \(\phi_\Delta\).
It has exactly \(d_\Delta=\sum_r|\Delta_r|\) defects and
\(\langle H\rangle=-720+3d_\Delta\). Hopping and magnetic expectations
vanish by Haar integration. The global source parity is
\((-1)^{\sum_r\Delta_r}\), since the anchor has 240 occupied modes.
These are fermionic or bosonic source parity labels, without an assigned
relativistic particle species.

## 2. A full-space sector lower bound

Reuse the [full-occupation electric Schur inequality](../particle_primary_occupation_bound/PROOF.md)
and [native-specific allocation](../particle_primary_branch_resolution/PROOF.md).
For retained onsite fractions \(r_j\), put

\[
 m_j=1-r_j,\quad
 d_* =4\sum_j d_j(\sqrt{m_j^2+1/4}-m_j),\quad
 w_j=(105C_j/2-d_*)^{-1}.
\]

Every square root is enclosed upward by a rational. The full nonconstant
link margin \(63/4-d_*\) is positive for the declared allocations.
The full two-site Fock inequalities
\(m_j(D_{j,x}+D_{j,y})-w_j\mathcal Q_j\ge-c_jI\) retain the
vacuum crossblock, all inactive rows and every occupation sector. The
charge-six block uses exact rational positive pivots on all 256 states;
the other blocks use the inherited occupation-sensitive inequalities.
No small compressed matrix substitutes for this operator bound.

Summing the twelve interior bonds, retaining the unallocated onsite
energy and bounding all six plaquettes gives

\[
 H\ge -720-12\sum_jc_j-6/7+3\sum_jr_jD_j.
\]

Using (1) on a fixed physical sector proves

\[
 H|_\Delta\ge L_\Delta I,\qquad
 L_\Delta=-720-12\sum_jc_j-6/7+3\sum_jr_j|\Delta_j|.        \tag{2}
\]

The retained fractions, in the stated native order, are

\[
\begin{array}{c|c}
(0,2,0,0,1)&(0,73/100,0,0,72/100)\\
(1,0,0,2,0)&(52/100,0,0,78/100,0)\\
(0,0,0,3,-1)&(0,0,0,81/100,66/100).
\end{array}
\]

They are choices within an operator inequality, not source-action changes
or measured-mass inputs. The exact rational certificate establishes
\(L_\Delta>U_0\) for each row, where \(U_0\) is regenerated from the
full-source variational producer. Consequently no full-source ground
state lies in any of these six sectors. The lower bound applies to
every energy in a sector, not only to the trial vectors.

The anchor-complement shortcut has a different scope. The
[retained-penalty obstruction](../particle_primary_branch_resolution/PROOF.md#7-an-exact-obstruction-to-retained-penalty-isolation)
shows that optimizing its native fractions cannot isolate the whole
anchor complement or place all A1u states above the authenticated ground
upper. The fixed occupation increments in (1) add actual sector
information absent from that comparison. Equation (2) does not evade
its conclusion for the A1u sector or prove a global vacuum gap.

## 3. Exact interacting variational uppers

The localized physical vector has a source-derived dressing in its same
occupation sector. Let \(h_{r,k}=P_{D=d_\Delta+k}T_r\phi_\Delta\), where
\(T_r\) is the actual hopping in native block r. Direct CAR action and
Haar orthogonality give only \(k=0,2\), with

\[
 \|h_{r,0}\|^2=\tfrac34|\Delta_r|,
 \qquad \|h_{r,2}\|^2=12d_r-\tfrac34|\Delta_r|.            \tag{3}
\]

The code evaluates every Wilson axis, native matrix entry, forward and
reverse term on the complete two-site vectors. Each entry has Haar
square norm \(1/d_r\). The three bonds touching the origin supply the
changed norms; the other nine retain their exact vacuum value \(d_r\).
Both signs of each increment are evaluated separately. The packets are
orthogonal by their native character, bond and defect count, and have
relative unperturbed energy
\(105C_r/2+3k\).

Hopping matrix elements between two such packets vanish because every
hop changes the parity of the total matter number on one bipartite
sublattice. The anchor has definite sublattice parity and every hopping
packet has the opposite parity. Magnetic matrix elements between hopping
packets vanish because a plaquette has at least two edges absent from
the two one-link packets; their single native characters have zero Haar
mean. The mixed hopping/loop terms vanish for the same unused-edge reason.

Include the five real native plaquette packets \(\ell_r=V_{B,r}\phi_\Delta\).
Their common norm squared is \(g=6/22050\), and their relative matrix is

\[
 A=\operatorname{diag}(210C_r)-\mathcal N/420,
\]

where \(\mathcal N\) is the exact cubic character integral from the
[variational calculation](../particle_primary_variational_readout/PROOF.md).
The matter witness factors out of these link-only inner products. All
packets lie in the actual degree-one product cutoff and physical carrier.
They supply the exact 13-dimensional compression, with its full-space
variational meaning; invariance of this small space is not assumed.

Choose the explicit rational trial

\[
 \psi_\Delta={\phi_\Delta-\sum t_{r,k}h_{r,k}-\sum u_r\ell_r
                 \over\sqrt{1+\sum n_{r,k}t_{r,k}^2+g\sum u_r^2}},
 \quad t_{r,k}=(105C_r/2+3k+4)^{-1},\quad u=(A+4I)^{-1}\mathbf1.
\]

Its exact Rayleigh quotient \(M_\Delta\) bounds the sector bottom from
above. The coefficients specify a trial, without changing the source or
claiming optimality. With the authenticated global ground floor \(L_0\),

\[
 L_\Delta-U_0\le\min\operatorname{spec}(H|_\Delta)-E_0
                         \le M_\Delta-L_0.              \tag{4}
\]

The upper quotient is evaluated directly, without a perturbative
extrapolation. Equality of the two signed sectors' exact spectra is not
assumed; their equal displayed trial quotients follow from independent
packet calculations.

The exact certificate gives these outward-rounded sector-bottom intervals
above the full ground, in source energy units:

| Occupation-sector pair | Lower threshold | Variational upper threshold |
| --- | ---: | ---: |
| `±(0,2,0,0,1)` | 1.54674 | 12.71314 |
| `±(1,0,0,2,0)` | 0.56642 | 12.72435 |
| `±(0,0,0,3,-1)` | 3.56682 | 15.71559 |

These are bounds on sector minima, with no assertion of equality to any
endpoint or of an isolated line. Decimal displays are rounded from exact
rational inequalities and are not measured inputs.

## 4. Complete sector preparation with a compensating memory

Two identical idle physical memories A and B start in their declared
filled-beta constant-link anchors. A joint algebraic Householder prepares
\(\psi_\Delta^A\otimes\psi_{-\Delta}^B\). The target and input are
orthogonal physical vectors with identical combined native occupations
and even total fermion parity. Extending this reflection by the identity
therefore preserves those charges and parity. It is an admitted finite
apparatus service, without an efficient compiler or autonomous controller
claim. The ready-memory allocation is at most 1392 qubits.

Fix the tensor order S,A,B,R, with S and A adjacent; no exchange across
the idle reference is implicit. Exchange S with A using the graded swap.
If the target parity is
p, the exchange contributes \((-1)^{p\,p(x)}\) to a source basis vector
x moved into A. Applying the even parity gate \(P_A^p\) removes this
phase. The exact operation, for arbitrary input/reference correlations,
is consequently

\[
 \rho_{SR}\otimes\Omega_A\otimes\Omega_B
   \longmapsto (\psi_\Delta)_S\otimes\rho_{AR}
                         \otimes(\psi_{-\Delta})_B.       \tag{5}
\]

Here symbols for pure vectors in a density expression denote their
rank-one projectors. Every original source/reference correlation is
retained in A; B retains the opposite occupation increment. The complete
operation preserves combined source-plus-memory native numbers and
fermion parity. Each carrier remains Gauss invariant. The source
Hamiltonian, reference and clock are unchanged, while its marginal state
changes by this explicit intervention. The exact three-state control
checks the entire unitary, its charges and arbitrary-reference action;
omitting the parity correction changes odd-sector coherences.

The source work is bounded by
\([M_\Delta-U_0-1/16,M_\Delta-L_0]\). This is source energy accounting
for a controlled intervention, not an autonomous energy-conserving
implementation of its gates. The full dense source vectors and clock
histograms are not numerically executed.

## 5. Energy records and the selected state's sector probability

On the prepared source, \(H-L_\Delta\ge0\) and its mean is
\(M_\Delta-L_\Delta\). Markov's inequality supplies the original
spectral weight

\[
 \Pr_{\psi_\Delta}\{H\in[L_\Delta,2M_\Delta-L_\Delta]\}
                            \ge1/2.                     \tag{6}
\]

One complete source clock has resolution \(\eta=1/64\), tail allowance
\(1/1024\) and implementation allowance \(1/1024\). Every clock output
and both physical memories and the reference are retained. In the interval
\([L_\Delta-\eta,2M_\Delta-L_\Delta+\eta]\), reported probability is
at least \(1/2-1/512=255/512\). Its lower endpoint exceeds the global
ground upper for all six sectors. No success-conditioned normalization
or initial selected-state stationarity assumption is used in this
one-clock experiment.

On the original selected state, (2) instead gives

\[
 \operatorname{Tr}\rho(P_\Delta+P_{-\Delta})
               \le\frac{1/16}{L_\Delta-U_0}.             \tag{7}
\]

The two projectors are orthogonal, and \(H-E_0\) is nonnegative on their
complement. Faithfulness and the explicit witnesses prove strictly
positive original weight, without equating it to (6). A prepared sector
experiment, a stationary thermal-sector population and a relativistic
particle residue are different observables.

These thresholds give exact finite sector information at fixed capacity.
They supply no isolated level within a sector, continuum survival,
Lorentz dispersion, scattering interpretation or GeV rest mass. Native
occupation protection under this Hamiltonian does not classify all
possible repair-protected topological defects or all apparatus actions.
