# Which port charge survives a seam write?

With electrically neutral memory, every conserved additive real charge of the sharp ternary loads assigns equal and opposite charges to a vacancy and an excess, relative to a unit-filled host. Allowing the seam memory to carry charge admits a second conserved quantity. It counts one for each defect and two for each memory recording contrast of magnitude two. Both conclusions use the full native write and its inverse. [Registered seam operator](../../../../candidate_laws/contrast_memory_spherical_v1/law.py) · [Exact integer checker](charge_identity.py)

## Conservation on every write column

Assign charge \(q_n\) to a port of load \(n\), the same at every port. Use a sharp blank as zero of memory charge and let \(m_{c,d}\) be the charge of a written coin and contrast label. For each of the 18 nonzero amplitudes in the native write, conservation requires
\[
q_a+q_b=q_{y_1}+q_{y_2}+m_{c,a-b}.
\]
Inputs of contrast zero or one in magnitude preserve their unordered pair of loads, giving \(m_{c,0}=m_{c,1}=m_{c,-1}=0\). Each of the two input orders 02 and 20 writes both coin values and balances the loads to 11. Hence
\[
m_{c,2}=m_{c,-2}=q_0+q_2-2q_1.
\]
These equations exhaust every written label and every column. If the memory is neutral, the right hand side vanishes. Subtracting the host charge leaves precisely \(q_n-q_1=A(n-1)\). With a charged memory, every port assignment extends uniquely to the sharp memory labels. Equivalently,
\[
q_n-q_1=A(n-1)+B(n-1)^2,\qquad m_{c,d}=2B\,\mathbf 1_{|d|=2}.
\]
On an initially blank excitation with one vacancy and one excess, total charge relative to the host is \(2B\). During a native 02 to 11 write the same \(2B\) moves into the seam record. Forcing that record to be neutral would violate all four signed-extreme input and coin combinations for \(B\ne0\).

Write \(Q_B,Q_W\) for the additive charge on the blank and written blocks. The column equations give \(Q_W W=W Q_B\). Taking adjoints shows that \(Q_W\) commutes with \(W W^*\), which proves commutation with the entire 99-state echo, including the written complement. The memory term commutes with contrast reversal and squared contrast. Thus the charge also commutes with the native seam generator for any contrast coefficient and with a sum of these contact generators when each physical port and memory is counted once. A replacement or clearing of charged memory requires accounting for its destination.

The blank seam at nonzero contrast also fixes the choice of sharp load basis for a common one-port Hermitian charge commuting with the generator. The sum of two identical one-port charges commutes with SWAP. For \(n'\ne n\) the contrast-square energies of \(|n,n\rangle\) and \(|n',n\rangle\) differ, so their commutator entry forces the single-port matrix element \(q_{n'n}\) to vanish. If one permits a separate diagonal assignment at each site, native blank swaps force their differences between load levels to agree across each connected component.

## Coherent memory charges

A general Hermitian memory charge with the blank as a sharp zero-charge state gives the same two parameters on all prepared blank inputs. For \(d=\pm1\) the two written coin values multiply distinct load words. The write forces the memory operator to annihilate each of those four sharp labels. At \(d=0,\pm2\) the two coin values multiply the same load word. The equal-coin superpositions must have eigenvalues \(0,2B,2B\), respectively. Hermiticity leaves only the three orthogonal coin-odd directions at these even contrasts free. They have zero amplitude in every written preparation.

At the native coefficient eight, the memory generator restricted to these directions has three distinct energies, 0, 32 and 33. Their vectors are the coin-odd state at zero contrast and the symmetric and antisymmetric combinations of the two coin-odd extreme contrasts. An arbitrary real charge on each of these three modes gives all additional Hermitian memory freedom compatible with both the write and generator. The checker verifies the three projectors on the full echo blocks. These choices leave the charge of the blank prepared pair equal to \(2B\).

## What identifies the electric source?

The signed ramp and its two separate material standards resolve the \(A\) charge. The unsigned ramp resolves the defect term on blank material. A read that a vacancy and an excess have equal and opposite *electric* charge, or that a written extreme-contrast memory has zero electric charge relative to blank, would set \(B=0\) for an independently identified additive electric source of this form. For a matched-port comparison, writing 02 gives two unit loads and a coin-even memory of contrast −2. An ordinary pair of unit loads with blank memory has the same terminal port census. The two conserved assignments with coefficients (A,B) = (1,0) and (1,1) predict respective charge differences zero and two between these preparations. [Material impulse](README.md) · [Remote signed and parity reads](../remote_charge/README.md)

The mobile-charge boundary flux in the written-field Laplacian uses a different seam generator. The finite Maxwell-action instrument reads the material load and enters its excess over one as the specified source density. The blank vacancy loop measures path-dependent permutations of the material register. Each tests a separate connection between load and observed records. The free bound pair requires an independently calibrated electric-flux response on its own law or a demonstrated transfer of that response to its source. [Written-field flux](../../../../candidate_laws/interval_echo/openwave/electron_field/README.md) · [Material source in the finite Maxwell action](../../../interval_echo_gauge/temporal_action/README.md) · [Vacancy loop measurements](../../../interval_echo_gauge/loop_transport/README.md)
