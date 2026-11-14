# Native records and source-transfer checks

For one integer load on a finite undirected graph with (M) seams,

\[
T=I-\frac{L}{2M},\qquad G=n+AF.
\]

Here (T) is the uniformly scheduled native transition, (L) is the graph
Laplacian, and (A) is the half-seam incidence matrix. The modular triangle
has six position coordinates, six records in \(\mathbb Z_6\), mass \(m=.23\),
and one coherent development interval \(t=.37\). Its uniform-fibre isometry
\(V\) satisfies

\[
HV=VH_0,\quad K_rV=V\widetilde K_r,
\quad P=VV^*,\quad
\operatorname{tr}\left[P\sum_f\Pi_f V|0\rangle\langle0|V^*\Pi_f\right]=\frac16.
\]

The \(K_r\) are complete coarse native outcome operators. The \(\Pi_f\)
measure one half-seam record. The first two equalities allow reduced repeated
native evolution; the last equality excludes that reduction after this record
measurement. Terminal effects can be compressed as \(V^*EV\) without implying
that their continuing measurement maps preserve the prepared subspace.

The interacting all-number check uses the specified six-orbital CAR source:

\[
H=H_m+\frac{.37}{2}\sum_h r_6(F_h)^2
       +.61(2I-W-W^*),\qquad B=3-N.
\]

The cycle shift (W) adds `(1,-1,1,-1,1,-1)` modulo six. The fixed Gauss
background is `(1,1,1,0,0,0)`. Direct solution of the Gauss equations gives
six record configurations for each of 64 occupation masks, or 384 states.
The reference is the unique ground projector of this complete Hamiltonian.
For a native channel \(\mathcal C\), the measured survival and work are
\(\operatorname{tr}(\rho_0\mathcal C(\rho_0))\) and
\(\operatorname{tr}[H(\mathcal C(\rho_0)-\rho_0)]\).

| Equation or observable | Implementation |
| --- | --- |
| Literal endpoint average and (T) | [native_scope.py](native_scope.py), `endpoint`, `one_load_kernel` |
| Gauss fibres, uniform code and measurement leakage | [native_scope.py](native_scope.py), `triangle`, `evaluate` |
| Complete all-number basis and interacting (H) | [joint_source_audit.py](joint_source_audit.py), `enumerate_basis`, `operators` |
| Native branch amplitudes, survival and work | [joint_source_audit.py](joint_source_audit.py), `commit_reference`, `run` |

The native endpoint average agrees with the repository's literal rule on
578 signed endpoint/coin branches. All 64 triangle cuts obey the modular
Gauss identity on all 36 basis states. One record measurement leaves
`1/6` in the uniform code; continued position probabilities change by
`0.0001328839493` at the stated time. Starting with definite records gives
transit probability `0.2379296315`; the native port-only chain assigns zero
to that event at every scheduling rate. This is a full-outcome distinction.
It does not negate a separately stated conditional port comparison.

The independent all-number construction imports no producer module. It
solves the basis constraints directly, counts intervening occupations for
fermionic signs, diagonalizes the whole Hamiltonian, and builds native branch
amplitudes using endpoint averaging. The complete basis and Hamiltonian
match the supplied archive, with spectrum error below `2.5e-14`. The native
C, E and equal-mixture ground survivals are `0.6121422563`, `0.3638641361`
and `0.4880031962`; the corresponding dimensionless work values are
`1.1994114010`, `2.3602917513` and `1.7798515762`.

Both scripts retain mutations that reject the relevant invalid transfer.
The interacting check hashes its inputs before and after execution. Its
receipt fixes the precise artifact bytes reviewed; a changed source needs
a separate comparison. These are development checks, with no held-out
acceptance configuration consumed.

From the repository root:

```sh
python3 code/openwave_audit/native_scope.py --output code/openwave_audit/native_scope_results.json
python3 code/openwave_audit/joint_source_audit.py --c10 ../openwave-tuning/next/C10 --output code/openwave_audit/joint_source_results.json
```

The second command requires the explicitly supplied C10 artifact directory.
The first command uses the repository native source for its arithmetic
comparison. Scientific results and version/source identities are retained in
[native_scope_results.json](native_scope_results.json) and
[joint_source_results.json](joint_source_results.json).

The native transition identity also simplifies the scope of clock claims.
For a connected graph, (L) has a unique zero mode and
\(0\leq\lambda(L)\leq2M\). Thus every nonconstant one-load state mode has a
real transition eigenvalue in \([0,1)\). It decays rather than furnishing a
nonzero frequency of persistent oscillation. This statement concerns
state-only observables of the specified memoryless one-load chain. Added
history registers, feedback, coherent continuations, and different preparations
require their own analysis. Ensemble diffusion does not imply loss or
fragmentation of the unit in an individual trajectory.

The [repair-to-quotient proof](../../model/repair_sector_bridge/PROOF.md)
supplies the corresponding finite-source reasoning for the two-load distance
diagnostic. The [clock phase theorem](../../model/openwave_clock_phase_resource/PROOF.md)
specifies the energy-compatibility hypotheses behind its conditional phase
obstruction. Neither result automatically transfers to a changed source,
state, clock, or measurement instrument.

The actual compact energy instrument has a separate
[archive comparison](spectral_archive_audit.py). It independently reconstructs
number-sector spectral projectors and checks all 5,916 serialized branches
for 17 source inputs, including full complex coefficient matrices and zero
outcomes. It also reconstructs the selected native protocol, the populated
doublet's continuing state, and both continuation controls. Maximum
eigen-residual is below `9.2e-15`; the other compared errors are below
`1e-15`. [Its receipt](spectral_archive_results.json) fixes the exact producer
archive and implementation hashes. This is an actual producer comparison;
an abstract proof of the Lüders rule alone would not validate its serialization.

```sh
python3 code/openwave_audit/spectral_archive_audit.py --c10 ../openwave-tuning/next/C10 --output code/openwave_audit/spectral_archive_results.json
```

The [conditional energy-range check](energy_scale.py) uses
\(\|H_m\|\leq\sqrt{m^2+12}\) on a unit-weight graph of maximum degree six.
The eigenvalue bound follows from \(H_m^2=\operatorname{diag}(m^2+BB^*,m^2+B^*B)\)
and \(\|BB^*\|\leq 12\). Under the declared anchor
\(E_0J=.51099895\) MeV, every same-source energy difference is at most
\(2J\sqrt{m^2+12}\). It is `2.1295` MeV for the lowest flat-level anchor
and `15.4265` MeV for the supplied mass `.23`. A physical electron identified
with this same one-load source cannot then cover the 104.5 GeV beam energy
reported in [CERN's LEP RF report](https://cds.cern.ch/record/507409/files/mpph123.pdf).
This is a conditional energy-identification obstruction. A different
many-body or interacting construction does not inherit it.

The separate `4.56665e-18` m spacing calculation requires physical momentum
to equal \(\hbar k\), with \(k\) restricted to the first Brillouin zone.
It is not an independent experimental electron-radius limit. The energy
bound avoids that additional momentum dictionary. The script includes a
small supplied-mass case that the energy-range test does not exclude;
passing this necessary test establishes no physical particle.

The [reservoir apparatus audit](reservoir_audit.py) constructs its reference
as a signed permutation of seven ordered CAR modes. A subprocess exports the
actual producer operators for comparison; none of the producer's algebra,
sector enumeration, references or spectral routines supplies the reference
algorithm. All four operators agree exactly on the 384-state source.
All eight archived reference branches and 2,784 number/energy records agree
within `2.3e-16`. The mean apparatus work is `1.160880350314914` in the supplied
units. Both the source background and the exact reservoir charge
\(n_R-B\) are conserved, with \(B=3-N\).

The fair reservoir preparation has weight `1/2` in each branch. Its energy
is `.23 n_R`, and extra boundary energy is zero by construction. The pulse
commutes with the declared field term, but no physical pulse duration or
complete apparatus cost is derived. Fixed-number reference branches cannot
detect a removed spectator-parity sign. A separate cross-number coherence
control detects that mutation with conditional trace distance one; this
control is an admitted mathematical state, not a physical preparation claim.
The [receipt](reservoir_results.json) binds the target code before and after
extraction and retains all outcomes, including unchanged source number.

```sh
python3 code/openwave_audit/reservoir_audit.py --c10 ../openwave-tuning/next/C10 --output code/openwave_audit/reservoir_results.json
```

**Scope**: finite source arithmetic, complete modular records, specified
native, energy and reservoir instruments, and a conditional one-load energy bound.
No physical electron, winding identification, autonomous rest-energy anchor,
self-starting clock, radiation, infinite-register limit, regulator equivalence,
or OpenWave row award is established. No paper-level claim or platform submission follows from the receipts
alone.
