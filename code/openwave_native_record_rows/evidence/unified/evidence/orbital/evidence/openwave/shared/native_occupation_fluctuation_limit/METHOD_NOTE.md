# Method note

P=I-L/(2W), Y_n=sqrt(M/(4Wn)) sum_t(e_Xt-1/M),
V_n=Pi sum_t(2c_t-1)/sqrt(4Wn), F_n=Y_n+V_n.
Cov(F_n) tends to L^+; its error on a P eigenvalue z<1 is
z(1-z^n)/(2Wn(1-z)^2). All higher fixed moments tend to Gaussian moments by
finite analytic transfer perturbation. Rooted differences have covariance
L_D^-1 in the limit.

Code map: inputs/native.py:native executes the endpoint rule. probe.py:run
enumerates720 one-unit branches and24 auxiliary zero-pair coin branches.
Current continuity rebuilds the states. probe.py:moments propagates exact
integer path counts of powers up to4; it sums all words by a proved Markov
recurrence, without claiming literal long-tape executions. probe.py:inv
independently computes the graph Green target and rooted identity. headline.py
regenerates and byte-compares the complete receipt.

The stationary preparation, independent events, window and normalization are
supplied apparatus. No nature value, Gaussian random draw, oscillator step or
physical Hamiltonian is input. The target Laplacian is a comparator derived
from the declared graph; the moment transfer is reconstructed from native
current records. HOLD fails the same variance comparison. Population moments
are exact rational values, with no sampling confidence interval.

Scope: finite source Gaussian closure, positive mathematical scalar limit;
no implemented real-time instrument, physical metric, interaction, atomic
level or completed OpenWave criterion. The n512 failure and n2048 resource
amendment are retained and described in TASK.md.
