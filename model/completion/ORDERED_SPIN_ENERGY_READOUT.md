# Certified ordered-spin energy readout and the remaining shell proof

The finite spectrum is calculated and certified. The sharp particle mass shell
of the infinite ordered ground is **not proved**. These are different results:
the finite positive selection field and the finite box cannot be treated as
physical rest-mass terms.

## Calculated relative energies

The Hamiltonian is the existing FA98 example,

\[
H_B=J\sum_{\langle xy\rangle}(\mathbf S_x\cdot\mathbf S_y-1/4)
 -B\sum_x\epsilon_x S_x^1,
\qquad B/J=1/4,
\]

on the open eight-site cube with twelve bonds. The reference energy is its
ground energy, and the relative unit is its first positive excitation gap:
\(\Delta_i=E_i-E_0\), \(r_i=\Delta_i/\Delta_1\). The first six separated
positive spectral enclosure groups are shown below. Enclosure midpoints are
rounded to nine decimal places. Grouping overlapping enclosures does not assert exact degeneracy or
identify particle species.

| Enclosure group | Excitation energy \(\Delta/J\) | Relative energy \(\Delta/\Delta_1\) |
| --- | ---: | ---: |
| 1 | 1.092284678 | 1.000000000 |
| 2 | 1.442323755 | 1.320465062 |
| 3 | 2.166719726 | 1.983658445 |
| 4 | 2.543094000 | 2.328233703 |
| 5 | 2.665699561 | 2.440480596 |
| 6 | 2.729285257 | 2.498694079 |

In particular, the first gap is rigorously enclosed by
\([1.092284677676,1.092284677698]J\). The
[complete machine-readable readout](ordered_spin_energy_readout.json) retains
all 256 ordered levels, exact rational bounds, outward decimal bounds and
the enclosure grouping. These are spectral calculations; no claim is made
that all levels have nonzero weight in a particular detector.

The [readout tool](../../tools/read_ordered_spin_energies.py) verifies the pinned
FA98 witness and rational certificate before computing the ratios. It uses
the same ground energy in numerator and denominator, bounds the ratio over
the allowed interval corners, and rejects a nonpositive denominator.
It introduces no fitted physical energy scale. The
[original construction](../native_ordered_ground/README.md) supplies the
Hamiltonian, certificate proof and actual finite native preparation.

## What the gapless result says about rest energy

The existing infinite-ground result gives nonzero response in a broad band
\(\alpha Jp\le E\le\beta Jp\) for arbitrarily small \(p>0\), with fixed
positive constants. Consequently that response has zero positive-energy
threshold. If a sharp branch is subsequently proved within this band, then
its rest-energy limit is zero: the upper bound tends to zero as \(p\to0\).
If that branch also has a relativistic dispersion
\(E(k)^2=v^2|k|^2+m^2v^4\) with fixed \(v>0\), taking the limit gives
\(m=0\). This conditional implication neither constructs the branch nor
identifies it with a measured massless particle. Ratios of two zero rest
masses are undefined, not the finite energy ratios in the table.

## The unproved ground-state cancellation

For the zero-field ordered ground \(\Omega\), the exact identity is
\[
K^2A_f\Omega=J^2 e(P)A_f\Omega+F_f\Omega,
\qquad e(k)=\sum_j(1-\cos k_j).
\]
The [shell diagnostic](../spin_shell_diagnostic/README.md) proves the current
and three-spin force identities and the spectral criterion. For one common
speed parameter \(c>0\), the estimate
\[
\left\|G(K/(Jp))
 [F_f+(1-2c^2)J^2e(P)A_f]\Omega\right\|^2=o(J^4p^5)
\]
would concentrate the entire selected response on the scaling cone. A weaker
alternative is to prove a positive retained fraction in shrinking cone tubes,
allowing a continuum in the complementary response. Neither conclusion has
been established for this ordered ground. The proved first/inverse moments
and Ward identity alone allow continuous spectral measures, so using them as
a proof of this cancellation would add an unsupported assumption.

The Hamiltonian and finite preparation are specified. The missing input is a
theorem controlling its infinite-ground correlations and spectral weight, not
an unspecified numerical coupling that can be filled in from the finite table.
