# Joint clock and transport records with a supplied spherical ruler

The interacting incidence candidate uses one cell–seam kinetic coefficient for charge and cage
and an onsite contact controlled by a native load pair. Both objects carry a mass-tag clock.
The clock, native swap and contact are active in the same free interval. A parent-consistent
spherical chart sets the link lengths, and a separate rest-clock read sets the timing unit.
[Incidence model](../relativity/README.md) · [Joint replay](experiment.py)

## A ruler and a reference clock on the native carriers

Associate each stored cell with the normalized centroid of its spherical triangle. The base
icosahedral chart is matched to the actual cell dual. At a refinement the central child and
the three corners are identified from their internal star and the neighboring parent seams.
The chart verifies the full fine edge set. Supply one base separation as the ruler unit.
It fixes the sphere radius. Every seam length $\ell_s$ is then its independently calculated
centroid arc. The central descendants of the tested base beacons keep their chart positions.

At an isolated massless link use the same incidence generator below. The three states are
the sending cell, the seam and the receiving cell. For coefficient $c$ the receiver has
$P(t)=(1-\cos(\omega t))^2/4$, where $\omega=\sqrt2 c/\ell_s$.
A separate rest oscillator reads the time from two quadratures. The link probability
and the measured time on the reference seam give $c$. With that value fixed, relay
pulses timed from the atlas reproduce full transfers on 629 held-out seams. Round trips
read the elapsed rest-clock times between the declared beacons. The relay advances
at $\sqrt2 c/\pi$ during the supplied pulse sequence. A detuned read is a control.
[Atlas and calibration replay](geometry.py)

Uniform hop lengths at each level fail to preserve distances on the native refinement.
Under the central-child map, coarse hop distances one and two become three and
four. The required scale ratios differ. The replay tests both bounds on both stored
refinements. Weighted routes are compared with the independent great-circle beacon
distance as an atlas audit. The chart, its ruler and its reference oscillator are supplied
standards. The clock and relay measurements calibrate lengths and times within this model.

## The interacting generator

For actual oriented seam $s$, the corresponding column of $B$ has entries
$1/\ell_s$ and $-1/\ell_s$ at its two incident cells. The charge and cage use the same
read value of $c$ in the cell–seam generator

$$
D_m=\begin{pmatrix}mI & cB\\ cB^* & -mI\end{pmatrix}.
$$

The cage cell mode denotes its intact selected ring excitation. The seam mode is an
auxiliary label for that excitation. There are two adjacent binary load configurations
on one actual complete-cell face. They differ by one native swap. Their internal-cell
contrast costs agree and their addressed external contrast-square values are one and zero.
The supplied pulse keeps that swap and parks the other native swaps. It removes the common
internal contrast energy and omits the external blank-energy bias, leaving the half-swap
generator $H_L$. The addressed external square enters through the contact below.

The full source is

$$
H = H_L + D_{m_q}^{(q)} + D_{m_c}^{(c)}
      + g A_s n^{(q)}_s n^{(c)}_s ,
\qquad A_s=|0\rangle\langle0| .
$$

The two number operators select the charge and cage seam modes. The onsite contact
$g=1$ preserves the fixed metric and all kinetic coefficients. It fails to commute
with the live native swap and transport. With the same preparation, the additive
control sets just $g$ to zero. The binary loads stay adjacent under this selected
pulse. Full-link experiments and addressed-seam reductions test the 20-, 80- and
320-cell carriers, with both excitations initially dressed from the sending cell.
The supplied one-body preparations have positive energy. The interacting evolution
retains both signs of one-body energy. [Source and basis](source.py)

## Read the clock during contact

The mass tag begins in a definite label. Its first $X$ read starts the relative phase
and records the outcome. Its two free arms are $H_\pm=H\pm\delta\beta/2$, with
$\beta=+1$ on the selected object's cell modes and $-1$ on its seam modes.
Every other term, including the interaction and native swap, is active in both arms.
The terminal $X$ and $Y$ counts recover
$Z_\delta=\langle e^{-itH_-}\psi,e^{-itH_+}\psi\rangle$.
All first and last outcomes and all terminal groups enter the inclusive phase. The
common rest preparation is the uniform state in the connected cell support with no
seam component. It has $\beta=1$, the contact annihilates it, and its measured rate
is exactly the supplied mass split on the same source.

Duhamel's derivative gives

$$
-\left.\partial_\delta\operatorname{Im} Z_\delta\right|_0
= \int_0^t \langle\beta\rangle_s\,ds .
$$

The replay checks this by independent sparse block propagation. When $x=|\delta|t\le1/2$,
Taylor bounds give $|\operatorname{Im}Z_\delta+\delta\int_0^t\langle\beta\rangle_s ds|
\le x^3/6$ and $\operatorname{Re}Z_\delta\ge1-x^2/2$. The normalized phase rate differs
from the exact average grading by at most $2x^2$. Both moving and rest clocks, and the
comparison with zero contact, are read from counts. At $t=3/4$ and $\delta=1/128$,
the contact-induced phase-rate shift exceeds the sum of the two phase error bounds
for each object on each chart. [Clock assay](clock.py)

## A common velocity bound in the interacting history

Let $F$ be the distance to the sending beacon in the supplied chart, with a seam's
value the midpoint of its endpoints. The onsite contact and native load swap commute
with each object's $F$. The radial velocity $V_F=i[H,F]$ is odd in that object's $\beta$.
Each cell–seam entry of $V_F$ has magnitude at most $c/2$ by the metric triangle
inequality. The cell row sum is at most $3c/2$ and the seam column sum at most $c$.
Consequently both objects obey the same bound $v_*=c\sqrt{3/2}$ for this read.
This bounds radial drift in the declared lengths and times.

Oddness and $\beta^2=I$ give $\langle\beta\rangle^2+
(\langle V_F\rangle/v_*)^2\le1$. The time average lies in the same unit disk,
since the disk is convex. The exact clock rate and finite-interval radial drift in
the actual interacting source therefore satisfy

$$
r_{\rm clock,0}^2+
\left(\frac{\langle F(t)\rangle-\langle F(0)\rangle}{t v_*}\right)^2 \le 1.
$$

The finite-gap clock has the phase error bound above. The same histories give an
obstruction to using the expectation of position as a proper-time trajectory.
The initial state and Hamiltonian are real in the detector basis. For any real diagonal
position read, the initial expectation of its current is zero because the current
is imaginary antisymmetric. A proper-time clock based on the centroid speed therefore
starts at the rest rate for any positive limiting velocity. The infinitesimal
mass-tag response is $1-2p$ for initial seam probability $p>0$. This obstruction
is specific to centroid-based trajectories.

The radial drift bound concerns a mean distance. A strict signal cone fails for
this finite continuous-time generator on sharp cell seeds. Across an enabled seam,
the first nonzero amplitude between its endpoint cells is of order $t^2$. Only the
cell–seam–cell path contributes to that coefficient. The product of its oriented
hops is $-c^2/\ell_s^2$, so the amplitude is $c^2t^2/(2\ell_s^2)+O(t^3)$ even with
the other registers and the contact present. This gives nonzero receiver
probability at arbitrarily short positive waits. It leaves open an effective
continuum cone with a bound on the escaping probability.

## Obstruction to carrying over the positive-shell clock

For each free $D_m$, let $E=\sqrt{D_m^2}$ and $\Pi_\pm$ be its energy projectors.
Within the respective bands the grading compresses to $\pm m/E$.
An onsite seam contact can populate the negative band. On the addressed seam,
the dressed sender has equal norm in the symmetric and antisymmetric components.
Denote the initial seam probabilities by $p_q,p_c$. The negative charge-band
weight of the seam ket is $1-2p_q$, so the squared leading amplitude for a
negative charge-band outcome is $g^2(1-2p_q)p_qp_c>0$.
The replay checks a strictly positive short-wait read
against a Taylor bound and the zero-contact control. Thus the exact stationary
positive-shell formula of the incidence model does not apply to the whole
interacting history.

There is a restriction on diagonal replacements for this contact. Fix a connected
cell graph and let $P_m$ be the free positive-band projector. Its cell block is
$(I+mA^{-1/2})/2$ for $A=m^2 I+c^2BB^*$. The integral formula for $A^{-1/2}$ in
terms of $e^{-uA}$ gives strictly positive entries in this block, since $-BB^*$
generates an irreducible finite cell walk. The cell-to-seam block is
$c A^{-1/2} B/2$, with a nonzero entry in each column, and every seam has positive
diagonal weight in $P_m$. If a Hermitian contact is diagonal in the detector basis
and preserves the product of the two free positive sectors, it commutes with the
product projector. For fixed load and spectator mode the entries just described
force its diagonal values to agree across all cells and seams of the other object.
Applying the condition in both factors makes the contact independent of both
objects' positions. An exact free positive sector thus rules out a nonconstant
position-diagonal contact on this incidence source.

The same history has a quantitative comparison. If $f$ is the norm in the
negative band, the difference between the grading read and the expectation of
$m/E$ is at most $2f^2+2f$. In the free interaction picture $f$ is $g$-Lipschitz,
and the expectation of $m/E$ is $2g$-Lipschitz. Midpoint enclosures give an
error bound on the integrated comparison. A physical proper-time identification
requires a trajectory or geometric momentum attachment with controlled
refinement, a physical ruler, and an interacting vacuum. The calibrated incidence
coefficient, scheduled relay velocity and radial envelope have distinct measured
meanings in this candidate.

## Replay

Use `model/requirements.txt` from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.scalar_channel.mobile_memory.calibrated
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m pytest -q model/scalar_channel/mobile_memory/calibrated
```

The executable runs the full-link experiment at all three levels. The unit-test
replay uses full-link charts through 80 cells and addressed reductions through
320 cells.
