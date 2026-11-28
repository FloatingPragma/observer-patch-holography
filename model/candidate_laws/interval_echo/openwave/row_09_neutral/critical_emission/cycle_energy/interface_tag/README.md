# Before-absorption charge records and the distal daughter tag

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

Opening the last link for one original dwell gives the already identified daughter either sign of local charge at the last dimer, with exactly opposite charge on its neighbor. Two independent pointers can read these signs before the absorber. The original guarded absorber and vacuum tag have exactly zero joint success on either charged branch. A full native swap after both records restores the balanced aperture and gives measurable calculated joint success on both prior signed branches. [Protocol](PROTOCOL.md) · [Complete joint record](results.json)

## Reading the two signs without erasing the neutral mode

For a data port of load $a$, the native diagonal contact costs to the two pointer branches differ by $V(a,2)-V(a,0)=65(1-a)/2$. Routing both rails of one dimer to the same refocused contact gives a phase proportional to its charge $q=a+b-2$. The original parity and sign settings alias $q=0$ with $q=\pm2$. The same pulse circuit with half the phase exposure gives separate magnitude and sign reads on the even branch. The original settings still decode $\pm1$ if the first parity read gives odd. The two observers use distinct pointer rails, 25–24 and 61–60, and return each data visit by a checked native transposition. The complete five-qutrit pulse columns give the expected charge for every ternary word. They preserve the relative phase of the two neutral orientations. All original parent and one-mode refusals remain in the bank. [Calibration and routing](measure.py)

## What the guarded absorber can report

On the 02/20 guide let $D$ select unlike loads at the two contacting ports and let $D_-$ and $D_+$ select their two possible orders. On this link a complete native contact of duration $t=1/260$ has the same diagonal phase as the original guarded pulse. It multiplies the no-exchange part by $c=\cos(t/2)$ and the swapped part by $i s=i\sin(t/2)$. A swap yields charge $\pm2$ on the receiver and $\mp2$ on its neighbor. After the charge read, the neutral branch carries $I-(1-c)D$ relative to the guarded continuation. Each charged branch carries $i s\,\mathrm{SWAP}D_\pm$. The original subsequent guarded calendar, local phase operation and final vacuum tag preserve each dimer's amount. The terminal vacuum belongs to the balanced support, so any continuation with those properties assigns exactly zero vacuum amplitude to both charged branches.

Use $\beta$ for the calibrated original one-mode vacuum-amplitude row, and $\psi$ for the enrolled lowest mode. Its neutral tag probability is $|\beta\psi-(1-c)\beta D\psi|^2$. The separate calibration arm performs a full native swap on the same link only after recording opposite nonzero charges. This turns the charged branch into $i s D_\pm\psi$ up to a common phase. The odd parts of $D_\pm$ under global orientation reversal have zero vacuum-tag matrix element between one-mode inputs. Each retained signed branch therefore gives tag probability $s^2|\beta D\psi|^2/4$. The bound from the actual bright leakage and the lowest-mode coefficient is strictly positive for both signs. An independent calculation on all 256 physical extreme-load words of four dimers checks the contractions, the charge branch and the native swap.

| Calculated probability per original launch | 23 dimers | 45 dimers |
| --- | ---: | ---: |
| Either specified nonzero charge sign at the receiver before absorption | $2.16573\times10^{-9}$ | $1.35471\times10^{-9}$ |
| Zero charge and direct vacuum tag | $4.24363\times10^{-4}$ | $1.28895\times10^{-4}$ |
| Either specified nonzero charge sign and direct vacuum tag | 0 exactly | 0 exactly |
| Either specified nonzero charge sign and vacuum tag with the later full swap | $2.26976\times10^{-11}$ | $7.14162\times10^{-12}$ |

The receipt keeps the tag complement for each sign and includes every launch. In $2^{42}$ seeded draws per arm the direct charged bins contain 9,314 and 9,509 refusals at 23 dimers, and 5,953 and 5,994 at 45. The charged-arm return gives 108 and 104 tagged draws at 23, and 39 and 26 at 45. The exactly computed zeros follow from the continuation's support, independent of these draws. The explicit return calibrates a detector after a charged intermediate record and supplies the compensation through the same link.

The two before-absorption charge readings establish zero sum on the enlarged interface even in the charged branches. All other dimer amounts in this preparation were fixed by the original guard. The original tag therefore has a charge-zero interpretation on its balanced continuation, and its success alone says nothing additional about amounts on the branches it excludes. A conserving inference about a released daughter on a more general aperture needs a before-absorption comparison with the unchanged parent and the remaining material, or a theorem identifying those records for that aperture. The present three-outcome charge read and guarded tag do not supply a vacuum decoder on continuously unguarded matter. [Row disposition and exact charge premise](../../../../../../../../openwave/interval_echo/preparation/row_09/science_gate/README.md)

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_09_neutral.critical_emission.cycle_energy.interface_tag.measure
~~~
