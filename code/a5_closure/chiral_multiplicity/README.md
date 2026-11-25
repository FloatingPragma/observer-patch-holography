# What do the anomaly equations count?

Inside the ten-component exterior matter menu, arbitrarily many copies of each component are allowed in this test. Cancellation of the four registered anomaly traces forces **the same signed multiplicity for all five matter types**. A saturated measurement of the net number of quark doublets would therefore determine the net number of charged-lepton singlets and lepton doublets. The common integer can be any value. The result extends the mask theorem, which selects one generation or its conjugate when a component can appear at most once and opposite orientations cannot coexist. [Imported exterior tables and mask theorem](../../../Lean/Screen/ExteriorSelection.lean) · [Integer multiplicity theorem](../../../Lean/Screen/ExteriorMultiplicity.lean)

## Which multiplicities are possible?

Choose the displayed orientation $G=(Q,u^c,d^c,L,e^c)$, with charges in sixths of hypercharge $(1,-4,2,-3,6)$. Write $q,u,d,l,e$ for the five component counts minus their respective conjugate counts. The ten source rows give the following exact, linear responses.

| Trace | Coefficients on $(q,u,d,l,e)$ |
| --- | --- |
| Gravity squared times hypercharge | $(6,-12,6,-6,6)$ |
| Colour squared times hypercharge | $(2,-4,2,0,0)$ |
| Weak isospin squared times hypercharge | $(3,0,0,-3,0)$ |
| Cubic hypercharge | $(6,-192,24,-54,216)$ |

If these traces vanish, the weak row gives $l=q$, the colour row gives $d=2u-q$, and the gravity row gives $e=l$. The cubic row then reads $144(q-u)=0$. Consequently $(q,u,d,l,e)=k(1,1,1,1,1)$. Every integer $k$ solves the four equations. The cubic colour anomaly is $2q-u-d=0$. Counting colour copies and both orientations gives $4k+6\bar n_Q+2\bar n_L$ weak doublets, an even number, so the usual mod-two Witten check for this doublet menu passes. This finite representation calculation uses the ordinary spin-manifold parity test, with the global topology and physical chiral realization specified separately.

With nonnegative component counts there is a unique inventory in a fixed orientation. Remove from each conjugate pair the smaller of its two counts as matched vectorlike pairs. The residue is $k$ complete $G$ for $k\ge0$ or $|k|$ complete conjugate generations for $k<0$. Choosing chirality in the original mask sense excludes all vectorlike pairs. The primitive nonzero chiral solution then has rank 15, and three such copies have rank 45. Gauge rank, faithfulness and anomaly freedom by themselves supply no reason to choose the third multiple. Neutral singlets, other gauge representations and strongly interacting emergence lie outside the ten-slot theorem.

## Which extra premises select three?

| Assumptions on this matter menu | Consequence |
| --- | --- |
| Full charged chiral packet and anomaly cancellation | Any positive number of complete copies in either chosen orientation. The same diagonal $\mathbb{Z}_6$ acts trivially for every such count. |
| One Higgs, intrinsic quark CP capability, positive one-loop weak coefficient | Counts $3,4,5$ all pass. With $n$ copies the coefficient is $(43-8n)/6$. An explicit three-dimensional Fourier mixing block has squared Jarlskog invariant $1/108$ and extends to four or five copies. |
| Multiplicity must be one complete faithful band of the actual twelve-port adjacency | Available band dimensions are $3,3,5$. The four-dimensional abstract $A_5$ representation is absent from this particular adjacency spectrum. |
| The family carrier occupies the lowest-cost admissible band | Its rank is exactly three, since $5-\sqrt5<6<5+\sqrt5$. Applying this count to particles requires attachment of that band to the measured chiral matter response. |

The finite gauge derivation supplies its conditional weak-doublet, colour and hypercharge types from the declared exterior matter module. The inherited spatial spinor and statistics require the spin and matter join. The simulator's orbital winding doublet has a different source map. A local rotation or exchange on a named mobile carrier, or a mass map on that carrier, needs its own observable evidence. [Conditional matter image and family selection](../../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) · [Family band theorem](../../../Lean/Screen/A5FamilyBand.lean) · [Generation-window certificate](../multiplicity_window_certificate.py)

## What could an observer measure?

On a complete single-orientation packet, calibrate the amplitude of a weak raising transition of one fundamental doublet to one. Then the inclusive transition strength $\mathrm{Tr}(T_-T_+)$ equals $4|k|$, counting three colours from $Q$ and one doublet from $L$ per family. Measuring twelve such units counts three generations only under the stated completeness and chirality assumptions. A normalized expectation per particle discards the extensive normalization. The code also enumerates the kernel of the central action and the electromagnetic charge-square trace on each finite representation.

Three oriented generations and two oriented generations plus one conjugate generation both give dimension 45, weak strength 12, and $\mathrm{Tr}(Q_{\rm em}^2)=16$. Their signed indices are 3 and 1. A separate example adds a vectorlike pair of charged weak singlets to three generations. The weak strength stays 12 while the dimension increases to 47 and the charge-square trace to 18. The mirrored examples obey the anomaly conditions and deliberately relax the original chiral selection. A measurement that resolves handedness or detects the additional singlet channels distinguishes these alternatives.

For each count from one to six, the control also tensors the charged subspace with a separately supplied spinor and verifies opposite orders of perpendicular half-turns differ by a minus. In units with $v/\sqrt2=1$, both the diagonal charged-lepton Yukawa assignments $(1,\ldots,n)$ and $(1,\ldots,1)$ commute, as mass operators on this subspace, with the same spatial turns and electric charge. These are allowed mass choices in a comparison model. The finite theorem fixes the chiral index relation; charged-lepton mass eigenvalues and distinct propagating species require an independently measured dynamics. [Exact checks and countermodels](certificate.py) · [Reproducible receipt](results.json)

Run from the repository root with a single numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/a5_closure/chiral_multiplicity/certificate.py
(cd Lean && lake build +ExteriorMultiplicity)
```
