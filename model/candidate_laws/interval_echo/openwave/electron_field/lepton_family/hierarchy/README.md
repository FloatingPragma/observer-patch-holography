# How much energy range can a single-charge guide provide?

A written guide with one configuration per distinct charge port, one enabled contact of weight four per active seam and no independent ring state has generator
\(H=4L_G+8\operatorname{diag}(6-d_v)\). Here \(d_v\) counts seams on the physical board and \(L_G\) uses the active edges. Either sign has this operator at load reference one with its compensator outside the guide. The neutral vacuum has energy zero. The bounds cover all connected scalar graphs of this form on the stored boards. More than one head configuration at a port changes the state space. A measured energy difference between two charged levels also has a different reference.

| Federation | Longest degree-six run bound | Ratio ceiling, simple cycle | Ratio ceiling, open path | Ratio ceiling, any connected scalar graph |
|---|---:|---:|---:|---:|
| 20 cells | 6 | 48.59 | 137.15 | 440.95 |
| 80 cells | 6 | 48.59 | 137.15 | 440.95 |
| 320 cells | 12 | 137.15 | 444.45 | 1272.73 |

Each ceiling bounds **any two eigenenergies above the same neutral vacuum** in a connected guide containing a degree-five port. The tables round up. Thus the tau/electron hierarchy lies outside every such guide, even with arbitrarily branched routing. The muon/electron hierarchy also lies outside every simple cycle. The degree-six induced graphs have maximum degree two and contain no cycles. A connected guide made wholly of those ports is a path of at most twelve sites with zero ground energy. The ratio of its highest to lowest positive standing-wave level is at most \(\cot^2(\pi/24)<58\). Comparing different neutral media or choosing a new energy zero for each level would change the question. [Geometry, rational certificates and spectra](results.json)

## Why the bounds hold

On a cycle, split the quadratic form into intervals between degree-five ports. Each endpoint contributes half its norm and half its onsite energy to each adjacent interval. Each intervening degree-six port contributes its full norm. On an open path add the terminal intervals with a free end. The census bounds the number of intervening sites, so there are finitely many local stiffness and norm matrices. Their lowest generalized Rayleigh quotient bounds that of the entire guide. The cycle inequality also holds with a boundary phase since replacing amplitudes by their moduli can only decrease edge energy.

For a branched guide, each connected cluster of degree-six sites is a subpath of one of the censused paths. Keep one available edge from each cluster to a neighbouring degree-five attachment and discard the other edges between clusters and degree-five sites. The discarded terms are positive. Give each cluster one fifth of its attachment's norm and onsite energy. At most five clusters can use one attachment. Its unused fraction has Rayleigh quotient eight. Test every cluster size and every position of its attachment. For each block, exact positive pivots of the integer matrix obtained from \(A-\lambda B\) certify the rational lower bound. The upper energy bounds are 24 on paths and cycles and 56 on general guides, by the degree bound on the weighted Laplacian. [Executable certificate](search.py)

## Searching with Koide balance

For ordered positive energies normalize their square roots to \(1,b,s\) with \(1\le b\le s\). At exact balance, \(3(1+b^2+s^2)=2(1+b+s)^2\). Substituting one for the middle root increases the left side minus the right by \((b-1)(3+4s-b)\). Hence \(s^2-8s-2\ge0\) and the largest energy is at least \((4+3\sqrt2)^2=67.9411255\) times the smallest. Equality is attained when both light levels coincide. The span is a necessary screen even when the circulant phase is unknown. [Lean proof](../../../../../../../Lean/ObserverPatchHolography/KoideCirculant.lean)

The direct spectral scan includes every simple physical path starting at its sole degree-five port and continuing exclusively through degree-six ports. The 20-, 80- and 320-cell boards admit respectively 984, 3,296 and 17,664 such paths. With equivalent path lengths counted once, the largest board has 12 spectra and 1,001 distinct index triples across them. Its closest value is \(Q=0.669634316\) at indices \((0,1,11)\) on a 13-port path, with native ratios \(1:8.912719:265.393338\). This was selected by closeness to balance after scanning. It demonstrates why a numerical near-hit needs a phase test and an independently propagating species label. [All template counts and best triple per board](results.json)
