# A signed charge winds around a separately located field excitation

A phase preparation on a path of written seam heads puts measurable ring energy at its two endpoint faces. A controlled material tour detects which side of a six-link boundary each endpoint occupies. This occurs on the 30 seams and 20 active triangular field terms of one actual cell of the positive-Laplacian repair. The initial field state, the phase preparation, addressed tour and interference reference are specified apparatus inputs. [Design](DESIGN.md) · [Calculation](search.py) · [All outcomes](results.json)

## Where the field energy goes

A circulating pair of sharp head words on face f has energy operator J(I − X) with J = 1/2. The uniform superposition of all 458 words connected to the specified initial word has zero ring energy. On a path between two faces of the dual graph, multiply each amplitude by the product of head parities on the crossed physical seams. A ring flip reverses this product precisely when its face is an endpoint. Hence each enabled undirected flip there contributes 4J/458 to the mean ring energy, and every other flip contributes zero. A different dual route gave the same state up to global phase for every pair. A matched diagonal field law with the same eligibility and rate produces zero excess energy on every face under this preparation.

All **190** unordered choices of two distinct faces gave nonzero energy at both selected ends. For the antipodal faces (0, 1, 5) and (2, 3, 6) the energies were 0.323144 and 0.441048 in native units. Both faces are vertex-disjoint from the physical tour. Under free evolution with all 20 ring terms open, the fraction of energy on these two faces fell from 1 to 0.704975 at time 1 and 0.267456 at time 2. The localized phase preparation is an evolving field excitation.

## What the outer charge detects

The physical perimeter (4, 9, 8, 7, 10, 11) encloses ten faces. At time π/8 an isolated enabled repaired material contact swaps the charge position and seam head with amplitude +1. The fixed forward and truncated reverse sweep makes the two circulating inputs wind in opposite directions and return with their boundary heads reversed. It acts on a signed material load, with a fixed remote compensator and the integer Gauss constraint intact. In the uniform ground, 24 of 458 words wind for each sign. Of the other words, 400 return unchanged for an excess and 198 for a vacancy. The remaining 34 or 236 finish away.

| Endpoint faces in the disk | Number of pairs | Signed amplitude within the winding branch | Inclusive bright probability, excess | Inclusive bright probability, vacancy |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 45 | +1 | 0.962882 | 0.742358 |
| 1 | 100 | −1 | 0.910480 | 0.689956 |
| 2 | 45 | +1 | 0.962882 | 0.742358 |

The inclusive probability uses an ideal equal-amplitude identity arm and counts every input, including nonreturn. A full forward sweep followed by its exact inverse gives bright probability 1 for every input. The minus sign measures mutual winding of a material load and a prepared field endpoint. Its origin is the odd intersection of the perimeter with the preparation path. The earlier product preparation on two separate circulating rings changes the inner ring independently and leaves the same outer tour insensitive to it. Here the field energy is local to the endpoints while the relative phase extends over the path joining them.

## Does the preparation affect unsteered motion?

Start the material at port 0 with the two antipodal ends, open every internal written contact and triangular field term, retain the physical port-degree potential 8(6 − d), and keep the compensator outside the cell. The charge-plus-field sectors have 1,954 and 5,790 states for excess and vacancy. At time 0.5 the excess port distribution differs by 0.310983 in total variation from a release with a uniform field and the same charge. The vacancy's largest difference on the time grid is 0.088742 at time 1.625. The receipt includes all port probabilities and all face energies from time 0 to 2. At time 0.5, 50.95% of the excess state's ring energy lies on faces incident to its charge, against 29.52% from multiplying the separately measured marginals. This energy-weighted correlation and the reference-arm signal establish a physical interaction in the specified extension. Identifying a spin-carrying mobile particle would additionally require a bound pair of internal states and their self-exchange.

A follow-on read asks whether the distant charge arrives together with a field pattern characteristic of a displaced endpoint. At each eligible receiving port we construct the uniform field there and the field with the original two endpoints. We project their span out of a third field template which places the near endpoint on a face incident to the receiving port and leaves the remote endpoint fixed. The normalized remainder defines a joint charge-position and whole-cell coherent read. Every receiving face is vertex-disjoint from the source face. The old and bare templates have exactly zero overlap with this read at the receiving port. Free evolution from a bare source can nevertheless scatter into it, so we retain that control.

At time 0.625 the excess reaches the port-8 read with a displaced field signature in **1.4957%** of string-input trials and **0.2188%** of uniform-field trials. The vacancy reaches the port-9 read in **0.6623%** and **0.0151%**, respectively. Each reported setting maximizes the string-input probability on the 17-time grid across the retained port/face filters, ten for the excess and fifteen for the vacancy. These are inclusive joint probabilities. They supply a charge-correlated displaced-field signal within one cell. The filtering requires coherent control over the field of that cell, and the observations do not determine the lifetime or exchange of a bound composite. [Full receiving templates and uniform-field controls](results.json)

From the repository root, run the following command into a fresh path with NumPy, SciPy and NetworkX available.

~~~sh
OPENBLAS_NUM_THREADS=1 python3 -B model/candidate_laws/interval_echo/openwave/electron_field/flux_string/search.py --output /tmp/echo_flux_string.json
~~~
