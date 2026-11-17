# Prospective robustness repair, v2

Candidate **A-native-incidence-rotation-filtered-v2** descends from the v1
candidate in this directory. Its new operation is three common optical-sector
filters before the main path is split. V1's failed phase-error diagnostic is
retained: a1e-5-radian preparation phase error gave two-turn contrast fidelity
about0.399 and visibility0.622. This motivates, but does not validate, v2.

Each new filter coherently compares a fair native full-triangle event against
a native auxiliary-vacuum event, erases its18 fine flags, then accepts the
auxiliary path-minus result. All36 (eraser,path) outputs are kept. Its success
Kraus map is (I-P_all)/(2sqrt3). The native endpoint determines
P_all=Q_constant+P_optical/2, so success is P_optical/(4sqrt3); the swap law
instead gives P_optical/(2sqrt3). This is an entire two-dimensional plane
filter, not a chosen ray or rotation phase. It commutes with all triangle
permutations and acts identically on both future interference arms. Therefore
it cannot itself impose the law-dependent full-turn sign. The same swap
control and reader remain mandatory.

The filter is implemented from its complete fine records, not substituted by
its analytic projector expression. Three ideal filters cost48^-3 under native
law and12^-3 under swap, in addition to every already reported cost. Arbitrary
input robustness remains the exact preparation theorem. The new diagnostic
also perturbs all preparation/filter path analyzers by a common frozen phase
error delta in {0,+/-1e-5,+/-1e-3,+/-1e-2,+/-1e-1}. This is an adversarial
systematic imperfection, not an installed rotation phase. Every accepted and
rejected output is retained at the instrument's documented compression level.

Frozen gates: all28 v1 ideal cases must keep their exact fringe and ray values,
and their success must be multiplied by48^-3 (native) or12^-3 (swap). At h4,
uniform input, one/two turns, errors with |delta|<=0.01 must keep native
visibility<-0.99 / >0.97 and rotated contrast fidelity>0.999. The0.1 case is
an explicitly ungated stress diagnostic and remains even if it fails. No
claim of tolerance beyond those gates, and no removal of failed v1 cases.

These gates were written before v2 execution. Ideal comparisons are predicted
reproductions; only the new complete filter composition and finite error
challenge are new tests. All source and flag probabilities remain charged.
