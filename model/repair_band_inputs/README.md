# Finite repair-band comparison inputs

These numerical tables describe twisted port Laplacians and coherent walks
on the declared two-carrier honeycomb cell. The band levels use the
relaxation clock; the walk table uses phase per walk step. The chosen
section, gluing, frame rule and address count are recorded in each producer.
They supply cutoff-scale comparison bands without a particle identification.

`bands/` contains the four producers and their full JSON tables. Run a
producer with its working directory set to `bands`; it uses NumPy and
imports only its sibling producers. The formulas and finite geometry are
specified in the source docstrings.

`cf1/runs/H2_coherence_48.json` records mean free paths for routing in a
link-flux gas with fixed level-10 comparison densities. These are separate
from the active universe census. `manifest.json` binds each retained input
by its source identifier, length and SHA-256.
