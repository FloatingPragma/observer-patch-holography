# Frozen finite-error amendment

This amendment tests periodic whole-plane filtering on the fine observer,
before executing its error cases. The v2 coarse process's tolerances are known;
no fine-frame finite-error result has yet been run. The ideal fine-frame pilot
is explicitly prior knowledge.

For M4/M8, h4, one/two turns, a uniform input, and delta in {0,+0.01,-0.01},
apply that SAME systematic phase error to the initial native contrast
preparation analyzer and EVERY subsequent native plane-filter analyzer.
There is no added relative phase between the main rotating/reference arms.
Require one-turn visibility<-0.99, two-turn visibility>0.97, and final rotating
source contrast fidelity>0.9999. Retain M8 two-turn +/-0.1 stress cases without
a success gate. Every preparation, marker, target, filter and final path
outcome is retained. No gate changes or error-case removal after execution.

The source instruments are those of the current frozen candidate; complex
analyzer phase is the only perturbation. At every operation all flag
probabilities are emitted. Incoming amplitude ensembles and the bound Kraus
catalog reconstruct every post-outcome state. Uniform-input warmup compression
to three orthogonal modes is checked against all18 initial native branches;
it is not a new choice of source state or hidden projection. Marker outcome
histories are the ideal candidate's independently decoded histories, since
this error affects neither marker endpoint nor classical controller.
