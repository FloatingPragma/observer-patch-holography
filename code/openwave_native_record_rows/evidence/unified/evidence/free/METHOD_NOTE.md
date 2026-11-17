# Method note: source records to common free waves

P=I−L/(2W), L=L_internal+a^-2 L_cubic4. Native branch counts fix L.
For normalized occupation sums, C_n=(1/(4W))[Π+2Σ_(j=1)^(n−1)(1−j/n)ΠP^jΠ]
and C_infinity=L^+−Π/(4W). The finite transfer recurrence gives all moments;
the native all-moment limit proves Gaussianity. A source field is never
installed in place of the event.

The internal source spectrum is0(1),5−√5(3),6(5),5+√5(3). Scalar band
symbols are λ+4a^-2Σ_(μ=0)^3 sin²(ap_μ/2). Slice energy is
E_a=2asinh(a sqrt(λ+q_a)/2)/a and
0<=λ+|p|²−E_a²<=a²[Σp_j^4+(λ+|p|²)²]/12.
The same bound yields the common limiting speed. Four constant-sector
channels give F_μν=D_μ A_ν−D_ν A_μ. The reflected physical Gram is positive
of rank two; scalar/vector internal projectors commute, so one joint source
supports both algebras. Full covariance/moment limits yield local causal
Maxwell and Klein–Gordon fields.

| Term or result | Exact implementation or proof |
| --- | --- |
| native event and full branch compression | common_cubic/probe.py:main, frozen inputs/native.py; original37908 branches and384 full-vector checks |
| current position, occupation and same-tape clock | common_cubic/probe.py tape loop and independent reconstruction loop; every wait and final open interval retained |
| source L, bands and projections | common_cubic/probe.py integer moves/Laplacian/Kronecker block and projector loop |
| finite C_n, contact, clock, EE/BB and IID control | common_cubic/probe.py population loop; explicitly compressed n=100000000 |
| higher moments | shared/native_occupation_fluctuation_limit/PROOF.md; fixed-graph transfer proof and retained numerical gates |
| two physical polarizations | record_maxwell/PROOF.md, direct reflected Gram and Gaussian lift |
| common masses/speed | slice_refinement/PROOF.md and continuum_reconstruction/PROOF.md |
| same source stopping-clock dictionary | continuum_reconstruction/SOURCE_CLOCK_JOIN.md and COMMON_NATIVE_FIELD_V1.md |

The source population checks pass their filed thresholds: eigen residual
2.42e-14<1e-10, covariance error6.764e-5<0.001, clock error4.18e-14<1e-10,
EE/BB error0.000461<0.001. All169364 native calls and exact artifact bytes
replay. The four short tapes have clock counts4,5,3,0; their unfinished
intervals and tail remain. Negative finite contact eigenvalues are recorded.

Scope: one declared native free-field source with a controlled full
reconstruction and actual finite source execution; no primitive field
update, charged matter/detector preparation, laboratory constant fit or
upstream acceptance is asserted.
