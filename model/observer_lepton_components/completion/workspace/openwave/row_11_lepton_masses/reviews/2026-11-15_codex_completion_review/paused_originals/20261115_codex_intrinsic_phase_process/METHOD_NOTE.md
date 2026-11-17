# Method note

    p=Pr(1[j_X=1] XOR 1[j_Y!=0])=2/9,
    J=r(1−2o), pi=(5/12,5/12,1/12,1/12),
    P=T D T,
    D=(I−epsilon P0)(I−epsilon p A1 P1)(I−epsilon p A2 P2),
    A1=1[J<=0], A2=1[J=-1],
    H_epsilon=−log(S)/g_epsilon,
    S=sqrt(rho) P/sqrt(rho), g_epsilon=−log(1−epsilon),
    H_epsilon -> P0+p(A1P1+A2P2),
    Pr(N=n)=epsilon(1−epsilon)^(n−1).

| Term | Source/code | Check and scope |
| --- | --- | --- |
| j, p, pi | pinned `native.py::native`; `headline.py::run` source-role and pair loops |720 literal double columns; two3600-pair count tables; original outcomes retained |
| epsilon, delta | `run` zero-coin role loops; epsilon2^-k, deltaepsilon^3 |240 literal coin columns; product-word compression explicit |
| P, T, D | `run::dkernel`, direct refresh/flip sum and exact multiplication | All1024 row transitions at k1,2,3,6; detailed balance and factor positivity |
| all characters | `run` eight-character loop | Four exact4-by4 blocks per character family; all32 source states |
| g and r=h/g | `loginterval`; reference covariance loop; PROOF post-run pole lemma | Exact rational logarithm enclosures and executed reference lags0..8. Analytical relative-pole attachment uses lambda_r=(1−delta)^2 and lambda_p=lambda_r(1−epsilon p), each with original weight5/6; no new check count |
| full algebra | matrix-unit loop; PROOF sections3,5 |4096 one-step words; positive density belongs to symmetric S representation |
| first mark | first-mark loop and PROOF section6 | Eight complete finite-duration laws plus timeout for every refinement; infinite geometric tail and resolvent proved |
| controls | HOLD source loop; full changed kernels |480 changed-law calls; AND decoder and count-one plant separately named |

Results:37 checks pass. The full generator error bound at epsilon1/64 is
0.003252541649632; the phase/reference ratio is approximately0.22086522317,
with controlled limit2/9. Full byte replay passes. The first-mark stop yields
an actual Euclidean resolvent. Under an additional coherent family reading,
its limiting phase harmonic is1/(1−iqp), rather than exp(iqp).

Scope: N+ positive-memory-history component with supplied preparation,
memory/routing and clock convention. No native selection of a fixed family
action window, no executed family analyzer, no charged-lepton identification
and no physical mass ratio. Same-group review passes no OpenWave gate.
