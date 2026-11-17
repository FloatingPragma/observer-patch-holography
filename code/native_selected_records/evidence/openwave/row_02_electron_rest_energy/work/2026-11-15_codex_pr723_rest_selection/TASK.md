# Prospective finite task

Candidate: A-native-discovered-bond-v1. Law, memory, state and decoder are fixed in
DESIGN.md and experiment.py. Deterministic exhaustive columns; no development or
held-out random seeds, no stochastic significance claims. Inputs: copied graph,
literal native.py (sha256 dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d),
all 12 one-unit starts, both signs, 30 bonds and both orientations. Search cap256,
return cap32, feedback reset cap8. These caps were fixed before execution.

Gates registered before the first claimed run:
1. Full graph240/630, compact f/g satisfy Lf=6f; 80 eigenvalues below6.
2. Actual fine-current unread channel has compact fidelity<1 and outside mass>0.
3. Every initial local signed unit has first movement probability1/12.
4. Every directed bond/sign gives stationary kernel [[3,1],[2,2]]/4 and
   stationary weights2/3,1/3. Disabled feedback gives uniform rows1/2,1/2.
5. Every displaced state returns by32 ticks with failure<=33/2^32; failures kept.
6. Separate native work pulse yields W/log2=-1,0,+1 with probabilities1/6,1/2,1/3
   in all120 contexts, from a reference derived independently of that pulse.
7. The matched exponential-work mean equals1; capped feedback-reset mean differs.
   This is a controlled information-work identity, not a mass criterion.
8. Independent paired epochs yield the literal product native kernel; product
   modular costs and relaxation eigenvalues add logarithmically. Exclusive
   scheduling differs. No physical-time equivalence is assumed.
9. Retain the negative physical gates: selected stationary transfer energy=0,
   no physical work calibration, no electron identification, no independent mass
   or frequency measurement. Existing compact pair exceeds1022keV.

Exact rational comparisons are used except full240 spectral census and displayed
logarithms (absolute tolerances1e-8 and1e-9 in the independent reader).
run: python3 experiment.py --write; replay: python3 headline.py;
separate reader: python3 independent_check.py. All are offline, no install or fetch.
Registration is local only; no retrospective platform preregistration is claimed.
