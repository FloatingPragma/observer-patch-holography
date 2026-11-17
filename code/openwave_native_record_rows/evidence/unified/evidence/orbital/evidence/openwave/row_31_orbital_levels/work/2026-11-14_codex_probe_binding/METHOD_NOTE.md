# Method and result

P_out=I-L_eff/(2W), L_eff=L_external+(1-h)L_internal; the conditional
move probability is(1-h)/2. h is derived from actual native first-hit
records. Q=(1/22)[[16,4,0],[2,16,2],[0,3,16]] gives h=(5/22,1/11,1/22).
The controller never reads this matrix, h table, a gamma or a potential.

`probe.py` calls the frozen endpoint function in `inputs/native.py` for
44100 census branches and every reset, probe, redirected hold and outer
attempt,305157 calls total. All records replay from their disclosed
starting positions. The exact native kernel, capped targets and finite
well eigenvalues are calculated from branch counts in that same script.
The exact finite matrix assumes almost-sure preparation; the executed
bounded variant has conditional history error at most25·2^-128 per attempt
as specified in SCOPE_ADDENDUM.md.

All registered comparisons pass. Feedback100 moves0.384521484375 against
0.3863636363631919; disabled movement0.509521484375 against0.5. Simultaneous
99% interval halfwidth0.03181913310349, plus the preparation allowance.
No preparation or probe timeouts occurred; every trial remains included.
The probe cap tail bound is2.532051878034809e-11. Three finite lowest
squared energies are2.2720093211074106,4.218897837046055,4.612330904148159,
strictly positive and below their unfitted native thresholds. Source and
receipt gzip bytes match on a full rerun.

Scope: local native feedback and a source-derived finite well; no executed
large-shell spectrum, long production-field tape, physical charge or
hydrogen identification.
