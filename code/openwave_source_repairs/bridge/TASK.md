# Frozen finite bridge task

Inputs: for p=0…11, r=(p+1)mod12, use2e_p−e_r,8e_p−7e_r and17e_p−16e_r
on the240-entry source. All36 inputs have total1 and span>1. Selection
seed72321000+i, discovery seed72322000+i, census seed72323000+i,
holding seed72324000+i. All seeds fixed before first run.

Selection:1024 body attempts, each uniform first-cell address sampled by
five blank native coin records,16 selector retries maximum. Native body
and matched swap/hold body use identical random records. Census cap32
per edge. Discovery cap256, same16-retry selector. Holding64 events.

Gates:36/36 native selections and0/36 control selections; zero native
post-entry exits;36 successful record-only source handoffs; zero actual/
observed handoff mismatches;36 native discoveries and64 retained holding
events each; native total-variation distance0 between opposite endpoint
one-step laws and deterministic HOLD/SWAP distances1 for all discovered
bonds; all records/selector rejections/timeouts retained. Each source event
and every blank reference event is a literal call to the archived native
function. Control body events are counted separately.

Forced checks: census cap1 reports inconclusive; forcing selector label31
for16 trials reports a selector failure with all80 native records; neither
can produce a successful handoff. A real preselection2/−1 source is read
as multi-load by the same census. A changed current record is an integrity
check only and earns no mechanism credit.

Prepare PREREG and ledger before execution. No earlier new bridge run
exists. Earlier selection and row7 discovery components were known. This
is prospective implementation/integration validation, not a blind new
physical prediction. Standalone headline and independent endpoint replay
target under60seconds each. No platform submission or award is asserted.
