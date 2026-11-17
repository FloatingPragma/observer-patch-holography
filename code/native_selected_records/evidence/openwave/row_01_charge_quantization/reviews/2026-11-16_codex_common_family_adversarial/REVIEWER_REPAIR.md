# Reviewer implementation repair

The first k=0 verification attempt failed at the first nonzero-time Ramsey check because this review script reused the family refinement variable k as the transport-word loop index. That reviewer-only variable shadowing made the reconstructed measurement stride incorrect. The original reviewer source is retained as audit_initial.py; the producer source and evidence were unchanged. The loop index is renamed primitive_id before rerunning the same checks. This is a verification implementation repair, not a changed scientific gate or producer finding.

The first selector-census reviewer attempted little-endian words and the outcome label reject. The declared selector uses big-endian words and the label retry. The source was read and the reviewer census corrected to that declared encoding. Original reviewer source is retained as ticket_audit_initial.py. This does not change any producer record or scientific gate.
