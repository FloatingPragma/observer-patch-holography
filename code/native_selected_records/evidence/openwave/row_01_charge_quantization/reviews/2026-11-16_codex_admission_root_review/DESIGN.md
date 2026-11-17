# Admission audit on independently generated mixed inputs

The reviewed source is the sealed record-only admission candidate with manifest
ba4fc9e36bf41592b272feee134c4bdd54d67c51416c8f992416535d536c9f2f.
Its existing 1752-case result and proof are known. This audit exercises its
actual controller against a separately implemented source, using floor/ceiling
of the endpoint sum rather than the producer's current formula.

For each public B=2,3,4 and M=4,8, the input menu includes eight structural
preparations and 32 hidden mixed vectors generated independently of the source
coin stream. Each input receives each of two fixed public coin streams. The
controller receives only the event closure and public graph. It receives no
initial vector, total, seed or case identifier. Source inspection belongs to
this audit and occurs outside the controller interface.

Every complete run must reconstruct the initial core exactly. Native source
events must preserve the total and obey the independently checked current
formula. Zero total must terminate without a unit. Every nonzero input must
extract one same-sign target, keep center0 empty at admission and retain all
excess in the reservoir. The preserved marker and blank/copy regions must
agree with the complete source state. Every reported strict decrease must
actually lower the quadratic potential. A capped timeout is retained and
cannot count as an extracted target. A forced constant coin must produce the
explicit probe-timeout outcome.

This is held-out verification of the stated all-input theorem, not a new
prospective physical discovery. The audit imports the controller being tested,
but no producer source law or preparation code. New case inputs, code and gates
are frozen before their first execution. Default replay compares saved evidence;
an explicit exclusive first-record destination is required to create it.
