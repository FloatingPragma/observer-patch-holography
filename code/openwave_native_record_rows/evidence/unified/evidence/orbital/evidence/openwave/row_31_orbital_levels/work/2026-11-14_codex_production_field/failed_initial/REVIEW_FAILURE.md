# Initial attempt failed independent band-decoder review

The producer printed all registered checks pass and its own replay was
byte-identical. Independent review then found an implementation error:
`u=proj[:,argmax(diag(proj))]` is a view, and in-place normalization rescales
one column of the already stored projector. Thus the retained band readouts
are not orthogonal internal-band projections. The initial band-readout
clause fails. Its code, preregistration, receipt, archive, logs and producer
pass claims are preserved here without alteration. The raw native events,
occupation counts, clocks and curl are separately checked on amendment.

The correction must copy the vector and add independent projector identity
checks, with preregistered amendment before rerun. This is not a tolerance
change and does not erase this failed attempt.
