# Numerical implementation correction, before a completed v2 error receipt

The first frozen `filtered_error.py` aborted during a preparation probability
audit: a trace's imaginary residual failed a relative tolerance at a nominally
zero-probability output. It produced no completed receipt; earlier cases may
have evaluated internally before the exception. No gate outcome was printed.
That original source is retained verbatim as
`filtered_error_failed_v1.py`, matching its original preregistration hash.

The corrected implementation propagates the18 warmup amplitude branches as
an ensemble. Every probability is a sum of squared norms and cannot suffer
the cancellation inherent in the former density trace. The complete outcomes,
input errors, native maps, parameters and frozen scientific gates are identical.
This amendment is hashed before the corrected execution. No failed scientific
gate was changed; the v1 scientific robustness failures remain separately.
