# Reader-only receipt schema correction

The original sealed verifier compares an in-memory Python object containing tuple
coordinates with its serialized JSON receipt, whose coordinates are lists. Python
container equality distinguishes them. This amendment changes only that comparison
to compare JSON values. All source calls, scientific checks, source-pin equality,
complete retained arrays and exact array comparisons remain unchanged.

The wrapper checks the original file hash, changes one exact line in memory, keeps
its original __file__ so source paths and pins remain identical, and runs without
--record. It writes no original packet file. A mismatch after canonicalization is
a failure and must be investigated; no numerical tolerance or gate is changed.

This is a new reader amendment, not retroactive success of the failed original
replay. The failed root output is retained. Original packet and receipt stay sealed.
