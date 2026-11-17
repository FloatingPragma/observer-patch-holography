# Same-group review of the compatible-current witness

Reviewer: Codex native_gravity_bridge. Written code/receipt review; no new source execution. This is not G6 or independent cross-group review.

Reviewed frozen/run artifacts:

- run.py: d8415fd77efaa79561ec5d460f6b98d650377b4bf6da08804c7ab43ea8670831
- verify.py: 566f57b8b7ebce0fb78c4fcdfc4bff4f06ba16e33f6679dcbf5138ca9e9c6527
- TASK.md: 842fa6fa68f0f3f50445a0136cdf8ec085e8188dbfcf31c108ad25eac89221ee
- PREREG.json: fd8568caebce6e0cbb1ae6ccc1b775485c6f63e0f88a2573dfd03126e78cc370
- RECEIPT.json: aa9b72fc6bbb20c5e9f0705c6aabb1303a79ce9c488b9177c2f1f544ef7ec3d6
- VERIFICATION.json: 501d0dfb1341b5e1c91a7f65a00e4d4795dd7d19e9c7ac7a44fbe1fc77d7eb9b
- RECORDS.jsonl.gz: 43bec8c495e9740e87947bbf3f0d0988c6fc084c4f1794ee5c6733241056a457

The producer enumerates both initial orientations and all4096 twelve-bit words, so the persistent two-port current history is complete. The gate/sign readout uses absolute native currents, not only empty-seam coin labels. Conditional fair current bits imply independent four-bit ternary coordinates with weights1,14,1. This implication remains valid despite the shared persistent unit: its location is tracked and its next absolute current is fair in either orientation.

For one orientation, the compatible sum-zero weight is14^3+6*14=2828. Both orientations give5656 accepts out of8192, hence707/1024; the nonzero compatible words total2*6*14=168. Each conditional coordinate's nonzero count is112 over both orientations, so the sqrt8-scaled diagonal covariance is8*112/5656=16/101=112/707. Sum-zero symmetry gives off-diagonal-8/101. The finite fourth moment128/101 follows from the same count. These independently derived integers agree with the receipt.

The hold mutation explicitly replaces the event by the identity and records each branch as current0. It therefore makes all ternary coordinates zero, accepts all8192 words, and destroys source flux variance. Its event-loop count is correctly separate from98304 calls to the pinned endpoint. Both preparations and all rejected words are disclosed. The mutation is a control implementation, not a second execution of native.py.

The producer's convolution powers and independent reader's direct multinomial coefficients implement different exact population calculations. The constraint indexk=3n-i-j correctly enforces z1+z2+z3=0. Factors8/n and64/n^2 in the second/fourth moments agree with F=sqrt(8/n)Z. The long-n values are explicitly compressed exact source-law populations, not literal long native tapes or exact finite Gaussians.

The independent current reader checks each archived event against the closed form j=left-coin and updates its tracked location. It reads the entire source/control archive and independently recomputes every population. This review inspected that code and the completed verification receipt; it did not repeat the reader or count its work as a new experiment.

No scoped blocker found. This is a small exact source-sensitive compatibility instrument and a valid finite witness for the conditional theorem. It is not an executed spacetime Maxwell field, curved stress measurement or charged-state preparation.
