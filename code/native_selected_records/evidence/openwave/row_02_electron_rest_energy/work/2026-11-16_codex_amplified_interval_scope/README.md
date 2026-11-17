# Impurity-bound scope correction

The sealed amplified interval certificate has 112 native cases and 144
rotation-control cases. Its native pair impurity is below 0.005486130010;
the valid bound across all 256 cases is 0.005486165872. Both are below the
unchanged registered eta=0.01 gate.

The narrower value in the sealed source PROOF and its FINAL_RECEIPT summary
must be read as **native-only**. Their numerical certificate, source code,
gates and every individual interval remain unchanged. The independent RER
integration review identified this prose-scope error. This packet corrects
it without replacing a sealed artifact or weakening a scientific gate.

Run `python3 -B check.py` for the exact-dyadic, standard-library-only check.
It reads the sealed certificate, writes nothing, and rejects failures even
under optimized Python. The result is a known-data scope verification,
not a new experiment or an external G6 review.
