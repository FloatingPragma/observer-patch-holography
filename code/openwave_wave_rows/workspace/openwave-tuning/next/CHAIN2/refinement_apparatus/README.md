# Native wave-packet records

A native preparation circuit produces a nonuniform two-mode scalar packet on the cubic LOCAL12 walk. Native C/E pointers read its terminal position, displacement and shared clock. The construction gives an error bound for the entire record law that tends to zero along a dyadic refinement family.

Run `OPENBLAS_NUM_THREADS=1 python3 test_apparatus.py` in this directory. It executes development graphs of sizes 2 and 4, verifies the native blocks and writes the finite macro tape, retained native records and result. `python3 generate.py` regenerates the size-two macro tape independently. The source walk uses certified block compression, with its error stated explicitly. Larger sizes in the result are formal bound evaluations only.

The preparation covers one named scalar Fourier-support family. It supplies no physical clock calibration, autonomous-field limit or particle-species identification. The [construction contract](CONTRACT.json), [proof](PROOF.md) and [executable result](RESULT.json) state the complete domain.
