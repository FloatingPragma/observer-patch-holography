# First independent review correction

This annex supersedes the dictionary-uniqueness wording in `PROOF.md` §2.
The original package is byte-frozen in manifest
`c4f6db87b1d2c7ddd4caa4c8ac5a95f09c33679f060d50845aa0bf42f5e54d12`.
Its executed source, topology, laws and examples do not change.

The abstract isomorphism Z^30/im(A) ≅ Z^11 does **not** select a unique
charge dictionary up to sign. Z^11 admits GL(11,Z) automorphisms, including
nontrivial shears. The correction executable exhibits one with determinant
one that is neither I nor -I.

The candidate instead declares a **local circle-target map on each actual
dual edge**, with the spatial incidence intertwiner B and fixed regional
boundaries. For that chosen circle target, its primitive fundamental-group
generator has the two orientation conventions +1 and -1. The exact lattice
calculation proves that this declared spatial map is faithful and complete;
it makes no uniqueness assertion about arbitrary nonlocal charge dictionaries.
The doubled and sign-only mutations test the declared local map and its
actual pointer dictionary. They are not a theorem ruling out every other
possible model dictionary. These architectural choices are explicit in
`CONTRACT.json` and are not empirical charge calibration.

`topology.raw_pointer_degree` is a numeric decoder **after complete tape
validation**, not a standalone admission function. For example, an empty
dictionary marked successful would have raw count zero; it is not a valid
zero-charge instrument. The public `reviewed.admit_and_decode` first calls
the frozen observer-only facade. Only an admitted complete native tape can
reach the raw decoder, and its result must agree with the admitted signed
charge. The added executable rejects an empty-success claim. Genuine zero
readout still performs the actual discovery and reset documented in
`VACUUM.json`.

The correction adds no physics premise and changes no scientific result.
The root reviewer identified both issues independently. The full-row verdict
remains pending the two reviews of the base package plus this annex.
