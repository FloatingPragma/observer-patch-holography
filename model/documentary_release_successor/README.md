# Documentary release successor

This is a documentary successor to the sealed r2064 manifest. It preserves
that original manifest, every original scientific authority, the physical
program class and all original test files. The current release keeps its
identity and date. Exactly three paper artifact rows and the canonical book
row may have their reviewed PDF hash and size updated; membership, paths,
all other fields and all other artifact rows remain unchanged.

The current manifest and every actual current PDF are checked by the
unchanged strict release validator. The reviewed authority also binds the
current source inputs and the actual build record. This does not replace a
paper build with a hash update, and it does not accept a physical gate.

The original current-record navigation checker sealed the former publication
manifest in addition to its scientific inputs. Its unchanged scientific
checks now run in a private view containing the current repository files,
with exactly `paper/paper_release_manifest.json` replaced by its archived
sealed bytes. The view is not claimed to be a coherent publication build:
its sole purpose is to check the original scientific custody contract with
its original documentary identity. Actual publication artifacts are validated
separately against the actual current manifest.

The private view preserves every other current file byte-for-byte. It
contains no historical checkout of scientific code. Every file is copied;
none is hardlinked to the actual source. The copied manifest is replaced by
archived bytes, so the main checkout is never modified. Every
input is checked before and after the scientific validation. Only the two
explicit current-record custody/navigation test files use this execution
view; their original test bodies and cases run unchanged. All case IDs,
counts, order and actual failure outcomes are retained. Per-case properties
and session output identify the documentary substitution. All other tests
keep their assigned current or separately declared historical domain.

`python -B -W error -m model.documentary_release_successor --check` checks the current
publication artifacts and executes the unchanged current scientific custody
validator in that one-substitution view. The source snapshot and actual test
journals remain distinct from the frozen original authorities.
