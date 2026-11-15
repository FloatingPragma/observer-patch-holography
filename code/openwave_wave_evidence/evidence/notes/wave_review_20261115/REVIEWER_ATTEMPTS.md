# Reviewer execution history

These are review-tool failures, not new scientific trials.

1. The original custody command in the Windows working tree failed on
   `evidence/boost/CELL.md`. `git ls-files --eol` showed LF Git objects and
   CRLF working files.
2. The reviewer's first `reproduce.py` attempt also failed custody, before any
   scientific replay. Investigation found that `git archive` inherits
   `core.autocrlf` and had exported CRLF. For example, the native source's
   committed SHA-256 is `dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d`,
   but that archive produced
   `a31b6fae4de2ea15bc0766e436f6b50732ed156d60a009cb8e9b572c12c9bbe7`.
   Normalizing that archive's CRLF recovered the exact committed blob.
3. The extractor now sets `core.autocrlf=false` and `core.eol=lf` only for the
   archive command and independently verifies **every extracted Git blob ID**
   against `git ls-tree`. It does not change repository configuration,
   original files, scientific source or manifests.

The successful extraction and subsequent scientific/verification outcomes
are recorded separately in `results/`.

4. The first WSL invocation stopped before running Python because WSL stripped
   backslashes from the path passed to `wslpath`. The reviewer launcher now
   passes a forward-slash Windows path. This is only an execution-path fix;
   source/manifest/acceptance logic is unchanged. Each environment retains
   the exact reviewer script it executed as `REPRODUCE_SOURCE.py`.
