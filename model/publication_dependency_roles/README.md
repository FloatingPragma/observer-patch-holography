# Scientific inputs and presentation provenance

The exact register classifies three receipt edges as presentation provenance.
Their producer read traces use the document bytes only for hash metadata or
required binding membership. The checker authenticates the current owner,
provider, presentation document and every operative input. It rejects an
extra document edge and any numerical use of these notes, including a read
hidden inside a binding function or hash comprehension.

Run `python3 -m model.publication_dependency_roles --write` to generate the
current role-check result and `--check` to replay it. [roles.json](roles.json)
contains the exact paths and digests; [check.py](check.py) checks their roles.
The result is a provenance check. It does not rerun scientific calculations
or transfer an independent component's physical claims to the official
common source.

The Yang--Mills presentation has one-way publication provenance. Its
mathematical criterion, source code, Lean proofs and data remain inputs of
their scientific owners. The paper binds those outputs without becoming an
input to their calculations. After a paper edit, run
`python3 -m model.publication_dependency_roles.paper_check --write`, then
`--check`. Repeated binding edges must agree and their exact union is hashed
once per invocation. Operative changes invalidate the check; there is no
cross-run cache or document exemption.
