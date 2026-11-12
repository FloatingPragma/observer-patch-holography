# Keeping the complete lepton archive in fewer qubits

The prepared lepton source needs at most two quantum archive qubits after
any of the tested one, two or four rounds. Every original return history
is retained in two classical bits per round, and an explicit isometry
reconstructs the complete source archive, including reference correlations.
The exact four-round ranks are certified without floating eigenvalue
thresholds, which incorrectly discard small directions in 136 histories.

[Proof](PROOF.md), [receipt](results.json) and [code](compression.py) give
the native compression circuits and their additional service costs.
The construction supplies coherent control and routing; the information
dimension is not a physical work or autonomous-reset calculation.

Run `python3 -m model.lepton_archive_compression --check` from the repository
root to replay the native services and every complete history.
