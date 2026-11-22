"""Read the frozen flight and return acquisitions against their native detector laws."""
from hashlib import sha256

import numpy as np

from . import shell_beats, shell_spectroscopy, shells
from .shell_source import check

PUBLISHED_COUNTS = {
    "shell": "9088a369367e3fd5fe4ec4c61a0ffc39b5ac1c82d95bd359aed43885a6f5fd56",
    "shell_spectroscopy": "2e7b438ee945f511729c9556c00dff0b8e4e83aa290bfc38031fd571e5042893",
    "shell_beat": "ba5ada96b04969e3a8a4ce1900a971e33901222a41951e64a3982aa719408cba",
}


def check_cohort(counts, probabilities, shots, error):
    check(counts.shape == probabilities.shape and counts.dtype == np.int32
          and np.all(counts >= 0) and counts.sum(dtype=np.int64) == shots,
          "incomplete archived detector cohort")
    check(np.isfinite(probabilities).all() and probabilities.min() >= 0
          and abs(probabilities.sum() - 1) < 1e-10,
          "incomplete native detector probabilities")
    check(np.max(np.abs(counts / shots - probabilities)) <= error,
          "archived counts outside simultaneous detector bound")


def replay(name, run, protocols, shots):
    path = shells.HERE / f"{name}_counts.npz"
    check(sha256(path.read_bytes()).hexdigest() == PUBLISHED_COUNTS[name],
          f"changed published acquisition: {name}")
    with np.load(path, allow_pickle=False) as archive:
        records = {key: archive[key] for key in archive.files}
    bins = sum(counts.size for counts in records.values())
    # A union bound across every bin of all three runs, plus CDF quantization.
    alpha = .01 / len(PUBLISHED_COUNTS)
    error = np.sqrt(np.log(2 * bins / alpha) / (2 * shots)) + 2.**-36 + 1e-10

    def observe(key, index, probabilities):
        counts = records[key][index]
        check_cohort(counts, probabilities, shots, error)
        return counts

    observed, summary = run(*protocols, observe=observe)
    shells.save(name, observed, summary, output=None)


def main():
    flight = shells.protocol("shell_protocol.json")
    spectra = shells.protocol("shell_spectroscopy_protocol.json")
    beats = shells.protocol("shell_beat_protocol.json")
    replay("shell", shells.run, (flight,), flight["shots"])
    replay("shell_spectroscopy", shell_spectroscopy.run, (spectra, flight), spectra["shots"])
    replay("shell_beat", shell_beats.run, (flight, spectra, beats), spectra["shots"])


if __name__ == "__main__":
    main()
