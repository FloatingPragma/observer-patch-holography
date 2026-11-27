"""Integer full-column checks for a pinned physical port on a registered seam."""

import numpy as np
from model.candidate_laws.interval_echo import candidate


def seam_parts():
    law = candidate.native_law
    pattern = law._write_pattern()
    reversal, square = law._label_ops()
    identity = np.eye(90, dtype=np.int64)
    twice_complement = 2 * identity - pattern @ pattern.T
    exchange = np.kron(np.eye(10, dtype=np.int64) - reversal, np.eye(9, dtype=np.int64))
    contrast = np.kron(square, np.eye(9, dtype=np.int64))
    parts = []
    for memory_term, multiplier in ((exchange, 1), (contrast, 2)):
        if np.any(pattern.T @ memory_term @ twice_complement):
            raise ValueError("Blank/written cross block is nonzero")
        scaled = np.zeros((99, 99), dtype=np.int64)
        scaled[:9, :9] = 2 * multiplier * pattern.T @ memory_term @ pattern
        scaled[9:, 9:] = multiplier * twice_complement @ memory_term @ twice_complement
        parts.append(scaled)
    return parts


def audit():
    law = candidate.native_law
    scaled_parts = seam_parts()
    for scaled, physical in zip(scaled_parts, law.local_term_parts(), strict=True):
        np.testing.assert_allclose(scaled / 8, physical, atol=2e-14, rtol=0)
    words = np.tile(np.array(law.SRC), (11, 1))
    records = []
    for port in (0, 1):
        for fixed in range(3):
            take = np.flatnonzero(words[:, port] == fixed)
            other = words[take, 1 - port]
            bad = other[:, None] != other[None, :]
            residuals = []
            for term in scaled_parts:
                entries = term[np.ix_(take, take)]
                residuals.append(int(abs(entries[bad]).max()))
            records.append((port, fixed, len(take), residuals))
    return records


if __name__ == "__main__":
    print(audit())
