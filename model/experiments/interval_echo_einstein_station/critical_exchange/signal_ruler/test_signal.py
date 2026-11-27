import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from . import measure


def test_receiver_excess_on_an_entangled_non_slater_source():
    n = 4
    h = np.diag(np.full(n - 1, -0.5), 1)
    h += h.T.copy()
    spin = np.zeros((2**n, 2**n))
    for word in range(2**n):
        for site in range(n - 1):
            mask = 3 << site
            if (word & mask).bit_count() == 1:
                spin[word ^ mask, word] = -0.5
    # A superposition of two disjoint occupied pairs has Slater rank two.
    psi = np.zeros(16, complex)
    psi[3], psi[12] = 0.6, 0.8j
    rho = np.outer(psi, psi.conj())
    remaining = rho[::2, ::2] + rho[1::2, 1::2]
    arms = []
    for bit in (0, 1):
        a = np.zeros_like(rho)
        a[bit::2, bit::2] = remaining
        arms.append(a)
    for t in (0, 1.7, 3):
        unitary, single = expm(-1j * t * spin), expm(-1j * t * h)
        one, full = [(unitary @ a @ unitary.conj().T).diagonal().real for a in arms]
        measured = [sum((full - one)[(np.arange(16) & (1 << j)) != 0]) for j in range(n)]
        np.testing.assert_allclose(measured, abs(single[:, 0])**2, atol=2e-14)


def test_native_contrast_control_and_peak_forecasts():
    local = candidate.native_law.seam_term(candidate.KAPPA)
    binary = [0, 1, 3, 4]
    expected = np.diag([0, 8.5, 8.5, 0])
    expected[1, 2] = expected[2, 1] = -0.5
    np.testing.assert_allclose(local[np.ix_(binary, binary)], expected)
    with np.load(measure.HERE / "counts.npz", allow_pickle=False) as archive:
        data, counts, _ = measure.run(archive)
    assert len(counts) == 3 and data["tested_arm_probabilities"] == 2916
    assert max(data["maximum_absolute_signal_error"].values()) < 2e-13
    assert data["maximum_single_arm_change_for_rotated_source"] > 0.004
    for row in data["peaks"]:
        assert all(row["peak_clock"] in obs["compatible_peak_ages_99_percent"] for obs in row["sources"].values())
        if row["separation"] in (1, 4):
            assert all(obs["compatible_peak_ages_99_percent"] == [row["peak_clock"]] for obs in row["sources"].values())
        if row["separation"] in (16, 51, 55):
            assert all(len(obs["compatible_peak_ages_99_percent"]) >= 3 for obs in row["sources"].values())
    short = data["four_site_at_time_three"]
    assert abs(short["open_quadratic"][0]["excess"] - short["open_quadratic"][1]["excess"]) < 2e-14
    assert short["open_native_blank_contacts"][0]["excess"] > 0
    assert short["open_native_blank_contacts"][1]["excess"] < -0.98
    assert abs(short["literal_spin_ring_closure"][0]["excess"] - short["literal_spin_ring_closure"][1]["excess"]) > 0.63
