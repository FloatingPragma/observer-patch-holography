"""Shared aperture, cyclic spectral equivalence, exact selection and blind controls."""

import json

import numpy as np
from ..apparatus import native_graph
from ..measure import ONE
from .instruments import extract, material, loan, BASES
from .measure import HERE
from .spectrum import exact_moments, symmetry


def test_real_material_emission_and_full_register_phase_preserving_extraction():
    rng = np.random.default_rng(309)
    state = rng.normal(size=(8, 3)) + 1j * rng.normal(size=(8, 3))
    for side in ("hub", "tip"):
        out, storage, rejected = extract(side)
        assert not (set(storage) & set(rejected))
        phase = np.diag(out[:, ONE])
        np.testing.assert_allclose(out @ state, phase[:, None] * state[ONE], atol=1e-12)
        for name, M in material(side).items():
            gram = M[ONE].conj().T @ BASES[name]
            np.testing.assert_allclose(abs(gram), np.eye(3) / np.sqrt(3), atol=1e-12)
    assert len(loan()) == 2


def test_both_conjugate_characters_have_identical_spectral_weights_and_distinct_uniform_moment():
    h = native_graph()
    g, _ = symmetry(h)
    assert g["cyclic_ranks"] == [1200, 1192, 1192]
    m = exact_moments(h)
    assert m["power_of_2H"] == 10 and m["difference"] == 384


def test_measured_qutrit_transfer_and_local_selection_rules():
    r = json.loads((HERE / "results.json").read_text())
    for b in r["observed"]:
        assert b["entanglement_fidelity_lower_99"] > 0.73
        assert all(sum(n) == r["starts_per_setting"] for _, _, _, n in b["counts"])
    assert r["counted_dephasing_control"]["entanglement_fidelity_lower_99"] < 0.2
    for section in r["joint_character_read"]:
        assert section["forbidden_count"] == 0
        assert section["allowed_net_onehot"] == 9
        assert section["original_accepted_allowed"] > 120000
        for i, a, b, _, n in section["rows"]:
            j = 0 if a == "all_unwound" else a
            assert sum(n) == r["starts_per_setting"]
            assert (n[0] > 0) == ((b + j - i) % 3 == 0)


def test_independent_phase_conventions_preserve_each_simulated_count():
    f = json.loads((HERE / "frames.json").read_text())
    assert f["settings"] == 90 and f["max_bank_probability_error"] < 3e-17
    same, transported, wrong = f["character_counts_shared_transported_untranslated"]
    assert same == transported == {"accepted": 245233, "forbidden": 0}
    assert wrong["forbidden"] > 200000
    for _, _, _, _, _, original, untranslated in f["rows"]:
        assert sum(original) == sum(untranslated) == f["starts"]
