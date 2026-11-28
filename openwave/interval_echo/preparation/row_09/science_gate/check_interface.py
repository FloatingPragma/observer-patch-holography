"""Independently reduce the signed-interface coincidence records."""

from hashlib import sha256
import json
from math import cos, fsum, isclose, sin
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
GUIDE = (
    ROOT
    / "model/candidate_laws/interval_echo/openwave/row_09_neutral/critical_emission"
)
INTERFACE = GUIDE / "cycle_energy/interface_tag"


def require(test, message):
    if not test:
        raise ValueError(message)


def compare(data, old):
    early = old["parent_refused_mode_refused_admitted"]
    beta = complex(*data["original_low_mode_amplitude_real_imaginary"])
    unlike = data["one_mode_matrix_element_of_total_unlike_projector"]
    mixed = complex(*data["coherent_end_contraction_real_imaginary"])
    t = old["guard_only_terminal_charge_control"]["dwell"]
    u = early[2]
    signed = u * sin(t / 2) ** 2 * unlike / 2
    direct = u * abs(beta - (1 - cos(t / 2)) * mixed) ** 2
    restored = u * sin(t / 2) ** 2 * abs(mixed) ** 2 / 4
    require(
        0
        < data["charged_tag_probability_per_original_lower_after_deliberate_return"]
        <= restored,
        "No certified positive charged-arm restoration",
    )
    for arm, branch in (("direct", 0.0), ("charge_conditioned_return", restored)):
        vector = [
            *early[:2],
            0.0,
            branch,
            signed - branch,
            direct,
            u - 2 * signed - direct,
            branch,
            signed - branch,
        ]
        record = data["guarded_and_return_arms"][arm]
        p, c = record["probabilities"], record["seeded_draws"]
        require(
            len(c) == 9
            and all(type(n) is int and n >= 0 for n in c)
            and sum(c) == 2**42,
            "An original outcome is missing",
        )
        require(
            len(p) == 9
            and min(p) >= 0
            and isclose(fsum(p), 1, abs_tol=8e-15, rel_tol=0),
            "Invalid complete joint law",
        )
        require(
            all(
                isclose(a, b, abs_tol=2e-17, rel_tol=5e-10)
                for a, b in zip(p, vector, strict=True)
            ),
            "A native charge or optical-tag amplitude does not reproduce",
        )
        if arm == "direct":
            require(
                p[3] == p[7] == 0 and c[3] == c[7] == 0,
                "Original guarded direct tag must reject both charged reports",
            )
        else:
            require(
                c[3] > 0 and c[7] > 0,
                "Charged return control must resolve both prior signs",
            )
    print(
        f"{data['dimers']} dimers: each signed read {signed:.8g}, direct tag zero, charged-arm restored tag {restored:.8g} per original"
    )


def main():
    archive = json.loads((INTERFACE / "results.json").read_bytes())
    for path, key in (
        (INTERFACE / "PROTOCOL.md", "protocol_sha256"),
        (GUIDE / "emitter.json", "emission_sha256"),
        (GUIDE / "cycle_energy/results.json", "cycle_energy_sha256"),
    ):
        require(
            sha256(path.read_bytes()).hexdigest() == archive[key],
            f"Input digest differs for {path.name}",
        )
    earlier = json.loads((GUIDE / "cycle_energy/results.json").read_bytes())
    panels = {p["dimers"]: p for p in earlier["panels"] if p["schedule"] == "native"}
    require(
        archive["receiving_and_neighboring_charge_per_valid_record"]
        == [[-2, 2], [0, 0], [2, -2]],
        "The independently reported compensating charge was lost",
    )
    for data in archive["cases"]:
        compare(data, panels.pop(data["dimers"]))
    require(not panels, "A held guide length is missing")
    print(
        "The direct tag selects the balanced continuation. The charged-arm calibration restores it after recording the sign."
    )


if __name__ == "__main__":
    main()
