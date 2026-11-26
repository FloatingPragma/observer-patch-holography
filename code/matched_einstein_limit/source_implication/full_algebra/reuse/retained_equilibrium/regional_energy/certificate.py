"""Check both sources using independently registered cell Hamiltonian energies."""

import json
from pathlib import Path

import sympy as sp
from mpmath import iv
from mpmath.libmp import round_ceiling, round_floor, to_int

from .....directional_stress import bridge

HERE = Path(__file__).parent
SOURCE = HERE.parents[3]


def enclosure(text):
    return iv.mpf(text.strip("[]").split(","))


def band(v):
    lo, hi = (v * 10**8)._mpi_
    return [to_int(lo, round_floor), to_int(hi, round_ceiling)]


def verify():
    iv.dps = 42
    record = json.loads((HERE / "effects.json").read_text())
    source = bridge.action.front.inputs()
    base = json.loads((SOURCE / "full_algebra/certificate.json").read_text())
    signs = base["signed_transfers"]
    entropy_o, entropy_c = map(
        enclosure, (signs["ordinary_entropy"], signs["pairing_entropy"])
    )
    energy_o = enclosure(signs["ordinary_energy"])
    eta = enclosure(base["calibrated_gain"])
    ratio = entropy_c / entropy_o
    charge = {}
    for q in ("-1", "1"):
        effects = [
            [
                (sp.Matrix(re) + sp.I * sp.Matrix(im)) / record["scale"]
                for re, im in sector
            ]
            for sector in record["rows"][q]
        ]
        errors = record["cell_error_trillionths"][q]
        assert sorted(record["marks"]) == [0, 3]

        def read(indices):
            entries = [sum((sect[i] for i in indices), sp.zeros(4)) for sect in effects]
            flats = [m.trace() / 4 for m in entries]
            # The scalar parts distinguish only the reversal masses. Both target
            # rotations and the correlation intervention preserve them exactly.
            traceless = [m - a * sp.eye(4) for m, a in zip(entries, flats, strict=True)]
            if any(m != sp.zeros(4) for m in traceless):
                ans = bridge.action.source_jets(traceless, source)
            else:
                ans = {"ordinary": iv.mpf(0), "pairing": iv.mpf(0)}
            s = flats[1] - flats[0] + int(0 in indices)
            ans["ordinary"] += energy_o * bridge.action.front.full.ri(s)
            eps = sum(errors[i] for i in indices) * iv.mpf("1e-12")
            ans["ordinary"] += iv.mpf([-999, 999]) * eps
            ans["pairing"] += iv.mpf([-299, 299]) * eps
            return ans

        all_cells = [read([i]) for i in range(80)]
        core = read([0, 3])
        total = read(list(range(80)))
        twelve = read(record["chart"])
        response = [v["pairing"] - ratio * v["ordinary"] for v in all_cells]
        leakage = sum(
            max(s.b, iv.mpf(0)) for i, s in enumerate(response) if i not in (0, 3)
        )
        free_o = [v["ordinary"] for i, v in enumerate(all_cells) if i not in (0, 3)]
        free_c = [v["pairing"] for i, v in enumerate(all_cells) if i not in (0, 3)]
        lower_o = core["ordinary"] + sum(min(z.a, iv.mpf(0)) for z in free_o)
        upper_o = core["ordinary"] + sum(max(z.b, iv.mpf(0)) for z in free_o)
        lower_c = core["pairing"] + sum(min(z.a, iv.mpf(0)) for z in free_c)
        upper_c = core["pairing"] + sum(max(z.b, iv.mpf(0)) for z in free_c)
        source_peak = (
            max(response[0].b, response[3].b)
            + max(iv.mpf(0), min(response[0].b, response[3].b))
            + leakage
        )
        pair_peak = (core["pairing"] - ratio * core["ordinary"]).b + leakage
        ratio_region = twelve["pairing"] / twelve["ordinary"]
        assert lower_o.a > 0
        charge[q] = dict(
            head_ordinary=[band(all_cells[i]["ordinary"]) for i in (0, 3)],
            head_collar=[band(all_cells[i]["pairing"]) for i in (0, 3)],
            head_isentropic=[band(response[i]) for i in (0, 3)],
            positive_exterior_allowance=band(leakage)[1],
            profile_isentropic_upper_one_mark_peaks=band(source_peak)[1],
            profile_isentropic_upper_both_marks_peak=band(pair_peak)[1],
            both_mark_profiles_ordinary=[band(lower_o)[0], band(upper_o)[1]],
            both_mark_profiles_collar=[band(lower_c)[0], band(upper_c)[1]],
            both_mark_profiles_abs_collar_over_ordinary_upper=band(
                max(abs(lower_c.a), abs(upper_c.b)) / lower_o.a
            )[1],
            twelve_cell_ordinary=band(twelve["ordinary"]),
            twelve_cell_collar=band(twelve["pairing"]),
            twelve_cell_pairing_over_ordinary=band(ratio_region),
            full_flight_ordinary=band(total["ordinary"]),
            full_flight_collar=band(total["pairing"]),
        )
    # Independent vacancy calibration holds its twelve-cell energy constant.
    # Its ordinary response is transverse on the constant-volume sheet.
    quotient = charge["-1"]["twelve_cell_pairing_over_ordinary"]
    k = iv.mpf([str(x) for x in quotient]) * iv.mpf("1e-8")
    other = charge["1"]

    def decoded(key):
        return iv.mpf([str(x) for x in other[key]]) * iv.mpf("1e-8")

    held = decoded("twelve_cell_collar") - k * decoded("twelve_cell_ordinary")
    entropy = entropy_c - k * entropy_o
    return dict(
        units="outward endpoints in hundred-millionths",
        tracers=charge,
        fixed_vacancy_region=dict(
            excess_energy=band(held),
            source_entropy=band(entropy),
            survey_area=band(-entropy / eta),
        ),
    )


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
