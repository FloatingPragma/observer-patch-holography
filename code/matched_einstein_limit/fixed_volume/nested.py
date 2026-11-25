"""Compare the fixed-radius volume response with the matched-volume area in fresh counter banks."""

import json
from pathlib import Path


def read(archive):
    panels = []
    for panel in archive["panels"]:
        reference = panel["flat_comparison"]["regions"]
        cases = []
        for case in panel["cases"]:
            records = case["record"]["regions"]
            volume_counts = (
                records["fixed_radius"]["held_interior_count"]
                - reference["fixed_radius"]["held_interior_count"]
            )
            area_counts = (
                records["volume_threshold"]["surface_count"]
                - reference["volume_threshold"]["surface_count"]
            )
            # The two counter densities are N/(4πr³/3) and N/(4πr²).
            ratio = volume_counts / area_counts if area_counts else None
            cases.append(
                dict(
                    curvature=case["prepared_sectional_curvature"],
                    volume_counts=volume_counts,
                    area_counts=area_counts,
                    three_volume_over_radius_area=ratio,
                )
            )
        panels.append(dict(radius=panel["flat_radius"], cases=cases))
    return panels


if __name__ == "__main__":
    data = json.loads(Path(__file__).with_name("measurements.json").read_text())
    print(json.dumps(read(data), indent=2))
