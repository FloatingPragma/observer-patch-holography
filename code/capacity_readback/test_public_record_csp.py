from itertools import product

from public_record_csp import public_global_sections_csp


def brute_sections(observers, interfaces):
    ids = sorted(observers)
    result = []
    for values in product(*(observers[observer] for observer in ids)):
        candidate = dict(zip(ids, values, strict=True))
        if all(
            interface["left_readout"][candidate[interface["left_observer"]]]
            == interface["right_readout"][candidate[interface["right_observer"]]]
            for interface in interfaces
        ):
            result.append(candidate)
    return result


def test_csp_matches_cartesian_enumerator_on_generic_noninjective_readouts():
    observers = {
        "a": ["a0", "a1", "a2"],
        "b": ["b0", "b1"],
        "c": ["c0", "c1"],
    }
    interfaces = [
        {
            "left_observer": "a",
            "right_observer": "b",
            "left_readout": {"a0": "x", "a1": "x", "a2": "y"},
            "right_readout": {"b0": "x", "b1": "y"},
        },
        {
            "left_observer": "b",
            "right_observer": "c",
            "left_readout": {"b0": "u", "b1": "v"},
            "right_readout": {"c0": "u", "c1": "v"},
        },
    ]
    assert public_global_sections_csp(observers, interfaces) == brute_sections(observers, interfaces)


def test_csp_handles_connected_twenty_four_atom_diagram_without_cartesian_blowup():
    atoms = [f"x{i}" for i in range(24)]
    observers = {f"o{i}": atoms for i in range(12)}
    interfaces = [
        {
            "left_observer": f"o{i}",
            "right_observer": f"o{i + 1}",
            "left_readout": {atom: atom for atom in atoms},
            "right_readout": {atom: atom for atom in atoms},
        }
        for i in range(11)
    ]
    sections = public_global_sections_csp(observers, interfaces)
    assert len(sections) == 24
    assert all(len(set(section.values())) == 1 for section in sections)

