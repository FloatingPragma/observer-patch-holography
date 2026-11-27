"""Exact finite observer distributions for a cycle of four repaired interfaces."""
import argparse
import json
import math
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

from .regional import PartitionRepair

HERE = Path(__file__).resolve().parent
EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
REGIONS = tuple(r for n in range(5) for r in combinations(range(4), n))
NAMES = ("one_face", "held_faces_and_pairing", "correlation_only", "volume_only")
MASS = 5**4 * 3**4
AMPLITUDE = F(1, 16)


def populations():
    """Enumerate admissible global records once; no Gibbs weights enter the repair."""
    populations = {r: {} for r in REGIONS}
    for bonds in product(range(5), repeat=4):
        a, b, c, _ = (5*int(i != 0)-4 for i in bonds)
        for volumes in product(range(2), repeat=4):
            # The probabilities of the two volume-bit outcomes are 2/3 and 1/3.
            m = 2 ** (4 - sum(volumes))
            dt = [a, b+c+a*c, a*c, 3*volumes[0]-1]
            if 1+AMPLITUDE*min(dt) <= 0:
                raise ArithmeticError("tilt leaves the positive law")
            for region, rows in populations.items():
                key = PartitionRepair.read_region(bonds, volumes, EDGES, region)
                entry = rows.setdefault(key, [0]*5)
                entry[0] += m
                for j, g in enumerate(dt, 1):
                    entry[j] += m*g
    return populations


def labels(region):
    output = [("edge", j, side) for j, e in enumerate(EDGES)
              for side, cell in enumerate(e) if cell in region]
    return output + [("volume", node) for node in sorted(region)]


def check_restrictions(tables):
    joins = 0
    for big, table in tables.items():
        slots = labels(big)
        for small in REGIONS:
            if not set(small).issubset(big) or small == big:
                continue
            indices = [slots.index(key) for key in labels(small)]
            seen = defaultdict(lambda: [0]*5)
            for record, count in table.items():
                observed = tuple(record.word[i] for i in indices)
                for i, value in enumerate(count):
                    seen[observed][i] += value
            assert dict(seen) == {r.word: val for r, val in tables[small].items()}
            joins += 1
    return joins


def v2(integer):
    return (integer & -integer).bit_length()-1


def inspect(region, rows):
    # Compare an independently observed log probability with the separately
    # counted aperture area and volume. Same odd factor leaves exact dyadic ratios.
    odd_parts = {ct[0] >> v2(ct[0]) for ct in rows.values()}
    assert len(odd_parts) == 1 and len({r.word for r in rows}) == len(rows)
    assert sum(ct[0] for ct in rows.values()) == MASS
    assert len({v2(ct[0])-r.open_boundary+r.volume for r,ct in rows.items()}) == 1
    arms = {}
    for j, name in enumerate(NAMES, 1):
        assert sum(ct[j] for ct in rows.values()) == 0
        ds = F(-sum(ct[j]*v2(ct[0]) for ct in rows.values()), MASS)
        da = F(sum(ct[j]*r.open_boundary for r,ct in rows.items()), MASS)
        dv = F(sum(ct[j]*r.volume for r,ct in rows.items()), MASS)
        variation = AMPLITUDE*sum(abs(F(ct[j], 2*MASS)) for ct in rows.values())
        kl = math.fsum(float((ct[0] + AMPLITUDE*ct[j])/MASS) *
                       math.log1p(float(AMPLITUDE*ct[j]/ct[0])) for ct in rows.values())
        assert ds + da - dv == 0
        arms[name] = dict(entropy_slope_in_bits=str(ds), aperture_slope=str(da), volume_slope=str(dv),
                          finite_deficit_nats=kl, pinsker_lower_nats=str(2*variation**2))
    return dict(cells=list(region), visible_words=len(rows), arms=arms)


def finite_repair():
    k = PartitionRepair.bond_kernel()
    ref = [F(1, 5)]*5
    assert all(sum(ref[i]*k[i][j] for i in range(5)) == ref[j] for j in range(5))
    assert all(k[i][j] == k[j][i] for i in range(5) for j in range(5))
    two = [[sum(k[i][m]*k[m][j] for m in range(5)) for j in range(5)] for i in range(5)]
    floor = min(min(row) for row in two)
    vol = [[F(sum(PartitionRepair.advance_volume(i, c) == j for c in range(4)), 4)
            for j in (0, 1)] for i in (0, 1)]
    assert F(2, 3)*vol[0][1] == F(1, 3)*vol[1][0]
    # With A = 4K the five-cycle satisfies A² − 3A + I = 11ᵀ.
    # On zero-sum modes this gives the two roots (3 ± sqrt(5))/8 for K.
    if not all(16*two[i][j] - 12*k[i][j] + int(i == j) == 1
               for i in range(5) for j in range(5)):
        raise ArithmeticError("Interface repair fails the five-cycle spectral identity")
    if sum(vol[i][i] for i in range(2)) - 1 != F(1, 4):
        raise ArithmeticError("Unexpected volume relaxation factor")
    p = [F(1), F(0), F(0), F(0), F(0)]
    pvol = F(0)
    blanks = []
    for tick in range(65):
        if tick in (0, 4, 16, 64):
            blanks.append(dict(sweeps=tick, bond_open_probability=str(1-p[0]), volume_one_probability=str(pvol),
                               bond_TV=float(sum(abs(x-r) for x,r in zip(p,ref,strict=True))/2)))
        p = [sum(p[i]*k[i][j] for i in range(5)) for j in range(5)]
        pvol = (1-pvol)*vol[0][1]+pvol*vol[1][1]
    return dict(bond_kernel=[[str(p) for p in row] for row in k], two_step_probability_floor=str(floor),
                bond_two_step_TV_contraction_upper=str(1-5*floor),
                volume_kernel=[[str(p) for p in row] for row in vol], from_all_blank=blanks)


def biased_control():
    # A modified acceptance on one word gives it twice the reversible weight.
    weights = (1, 2, 1, 1, 1)
    proposal = PartitionRepair.bond_kernel()
    kernel = []
    for i, row in enumerate(proposal):
        accepted = [p*min(F(1), F(weights[j],weights[i])) if i != j else F(0)
                    for j,p in enumerate(row)]
        accepted[i] = 1-sum(accepted)
        kernel.append(accepted)
    assert all(weights[i]*kernel[i][j] == weights[j]*kernel[j][i]
               for i in range(5) for j in range(5))
    left = [F(sum(w for (a,_),w in zip(PartitionRepair.bond_words,weights,strict=True) if a==v), sum(weights))
            for v in range(3)]
    return dict(endpoint_frequency=list(map(str,left)), aperture_slope=0, volume_slope=0,
                biased_pair_fixed_area_entropy_slope_nats=math.log(float(left[2]/left[1])))


def refining_faces():
    # Count-only example: d or d+1 legal endpoint partners at each face.
    cases = []
    for k in (4, 8, 16, 32):
        ell, area = F(1,k), F(1,k**4)
        d, faces = k**4, k**2
        p = F((d+1)**2, d*d+(d+1)**2)
        shift = ell**2/8
        q = p+shift
        da = faces*area*shift
        divergence = faces*(float(q)*math.log1p(float(shift/p))
                            +float(1-q)*math.log1p(float(-shift/(1-p))))
        lower = 2*da**2/(faces*area**2)
        assert da == ell**4/8 and lower/ell**4 == F(k*k, 32)
        cases.append(dict(ell=str(ell), face_area=str(area), boundary_faces=faces,
                          gain=math.log1p(1/d)/float(area), area_contrast=str(da),
                          deficit_over_ell4=divergence/float(ell**4), lower_over_ell4=str(lower/ell**4)))
    return cases


def unequal_footprint_control(calibration):
    ds = F(calibration["entropy_slope_in_bits"])
    da = F(calibration["aperture_slope"])
    # An otherwise identical aperture with twice the measured footprint.
    return [str(-ds/(a*da)) for a in (1, 2)]


def check_saved(actual, saved, path="receipt"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or actual.keys() != saved.keys():
            raise ValueError(f"{path}: receipt keys differ")
        for key, value in actual.items():
            check_saved(value, saved[key], f"{path}.{key}")
    elif isinstance(actual, (tuple, list)):
        if not isinstance(saved, list) or len(actual) != len(saved):
            raise ValueError(f"{path}: list length differs")
        for i, (left, right) in enumerate(zip(actual, saved, strict=True)):
            check_saved(left, right, f"{path}[{i}]")
    elif isinstance(actual, float):
        if not math.isclose(actual, saved, rel_tol=1e-12, abs_tol=1e-14):
            raise ValueError(f"{path}: measured {actual!r}, archived {saved!r}")
    elif actual != saved:
        raise ValueError(f"{path}: measured {actual!r}, archived {saved!r}")


def region_receipt():
    tables = populations()
    reads = [inspect(r, entries) for r,entries in tables.items()]
    calibration = next(row for row in reads if row["cells"] == [0])["arms"]["one_face"]
    gain = -F(calibration["entropy_slope_in_bits"])/F(calibration["aperture_slope"])
    assert gain == 1
    return dict(finite_repair=finite_repair(), global_words=5**4*2**4,
                proper_nonempty_subsets=len(reads)-2, readbacks=reads,
                strict_restrictions_verified=check_restrictions(tables),
                tilt_amplitude=str(AMPLITUDE), calibrated_gain_bits_per_area_unit=str(gain),
                unequal_area_control_bits_per_area=unequal_footprint_control(calibration),
                changed_law=biased_control(), finite_area_refinement=refining_faces())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    output = region_receipt()
    path = HERE / "regional_results.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(output, indent=2)+"\n")
    else:
        saved = json.loads(path.read_text())
        check_saved(output, saved)
        print("Verified all overlapping regions, repair kernels and signed transfers.")
