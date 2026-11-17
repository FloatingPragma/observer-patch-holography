"""Independent exact verification of the frozen reduced-source certificates.

The producer's CAR compression matrices are consumed, not independently
reconstructed here. Inertia, dimensional completeness and energy subtraction
are independently checked with the root's fraction-free verifier.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
from nuclear_threshold import inertia

HERE = Path(__file__).resolve().parent


def dimension(highest):
    result = F(1)
    for i in range(len(highest)):
        for j in range(i+1, len(highest)):
            result *= F(highest[i]-highest[j]+j-i, j-i)
    assert result.denominator == 1
    return int(result)


def main():
    source = HERE/'hadrons/SPIN_FLAVOUR_SPECTRUM.json'
    manifest = HERE/'hadrons/SPIN_FLAVOUR_MANIFEST.json'
    result = json.loads(source.read_text())
    paths = [source, manifest, HERE/'nuclear_threshold.py']
    mass, actual = F(result['computational_mass']), F(result['source_mass'])
    error = F(result['mass_operator_error_bound'])
    assert error == 24*abs(mass-actual)
    t = F(result['t'])
    blocks, dimensions = {}, {}
    for record in result['results'][0]['all_competing_irrep_records']:
        weight = tuple(record['U4_highest_weight'])
        path = HERE/'hadrons/spin_flavour_cache'/('block_'+'_'.join(map(str, weight))+'.json')
        paths.append(path)
        b = json.loads(path.read_text())
        assert b['multiplicity'] == record['multiplicity'] == len(b['G'])
        assert dimension(weight) == record['SU4_dimension'] == b['SU4_dimension']
        G = b['G']
        assert all(G[i][j] == G[j][i] for i in range(len(G)) for j in range(len(G)))
        assert inertia([[F(-v) for v in row] for row in G]) == len(G)
        for M in b['components'].values():
            assert all(M[i][j] == M[j][i] for i in range(len(G)) for j in range(len(G)))
        blocks[weight] = b
        N = sum(weight)
        dimensions[N] = dimensions.get(N, 0)+dimension(weight)*len(G)
    assert len(blocks) == 43
    for N, d in dimensions.items():
        k = N//3
        assert d == dimension([3]*k+[0]*(8-k))
    pins = {str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    checks, conclusions = [], []
    for member in result['results']:
        kap = F(member['kappa'])
        sectors = {}
        for record in member['all_competing_irrep_records']:
            weight = tuple(record['U4_highest_weight'])
            b = blocks[weight]
            G, M, K, E = b['G'], b['components']['M'], b['components']['K'], b['components']['E']
            n = len(G)
            H = [[mass*M[i][j]+t*K[i][j]+kap*F(E[i][j], 6) for j in range(n)] for i in range(n)]
            lo, hi = map(F, record['computational_inertia_interval_exact'])
            physical = tuple(map(F, record['lowest_energy_exact_interval']))
            if record['exact_scalar_zero']:
                assert n == 1 and H == [[0]] and lo == hi == 0 and physical == (0, 0)
                counts = None
            else:
                counts = [inertia([[H[i][j]-z*G[i][j] for j in range(n)] for i in range(n)]) for z in [lo, hi]]
                assert counts == [0, record['certified_eigenvalue_count_in_ground_interval']]
                assert physical == (lo-error, hi+error)
            checks.append(dict(kappa=str(kap), highest_weight=list(weight), dimension=n, inertia=counts))
            N = sum(weight)
            if N not in sectors:
                sectors[N] = physical
            else:
                sectors[N] = min(sectors[N][0], physical[0]), min(sectors[N][1], physical[1])
        vacuum = sectors[12]
        assert all(vacuum[1] < e[0] for N,e in sectors.items() if N != 12)
        gap = sectors[9][0]-vacuum[1], sectors[9][1]-vacuum[0]
        binding = 2*sectors[9][0]-sectors[6][1]-vacuum[1], 2*sectors[9][1]-sectors[6][0]-vacuum[0]
        assert gap == tuple(map(F, member['baryon_gap_exact_interval']))
        assert binding == tuple(map(F, member['two_baryon_binding_exact_interval']))
        assert binding[0] > 0
        conclusions.append(dict(kappa=str(kap), vacuum_N=12,
                                exact_baryon_gap=list(map(str, gap)), exact_two_baryon_binding=list(map(str, binding))))
        print('verified', str(kap), 'all43 blocks and binding', flush=True)
    assert pins == {str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    out = dict(schema='oph.obstruction_repair.spin_flavour_certificate_review.v1',
               CAR_matrices_consumed_not_independently_reconstructed=True,
               independent_fraction_free_inertia=True, all43_Gram_forms_positive=True,
               Weyl_dimension_census_verified=dimensions, block_member_checks=len(checks),
               checks=checks, conclusions=conclusions,
               upper_number_particle_hole_transfer_consumed_not_reproved=True,
               finite_two_baryon_binding_component_verified=True,
               physical_nucleon_or_nuclear_row_solved=False, pins=pins)
    (HERE/'SPIN_FLAVOUR_CERTIFICATE_REVIEW.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
