"""Exact sparse local-singlet screening trial; no eigenvalue approximation."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib
import json
from static_colour_probe import current, fermion_hop

HERE = Path(__file__).resolve().parent


def apply_current(vector, sites, has_q, has_a, a, b):
    out = {}
    for state, coefficient in vector.items():
        for z, c in current(state, sites, has_q, has_a, a, b).items():
            out[z] = out.get(z, 0)+coefficient*c
    return {z: c for z, c in out.items() if c}


def trial(d):
    b = min((x for x in range(4) if x not in [0, d]), key=lambda x: ((-1)**x, x))
    vector = {}
    for q, anti in product(range(3), repeat=2):
        # epsilon(q,b,c) two-quark antitriplet at0, fundamental quark
        # contracted with external antiquark atd, and one local baryon.
        modes = [c for c in range(3) if c != q]+[3*d+anti]+[3*b+c for c in range(3)]
        inversions = sum(modes[i] > modes[j] for i in range(6) for j in range(i+1, 6))
        state = sum(1 << m for m in modes), q, anti
        vector[state] = (-1)**(q+inversions)
    norm = sum(c*c for c in vector.values())
    assert norm == 9
    # All prefixes are exact local SU3 singlets, including the probes.
    for cut in [1, 2, 3, 4]:
        trace = 2+int(d < cut)+3*int(b < cut)+1-int(d < cut)
        assert trace % 3 == 0
        for a, c in product(range(3), repeat=2):
            value = apply_current(vector, range(cut), True, d < cut, a, c)
            expected = {z: trace//3*v for z, v in vector.items()} if a == c else {}
            assert value == expected
    hopped = {}
    for (word, q, anti), amplitude in vector.items():
        for x, c in product(range(3), range(3)):
            for u, v in [(x, x+1), (x+1, x)]:
                w, coefficient = fermion_hop(word, 3*u+c, 3*v+c)
                if w is not None:
                    z = w, q, anti
                    hopped[z] = hopped.get(z, 0)-amplitude*coefficient
    assert sum(amplitude*hopped.get(z, 0) for z, amplitude in vector.items()) == 0
    mass_values = {sum((-1)**x*((word >> (3*x)) & 7).bit_count() for x in range(4))
                   for word, _, _ in vector}
    assert mass_values == {2+(-1)**d+3*(-1)**b}
    return dict(separation=d, spectator_baryon_cell=b, exact_norm_squared=norm,
                all_prefix_colour_fields_zero=True, hopping_expectation_exact='0',
                mass_energy_coefficient=mass_values.pop(),
                complete_trial_coefficients=[dict(word=w, probe=q, antiprobe=a, coefficient=c)
                                             for (w,q,a),c in sorted(vector.items())])


def main():
    files = [HERE/'hadrons/STAGGERED_CHAIN4_RESULT.json', HERE/'STATIC_COLOUR_RESULT.json',
             HERE/'static_colour_probe.py']
    pins = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    ref = json.loads(files[0].read_text())
    mass = F(ref['actual_mass_exact'])
    trials = [trial(d) for d in [1,2,3]]
    bounds = []
    for row in ref['results']:
        v = next(e for e in row['complete_sector_records'] if e['N'] == 6)
        vlo = F(v['ground_interval_exact'][0])
        for tr in trials:
            # E_with_probes <= trial energy, with exact vacuum lower bound.
            upper = tr['mass_energy_coefficient']*mass-vlo
            bounds.append(dict(kappa=row['kappa'], separation=tr['separation'],
                               potential_upper_bound_exact=str(upper),
                               potential_upper_bound=float(upper),
                               independent_of_bare_string_length_cost=True))
    assert pins == {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    out = dict(schema='oph.obstruction_repair.colour_screening_trial.v1',
               exact_local_Gauss_and_prefix_flux_checks=True, fixed_fermion_number=6,
               trial_states=trials, potential_upper_bounds=bounds,
               trial_not_claimed_ground_state=True,
               physical_intermediate_string_tension_solved=False, pins=pins)
    (HERE/'COLOUR_SCREENING_BOUND.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'exact_screened_trials':len(trials),'potential_upper_bounds':len(bounds)}))


if __name__ == '__main__':
    main()
