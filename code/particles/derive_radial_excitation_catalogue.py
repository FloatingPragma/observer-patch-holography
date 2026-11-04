#!/usr/bin/env python3
"""Exact source-qualified catalogue; no physical species identification."""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
import hashlib
import importlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools import strict_json
from model.radial_interacting_core.core import core_bounds, encode
from model.source_mass_shell import units

OUTPUT = ROOT/'code/particles/runs/catalogue/radial_excitation_catalogue.json'
MARKDOWN = ROOT/'code/particles/RADIAL_EXCITATION_CATALOGUE.md'
PARENTS = ('bounded_defect_index', 'radial_defect_gap', 'radial_interacting_core', 'source_mass_shell')
OWN_SOURCES = ('code/particles/derive_radial_excitation_catalogue.py',
               'code/particles/test_radial_excitation_catalogue.py')


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def decimal_bound(value, *, lower=True, precision=14):
    value = F(value)
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_FLOOR if lower else ROUND_CEILING
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def display_interval(values):
    return '['+decimal_bound(values[0])+', '+decimal_bound(values[1], lower=False)+']'


def parent_receipts():
    receipts, pins = {}, {}
    for name in PARENTS:
        relative = f'model/{name}/results.json'
        path = ROOT/relative
        if not path.is_file():
            raise ValueError('complete final source receipt required: '+relative)
        value = strict_json.load(path)
        expected = importlib.import_module(f'model.{name}.packet').source_bindings()
        if (value.get('schema') != 'oph.'+name+'.v1'
                or value.get('all_pass') is not True
                or value.get('source_bindings') != expected):
            raise ValueError('source receipt is incomplete, stale or failing: '+name)
        for source, digest in expected.items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != digest:
                raise ValueError('stale parent source: '+source)
            if source in pins and pins[source] != digest:
                raise ValueError('parent manifests disagree: '+source)
            pins[source] = digest
        pins[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        receipts[name] = value
    for source in OWN_SOURCES:
        pins[source] = hashlib.sha256((ROOT/source).read_bytes()).hexdigest()
    return receipts, dict(sorted(pins.items()))


def spectral_probe_bounds(core):
    """Normalize each low-pole residue separately, retaining its full trace."""
    if (core['radius'] != 1 or core['dimensions'] != 3 or core['sites'] != 27
            or core['each_pole_matrix_rank'] != 2):
        raise ValueError('rank-two pole on the complete radius-one cubic source required')
    t, alpha, a = (F(core[key]) for key in
                    ('projection_angle_bound', 'origin_weight_lower', 'spacing'))
    if t < 0 or not 0 < alpha <= 1 or a <= 0:
        raise ValueError('invalid source angle, origin weight or spacing')
    sigma = F(101, 1000)
    if 26*(1-alpha) >= sigma*sigma or alpha <= sigma or 6*t >= 1:
        raise ValueError('the exact Fourier and trace enclosures do not apply')
    r, trace_upper = alpha-sigma, 1+6*t
    center = (alpha-6*t)/trace_upper
    fourier = (r*r/27-6*t)/trace_upper
    block = r*r/54-3*t
    if min(center, fourier, block) <= 0 or 27*fourier > 1 or center > 1:
        raise ValueError('positive normalized source visibility was not certified')
    return encode({
        'position_definition': 'P_x=Tr(M_xx)/Tr(M), both internal source components',
        'Fourier_definition': 'P_k=Tr(F_k^* M F_k)/Tr(M), 27 orthogonal Fourier probes',
        'normalization': 'each addition or removal low-pole residue M separately',
        'reference': 'actual interacting ground mixture G/2',
        'each_residue_rank': 2, 'residue_operator_error': 3*t,
        'total_pole_weight_interval': [1-6*t, trace_upper],
        'center_probability_lower': center,
        'outside_center_probability_upper': 1-center,
        'RMS_position_squared_about_origin_upper_model_length_squared': 3*a*a*(1-center),
        'each_Fourier_probe_probability_lower': fourier,
        'each_Fourier_probe_probability_upper': 1-26*fourier,
        'all_Fourier_block_minimum_singular_value_lower': block,
        'Fourier_lower_auxiliary_r': r,
        'Fourier_probe_wave_numbers': 'p_j=2*pi*k_j/(3*a), k_j in {-1,0,1}',
        'Fourier_probe_count': 27, 'internal_rows_per_probe': 2,
        'off_diagonal_Fourier_blocks_retained': True,
        'position_trajectory_inferred': False,
        'globally_conserved_momentum_inferred': False,
    })


def scientific_catalogue():
    core = core_bounds()
    low, high = core['actual_charged_gap_lower'], core['actual_charged_gap_upper']
    probes, si = spectral_probe_bounds(core), units.certificate()
    nulls = {'observed_species': None, 'physical_rest_mass_GeV': None,
             'moving_defect_rest_energy': None, 'conserved_global_momentum': None}
    rows = [{
        'id': 'one_body_zero_doublet', 'name': 'Protected one-body zero doublet',
        'domain': 'one-body radial source, infinite lattice and symmetric finite cubes in the strong-mass regime',
        'energy_kind': 'one-body eigenvalue in the original source energy origin',
        'energy_model_interval': ['0', '0'], 'multiplicity': 2,
        'origin_probability_lower': str(core['origin_weight_lower']),
        'one_body_nonzero_spectral_magnitude_lower': str(core['one_body_gap_lower']),
        'energy_above_interacting_ground': None,
        'source_charge_relative_to_interacting_ground': None,
        'proofs': ['model/bounded_defect_index/PROOF.md', 'model/radial_defect_gap/PROOF.md'],
        **nulls,
    }]
    for number, sign, label in ((26, -1, 'removal'), (28, 1, 'addition')):
        rows.append({
            'id': 'interacting_'+label, 'name': 'Actual interacting charged '+label+' state',
            'domain': 'complete 27-site Fock source at fixed spacing and contact',
            'particle_number': number, 'source_charge_relative_to_interacting_ground': sign,
            'charge_unit': 'conserved source fermion number; physical electric charge not identified',
            'multiplicity': 1, 'residue_rank': 2,
            'energy_kind': 'positive excitation above the actual N=27 ground doublet',
            'energy_model_interval': [str(low), str(high)],
            'signed_retarded_pole_interval': [str(-high), str(-low)] if sign < 0 else [str(low), str(high)],
            'other_charged_excitation_energy_lower': str(core['other_charged_excitations_threshold_lower']),
            'spectral_probe_bounds': probes,
            'momentum_role': 'all Fourier probes see the same localized pole; nonzero off-diagonal mixing',
            'physical_units': si['radial_charged_gap_scale_family'],
            'proofs': ['model/radial_interacting_core/PROOF.md', 'model/defect_core_response/PROOF.md'],
            **nulls,
        })
    for sign, label in ((-1, 'negative'), (1, 'positive')):
        rows.append({
            'id': 'exterior_'+label, 'name': 'Exterior one-body '+label+' branch',
            'domain': 'infinite radial one-body exterior packets and their fixed-v,m continuum scaling',
            'energy_kind': 'signed continuum-shell energy, separate from the interacting-ground Lehmann measure',
            'energy_sign': sign,
            'dispersion_model': ('-' if sign < 0 else '+')+'sqrt(|p|^2+64^2)',
            'rest_energy_magnitude_model': '64',
            'signed_zero_wave_number_energy_model': str(sign*64),
            'multiplicity_at_each_signed_energy': 8,
            'full_folded_fiber_rows': 16, 'microscopic_valleys_retained': 8,
            'four_component_Dirac_copies_full_fiber': 4,
            'momentum_role': 'asymptotic packet angular wave number, not a global conserved radial momentum',
            'packet_position': 'anchor (0,0,R), support n=(0,0,R)+2r+s; R/L tends to infinity',
            'finite_spacing_dispersion': ('-' if sign < 0 else '+')+'sqrt((v/a)^2*sum_j sin(a*p_j)^2+m^2)',
            'continuum_mass_shell': 'E^2-v^2*|p|^2=m^2',
            'supplied_mass_parameter': '64',
            'interaction_included_in_shell': False,
            'same_finite_interacting_reference_used': False,
            'physical_units': si['radial_exterior_one_body_rest_energy_scale_family'],
            'proofs': ['model/source_mass_shell/PROOF.md', 'model/source_mass_shell/CALIBRATION_PROOF.md'],
            **nulls,
        })
    return {
        'source': encode({
            'hopping': 'h_(n,n+e_j)=i*v*sigma_j/(2*a)',
            'texture': 'B_n=m*sum_j (-1)^n_j*(n_j/|n|)*sigma_j for n!=0; B_0=0',
            'interaction': 'U*sum_x (n_up-1/2)*(n_down-1/2)',
            'spacing': core['spacing'], 'velocity': core['velocity'], 'mass': core['mass'],
            'contact': core['centered_contact_coupling'], 'finite_cube_radius': 1,
            'finite_sites': 27, 'one_body_rows': 54, 'full_Fock_dimension': core['full_Fock_dimension'],
            'source_time': 'exp(-i H_model tau)', 'supplied_profile_and_parameters': True,
            'defect_origin_model_coordinates': [0, 0, 0],
            'profile_center_is_a_dynamical_coordinate': False}),
        'interacting_reference': {'particle_number': 27, 'ground_multiplicity': 2,
            'density': 'G/2, normalized actual complete-Fock ground projector',
            'absolute_ground_energy': None, 'full_Fock_matrix_numerically_diagonalized': False},
        'topology': {'proved_index': 2,
            'operator': 'P_- H P_+: inversion-even to inversion-odd states',
            'grading': 'nonlocal spatial inversion (R psi)(n)=psi(-n)',
            'protection': 'inversion-odd compact self-adjoint perturbations preserve the Fredholm index',
            'strong_mass_exact_kernel_dimension': 2,
            'electric_charge_equals_texture_winding': False,
            'internal_chiral_index_equals_reflection_index': False,
            'arbitrary_symmetry_breaking_perturbations_pin_zero': False},
        'rows': rows, 'physical_units': si,
        'scope': {'complete_low_energy_finite_subspace_classified': True,
            'complete_physical_particle_catalogue': False,
            'rows_are_distinct_observed_particle_species': False,
            'zero_modes_and_their_interacting_excitations_counted_as_separate_species': False,
            'full_bound_spectrum_outside_central_gap_classified': False,
            'bulk_rest_energy_equals_localized_charged_gap': False,
            'shared_finite_infinite_or_interacting_free_vacuum': False,
            'moving_localized_defect_dispersion': False,
            'conserved_global_radial_momentum': False,
            'absolute_physical_mass_prediction': False,
            'source_selects_texture_and_parameters': False,
            'native_repair_protection_of_poles': False},
    }


def report():
    receipts, pins = parent_receipts()
    science = scientific_catalogue()
    if receipts['radial_interacting_core']['core']['retained_27_site_mass64_source'] != encode(core_bounds()):
        raise ValueError('interacting source bounds differ from their retained receipt')
    if receipts['source_mass_shell']['physical_units'] != science['physical_units']:
        raise ValueError('physical unit contract differs from its retained receipt')
    return {'schema': 'oph.radial_excitation_catalogue.v1', 'all_pass': True,
            'physical_particle_identification_count': 0,
            'catalogue': science, 'source_bindings': pins}


def render_markdown(value):
    data = value['catalogue']
    charged = data['rows'][1]
    probes = charged['spectral_probe_bounds']
    gap = display_interval(charged['energy_model_interval'])
    center = decimal_bound(F(probes['center_probability_lower'])*100)
    rms2 = decimal_bound(probes['RMS_position_squared_about_origin_upper_model_length_squared'], lower=False)
    fourier = decimal_bound(probes['each_Fourier_probe_probability_lower'])
    scale = data['physical_units']['radial_charged_gap_scale_family']
    gev = display_interval(scale['energy_GeV_times_seconds_per_source_time_interval'])
    return f'''# Radial defect excitation catalogue

This catalogue calculates excitations of the supplied radial lattice model.
It associates a protected zero-mode doublet and two actual interacting
charged states with a defect fixed at the origin. Observed particle species
and absolute physical masses are unassigned. The rows describe distinct
source domains; they do not count five observed species.

The finite source has 27 sites, two fermion components per site,
`a=v=1`, radial mass parameter `m=64`, and centered contact `U=1/65536`.
Its interacting reference is the actual `N=27` ground doublet, averaged
with equal weights. The complete Fock complement is retained in the proof;
the `2^54` matrix is not numerically diagonalized.

| Excitation | Energy in model units | Multiplicity | Position and momentum |
| --- | --- | --- | --- |
| Protected one-body zero doublet | 0 | 2 one-body modes | Localized at the supplied defect origin |
| Interacting removal state, `N=26` | {gap} above ground | 1 state; pole residue rank 2 | Central spectral weight at least {center}%; every Fourier probe has positive weight |
| Interacting addition state, `N=28` | {gap} above ground | 1 state; pole residue rank 2 | Same localization bound; signed pole opposite to removal |
| Negative exterior one-body branch | `-sqrt(p²+64²)` in the continuum | 8 rows per wave number | Exterior wave packets, with explicit residual bound |
| Positive exterior one-body branch | `+sqrt(p²+64²)` in the continuum | 8 rows per wave number | Exterior wave packets, with explicit residual bound |

The charged states have equal positive excitation gaps `g`; their retarded
poles are `-g` and `+g`. Every other charged excitation is at least
`8126437/131072 = 61.999794006347656` above the ground.
The two zero modes generate the low-energy interacting core; they are not
additional particle species beside its charged states. The index protects
the one-body zero sector under its stated symmetries, not the value of `g`
under arbitrary changes of the interaction.

## Defect association and localization

The [bounded radial index proof](../../model/bounded_defect_index/PROOF.md)
gives Fredholm index two for the inversion-even to inversion-odd block.
The grading is spatial inversion. Compact self-adjoint perturbations that
anticommute with this grading preserve the index. This is not an
identification of texture winding with electric charge. The stronger
[radial gap theorem](../../model/radial_defect_gap/PROOF.md) proves exactly
two zero modes in its strong-mass domain. The
[interacting core proof](../../model/radial_interacting_core/PROOF.md)
derives the charged states from this same core, including the virtual
complement. Their charges `-1,+1` are conserved source fermion-number
differences; physical electric charge is unassigned.

For either actual low-pole residue `M`, define the spatial spectral-probe
distribution `P_x=Tr(M_xx)/Tr(M)`. This measures where that pole can be
excited by a local fermion probe. It is not a particle trajectory.
The theorem gives `||M-Z/2||<=3t`, where `Z` is the rank-two one-body zero
projector and `t=27/16252874`. Both `M` and `Z` have rank two, so eigenvalue
perturbation gives `1-6t<=Tr(M)<=1+6t`. At the origin,
`Tr(M_00)>=alpha-6t`, where `alpha>=7688/7691`. Consequently

```
P_origin >= (alpha_lower-6t)/(1+6t)
         = {probes['center_probability_lower']}.
sum_x |x|² P_x <= 3a²(1-P_origin) <= {rms2} a².
```

For the 27 orthogonal Fourier probes at `p_j=2*pi*k_j/(3a)`,
`k_j in {{-1,0,1}}`, put `r=alpha_lower-101/1000`.
The quaternionic zero frame satisfies `W_0=sqrt(alpha) I2` and
`W_x*W_x=rho_x I2`. Cauchy-Schwarz bounds the sum of the other 26 frames:
`s_min(F_p*W)>=(sqrt(alpha)-sqrt(26(1-alpha)))/sqrt(27)>=r/sqrt(27)`.
The exact inequality `26(1-alpha_lower)<(101/1000)^2` justifies the
rational choice of `r`. Thus each diagonal low-pole block has eigenvalues
at least `r²/54-3t`, and its normalized two-component probability is

```
P_p >= (r²/27-6t)/(1+6t) >= {fourier}.
```

The same singular-value argument for two different Fourier labels proves
nonzero off-diagonal momentum blocks. The full radial Hamiltonian has no
translation symmetry: translating it moves its unique zero onsite mass
block. Each Fourier probe sees the same localized gap `g`; its probe label
does not define a moving-defect dispersion or a conserved momentum.

## Bulk rest energies and units

The [exterior packet construction](../../model/source_mass_shell/PROOF.md)
uses the same infinite radial one-body operator. Its finite-spacing energy
is `E=+/-sqrt((v/a)^2 sum_j sin(a*p_j)^2+m²)`. Compact normalized packets
at `(0,0,R)` with `L` cells per axis have residual at most
`3*pi*v/[a(L+1)]+4*sqrt(3)*m*L/R`. Their continuum limit has
`E²-v²p²=m²`. All sixteen folded rows are retained, corresponding to four
Dirac copies in the pure-spin bulk family. These multiplicities have no
observed-species assignment. The bulk rest energy 64 is supplied through
the source mass parameter. It is distinct from the finite interacting gap.
The fixed-mass continuum refinement leaves the finite interacting core's
strong-mass proof domain; no interacting relativistic pole follows from it.

With the explicit identification `v=1` with physical light speed,
`length_unit=299792458*T` metres, where `T` is seconds per Hamiltonian time
unit. The [unit conversion](../../model/source_mass_shell/CALIBRATION_PROOF.md)
gives `energy_unit=hbar/T` and the charged energy interval
`{gev}/T GeV`. A positive clock gap identified with cesium in the same
generator would determine `T`. No such input is supplied here. The
catalogue therefore leaves absolute GeV masses and the localized defect's
rest energy null.

## Reproduction

Run from the repository root:

```sh
.venv/bin/python code/particles/derive_radial_excitation_catalogue.py --write
.venv/bin/python code/particles/derive_radial_excitation_catalogue.py --check
.venv/bin/python -m pytest -q code/particles/test_radial_excitation_catalogue.py
```

The [exact JSON catalogue](runs/catalogue/radial_excitation_catalogue.json)
retains rational endpoints, source hashes, all five entries, and their
domain restrictions. Its checker recomputes the catalogue and verifies the
parent source manifests. The owning source packages provide independent
proof controls and cold replays.
'''


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('catalogue differs from the source-bound calculation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument('--write', action='store_true')
    operation.add_argument('--check', action='store_true')
    args = parser.parse_args()
    expected = report()
    markdown = render_markdown(expected)
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(canonical(expected))
        MARKDOWN.write_text(markdown)
        print('wrote exact radial excitation catalogue and Markdown')
    else:
        if canonical(strict_json.load(OUTPUT)) != canonical(expected):
            raise ValueError('retained catalogue differs from its exact source calculation')
        if MARKDOWN.read_text() != markdown:
            raise ValueError('catalogue Markdown differs from its generated view')
        print('radial excitation catalogue check passed')


if __name__ == '__main__':
    main()
