"""Exact quadratic consequences of a supplied renormalized SM benchmark.

No empirical masses, source vacuum, continuum regulator or matching correction
are inferred by this module. Hypercharge is Q=T3+Y, with coupling g_prime.
"""
import sympy as s

SOURCE_INPUTS = (
    'model/standard_model_source_map/FERMION_ATTACHMENT.md',
    'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
    'model/charged_lepton_response/RELATIVISTIC_MASS.md',
    'model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md',
    'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',
    'model/particle_candidate_atlas/DEFECT_FIELD_JOIN.md',
    'model/compact_defect_pair/spectral.py',
    'model/compact_defect_pair/topology.py',
)
PRIMARY_REVIEWS = {
    'electroweak': 'https://pdg.lbl.gov/2026/reviews/rpp2026-rev-standard-model.pdf',
    'qcd': 'https://pdg.lbl.gov/2026/reviews/rpp2026-rev-qcd.pdf',
    'neutrinos': 'https://pdg.lbl.gov/2026/reviews/rpp2026-rev-neutrino-mixing.pdf',
    'higgs': 'https://pdg.lbl.gov/2026/reviews/rpp2026-rev-higgs-boson.pdf',
}


def representations():
    """One family of left-handed Weyl representations, with no spin factor."""
    return (
        dict(name='Q', color='3', color_dimension=3, weak_dimension=2,
             hypercharge='1/6', electric_weights=['2/3', '-1/3']),
        dict(name='u_c', color='bar3', color_dimension=3, weak_dimension=1,
             hypercharge='-2/3', electric_weights=['-2/3']),
        dict(name='d_c', color='bar3', color_dimension=3, weak_dimension=1,
             hypercharge='1/3', electric_weights=['1/3']),
        dict(name='L', color='1', color_dimension=1, weak_dimension=2,
             hypercharge='-1/2', electric_weights=['0', '-1']),
        dict(name='e_c', color='1', color_dimension=1, weak_dimension=1,
             hypercharge='1', electric_weights=['1']),
    )


def anomaly_coefficients():
    """Local anomaly sums and the weak-doublet parity of the benchmark."""
    rows = representations()
    y = {r['name']: s.Rational(r['hypercharge']) for r in rows}
    linear = sum(r['color_dimension'] * r['weak_dimension'] *
                 y[r['name']] for r in rows)
    cubic = sum(r['color_dimension'] * r['weak_dimension'] *
                y[r['name']]**3 for r in rows)
    return {
        'mixed_gravity_U1': linear, 'U1_cubed': cubic,
        'SU3_squared_U1': (2*y['Q']+y['u_c']+y['d_c'])/2,
        'SU2_squared_U1': (3*y['Q']+y['L'])/2,
        'SU3_cubed': s.Integer(2-1-1),
        'weak_doublets_per_family': 4, 'weak_doublets_three_families': 12,
        'up_Yukawa_hypercharge': y['Q']+s.Rational(1, 2)+y['u_c'],
        'down_Yukawa_hypercharge': y['Q']-s.Rational(1, 2)+y['d_c'],
        'charged_Yukawa_hypercharge': y['L']-s.Rational(1, 2)+y['e_c'],
        'Dirac_neutrino_Yukawa_hypercharge': y['L']+s.Rational(1, 2),
    }


def quadratic_action():
    """Differentiate |D H0|^2 and the scalar potential, without mass inputs."""
    g, gp, v, lam = s.symbols('g g_prime v lambda_H', positive=True)
    w1, w2, w3, b = fields = s.symbols('W1 W2 W3 B', real=True)
    sigma = (s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    h0 = s.Matrix([0, v/s.sqrt(2)])
    covariant = (g*sum((w*t/2 for w, t in zip(fields[:3], sigma)),
                      s.zeros(2)) + gp*b*s.eye(2)/2)*h0
    gauge_density = s.expand((covariant.H*covariant)[0])
    gauge_mass_squared = s.hessian(gauge_density, fields)
    x1, x2, x3, radial = scalars = s.symbols('x1 x2 x3 h', real=True)
    doublet = s.Matrix([x1+s.I*x2, v+radial+s.I*x3])/s.sqrt(2)
    modulus = s.expand((doublet.H*doublet)[0])
    potential = -lam*v*v*modulus+lam*modulus**2
    origin = dict.fromkeys(scalars, 0)
    scalar_mass_squared = s.hessian(potential, scalars).subs(origin)
    photon = s.Matrix([0, 0, gp, g])
    z = s.Matrix([0, 0, g, -gp])
    if gauge_mass_squared*photon != s.zeros(4, 1):
        raise ArithmeticError('unbroken electromagnetic null vector failed')
    if s.simplify(gauge_mass_squared*z-v*v*(g*g+gp*gp)*z/4) != s.zeros(4, 1):
        raise ArithmeticError('neutral massive vector failed')
    return {
        'symbols': (g, gp, v, lam), 'gauge_density': gauge_density,
        'gauge_mass_squared': s.ImmutableMatrix(gauge_mass_squared),
        'scalar_mass_squared': s.ImmutableMatrix(scalar_mass_squared),
        'photon_vector_unnormalized': s.ImmutableMatrix(photon),
        'Z_vector_unnormalized': s.ImmutableMatrix(z),
        'm_W_squared': g*g*v*v/4,
        'm_Z_squared': (g*g+gp*gp)*v*v/4,
        'm_H_squared': 2*lam*v*v,
        'electric_coupling': g*gp/s.sqrt(g*g+gp*gp),
        'Higgs_electric_charge': s.ImmutableMatrix((sigma[2]/2+s.eye(2)/2)*h0),
    }


def field_census(neutrinos='majorana_eft'):
    """On-shell quadratic field counts, including antiparticles where distinct."""
    if neutrinos not in ('minimal_massless', 'majorana_eft', 'dirac'):
        raise ValueError('an explicit neutrino completion is required')
    extra = 3 if neutrinos == 'dirac' else 0
    return {
        'neutrino_completion': neutrinos, 'family_copies': 3,
        'left_Weyl_components': 45+extra,
        'two_component_CAR_coordinates': 2*(45+extra),
        'charged_Dirac_fields_counting_color': 3+6*3,
        'charged_fermion_on_shell_degrees': 4*(3+6*3),
        'neutrino_on_shell_degrees': 12 if extra else 6,
        'boson_on_shell_degrees_unbroken': 2*(8+3+1)+4,
        'boson_on_shell_degrees_broken': 2*8+2+2*3+3+1,
        'neutrino_mass_values': None,
        'minimal_massless_fits_oscillations': False,
        'counts_are_free_field_counts': True,
    }


def benchmark_rows():
    """The eighteen atlas target labels; numerical reference values are separate."""
    rows = []
    for i, name in enumerate(('electron', 'muon', 'tau'), 1):
        rows.append(dict(target=name, spin='1/2', charge='-1', color='1',
                         antiparticle='positron' if i == 1 else 'anti_'+name,
                         tree_mass=f'v*sigma_{i}(Y_e)/sqrt(2)',
                         observable='charged infrared threshold' if i == 1
                         else 'inclusive resonance mass and width'))
    for names, sector, q in ((('up', 'charm', 'top'), 'u', '2/3'),
                             (('down', 'strange', 'bottom'), 'd', '-1/3')):
        for i, name in enumerate(names, 1):
            rows.append(dict(target=name, spin='1/2', charge=q, color='3',
                             antiparticle='anti_'+name,
                             tree_mass=f'v*sigma_{i}(Y_{sector})/sqrt(2)',
                             observable='renormalized mass in a stated scheme and scale; no isolated colored asymptotic particle'))
    for name in ('electron_neutrino', 'muon_neutrino', 'tau_neutrino'):
        rows.append(dict(target=name, spin='1/2', charge='0', color='1',
                         antiparticle='Dirac conjugate or Majorana helicity channel, according to the chosen extension',
                         tree_mass=None,
                         observable='flavor production and detection; superposition of nu_1, nu_2, nu_3'))
    rows += [
        dict(target='photon', spin='1', charge='0', color='1', antiparticle='self',
             tree_mass='0', observable='massless unbroken electromagnetic gauge field; two helicities'),
        dict(target='gluon', spin='1', charge='0', color='8_adjoint', antiparticle='real adjoint multiplet',
             tree_mass='0', observable='eight massless Lagrangian gauge fields; no isolated colored asymptotic particle'),
        dict(target='W_plus', spin='1', charge='1', color='1', antiparticle='W_minus',
             tree_mass='g*v/2', observable='massive vector resonance; three quadratic polarizations'),
        dict(target='W_minus', spin='1', charge='-1', color='1', antiparticle='W_plus',
             tree_mass='g*v/2', observable='massive vector resonance; three quadratic polarizations'),
        dict(target='Z', spin='1', charge='0', color='1', antiparticle='self',
             tree_mass='v*sqrt(g**2+g_prime**2)/2', observable='neutral massive vector resonance'),
        dict(target='Higgs', spin='0', charge='0', color='1', antiparticle='self',
             tree_mass='sqrt(2*lambda_H)*v', observable='neutral scalar resonance'),
    ]
    for row in rows:
        row.update(conditional_benchmark=True, source_predicted_physical_mass=None,
                   numerical_mass_input=None)
    return rows


def certificate():
    """JSON-safe conditional algebra; no full source producer is replayed."""
    a = anomaly_coefficients()
    if any(v != 0 for k, v in a.items() if 'doublets' not in k):
        raise ArithmeticError('benchmark anomaly or Yukawa selection check failed')
    q = quadratic_action()
    matrices = ('gauge_mass_squared', 'scalar_mass_squared',
                'photon_vector_unnormalized', 'Z_vector_unnormalized', 'Higgs_electric_charge')
    quadratic = {k: [[str(v) for v in row] for row in q[k].tolist()] for k in matrices}
    quadratic.update({k: str(q[k]) for k in
                      ('m_W_squared', 'm_Z_squared', 'm_H_squared', 'electric_coupling')})
    return {
        'schema': 'conditional_standard_model_v1',
        'status': 'conditional_renormalized_SM_benchmark',
        'hypercharge_convention': 'Q=T3+Y; g_prime is not GUT-rescaled',
        'representations_one_family': list(representations()),
        'anomalies_and_Yukawa_charges': {k: str(v) for k, v in a.items()},
        'quadratic_action': quadratic,
        'field_censuses': [field_census(x) for x in
                          ('minimal_massless', 'majorana_eft', 'dirac')],
        'target_rows': benchmark_rows(),
        'neutrino_mass_fields': [dict(field=f'nu_{i}', mass=f'm_nu_{i} >= 0',
                                      numerical_mass=None) for i in (1, 2, 3)],
        'neutrino_flavor_relation': 'nu_alpha_L=sum_i U_PMNS[alpha,i]*nu_i_L',
        'Dirac_neutrino_mass_convention': 'M_nu=v*Y_nu/sqrt(2)',
        'Weinberg_mass_convention': 'M_nu=v**2*kappa/(2*Lambda), kappa symmetric',
        'CPT_mass_degeneracy_requires': 'local Lorentz-covariant quantum theory with positive-energy unitary vacuum representation and the usual spin-statistics/CPT hypotheses; not charge conjugacy alone',
        'empirical_reference_values_recomputed': False,
        'all_reference_masses_inverted_as_tree_parameters': False,
        'physical_OPH_assignments': 0,
        'source_selected_renormalized_parameters': False,
        'source_selected_vacuum_and_clock': False,
        'naive_eight_valleys_removed': False,
        'interacting_continuum_existence_proved': False,
        'finite_scalar_identified_with_continuum_Higgs': False,
        'primary_reviews': dict(PRIMARY_REVIEWS),
    }
