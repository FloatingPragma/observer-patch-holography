"""Field meanings and exact fixed-action algebra for the horizontal atlas.

This module reads the conditional Standard Model owner. Its parameters are
declared action inputs; no numerical laboratory mass enters the calculation.
The finite source-current check and the continuum-action check have separate
domains. Neither identifies the finite 307-state scalar block with a Higgs
resonance in a laboratory.
"""
from functools import lru_cache
import json
from pathlib import Path

import sympy as s

from model.particle_candidate_atlas import conditional_standard_model as sm
from model.charged_lepton_response import relativistic_mass
from model.common_source import native_matter_attachment as matter


SOURCE_INPUTS = (
    'model/particle_candidate_atlas/conditional_standard_model.py',
    'model/particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md',
    'model/common_source/native_matter_attachment.py',
    'model/standard_model_source_map/PROOF.md',
    'model/charged_lepton_response/relativistic_mass.py',
    'model/charged_lepton_response/RELATIVISTIC_MASS.md',
    'model/maxwell_local_selection/PROOF.md',
    'model/maxwell_weak_field/PROOF.md',
    'model/particle_observer_outputs/HIGGS_RECORDS.md',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
)


def _clean(x):
    if isinstance(x, s.MatrixBase):
        return s.ImmutableMatrix(x.applyfunc(s.simplify))
    return s.simplify(x)


def _matrix_strings(x):
    return [[str(a) for a in row] for row in x.tolist()]


@lru_cache(maxsize=1)
def electroweak_algebra():
    """Use the owning Hessians; diagonalize without measured masses."""
    q = sm.quadratic_action()
    g, gp, v, lam = q['symbols']
    mass = q['gauge_mass_squared']
    neutral = mass[2:4, 2:4]
    photon = q['photon_vector_unnormalized']/s.sqrt(g*g+gp*gp)
    z = q['Z_vector_unnormalized']/s.sqrt(g*g+gp*gp)
    if neutral.rank() != 1 or _clean(neutral.det()) != 0:
        raise ArithmeticError('the nonzero-v neutral mass form must have rank one')
    if (_clean(mass*photon) != s.zeros(4, 1)
            or _clean(mass*z-q['m_Z_squared']*z) != s.zeros(4, 1)
            or _clean(photon.H*z) != s.zeros(1, 1)
            or _clean(photon.H*photon)[0] != 1
            or _clean(z.H*z)[0] != 1):
        raise ArithmeticError('neutral field directions failed')
    # W+ = (W1-i W2)/sqrt(2); W- is its conjugate field.
    charged_basis = s.Matrix([[1, 1], [s.I, -s.I]])/s.sqrt(2)
    charged = _clean(charged_basis.H*mass[:2, :2]*charged_basis)
    if charged != q['m_W_squared']*s.eye(2):
        raise ArithmeticError('conjugate charged vector masses differ')
    scalar_expected = s.diag(0, 0, 0, 2*lam*v*v)
    if q['scalar_mass_squared'] != scalar_expected:
        raise ArithmeticError('three orbit directions and radial curvature differ')
    if q['Higgs_electric_charge'] != s.zeros(2, 1):
        raise ArithmeticError('the chosen neutral background is electrically charged')
    t3 = s.diag(1, -1)/2
    raising = s.Matrix([[0, 1], [0, 0]])
    lowering = raising.H
    if (t3*raising-raising*t3 != raising
            or t3*lowering-lowering*t3 != -lowering):
        raise ArithmeticError('weak charged-current weights differ')
    # The quartic energy supplies a cubic and a quartic radial interaction.
    h = s.symbols('h', real=True)
    potential = s.expand(-lam*v*v*(v+h)**2/2+lam*(v+h)**4/4)
    radial = _clean(potential-potential.subs(h, 0))
    if _clean(radial-(lam*v*v*h*h+lam*v*h**3+lam*h**4/4)) != 0:
        raise ArithmeticError('radial self-interaction expansion differs')
    return dict(owner=q, neutral_mass_squared=s.ImmutableMatrix(neutral),
                photon_normalized=_clean(photon), Z_normalized=_clean(z),
                charged_mass_squared=charged, radial_potential=radial)


@lru_cache(maxsize=1)
def source_color_algebra():
    """The eight actual color currents annihilate the declared Higgs doublet."""
    currents = matter.canonical_currents()[:8]
    lifts = tuple(matter.spin_lift(x) for x in currents)
    rank = s.Matrix.hstack(*(s.Matrix(x).reshape(25, 1) for x in lifts)).rank()
    v = s.symbols('v', positive=True)
    h0 = s.Matrix([0, 0, 0, 0, v/s.sqrt(2)])
    images = tuple(_clean(x*h0) for x in lifts)
    if rank != 8 or any(x != s.zeros(5, 1) for x in images):
        raise ArithmeticError('native color currents do not preserve the Higgs direction')
    noncommuting = sum(lifts[i]*lifts[j] != lifts[j]*lifts[i]
                      for i in range(8) for j in range(i+1, 8))
    if noncommuting == 0:
        raise ArithmeticError('the color currents became Abelian')
    return dict(independent_color_currents=rank,
                noncommuting_unordered_pairs=noncommuting,
                color_action_on_Higgs_zero=True,
                Higgs_generated_color_mass_squared=s.zeros(8))


def dirac_algebra():
    """Same Clifford convention as the free source-field owner, energy units."""
    owner = relativistic_mass.clifford_certificate()
    c, mu = s.symbols('c mu', positive=True)
    px, py, pz, t = s.symbols('p_x p_y p_z t', real=True)
    pauli = (s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    beta = s.kronecker_product(pauli[0], s.eye(2))
    alpha = tuple(s.kronecker_product(pauli[2], p) for p in pauli)
    h = c*sum((p*a for p, a in zip((px, py, pz), alpha)), s.zeros(4))+mu*beta
    energy_squared = c*c*(px*px+py*py+pz*pz)+mu*mu
    if _clean(h*h-energy_squared*s.eye(4)) != s.zeros(4):
        raise ArithmeticError('the Dirac rest energy does not close the mass shell')
    left = s.diag(1, 1, 0, 0)
    right = s.eye(4)-left
    if left*beta*left != s.zeros(4) or right*beta*right != s.zeros(4):
        raise ArithmeticError('the mass term is not off diagonal in chirality')
    u0 = s.cos(mu*t)*s.eye(4)-s.I*s.sin(mu*t)*beta
    if _clean(u0.H*u0) != s.eye(4):
        raise ArithmeticError('rest evolution is not unitary')
    chirality_probability = _clean((u0[:, 0].H*right*u0[:, 0])[0])
    if chirality_probability != s.sin(mu*t)**2:
        raise ArithmeticError('rest chirality transition differs')
    return dict(owner=owner, hamiltonian=s.ImmutableMatrix(h),
                energy_squared=energy_squared, rest_energy=mu,
                inertial_mass=mu/c**2,
                rest_chirality_transition=chirality_probability)


def speed_algebra():
    """Continuum and lattice speed statements have explicitly different domains."""
    c, a, p = s.symbols('c a p', positive=True)
    continuum = c*p
    # One ray inside the positive sine branch: 0<a*p<pi/2.
    lattice = c*s.sin(a*p)/a
    return dict(continuum_massless_energy=continuum,
                continuum_group_speed=s.diff(continuum, p),
                lattice_ray_energy=lattice,
                lattice_ray_group_speed=s.diff(lattice, p),
                lattice_ray_domain='0 < a*p < pi/2',
                small_spacing_speed_limit=s.limit(s.diff(lattice, p), a, 0))


def link(label, path, anchor=''):
    return dict(label=label, path=path, anchor=anchor)


ACTION_LINKS = (
    link('Chosen action, exact quadratic spectrum',
         'model/particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md',
         'supplied-action-and-its-exact-quadratic-spectrum'),
    link('Particle paper: conditional field completion',
         'paper/deriving_the_particle_zoo_from_observer_consistency.tex', 'L931'),
    link('Owning symbolic Hessian',
         'model/particle_candidate_atlas/conditional_standard_model.py', 'L63'),
)


def family_rows():
    """Every v1 target has a named action meaning and an inclusive readout."""
    templates = {
        'charged_leptons': dict(
            title='Electron, muon and tau: three charged matter modes',
            source_operation='Encode the lepton doublet and charged singlet-conjugate currents; the chosen family register supplies three coefficients.',
            field_nature='Each mode has left and right chiral components. A Yukawa interaction with the neutral Higgs background couples those components. Diagonalizing that coupling gives three rest frequencies.',
            mass_reason='The mass term mixes chirality even at zero momentum. In the free Lorentz action, E²=c²p²+μ² and μ=v σ(Y_e)/√2; inertial mass is μ/c².',
            readout='Prepare a charged excitation, couple to a complete detector, and infer energy from its phase or inclusive spectrum. Electron thresholds and muon/tau decay resonances use different interacting observables.',
            diagram='Two chiral arrows exchange through a Higgs node, then three normal-mode tracks reach charged detector records.',
            caveat='The three family labels and numerical Yukawa matrix are program inputs unless computed from a declared source controller. Free action masses do not supply interacting measured masses.',
            links=[link('Source weights and family attachment', 'model/standard_model_source_map/PROOF.md'),
                   link('Same-coefficient Dirac dispersion', 'model/charged_lepton_response/RELATIVISTIC_MASS.md', 'exact-mass-shell-of-the-native-quadratic-action')]),
        'quarks': dict(
            title='Quarks: matter modes that also carry color',
            source_operation='The same internal current register supplies a color-triplet weak doublet and color-antitriplet weak singlets; attach their full fermion operators and gauge links.',
            field_nature='Up-type and down-type Yukawa matrices couple the chiral color-carrying fields to the Higgs background. Their singular values give the two sets of three action masses; their left-frame mismatch gives weak mixing.',
            mass_reason='M_u=v Y_u/√2 and M_d=v Y_d/√2. The color field interacts with these modes, so a running quark mass requires a renormalization scheme and scale.',
            readout='Color-singlet hadrons, jets and decay records constrain quark properties. An isolated quark is not a laboratory asymptotic particle.',
            diagram='Three color strands enter a Higgs coupling, then feed a jet/hadron detector rather than a free-quark click.',
            caveat='Matching the internal color and weak weights supplies candidates. The interacting continuum, confinement observables and numerical Yukawa matrices require their declared realization.',
            links=[link('Full-current gauge attachment', 'model/standard_model_source_map/GAUGE_ATTACHMENT.md', 'one-representation-and-a-full-fock-gauss-identity'), *ACTION_LINKS]),
        'neutrinos': dict(
            title='Neutrinos: neutral weak-production and detection channels',
            source_operation='The neutral member of each lepton doublet has Q=0. A chosen Dirac-singlet or Majorana effective interaction supplies a neutrino mass term.',
            field_nature='The fields ν1, ν2 and ν3 propagate with their mass energies. Electron-, muon- and tau-neutrino labels describe weak production/detection combinations; different phases change the flavor probabilities.',
            mass_reason='Dirac: M_ν=v Y_ν/√2. Majorana effective theory: M_ν=v²κ/(2Λ), diagonalized by Takagi factorization. Absolute masses and the extension are separate inputs.',
            readout='A source prepares a flavor channel and a weak detector records a charged lepton after propagation. Keep all flavor, no-detection and failure outcomes.',
            diagram='A weak source splits into three mass-phase tracks that recombine at electron/muon/tau detector branches.',
            caveat='The minimal fifteen-state matter representation has no neutrino mass interaction. Assigning all three masses zero fails the oscillation target; flavor labels are not individual mass eigenstates.',
            links=[link('Neutrino completion and counts', 'model/particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md', 'neutrino-extension-and-the-complete-elementary-channels')]),
        'photon': dict(
            title='Photon: the unbroken electromagnetic wave',
            source_operation='Read loop curvature from seam potentials. Gauge shifts change local coordinates while preserving that curvature; the finite Maxwell action and evolution govern the declared field response.',
            field_nature='In the chosen electroweak phase, one mixture of the neutral weak field and hypercharge field leaves the Higgs background unchanged. This is the electromagnetic field A.',
            mass_reason='A=(g′W3+gB)/√(g²+g′²) is the exact null direction of the Higgs gauge-mass form. It has no Higgs-generated rest-energy term; a moving photon carries E=c|p| in the Lorentz vacuum regime.',
            readout='Prepare an electromagnetic field excitation and couple its energy to an absorbing detector. The detector location is a measurement outcome; a universal sharp photon position operator is not assumed.',
            diagram='Two neutral gauge arrows combine into a background-preserving direction, then a transverse wave reaches an absorber.',
            caveat='The finite carrier has its own oscillator spectrum. Identifying its detector response with laboratory photons requires a spatial continuum, state, clock and electromagnetic coupling map. Finite-lattice group velocity need not equal c.',
            links=[*ACTION_LINKS,
                   link('Finite face curvature and gauge invariance', 'Lean/Screen/LocalFaceMaxwellAction.lean', 'L188'),
                   link('Finite local action classification', 'model/maxwell_local_selection/PROOF.md', 'exhaustive-classification-of-the-local-quadratic-action'),
                   link('Finite wave law', 'Lean/Screen/TemporalMaxwellEvolution.lean', 'L430')]),
        'gluon': dict(
            title='Gluons: the interacting color connection',
            source_operation='The native internal interface supplies eight independent noncommuting color currents. A declared gauge-link field transports those color coordinates between sites.',
            field_nature='The color gauge field responds to color currents and to itself because its generators do not commute. The Higgs doublet is a color singlet, so its background gives no color gauge-boson mass term.',
            mass_reason='All eight Higgs-generated color mass entries are zero. Yang-Mills self-interaction and confinement can generate a physical energy scale without inserting a gluon mass into the gauge action.',
            readout='Color-singlet hadrons and jets carry the observable evidence. A zero gluon Lagrangian mass does not imply an isolated massless colored particle in the spectrum.',
            diagram='Noncommuting color-link transports meet at a three-way interaction node and terminate in color-singlet detector products.',
            caveat='Eight current generators do not determine the strong coupling, confinement scale or complete scattering spectrum.',
            links=[link('Actual native current and Higgs direction', 'model/standard_model_source_map/PROOF.md'),
                   link('Dynamical finite gauge links', 'model/standard_model_source_map/GAUGE_ATTACHMENT.md'), *ACTION_LINKS]),
        'W_plus': dict(
            title='W+ and W−: charged weak-field excitations',
            source_operation='The source weak raising and lowering currents obey [Q,T±]=±T±. A declared weak gauge field couples to those currents.',
            field_nature='The W directions rotate the neutral Higgs background into a different weak direction. The background contributes a quadratic energy cost, giving the charged vector fields a rest frequency.',
            mass_reason='W±=(W1∓iW2)/√2 and m_W²=g²v²/4 in natural units. Three scalar orbit directions supply the longitudinal components of W± and Z.',
            readout='Weak decay and scattering products reconstruct a charged vector resonance, with all competing decay and no-detection channels retained.',
            diagram='A weak-current arrow tilts the Higgs background, then a charged vector resonance branches into decay records.',
            caveat='The exact quadratic mass follows from the declared action and background. Loop corrections and resonance widths require the interacting calculation.',
            links=[*ACTION_LINKS, link('Source charged-current commutators', 'model/standard_model_source_map/PROOF.md')]),
        'Z': dict(
            title='Z: the neutral weak-field direction that feels the background',
            source_operation='Combine the neutral weak current and hypercharge current using the declared couplings.',
            field_nature='The neutral combination orthogonal to the photon changes the Higgs background. It acquires a rest frequency while the photon direction stays massless in the quadratic action.',
            mass_reason='Z=(gW3−g′B)/√(g²+g′²), with m_Z²=(g²+g′²)v²/4. Its charge is zero and it has three massive-vector polarizations.',
            readout='Neutral-current scattering and the complete decay spectrum reconstruct the Z resonance.',
            diagram='Two neutral arrows split into a flat photon direction and a curved Z direction, with the Z branch ending in decay records.',
            caveat='The photon/Z mixing angle depends on the chosen gauge couplings; the action Hessian does not determine those couplings.',
            links=list(ACTION_LINKS)),
        'Higgs': dict(
            title='Higgs: the radial fluctuation of the symmetry-breaking field',
            source_operation='A declared complex scalar doublet couples to the source matter currents through gauge-invariant pair creation and annihilation. In the continuum completion, choose its potential and phase.',
            field_nature='The neutral background is a reference configuration. A Higgs excitation changes its radial size. The potential curves along that direction, producing a restoring frequency; its cubic and quartic terms also make it interact.',
            mass_reason='V=−λv² H†H+λ(H†H)² gives m_h²=2λv². The same background enters |D H|² and Yukawa terms, producing W/Z and fermion masses. It does not act as friction.',
            readout='A scalar excitation is inferred from its complete decay products and resonance energy distribution. Its spin-zero transformation follows from the declared scalar field action.',
            diagram='A radial pulse moves in the curved Higgs potential, with side branches showing W/Z and left-right matter couplings.',
            caveat='The 307-state finite quantum-Higgs block has its own computed transition gap and complete detector. Equating that gap to the continuum radial Higgs mass requires an explicit field/state/scaling map.',
            links=[*ACTION_LINKS,
                   link('Finite interacting scalar records', 'model/particle_observer_outputs/HIGGS_RECORDS.md', 'one-actual-source-carrier-and-hamiltonian'),
                   link('Complete finite apparatus and common program clock', 'model/particle_observer_outputs/READOUT_CONTRACT.md')]),
    }
    result = []
    for base in sm.benchmark_rows():
        name = base['target']
        if name in ('electron', 'muon', 'tau'):
            family = 'charged_leptons'
        elif name in ('up', 'down', 'charm', 'strange', 'top', 'bottom'):
            family = 'quarks'
        elif name.endswith('_neutrino'):
            family = 'neutrinos'
        else:
            family = 'W_plus' if name == 'W_minus' else name
        row = dict(base, **templates[family])
        row['family'] = family
        row['status'] = 'exact quadratic consequences of a declared action; physical identification requires its observer map'
        row['antiparticle_meaning'] = (
            'The neutral real field is self-conjugate.' if base['antiparticle'] == 'self'
            else 'The color adjoint is a real representation; individual color labels depend on the gauge basis.' if name == 'gluon'
            else 'Dirac neutrinos have distinct conjugate states. Majorana fields are self-conjugate, with distinct weak detection helicity channels.' if family == 'neutrinos'
            else 'The conjugate excitation carries opposite additive gauge charges and the same mass in the declared CPT-compatible relativistic action.')
        result.append(row)
    return result


def certificate():
    ew, color, d, speed = (electroweak_algebra(), source_color_algebra(),
                            dirac_algebra(), speed_algebra())
    q = ew['owner']
    return {
        'schema': 'oph.horizontal_atlas.field_nature.v1',
        'scope': 'A declared field action determines its quadratic modes. The native finite current and complete finite scalar records are separate computed source results.',
        'inputs': ['g and g_prime', 'Higgs norm v', 'positive lambda_H',
                   'Y_e, Y_u and Y_d', 'neutrino extension',
                   'canonical Lorentz kinetic terms and field domain'],
        'units': 'Quadratic Standard Model formulas use hbar=c=1. The Dirac panel writes c explicitly and mu is rest energy, so inertial mass=mu/c^2.',
        'empirical_mass_inputs': [],
        'quadratic': {
            'gauge_basis': ['W1', 'W2', 'W3', 'B'],
            'gauge_mass_squared': _matrix_strings(q['gauge_mass_squared']),
            'neutral_mass_squared': _matrix_strings(ew['neutral_mass_squared']),
            'neutral_rank_for_positive_parameters': 1,
            'photon_normalized': _matrix_strings(ew['photon_normalized']),
            'Z_normalized': _matrix_strings(ew['Z_normalized']),
            'charged_mass_squared': _matrix_strings(ew['charged_mass_squared']),
            'scalar_basis': ['x1', 'x2', 'x3', 'h_radial'],
            'scalar_mass_squared': _matrix_strings(q['scalar_mass_squared']),
            'radial_potential_above_minimum': str(ew['radial_potential']),
            **{k: str(q[k]) for k in ('m_W_squared', 'm_Z_squared', 'm_H_squared', 'electric_coupling')},
        },
        'native_color': {k: (_matrix_strings(value) if isinstance(value, s.MatrixBase) else value)
                         for k, value in color.items()},
        'dirac': {k: (_matrix_strings(value) if isinstance(value, s.MatrixBase)
                      else str(value) if isinstance(value, s.Basic) else value)
                  for k, value in d.items()},
        'speed': {k: str(value) for k, value in speed.items()},
        'panels': [
            {'id': 'what_a_particle_is', 'title': 'A particle is a mode seen through its records',
             'text': 'A source register carries currents and interactions. With a declared spatial field action, its patterns propagate. A particle label names a stable or resonant excitation distinguished by its energy, charge, transformation and detector response.',
             'boundary': 'A protected defect can seed or carry such an excitation. Its topology alone does not fix a particle name or mass.'},
            {'id': 'why_mass', 'title': 'Mass is a rest frequency, not drag',
             'text': 'A massive mode evolves even when its momentum is zero. In a Lorentz field its rest energy μ gives E²=c²p²+μ². The Higgs background supplies some quadratic terms that create this frequency; interactions change the measurable spectral response.',
             'boundary': 'A massless wave has no rest gap and can carry nonzero energy. Generic repair effort becomes a particle mass only through a proved energy and spectral readout map.'},
            {'id': 'why_c', 'title': 'The photon follows the light cone of its field action',
             'text': 'For the massless electromagnetic mode in a Lorentz-invariant vacuum, E=c|p| and a wave packet travels at c. The parameter c belongs to the common spatial and temporal units of that action.',
             'boundary': 'Lattice wavelengths and propagation through matter can have different group velocities. A finite carrier wave equation alone does not establish the laboratory speed of light.'},
            {'id': 'what_is_chosen', 'title': 'Choose the program, calculate the response',
             'text': 'The architecture may declare its fields, couplings, state preparation and apparatus. Its record probabilities and spectrum must follow from that program. Target masses used to choose a coupling are calibration inputs.',
             'boundary': 'The quadratic field census is complete for the chosen Standard Model action and neutrino extension. A complete source-to-laboratory particle identification requires the same interacting observer output map.'},
        ],
        'families': family_rows(),
        'finite_Higgs_vs_continuum': {
            'finite': 'H in a conserved 307-state sector has two bright levels and 305 retained dark levels. Delta=sqrt(d^2+4g_pair^2), g_pair^2=4 Tr(Y^dagger Y); the gap is inferred from complete sequential records.',
            'continuum': 'A chosen color- and charge-preserving Higgs phase has radial Hessian 2 lambda_H v^2. This is the tree mass-squared coefficient in its canonical continuum action.',
            'distinct_symbols': 'g_pair is the finite pair-vertex strength; electroweak g is a gauge coupling. These are different operators and domains.',
            'required_join': 'An explicit embedding, state and scaling map must transport the finite complete correlators and instruments into the continuum radial field before the gaps can be equated.',
            'identified_with_measured_Higgs': False,
        },
        'bosonic_degrees': sm.field_census()['boson_on_shell_degrees_broken'],
        'physical_particle_assignments': 0,
        'source_inputs': list(SOURCE_INPUTS),
    }


if __name__ == '__main__':
    path = Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(certificate(), ensure_ascii=False, indent=2)+'\n')
    print(path)
