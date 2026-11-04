"""Independent numerical CAR and native radial checks."""
import numpy as np
from scipy.linalg import expm

from model.defect_current_bridge import construction
from model.mcpr_native_memory_bridge.construction import design,direct_memory
from model.common_source import native_family_record_memory as memory
from . import algebra


def numpy_car(modes):
    lower = np.array([[0.,1.],[0.,0.]])
    parity = np.diag([1.,-1.])
    result = []
    for index in range(modes):
        value = np.ones((1,1))
        for j in range(modes):
            value = np.kron(value, parity if j < index else lower if j == index else np.eye(2))
        result.append(value)
    return result


def numerical_block(mass,coupling):
    if not np.isfinite([mass,coupling]).all() or mass <= 0 or coupling < 0:
        raise ValueError('positive mass and nonnegative finite coupling required')
    c = numpy_car(4)
    n = sum((a.T@a for a in c), np.zeros((16,16)))
    pairing = sum((c[i].T@c[j].T+c[j]@c[i] for i,j in ((0,1),(2,3))), np.zeros((16,16)))
    v = (n-2*np.eye(16))@(n-2*np.eye(16))
    return mass*pairing+coupling*v, v, c


def spectrum_control(mass,coupling):
    h, v, c = numerical_block(mass,coupling)
    eig, vectors = np.linalg.eigh(h)
    r = np.hypot(mass,coupling)
    expected = np.sort([0.]*5+[4*coupling,2*coupling-2*r,2*coupling+2*r]
                       +[coupling-mass]*4+[coupling+mass]*4)
    ground = vectors[:,0]
    n0,n2 = c[0].T@c[0],c[2].T@c[2]
    connected = float(ground@(n0@n2)@ground-(ground@n0@ground)*(ground@n2@ground))
    if abs(connected+coupling/(4*r)) > 1e-9:
        raise ArithmeticError('connected four-point density identity failed')
    weights = np.abs(vectors.T@c[0]@ground)**2+np.abs(vectors.T@c[0].T@ground)**2
    gaps = [2*r-coupling-mass,2*r-coupling+mass]
    residues = [sum(weights[abs(eig-eig[0]-gap) < 1e-10]) for gap in gaps]
    expected_weights = [(1+mass/r)/2,(1-mass/r)/2]
    error = max(np.max(abs(eig-expected)),np.max(abs(np.array(residues)-expected_weights)))
    if error > 1e-9:
        raise ArithmeticError('independent CAR spectrum or residue control failed')
    return {'mass': format(mass,'.15g'), 'coupling': format(coupling,'.15g'),
            'single_fermion_gaps': [format(x,'.15g') for x in gaps],
            'addition_plus_removal_weights': [format(x,'.15g') for x in expected_weights],
            'connected_cross_pair_density': format(-coupling/(4*r),'.15g'),
            'full_Fock_rows': 16, 'maximum_error_below': '1e-9',
            'number_and_pairing_do_not_commute': bool(np.linalg.norm(v@h-h@v)>0)}


def onsite_control():
    _, w, _, row = construction.radial()
    w = w.toarray()
    rho = row['density']
    # Two distinct internal flavors make the antisymmetric CAR sector
    # isometric to a tensor product of their spatial-spin coordinates.
    inclusion = np.kron(w,w)
    points = np.repeat(np.arange(len(rho)),2)
    onsite = (points[:,None] == points[None,:]).reshape(-1).astype(float)
    acted = onsite[:,None]*inclusion
    compressed = inclusion.conj().T@acted
    ipr = float(rho@rho)
    residual = acted-inclusion@compressed
    gram = residual.conj().T@residual
    if np.linalg.norm(compressed-ipr*np.eye(4)) > 1e-10 or np.linalg.norm(gram-ipr*(1-ipr)*np.eye(4)) > 1e-10:
        raise ArithmeticError('onsite two-particle leakage identity failed')
    if not 0 < ipr < 1:
        raise ArithmeticError('source must occupy more than one site')
    return {'all_two_flavor_spatial_rows': len(onsite), 'core_two_flavor_rows': 4,
            'ipr': format(ipr,'.15g'), 'unit_onsite_leakage': format(np.sqrt(ipr*(1-ipr)),'.15g'),
            'exact_compressed_operator': 'c*I, c=sum_x rho_x^2',
            'exact_leakage_Gram': 'c*(1-c)*I',
            'weighted_onsite_leakage_Gram': '(sum_x g_x^2*rho_x^2-(sum_x g_x*rho_x^2)^2)*I',
            'nonzero_onsite_density_interaction_preserves_core': False}


def source_controls():
    row = design('forward',0)
    return {'P_branch': 'forward', 'inverse_memory_branch': 0,
            'declared_coupling_over_v': '1/128',
            'lepton_blocks': [spectrum_control(m,1/128) for m in row['masses_over_v']],
            'blocks_are_invariant_slices_not_full_core_vacuum': True,
            'full_core_modes_retained_in_analytic_construction': 90}


def memory_control():
    """All source outcomes on a16-state nonlinear auxiliary probe.

    The probe reads one family matrix entry. It is a control of the general
    averaging theorem, not a reducing four-mode sector of every family pulse.
    """
    source = design('forward',0)
    law = direct_memory(source)
    transition,enabled = law['transition'],law['enabled_transition']
    projections = [float(np.asarray(p,float)[0,0]) for p in memory._projectors(False)]
    gain = source['gain_over_v']
    fixed,v,_ = numerical_block(.25,.5)
    pairing,_,_ = numerical_block(1.,0.)
    mean = fixed+gain*law['coefficient_matrix'][0,0]*pairing
    initial = np.zeros((16,16),complex)
    initial[0,0] = 1
    target_unitary = expm(-1j*mean)
    target = target_unitary@initial@target_unitary.conj().T
    errors = []
    for steps in (32,64,128):
        disabled_unitary = expm(-1j*fixed/steps)
        enabled_unitaries = np.array([expm(-1j*(fixed+gain*x*pairing)/steps) for x in projections])
        blocks = law['stationary'][:,None,None]*initial
        for _ in range(steps):
            disabled_blocks = np.einsum('jk,kab->jab',transition-enabled,blocks)
            enabled_blocks = np.einsum('jk,kab->jab',enabled,blocks)
            blocks = (disabled_unitary@disabled_blocks@disabled_unitary.conj().T
                      +enabled_unitaries@enabled_blocks@enabled_unitaries.conj().transpose(0,2,1))
        errors.append(float(np.sum(np.abs(np.linalg.eigvalsh(blocks.sum(axis=0)-target)))))
    if not errors[0]>errors[1]>errors[2]>0 or errors[0]/errors[-1] < 3.5:
        raise ArithmeticError('complete-memory nonlinear control failed its weak-step trend')
    return {'probe_Fock_rows':16, 'old_memory_inputs':12,
            'retained_fine_records_per_input':law['retained_records_per_input'],
            'steps':[32,64,128], 'trace_norm_errors':[format(x,'.10g') for x in errors],
            'stationary_product_initial_memory':True,
            'probe_is_full_physical_core':False,
            'noncommuting_quartic_and_pairing_pulses':True}
