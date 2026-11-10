"""Complete local three-family quantum memory on the native CAR and scalar operators."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import mpmath as mp
import numpy as np
import sympy as s
from model.coherent_defect_wave.chiral_joint.model import active_terms, sequence
from model.coherent_defect_wave.chiral_higgs.model import quantum_lepton_block,tensor_checks
from model.coherent_defect_wave.physical_attachment.chiral import representation

_, TENSORS, FIELDS, _ = active_terms()
Y = s.symbols("y0 y1 y2", positive=True)
OMEGA = s.Symbol("omega", positive=True)


def mode(family, spin, internal):
    return 6 * family + 3 * spin + internal


def column(key, ys=Y, omega=OMEGA):
    mask, scalar = key
    out = defaultdict(lambda: s.Integer(0))
    if scalar:
        out[key] += omega
    for f, coefficient in enumerate(ys):
        if coefficient == 0:
            continue
        for k in range(2):
            for i, j in product(range(3), repeat=2):
                tij = TENSORS[k][i, j]
                if not tij:
                    continue
                ann = ((mode(f, 0, i), False), (mode(f, 1, j), False))
                cre = ((mode(f, 1, j), True), (mode(f, 0, i), True))
                for scalar_op, word in ((FIELDS[k].H, ann), (FIELDS[k], cre)):
                    step = sequence(mask, word)
                    if step is None:
                        continue
                    for dest in range(5):
                        if scalar_op[dest, scalar]:
                            out[step[0], dest] += coefficient * tij * step[1] * scalar_op[dest, scalar]
    return {key: s.expand(value) for key, value in out.items() if value != 0}


def component(seed):
    seen, front, edges = {seed}, [seed], {}
    while front:
        key = front.pop()
        edges[key] = column(key)
        for target in edges[key]:
            if target not in seen:
                seen.add(target)
                front.append(target)
    keys = [seed] + sorted(seen - {seed})
    index = {key: i for i, key in enumerate(keys)}
    matrix = s.zeros(len(keys))
    for source, col in edges.items():
        for target, value in col.items():
            matrix[index[target], index[source]] = value
    assert matrix == matrix.H
    assert matrix[0, 0] == 0
    assert matrix[1:, 1:] == OMEGA * s.eye(len(keys)-1)
    weight = s.expand((matrix[0, 1:] * matrix[1:, 0])[0])
    return keys, matrix, weight


def creation(source, target, which_mode):
    target_index = {key: i for i, key in enumerate(target)}
    matrix = s.zeros(len(target), len(source))
    for j, (mask, scalar) in enumerate(source):
        step = sequence(mask, ((which_mode, True),))
        if step is not None:
            key = (step[0], scalar)
            assert key in target_index, ("operator image escapes retained support", key)
            matrix[target_index[key], j] = step[1]
    return matrix


def higgs_composite(source, target, family):
    """Independent CAR formula [H_Y,c_0,e^dag]=-y sum C_ej phi^dag c_1,j."""
    index = {key:i for i,key in enumerate(target)}
    matrix = s.zeros(len(target),len(source))
    for column_index,(mask,scalar) in enumerate(source):
        for k in range(2):
            for j in range(3):
                coefficient = -Y[family]*TENSORS[k][0,j]
                if not coefficient:
                    continue
                step = sequence(mask,((mode(family,1,j),False),))
                if step is None:
                    continue
                for dest in range(5):
                    value = FIELDS[k].H[dest,scalar]
                    if value:
                        assert (step[0],dest) in index
                        matrix[index[step[0],dest],column_index] += coefficient*step[1]*value
    return matrix


def check_native_one_family():
    native = quantum_lepton_block()[0]
    assert native.shape == (320, 320)
    for mask in range(64):
        for scalar in range(5):
            expected = {5*tmask + tscalar: value
                        for (tmask, tscalar), value in column((mask, scalar), (1, 0, 0), 1).items()}
            actual = {row: native[row, 5*mask+scalar] for row in range(320)
                      if native[row, 5*mask+scalar] != 0}
            assert actual == expected
    ward=tensor_checks()
    return {"complete_native_columns_compared": 320, "mismatches": 0,
            "actual_native_pair_Ward_checks":ward['full_twelve_current_pair_Ward_checks'],
            "actual_native_scalar_commutator_checks":ward['quantum_scalar_commutator_checks']}


def numeric_transfer(h0, hf, T, ys, omega):
    substitutions = dict(zip(Y, ys)) | {OMEGA: omega}
    n0 = np.asarray(h0.subs(substitutions).evalf(18).tolist(), dtype=float)
    nf = np.asarray(hf.subs(substitutions).evalf(18).tolist(), dtype=float)
    tn = np.asarray(T.tolist(), dtype=float)
    e0, v0 = np.linalg.eigh(n0)
    ef, vf = np.linalg.eigh(nf)
    g = v0[:, 0]
    amplitudes = vf.T @ tn @ g
    jn = nf @ tn - tn @ n0
    composite = vf.T @ jn @ g
    gaps = ef - e0[0]
    defect = np.max(np.abs(composite - gaps * amplitudes))
    assert defect < 2e-13
    assert np.max(np.abs(amplitudes[1:-1])) < 1e-12
    assert np.max(np.abs(composite[1:-1])) < 1e-12
    assert abs(sum(amplitudes**2) - np.linalg.norm(tn@g)**2) < 2e-13
    no_operation=1-sum(amplitudes**2)
    assert no_operation>=-2e-13
    S0 = float(sum((h0[0,j]**2).subs(substitutions) for j in range(1,h0.rows)))
    Sf = float(sum((hf[0,j]**2).subs(substitutions) for j in range(1,hf.rows)))
    w = float(omega)
    analytic = []
    N0sq = 1/(1+S0/(e0[0]-w)**2)
    for j in (0,len(ef)-1):
        Nfsq = 1/(1+Sf/(ef[j]-w)**2)
        residue = N0sq*Nfsq*(1+Sf/((e0[0]-w)*(ef[j]-w)))**2
        assert abs(residue-amplitudes[j]**2) < 2e-13
        analytic.append(residue)
    return {"source_reference_energy": float(e0[0]),
            "complete_target_dimension": len(ef),
            "dark_energy": float(ef[1]), "dark_multiplicity": len(ef)-2,
            "dark_total_elementary_residue": float(sum(amplitudes[1:-1]**2)),
            "lines": [{"energy_transfer": float(gaps[j]),
                       "elementary_residue": float(amplitudes[j]**2),
                       "composite_residue": float(composite[j]**2)} for j in (0, len(ef)-1)],
            "residue_sum": float(sum(amplitudes**2)),
            "elementary_preparation_no_operation_probability":float(no_operation),
            "complete_elementary_record_probability":float(sum(amplitudes**2)+no_operation),
            "closed_form_residues": analytic,
            "same_H_amplitude_identity_max_error": float(defect)}


def masses(ys, omega):
    S0 = 4 * sum(y*y for y in ys)
    D0 = mp.sqrt(omega*omega + 4*S0)
    # Rationalization preserves tiny charged lines at small y/omega.
    return [4*y*y/(D0 + mp.sqrt(D0*D0 - 8*y*y)) for y in ys]


def preparation_check(h0):
    """Exact algebraic stationary state and an explicit Householder preparation."""
    parameters=dict(zip(Y,(s.Rational(1,7),s.Rational(2,7),s.Rational(3,7))))|{OMEGA:s.Integer(1)}
    h=h0.subs(parameters)
    weight=(h[0,1:]*h[1:,0])[0]
    energy=(1-s.sqrt(1+4*weight))/2
    raw=s.eye(h.rows)[:,0]+h[:,0]/(energy-1)
    norm=s.simplify((raw.T*raw)[0])
    g=raw/s.sqrt(norm)
    assert s.simplify((g.T*g)[0]-1)==0
    assert (h*g-energy*g).applyfunc(s.simplify)==s.zeros(h.rows,1)
    blank=s.eye(h.rows)[:,0]
    vector=blank-g
    denominator=s.simplify((vector.T*vector)[0])
    # Full Householder I-2|v><v|/<v,v| is defined on the source carrier;
    # its rank-one change is confined to this graph-closed support.
    prepared=blank-2*vector*(vector.T*blank)[0]/denominator
    assert (prepared-g).applyfunc(s.simplify)==s.zeros(h.rows,1)
    return {'energy':str(energy),'norm_squared':'1','stationary_under_actual_H':True,
            'preparation':'Householder on blank-minus-ground, identity on orthogonal complement',
            'exact_Householder_maps_blank_to_ground':True,'preparation_success_probability':'1',
            'idle_reference_dimension':1,'full_source_global_ground_claim':False}


def compensated_current_check():
    currents=representation()[2]
    source_unit=s.eye(15)[:,9]
    for i,current in enumerate(currents):
        expected=source_unit if i==11 else s.zeros(15,1)
        assert current[:,9]==expected
    apparatus_charge=s.diag(0,-1)
    transition=s.Matrix([[0,0],[1,0]])
    assert apparatus_charge*transition-transition*apparatus_charge==-transition
    return {'actual_native_current_columns_checked':12,'source_creator_hypercharge':'1',
            'source_creator_weak_and_colour_singlet':True,'apparatus_dimension':2,
            'apparatus_success_charge':'-1','apparatus_Hamiltonian':'0',
            'all_twelve_compensated_Ward_identities':True,
            'elementary_instrument_effects':'T^dag T and I-T^dag T, both retained'}


@lru_cache(maxsize=1)
def calculate():
    mp.mp.dps = 70
    native = check_native_one_family()
    keys0, h0, S0 = component((0, 0))
    assert S0 == 4 * sum(y*y for y in Y)
    assert len(keys0) == 13
    preparation=preparation_check(h0)
    outputs = []
    for f in range(3):
        keysf, hf, Sf = component((1 << mode(f, 0, 0), 0))
        assert len(keysf) == 11
        assert s.expand(S0-Sf) == 2*Y[f]**2
        T = creation(keys0, keysf, mode(f, 0, 0))
        no_operation=s.eye(13)-T.H*T
        assert no_operation*no_operation==no_operation
        # Scalar number commutes with the elementary creator. The identity
        # compares actual rectangular source/target blocks on every input.
        ny0 = OMEGA*s.diag(0, *([1]*12))
        nyf = OMEGA*s.diag(0, *([1]*10))
        assert nyf*T - T*ny0 == s.zeros(11,13)
        J = (hf-nyf)*T - T*(h0-ny0)
        J_independent = higgs_composite(keys0,keysf,f)
        assert (J-J_independent).applyfunc(s.expand) == s.zeros(11,13)
        assert (hf*T - T*h0 - J).applyfunc(s.expand) == s.zeros(11,13)
        # c^dag maps the source bright leaf vector to the target bright leaf
        # vector. This proves the complete dark spectral residue is zero.
        e0 = s.eye(13)[:,0]
        ef = s.eye(11)[:,0]
        assert T*e0 == ef
        assert T*h0*e0 == hf*ef
        outputs.append({"family": f, "source_dimension": 13, "target_dimension": 11,
                        "S_source": str(S0), "S_target": str(Sf),
                        "same_H_rectangular_identity_exact": True,
                        "composite_independently_compiled_from_native_CAR_and_scalar": True,
                        "complete_dark_spectral_weight_exact": "0",
                        "diagnostic": numeric_transfer(h0,hf,T,[s.Rational(1,7),s.Rational(2,7),s.Rational(3,7)],s.Integer(1))})
    # The local mass response is a theorem for every positive family scale,
    # not only an expansion at y=0. Its sign excludes the requested negative
    # correction to the muon/electron ratio within this declared source.
    theta = mp.mpf(2)/9
    lambdas = [1/mp.sqrt(6)+1/mp.sqrt(3)*mp.cos(theta+2*mp.pi*k/3) for k in range(3)]
    rows = []
    baseline = sorted(l*l for l in lambdas)
    for scale in (mp.mpf("0.000001"),mp.mpf("0.001"),mp.mpf("0.01"),mp.mpf("0.1"),mp.mpf(1),mp.mpf(100)):
        mm = sorted(masses([scale*l for l in lambdas],mp.mpf(1)))
        ratios = [mm[1]/mm[0], mm[2]/mm[0]]
        assert ratios[0] > baseline[1]/baseline[0]
        assert ratios[1] > baseline[2]/baseline[0]
        rows.append({"y_over_omega": str(scale), "mu/e": str(ratios[0]), "tau/e": str(ratios[1]),
                     "Q": str(sum(mm)/sum(mp.sqrt(x) for x in mm)**2)})
    out = {"domain": "local three-family quantum-Higgs memory, no spatial hopping",
           "source_dimension": 5*2**18, "clock": "exp(-itH), hbar=1",
           "reference": "normalized negative-energy eigenvector of the vacuum-generated invariant star",
           "reference_stationary": True, "reference_claimed_global_ground": False,
           "explicit_preparation":preparation,
           "gauge_compensation": compensated_current_check(),
           "architectural_inputs": ["three copies of native lepton CAR modes", "one shared five-state scalar memory",
                                    "diagonal family Yukawa coefficients", "scalar penalty omega>0", "no hopping"],
           "same_native_one_family": native, "families": outputs,
           "symbolic_spectrum": {"S": "4 sum_f y_f^2", "D0": "sqrt(omega^2+4S)",
                                 "D_f": "sqrt(omega^2+4S-8y_f^2)",
                                 "low_transfer": "(D0-D_f)/2 = 4y_f^2/(D0+D_f)",
                                 "high_transfer": "(D0+D_f)/2", "dark_transfer": "(omega+D0)/2",
                                 "target_dark_multiplicity": 9,
                                 "low_transfer_leading_term": "2y_f^2/omega",
                                 "exact_heavy_light_ratio": "(y_h/y_l)^2 (D0+D_l)/(D0+D_h) > (y_h/y_l)^2",
                                 "composite_residue": "transfer^2 * elementary_residue"},
           "finite_memory_mass_rows": rows,
           "baseline_squared_circulant_ratios": [str(baseline[1]/baseline[0]),str(baseline[2]/baseline[0])],
           "physical_rest_mass_proved": False, "primary_common_source_transfer": False}
    out["empirical_inputs"] = []
    return out
