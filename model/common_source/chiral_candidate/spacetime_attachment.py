"""Exact principal-symbol attachment to the retained native metric chart.

This is a change of Clifford basis and a frozen-chart tetrad pullback. It
does not assert a shared full spatial Dirac history, matching boundary
domains, state preparation, global Spin transport, or physical particles.
The actual incident-clock collective response also supplies a sharp test of
the proposed history attachment. Its conditioned response is a strict
contraction; it cannot intertwine with an unmeasured closed quadratic Nambu
propagator. An effective field with detector/environment channels is outside
that obstruction. All retained packets are read-only. No measured masses enter.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path
from types import MappingProxyType

import sympy as s

from . import action
from ..hierarchy_diagnosis import collective_response as collective


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "model/observer_particles/results"
SOURCE_FILES = (
    "model/observer_particles/results/geometric_matter_join.json",
    "model/observer_particles/results/shared_spatial_matter.json",
    "model/observer_particles/geometric_matter_join.py",
    "model/observer_particles/shared_spatial_matter.py",
    "model/observer_particles/spacetime_correspondence.py",
    "model/relative_observer_geometry/intrinsic_geometry.py",
    "Lean/Screen/ObserverFieldInterfaces.lean",
    "model/common_source/hierarchy_diagnosis/collective_response.py",
    "model/repair_scheduler/coins.py",
)


def _clean(matrix):
    return s.ImmutableMatrix(matrix).applyfunc(s.simplify)


def _zero(matrix):
    return all(s.simplify(x)==0 for x in matrix)


def _entries(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def _packets():
    geometry = json.loads((RESULTS/"geometric_matter_join.json").read_text())
    spatial = json.loads((RESULTS/"shared_spatial_matter.json").read_text())
    # Check the complete existing provenance rather than running its producers.
    for packet in (geometry,spatial):
        for name,expected in packet["source_sha256"].items():
            if _hash(ROOT/name)!=expected:
                raise ArithmeticError(f"retained spacetime source changed: {name}")
        for name,expected in packet["artifacts"].items():
            if _hash(RESULTS/name)!=expected:
                raise ArithmeticError(f"retained spacetime artifact changed: {name}")
    return geometry,spatial


@lru_cache(maxsize=1)
def frozen_geometry():
    """Exact immutable A,L,N,g and coframe of the retained comparison chart."""
    packet = _packets()[0]["same_source_metric_and_clock"]
    parse = lambda rows:s.ImmutableMatrix([[s.sympify(x) for x in row] for row in rows])
    a,g,coframe = (parse(packet[name]) for name in
                   ("spatial_principal_tensor","metric","orthonormal_coframe"))
    lapse = s.simplify(s.sympify(packet["lapse"]))
    lower = _clean(a.cholesky())
    eta = s.diag(-1,1,1,1)
    if (not _zero(lower*lower.T-a) or s.simplify(lapse*lapse+g[0,0])!=0
            or not _zero(coframe.T*eta*coframe-g)
            or not _zero(g-s.diag(-lapse**2,lapse**2*a.inv()))):
        raise ArithmeticError("retained metric/tetrad identities failed")
    return MappingProxyType({"A":a,"L":lower,"N":lapse,"g":g,"coframe":coframe})


def pauli_matrices():
    return (s.ImmutableMatrix([[0,1],[1,0]]),s.ImmutableMatrix([[0,-s.I],[s.I,0]]),
            s.ImmutableMatrix(s.diag(1,-1)))


def clifford_intertwiner():
    """W maps retained Dirac alpha_j to diag(-sigma_j,+sigma_j)."""
    identity = s.eye(2)
    return s.ImmutableMatrix(s.BlockMatrix([[identity,-identity],[identity,identity]]).as_explicit()/s.sqrt(2))


def _momentum(momentum):
    values = tuple(momentum)
    if len(values)!=3 or any(isinstance(x,bool) for x in values):
        raise ValueError("three exact real momentum components required")
    values = tuple(s.sympify(x) for x in values)
    if any(x.is_real is not True or x.atoms(s.Float) for x in values):
        raise ValueError("use exact real momenta, not floating point inputs")
    return s.ImmutableMatrix(values)


def weyl_pullback(momentum):
    """q=L^T p in the same native-chart scalar principal tensor A=L L^T."""
    return _clean(frozen_geometry()["L"].T*_momentum(momentum))


def weyl_symbol(momentum):
    """Coordinate-time left-Weyl kinetic symbol -sigma.(L^T p)."""
    q = weyl_pullback(momentum)
    return _clean(-sum((q[j]*p for j,p in enumerate(pauli_matrices())),s.zeros(2)))


def coordinate_mass(proper_kernel):
    """Frozen static lapse scaling K_coordinate=N K_proper.

    The neutral Higgs component norm and every mass frequency transform by
    this factor. This does not assert that the new source clock has already
    been identified with native coordinate time or proper time.
    """
    if isinstance(proper_kernel,s.MatrixBase):
        return _clean(frozen_geometry()["N"]*proper_kernel)
    value = s.sympify(proper_kernel)
    if value.is_real is not True or value.atoms(s.Float):
        raise ValueError("an exact real proper mass parameter required")
    return s.simplify(frozen_geometry()["N"]*value)


def native_incident_response(law, deadline=2):
    """Actual collective response conditioned on one recorded incident tick.

    The source seam selector and finite proposal timeout remain unchanged.
    Success is a state-independent event: it depends on selected seams, not
    on currents or the work state. The complete outcome masses are returned
    separately from the normalized successful response.
    """
    if type(deadline) is not int or not 1 <= deadline <= 8:
        raise ValueError("a source proposal deadline from one through eight is required")
    data = collective.transfer(law)
    weights = collective.preparation.selectors(1)["completion_accepted_seam_weights"]
    cap = set(collective.preparation.CAPS[0])
    incident = sum((s.Rational(weight) for edge,weight in weights.items()
                    if cap.intersection(edge)),s.Rational(0))
    remote = sum((s.Rational(weight) for edge,weight in weights.items()
                  if not cap.intersection(edge)),s.Rational(0))
    timeout = 1-incident-remote
    left,right = data["incident"],data["remote"]
    one = s.ones(10,1)
    killed = data["killing"]
    if (any(value < 0 for value in (*left,*right))
            or left*one != (incident-killed)*one
            or one.T*left != (incident-killed)*one.T
            or right*one != remote*one or one.T*right != remote*one.T):
        raise ArithmeticError("native incident/remote row or column identities failed")
    waiting = sum((right**n for n in range(deadline)),s.zeros(10))
    clock_waiting = sum(remote**n for n in range(deadline))
    success = incident*clock_waiting
    terminal = {"success":success,"selector_timeout":timeout*clock_waiting,
                "deadline_timeout":remote**deadline}
    if sum(terminal.values()) != 1:
        raise ArithmeticError("native stopped clock lost a terminal flag")
    response = _clean(left*waiting/success)
    factor = 1-killed/incident
    if (not 0 < factor < 1 or response*one != factor*one
            or one.T*response != factor*one.T or any(value < 0 for value in response)):
        raise ArithmeticError("conditioned native response is not the claimed contraction")
    # Positive equal row/column sums give ||R||_2 <= sqrt(||R||_1||R||_inf)
    # =factor. The uniform vector attains equality, so this is the exact
    # operator norm and spectral radius, with no numerical rank decision.
    outer = tuple(i for i,(site,_,_) in enumerate(data["fields"])
                  if site not in collective.preparation.CORE)
    if len(outer)!=9 or any(response[i,j] or response[j,i]
                          for i in outer for j in range(10) if j not in outer):
        raise ArithmeticError("the charged outer response ceased to be invariant")
    return MappingProxyType({"response":response,"outer_indices":outer,
        "incident_probability":incident,"remote_probability":remote,
        "selector_timeout_probability":timeout,"coherence_killing":killed,
        "terminal_probabilities":MappingProxyType(terminal),
        "operator_norm_and_spectral_radius":factor,
        "relative_unitary_intertwiner_defect_lower":1-factor})


def record_response_control(law):
    """Sum actual complete Ramsey outputs, keeping all timeout probabilities."""
    source = collective.clock_experiment(law)
    data = native_incident_response(law,2)
    probabilities = {flag:s.Rational(value) for flag,value in
                     source["terminal_probabilities_exact"].items()}
    if probabilities != dict(data["terminal_probabilities"]):
        raise ArithmeticError("complete clock execution disagrees with incident normalization")
    contrasts = {flag:sum((s.Rational(value) for value in row.values()),s.Rational(0))
                 for flag,row in source["terminal_signed_Ramsey_response_exact"].items()}
    if contrasts["success"]/probabilities["success"] != data["operator_norm_and_spectral_radius"]:
        raise ArithmeticError("actual recorded Ramsey contrast is not the predicted one-tick loss")
    return {"native_law":law,"clock_cap":source["clock_cap"],
        "work_proposal_deadline":2,"successful_incident_ticks":1,
        "complete_terminal_probabilities":{k:str(v) for k,v in probabilities.items()},
        "raw_signed_Ramsey_contrast_by_terminal_flag":{k:str(v) for k,v in contrasts.items()},
        "successful_conditional_signed_Ramsey_contrast":str(
            contrasts["success"]/probabilities["success"]),
        "input_signed_Ramsey_contrast":"1",
        "enumerated_terminal_paths":source["enumerated_terminal_paths"],
        "complete_record_and_detector_sha256":source["complete_record_and_detector_sha256"],
        "state_preparation_and_coherent_Ramsey_detector_are_declared":True,
        "native_pulse_compilation_reexecuted":False,
        "conditional_normalization_replaces_complete_instrument":False}


def open_response_control(law):
    """A complete dephasing dilation defeats an overbroad unitary no-go.

    This mathematical two-Kraus control matches one uniform Ramsey response,
    not the full native coordinate channel or a selected physical apparatus.
    Its environment must be included in an effective-field history.
    """
    factor = native_incident_response(law)["operator_norm_and_spectral_radius"]
    probability = (1-factor)/2
    identity,z = s.eye(2),s.diag(1,-1)
    operators = (s.sqrt(1-probability)*identity,s.sqrt(probability)*z)
    field = s.Matrix([[0,1],[0,0]])
    charge = s.diag(1,0)
    channel = lambda x:sum((k*x*k.H for k in operators),s.zeros(2))
    if (sum((k.H*k for k in operators),s.zeros(2)) != identity
            or channel(field) != factor*field or channel(charge)!=charge):
        raise ArithmeticError("complete dephasing dilation control failed")
    return {"phase_flip_probability":str(probability),
        "K0":"sqrt(1-p)*I","K1":"sqrt(p)*Z",
        "trace_preserving_exact":True,"charge_populations_preserved_exact":True,
        "uniform_Ramsey_response_factor":str(factor),
        "full_native_spatial_channel_or_physical_detector_identified":False,
        "source_pulse_dilation_compiled":False}


def incident_history_attachment():
    """A source-measured limit on a closed quadratic field-history bridge."""
    laws = {}
    for law in ("C","E"):
        data = native_incident_response(law)
        deadlines = {}
        for deadline in (1,2,3,8):
            row = native_incident_response(law,deadline)
            if row["operator_norm_and_spectral_radius"] != data["operator_norm_and_spectral_radius"]:
                raise ArithmeticError("conditioned incident contraction depends on timeout deadline")
            deadlines[str(deadline)] = {k:str(v) for k,v in row["terminal_probabilities"].items()}
        sites = [site for site,_,_ in collective.fields()]
        before,after = sites.index(2),sites.index(7)
        one_tick_short = native_incident_response(law,1)["response"][after,before]
        one_tick_long = data["response"][after,before]
        if one_tick_short != 0 or one_tick_long <= 0:
            raise ArithmeticError("same-count conditional spatial-history control failed")
        laws[law] = {"incident_probability_per_proposal":str(data["incident_probability"]),
            "remote_probability_per_proposal":str(data["remote_probability"]),
            "coherence_killing_per_proposal":str(data["coherence_killing"]),
            "conditional_response_matrix_at_deadline_two":_entries(data["response"]),
            "exact_operator_norm_and_spectral_radius":str(data["operator_norm_and_spectral_radius"]),
            "relative_closed_unitary_intertwiner_defect_lower":str(data["relative_unitary_intertwiner_defect_lower"]),
            "uniform_outer_charged_eigenvector_dimension":9,
            "complete_terminal_probability_controls":deadlines,
            "same_incident_count_different_deadline_response":{
                "source_site":2,"detector_site":7,"incident_count":1,
                "deadline_one_conditional_response":str(one_tick_short),
                "deadline_two_conditional_response":str(one_tick_long),
                "full_native_waiting_record_can_be_replaced_by_tick_count_alone":False},
            "actual_detector_control":record_response_control(law),
            "open_field_environment_countercontrol":open_response_control(law)}
    return {"clock":"one actual cap-incident native work event; bounded proposal timeout",
        "source_field":"X_j=|0,j><5,j| on the fixed binary total-two source",
        "response_domain":"complete native records, coarse-grained to the explicitly conditioned successful incident response",
        "target_domain":"unmeasured closed quadratic 180-dimensional Nambu propagator exp(-i H Delta t)",
        "norm_proof":"nonnegative equal row and column sums lambda imply ||R||_2=lambda; the uniform vector attains it",
        "intertwiner_bound":"||U J-J R||_HS >= (1-lambda)||J||_HS for every unitary U and every linear J",
        "nonzero_exact_closed_quadratic_response_intertwiner_exists":False,
        "constant_tetrad_basis_change_or_clock_regraduation_removes_bound":False,
        "effective_field_with_detector_or_environment_channel_excluded":False,
        "interacting_thermal_or_continuum_particle_theory_excluded":False,
        "loss_factor_is_physical_rest_mass":False,
        "laws":laws}


def report():
    """JSON-ready exact frozen-chart bridge and actual finite-domain boundary."""
    geometry,spatial = _packets()
    data = frozen_geometry()
    w,identity = clifford_intertwiner(),s.eye(2)
    beta = s.diag(1,1,-1,-1)
    alpha = tuple(s.BlockMatrix([[s.zeros(2),p],[p,s.zeros(2)]]).as_explicit()
                  for p in pauli_matrices())
    if not _zero(w*w.H-s.eye(4)):
        raise ArithmeticError("Clifford basis change is not unitary")
    for native,p in zip(alpha,pauli_matrices(),strict=True):
        if not _zero(w*native*w.H-s.diag(-p,p)):
            raise ArithmeticError("native/chiral principal Clifford intertwiner failed")
    if not _zero(w*beta*w.H-s.BlockMatrix([[s.zeros(2),identity],[identity,s.zeros(2)]]).as_explicit()):
        raise ArithmeticError("native/chiral scalar mass intertwiner failed")

    p=s.ImmutableMatrix(s.symbols("p_r p_x p_y",real=True))
    energy,mass=s.symbols("E m_proper",real=True)
    h=weyl_symbol(p)
    scalar=(p.T*data["A"]*p)[0]
    if s.simplify((energy*s.eye(2)-h).det()-(energy**2-scalar))!=0:
        raise ArithmeticError("scalar/Weyl characteristic cones differ")
    q=weyl_pullback(p)
    old=sum((q[j]*alpha[j] for j in range(3)),s.zeros(4))+coordinate_mass(mass)*beta
    new=s.BlockMatrix([[h,coordinate_mass(mass)*identity],
                       [coordinate_mass(mass)*identity,-h]]).as_explicit()
    if not _zero(w*old*w.H-new):
        raise ArithmeticError("massive frozen-chart Dirac/Nambu symbol bridge failed")
    square=scalar+(data["N"]*mass)**2
    if not _zero(new*new-square*s.eye(4)):
        raise ArithmeticError("lapse-scaled mass shell failed")
    covector=s.ImmutableMatrix([energy,*p])
    if s.simplify((covector.T*data["g"].inv()*covector)[0]+mass**2
                  -(square-energy**2)/data["N"]**2)!=0:
        raise ArithmeticError("coordinate and proper Lorentz mass shells differ")

    kappa,c=s.symbols("kappa c",positive=True)
    pnew=s.ImmutableMatrix(s.symbols("pnew_r pnew_x pnew_y",real=True))
    native_p=c/kappa*data["L"].T.inv()*pnew
    if not _zero(kappa*data["L"].T*native_p-c*pnew):
        raise ArithmeticError("clock/chart regraduation failed")

    extra=spatial["extra_modes_and_vacuum"]
    if (extra["node_modes"],extra["edge_modes"],extra["zero_energy_modes"])!=(720,2850,2130):
        raise ArithmeticError("retained incidence-domain dimensions changed")
    boundary=geometry["radial_first_order_boundary_control"]
    bindings={name:_hash(ROOT/name) for name in (*SOURCE_FILES,*collective.preparation.SOURCE_PATHS)}
    bindings[Path(__file__).relative_to(ROOT).as_posix()]=_hash(Path(__file__))
    bindings["model/common_source/chiral_candidate/action.py"]=_hash(Path(action.__file__))
    for packet in (geometry,spatial):
        bindings.update(packet["source_sha256"])
        bindings.update({(RESULTS/name).relative_to(ROOT).as_posix():digest for name,digest in packet["artifacts"].items()})
    return {
        "schema":"oph.chiral_candidate.frozen_spacetime_attachment.v1",
        "status":"EXACT_COMMON_FROZEN_CHART_PRINCIPAL_SYMBOL_AND_CLIFFORD_BRIDGE",
        "frozen_geometry":{"spatial_principal_tensor_A":_entries(data["A"]),
                           "Cholesky_L":_entries(data["L"]),"lapse_N":str(data["N"]),
                           "Lorentz_metric_g":_entries(data["g"]),"coframe":_entries(data["coframe"]),
                           "signature":[-1,1,1,1],"finite_increment_is_continuum_tensor":False,
                           "measure_density_one_is_comparison_convention":True,
                           "native_readout":geometry["same_source_metric_and_clock"]["source_increment_readout"]},
        "Clifford_bridge":{"W4":_entries(w),"unitarity_exact":True,
                           "three_alpha_intertwiners_exact":True,"beta_intertwiner_exact":True,
                           "retained_native_Clifford_word_errors_diagnostic":geometry["native_Clifford_and_shared_symbols"]["native_basis_word_errors"],
                           "native_word_errors_newly_interval_certified":False,
                           "Nambu_redundancy_and_conjugated_gauge_action_must_be_retained":True},
        "same_metric_symbols":{"q_definition":"q=L^T p","left_Weyl":"h=-sigma.q",
                               "scalar_principal_square":str(scalar),
                               "Weyl_characteristic_polynomial":"E^2-p^T A p",
                               "same_scalar_characteristic_cone_exact":True,
                               "massive_coordinate_energy_square":"p^T A p + N^2 m_proper^2",
                               "proper_mass_shell_equals_coordinate_mass_shell_exact":True,
                               "full_curved_Dirac_square_equals_scalar_B":False},
        "mass_and_clock":{"lapse_rule":"K_coordinate=N*K_proper",
                          "clock_rule":"t_native=kappa*s_new; H_new=kappa*H_native",
                          "same_reference_thermal_rule":"beta_new=beta_native/kappa; beta*H is invariant",
                          "fixed_beta_one_references_require_kappa_one_in_this_dictionary":True,
                          "linear_chart":"x_native=(kappa/c)*L*x_new; p_native=(c/kappa)*L^(-T)*p_new",
                          "new_isotropic_principal_speed_exact":str(action.source_coefficients()["Omega"]/2),
                          "matching_principal_symbols_after_declared_regraduation_exact":True,
                          "mass_in_new_parameter":"K_new=kappa*N*K_proper",
                          "clock_rate_in_proper_time":"Omega_proper=Omega_new/(kappa*N)",
                          "frequency_ratio_when_both_rates_regraduated":"K_new/Omega_new = K_proper/Omega_proper",
                          "actual_native_incident_and_new_clock_histories_identified":False},
        "finite_domain_obstruction":{"incidence_node_modes":720,"incidence_edge_modes":2850,
                                     "complete_incidence_dimension":3570,"retained_edge_cycle_zero_modes":2130,
                                     "two_component_Weyl_on_same_scalar_nodes_dimension":1440,
                                     "complete_unitary_field_identification_possible":False,
                                     "discarding_edge_kernel_is_an_allowed_identification":False,
                                     "reason":"Different dimensions and retained observable persistent cycles; incidence form indices are not Weyl spin indices."},
        "existing_spatial_and_clock_contract":{
            "local_scalar_factorization":spatial["geometry"]["exact_local_factorization"],
            "ideal_joint_finite_clock_factorization":spatial["geometry"]["finite_angular_ideal_joint_factorization"],
            "Bwave_is_exact_Floquet_square_of_executed_word":spatial["geometry"]["Bwave_is_exact_Floquet_square_of_executed_two_axis_word"],
            "old_incident_clock":spatial["conditional_joint_action_state_clock"]["clock"],
            "old_declared_joint_reference":spatial["conditional_joint_action_state_clock"]["joint_reference"],
            "old_full_joint_reference_preparation_executed":spatial["conditional_joint_action_state_clock"]["full_joint_Gaussian_reference_preparation_newly_executed"]},
        "actual_incident_response_attachment":incident_history_attachment(),
        "radial_boundary":boundary,
        "remaining_operator_premises":geometry["full_operator_contract"],
        "scope":{"principal_and_frozen_chart_attachment_proved":True,
                 "new_full_native_spatial_Dirac_history_executed":False,
                 "same_complete_scalar_and_chiral_field_history_proved":False,
                 "same_modular_reference_and_vacuum_prepared":False,
                 "field_density_and_state_intertwiner_completed":False,
                 "boundary_spin_connection_and_curvature_controls_completed":False,
                 "physical_particles_or_Standard_Model_masses_identified":False,
                 "new_empirical_mass_inputs":[]},
        "source_bindings":dict(sorted(bindings.items()))}


if __name__=="__main__":
    print(json.dumps(report(),indent=2,allow_nan=False))
