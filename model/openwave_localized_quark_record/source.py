"""Complete relational string instruments, with native CAR and Gauss checks."""
from fractions import Fraction as F
from itertools import product
import sympy as s
from model.openwave_quark_flux_records.flux import native_holes,complete_native_preparation
from model.openwave_quark_strings.strings import reduced_string
from model.openwave_reference_binding.strings import gauss_image
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES

ELECTRIC=F(1,16)
COLOUR=F(1,8)


def path(length):
    if not isinstance(length,int) or length<1:raise ValueError('positive path length required')
    return tuple((edge,edge,edge+1,True)for edge in range(length))


def energy(q,length):
    """Full rotor Casimir and flux energy on the normalized string state."""
    return F(2)+length*(ELECTRIC*q*q/2+COLOUR*F(4,3)/2)


def exact_gauss(family,length):
    p=path(length);state=reduced_string(family,p)
    for vertex in range(length+1):
        for current in range(12):
            assert not gauss_image(family,p,state,vertex,current)
    return 12*(length+1)


def profile(vertices,endpoint,q):
    if not 0<endpoint<vertices:raise ValueError('endpoint must lie inside the carrier')
    charge=[F(0)for _ in range(vertices)];charge[endpoint]=F(q,3)
    reference=[F(0)for _ in range(vertices)];reference[0]=-F(q,3)
    flux=[q if e<endpoint else 0 for e in range(vertices-1)]
    casimir=[F(4,3)if e<endpoint else F(0)for e in range(vertices-1)]
    # Orientation is from the compensating reference towards the source hole.
    divergence=[(flux[x-1]if x else 0)-(flux[x]if x<vertices-1 else 0)for x in range(vertices)]
    assert all(F(divergence[x],3)==charge[x]+reference[x]for x in range(vertices))
    return {'vertices':vertices,'endpoint':endpoint,'source_charge_exact':list(map(str,charge)),
            'reference_charge_exact':list(map(str,reference)),'integer_flux_on_all_edges':flux,
            'colour_Casimir_on_all_edges_exact':list(map(str,casimir)),
            'complete_idle_energy_exact':str(energy(q,endpoint)),
            'outside_path_charge_and_flux_exact':'0'}


def records(q,length,phase,branch_weight):
    """A physical vacuum/string two-plane; its complete operator acts on full Fock."""
    w=s.Rational(branch_weight);phi=s.sympify(phase);E=s.Rational(energy(q,length))
    assert 0<=w<=1
    rows=[]
    # Exact times stated as fractions of the source's own coherence period.
    for phase_time in (s.Integer(0),s.pi/2,s.pi):
        coherence=2*s.sqrt(w*(1-w))*s.exp(s.I*(phi-phase_time))
        probabilities=[s.simplify((1+sign*component)/4)
                       for component in(s.re(coherence),s.im(coherence))for sign in(1,-1)]
        assert all(p>=0 for p in probabilities)and sum(probabilities)==1
        rows.append({'source_time_exact':str(phase_time/E),'X_plus_X_minus_Y_plus_Y_minus_exact':list(map(str,probabilities)),
                     'relational_string_mean_exact':str(s.simplify(coherence/2))})
    # A terminal alternative measures the entire diagonal charge/flux profile.
    return {'native_integer_charge':q,'path_length':length,'preparation_string_weight_exact':str(w),
            'preparation_phase_exact':str(phi),'all_time_records':rows,
            'alternative_charge_flux_records':{'vacuum_exact':str(1-w),'entire_string_profile_exact':str(w),'all_other_Fock_flux_outcomes_exact':'0'},
            'preparation_mean_work_exact':str(w*E),
            'pointer_nonselective_mean_work_exact':str((s.Rational(1,2)-w)*E),
            'diagonal_profile_nonselective_mean_work_exact':'0'}


def run():
    native=native_holes();pump=complete_native_preparation();gauges=representation()[2]
    checks=0;families=[]
    for row in native:
        family=row['family'];q=int(3*F(row['native_charge_exact']));ids=FAMILIES[family]
        # The actual current, not an assigned charge label, has nonzero U1 weight.
        assert gauges[10].extract(ids,ids)==s.zeros(3)
        assert -gauges[11].extract(ids,ids)==s.Rational(q,3)*s.eye(3)
        for length in(1,2,3):checks+=exact_gauss(family,length)
        profiles=[profile(L,x,q)for L,x in((3,1),(5,1),(9,1),(5,2),(9,2),(9,3))]
        for endpoint in(1,2):
            selected=[p for p in profiles if p['endpoint']==endpoint]
            assert len({p['complete_idle_energy_exact']for p in selected})==1
        families.append({'family':family,'native_charge_exact':str(F(q,3)),
             'native_hole_and_reference':row,'profiles':profiles,
             'complete_records':[records(q,l,phi,w)for l in(1,2,3)for phi in(0,s.pi/2)for w in(F(0),F(1,2),F(1))]})
    return {'schema':'oph.localized_native_quark_string_record.v1',
      'verdict':'PASS_NATIVE_LOCALIZED_RELATIONAL_STRING_RECORDS_AND_ONE_POINT_OBSTRUCTION',
      'declared_source':{'matter':'Complete thirty-CAR native source at each open-chain vertex, with its full thirty-CAR reference at vertex zero.',
          'links':'Full L2(SU3 x U1) on every edge; vacuum and fundamental words are states, not source cutoffs.',
          'idle_H':'Hole number plus reference number plus (1/32) sum E^2 +(1/16) sum C2; complete charge constraints retained.',
          'state':'Filled source, reference vacuum and Haar-constant links; gauge-invariant vacuum/string coherent preparation.',
          'external_native_hopping_or_magnetic_loop_added':False,
          'same_PW119_autonomous_source_claimed':False},
      'one_point_obstruction':{'exact_statement':'For a Gauss-invariant state and any local operator A of nonzero gauge weight, Tr(rho A)=exp(i theta q) Tr(rho A) for every theta, hence Tr(rho A)=0.',
          'native_hole_and_single_open_link_one_point_means_exact':'0',
          'relational_replacement':'B_p=P_ref1 (sum_ij r_i* [U_path]ij c_xj /sqrt(3)) P_refvac, with the native charge power and complete reference projectors.',
          'gauge_rule':'U_path transforms at its two endpoints; the native hole and quantum reference cancel both factors. Every internal factor telescopes.',
          'partial_texture_winding_inferred_from_charge':False},
      'full_instrument':{'native_complete_pump_check':pump,
          'B_squared_exact':'0','B_adjoint_B_bound':'0 <= B*B <= P_refvac',
          'quadratures':'X=B+B*, Y=i(B-B*); their norms are <=1 on the entire source.',
          'outcomes':'First test P_Gauss. The four physical pointer effects are P_Gauss(I+X)P_Gauss/4, P_Gauss(I-X)P_Gauss/4, P_Gauss(I+Y)P_Gauss/4, P_Gauss(I-Y)P_Gauss/4.',
          'exact_effect_sum':'The four physical effects sum to P_Gauss; together with I-P_Gauss they sum to I on full Fock and full rotors.',
          'seed_filter':'After the Gauss test use P_seed and P_Gauss-P_seed; stated seed succeeds with probability 1, and each failed preparation is an absorbing record.',
          'constraint_failure':'I-P_Gauss is the preceding test failure; it is not appended to an unprojected four-effect identity instrument.',
          'prepared_two_plane':'B|seed>=|path>, B*|path>=|seed>; all other states remain in the instrument.',
          'complement_work':'For general input Tr[H(Phi(rho)-rho)] is retained; displayed work is evaluated only on the specified two-plane.'},
      'exact_native_Gauss_generator_checks':checks,'families':families,
      'finite_theorem_shortcut':'The reference-vacuum CAR contraction, exact gauge telescoping and orthogonal-record instrument give a full-source observable without enumerating the infinite Peter-Weyl carrier.',
      'gauge_invariant_localized_order_parameter_record_constructed':True,
      'source_spatial_profile_resolved':True,'self_sustained_fractional_texture_winding_proved':False,
      'autonomous_moving_quark_or_confinement_proved':False,'physical_quark_rest_mass_GeV_derived':False,
      'new_empirical_calibration':False,'new_match_with_nature':False,
      'scope':'A localized, reference-supported open-string coherence and complete native fractional-charge/flux records under declared idle electric dynamics. Its endpoint is prepared by a controlled service. This is not a low-energy quark pole, a spontaneous soliton or a physical QCD identification.'}
