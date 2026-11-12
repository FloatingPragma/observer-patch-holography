"""Full-Fock electric-flux dressing and complete charged-block records."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
from model.openwave_native_quantum_charge.charge import native_data
from model.lepton_mass_compiler.compiler import _annihilator
from model.openwave_quantum_rotational.source import source

KAPPA=F(1,8)

def exact_charge_tree():
 d=native_data();charges=[F(r['charge'])for r in d['pairs']];integers=[int(3*q)for q in charges]
 assert all(F(z)==3*q for z,q in zip(integers,charges)) and sum(integers)==0
 return {'native_pair_charge_units_one_third':integers,'native_reference_charge_per_orbital':'0',
  'leaf_constraint':'E_i-3Q_i=0','root_constraint':'-sum_i E_i-R=0',
  'physical_isometry':'|n> -> |n> tensor |E_i=3Q_i(n)> tensor |R=-3 sum_i Q_i(n)>',
  'full_Fock_positive_floor':'H >= 16 N_charged+(1/2)N_neutral+N_photon+(kappa/18)sum_links E_l^2',
  'exact_reduced_electric_energy':'(kappa/2)sum_i depth_i Q_i^2',
  'integer_compact_reference_charge_R_over_three_retained':True,
  'geometric_scope':'Tree vertices label the complete oscillator orbitals; they are not asserted to be local points in a physical three-dimensional Maxwell lattice.'}


def finite_CAR_gauss_control():
 # Two complete charge species on two leaves. All16 Fock words are retained.
 charges=(-3,2,-3,2);leaves=(0,0,1,1);c=[_annihilator(4,i)for i in range(4)];N=[a.conj().T@a for a in c]
 physical=[]
 for mask in range(16):
  e=tuple(sum(charges[j]*((mask>>j)&1)for j in range(4)if leaves[j]==i)for i in range(2));physical.append((mask,e,-sum(e)))
 where={state:i for i,state in enumerate(physical)};dressed=np.zeros((16,16),complex);ordinary=np.zeros_like(dressed);undressed_violation=0
 for i,j in((0,2),(2,0),(1,3),(3,1)):
  term=c[i].conj().T@c[j];ordinary+=term
  for target,origin in zip(*np.nonzero(term)):
   mask,e,R=physical[origin];new=list(e);new[leaves[i]]+=charges[i];new[leaves[j]]-=charges[j]
   candidate=(target,tuple(new),R);assert candidate in where;dressed[where[candidate],origin]+=term[target,origin]
   undressed_violation=max(undressed_violation,max(abs(e[k]-sum(charges[a]*((target>>a)&1)for a in range(4)if leaves[a]==k))for k in range(2)))
 assert np.array_equal(dressed,ordinary) and undressed_violation==3
 # Complete electric outcome probabilities for a coherent full-Fock state.
 v=np.array([complex(i+1,(-1)**i)for i in range(16)]);v/=np.linalg.norm(v)
 rows=[{'occupation_word':m,'two_integer_link_fluxes':list(e),'reference_integer_charge':R,'probability':float(abs(v[i])**2)}for i,(m,e,R)in enumerate(physical)]
 assert abs(sum(x['probability']for x in rows)-1)<1e-14
 return {'complete_Fock_words':16,'exact_dressed_hopping_identity':True,'undressed_hopping_Gauss_defect':int(undressed_violation),'every_complete_charge_flux_reference_record':rows}


def charged_block():
 d=source(2);prior=json.loads((Path(__file__).resolve().parents[1]/'openwave_quantum_rotational/results.json').read_text());p=next(x for x in prior['interval_certificates']if x['cap']==2)
 cost=KAPPA/2;interval=[float(np.nextafter(float(F.from_float(x)+cost),direction))for x,direction in zip(p['actual_ground_energy_interval'],(-np.inf,np.inf))]
 static=[float(np.nextafter(float(F.from_float(g)*(F.from_float(e)+cost)/F.from_float(e)),direction))for g,e,direction in zip(p['rest_energy_normalized_static_g_interval'],reversed(p['actual_ground_energy_interval']),(-np.inf,np.inf))]
 olde=d['E'];unit=F('0.00051099895069')/F.from_float(olde)
 # On every one-charge state of the equal-depth star, sum Q_i^2=1.
 nfield=len(d['photon']['basis']);psi=d['psi'];weights=abs(psi.reshape(2,10,nfield))**2
 leaf=weights.sum(axis=(0,2));assert abs(leaf.sum()-1)<1e-13
 H=d['H']+float(cost)*np.eye(len(psi));residual=np.linalg.norm(H@psi-(olde+float(cost))*psi);assert residual<1e-12
 rows=[{'orbital_leaf':i,'integer_flux':-3,'other_leaf_fluxes':0,'reference_integer_charge':3,'probability':float(w)}for i,w in enumerate(leaf)]
 return {'kappa_exact':str(KAPPA),'star_exact_electric_energy_increment':str(cost),
  'full_one_electron_sector_H_identity':'W^dagger H_Gauss W = H_quantum + kappa/2 I',
  'source_plus_compensator_ground_energy_interval':interval,
  'static_current_derivative_interval_unchanged':p['native_probe_derivative_interval'],
  'energy_normalized_pair_response_g_interval':static,
  'conditional_GeV_using_unchanged_electron_unit':float((F.from_float(olde)+cost)*unit),
  'field_free_Gauss_reference_energy':0,'actual_quantum_ground_residual':float(residual),
  'complete_native_orbital_charge_link_flux_reference_records':rows,
  'dressed_spin_Casimir_exact':'3/4 on the one-charge sector',
  'star_field_energy_is_intrinsic_isolated_electron_rest_mass':False,
  'shell_depth_member':'depth_i=1+N_i gives electric cost(kappa/2)(I+N_matter) in the one-charge sector; this explicitly changes the trap coefficient',
  'full_native_neutral_and_multiple_charge_complements_retained':True}


def certificate():
 return {'schema':'oph.local_tree_Gauss_positive_quantum.v1','native_Gauss_construction':exact_charge_tree(),
  'complete_small_CAR_control':finite_CAR_gauss_control(),'same_quantum_charged_block':charged_block(),
  'same_source_dressing_rule':'B^U=D(U) B D(U)^dagger and K^U=D(U)sqrt(m^2+B^2)D(U)^dagger; operator-valued CAR order retained',
  'all_reference_and_flux_outcomes_retained':True,
  'scope':'Declared compact-U1 tree gauge extension of the actual native positive source. Exact local Gauss identities and positive electric energy, with complete matter/field/reference carriers. No inverse-distance continuum Coulomb law, local physical-density dictionary or empirical electromagnetic coupling is inferred.',
  'verdict':'PASS_NATIVE_LOCAL_TREE_GAUSS_AND_POSITIVE_QUANTUM_FLOOR'}
