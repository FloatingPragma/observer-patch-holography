"""Exact full-source cubic transport; no invariant-compression assumption."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_sector_threshold import model as sector
from model.particle_primary_composite_spin import model as spin
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
 'model/particle_primary_composite_transport/__init__.py',
 'model/particle_primary_composite_transport/__main__.py',
 'model/particle_primary_composite_transport/model.py',
 'model/particle_primary_composite_transport/test_model.py',
 'model/particle_primary_composite_transport/PROOF.md',
 'model/particle_primary_composite_transport/README.md',
 'model/particle_primary_sector_threshold/PROOF.md',
 'model/particle_primary_spectral_scale/PROOF.md',
 'model/particle_primary_locality/PROOF.md',
 'model/particle_common_contract/PROOF.md',
 *spin.SOURCE_INPUTS,
)


@lru_cache(maxsize=1)
def spin_blocks():
 pairs=sector.branch.response.beta_pairs()
 chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
 gamma=sector.branch.response.density.gamma_data();beta=gamma[3]
 answer=tuple(sector.branch.apparatus.clean(chart.H*((-beta+beta*g)/2)*chart) for g in gamma[:3])
 for h in answer:
  if h[:2,:2]!=-s.eye(2)/2 or h[2:,2:]!=s.eye(2)/2:
   raise ArithmeticError('actual Wilson beta-preserving spin block changed')
 return answer


def cubic_between_vectors(source_vector,target_vector,axis=0,link=None,*,electric_denominators=False):
 """Unnormalized exact 120-mode CAR matrix element for minimal transfer.

 All three required creation/annihilation pairs are taken from the literal
 Wilson matrix. At minimum order no other pair can connect these masks.
 The proof justifies identity-link evaluation and the actual link cutoff.
 """
 if axis not in (0,1,2):raise ValueError('three source axes required')
 R=s.eye(15) if link is None else s.Matrix(link)
 if R.shape!=(15,15):raise ValueError('actual fifteen-dimensional native link required')
 h=spin_blocks()[axis];omega=sum(1<<i for i in range(30,60));total=s.S.Zero
 for source,sv in source_vector.items():
  initial=source|(omega<<60)
  for target,tv in target_vector.items():
   final=omega|(target<<60)
   removed=tuple(i for i in range(120) if (initial&~final)>>i&1)
   added=tuple(i for i in range(120) if (final&~initial)>>i&1)
   if len(removed)!=3 or len(added)!=3:
    raise ValueError('the domain is a three-defect particle or hole transfer')
   terms=[]
   for i,j in product(added,removed):
    si,ii=divmod(i,60);sj,jj=divmod(j,60)
    a,u=divmod(ii,15);b,v=divmod(jj,15)
    if (si,sj)==(1,0):value=s.conjugate(h[b,a])*s.conjugate(R[v,u])
    elif (si,sj)==(0,1):value=h[a,b]*R[u,v]
    else:continue
    if value:terms.append((i,j,value))
   if electric_denominators:
    values={i:F(105,2)*C for (_,ix),(_,C,_) in zip(sector.native_blocks(),sector.refinement.native_control()) for i in ix}
    weighted=[(i,j,c/s.Rational(values[i%15].numerator,values[i%15].denominator)) for i,j,c in terms]
    vector=sector._act({initial:s.S.One},weighted)
    vector=sector._act(vector,terms);vector=sector._act(vector,weighted)
   else:
    vector={initial:s.S.One}
    for _ in range(3):vector=sector._act(vector,terms)
   total+=s.conjugate(tv)*sv*vector.get(final,0)
 return s.simplify(s.expand(total))


@lru_cache(maxsize=1)
def native_chart():
 density=sector.branch.response.density
 raw=json.loads((density.ROOT/density.SOURCE_INPUTS[-1]).read_text())
 source_masks=tuple(sum(1<<j for j in raw['basis_subsets'][i]) for i in raw['even_indices'])
 native_masks=tuple(i for i in range(32) if i.bit_count() in (2,4))
 permutation=tuple(source_masks.index(mask) for mask in native_masks)
 charges=density.source_data.load()['charges']
 expected={i:q for q,indices in sector.native_blocks() for i in indices}
 if any(charges[permutation[i]]!=expected[i] for i in range(15)):
  raise ArithmeticError('source/native exterior chart does not intertwine the central generator')
 return permutation


@lru_cache(maxsize=1)
def nontrivial_link():
 """Exact nonreal, nondiagonal S(U3 x U2) element in the native chart."""
 z=(s.Integer(3)+4*s.I)/5
 u=s.Matrix([[s.Rational(3,5),4*s.I/5],[4*s.I/5,s.Rational(3,5)]])
 v=s.Matrix([[s.Rational(5,13),s.Rational(12,13)],[-s.Rational(12,13),s.Rational(5,13)]])
 g=s.diag(z**-2*s.diag(u,s.ones(1)),z**3*(u*v))
 source_R=sector.branch.response.density.native_representation(g)
 permutation=native_chart();R=source_R.extract(permutation,permutation)
 masks=tuple(i for i in range(32) if i.bit_count() in (2,4))
 subsets=tuple(tuple(j for j in range(5) if mask>>j&1) for mask in masks)
 direct=s.Matrix(15,15,lambda i,j:g.extract(subsets[i],subsets[j]).det() if len(subsets[i])==len(subsets[j]) else 0)
 if sector.branch.apparatus.clean(R-direct)!=s.zeros(15):
  raise ArithmeticError('the actual exterior source chart conversion failed')
 if s.simplify(g.det())!=1 or sector.branch.apparatus.clean(g.H*g)!=s.eye(5):
  raise ArithmeticError('exact gauge sample is not in S(U3 x U2)')
 return R


@lru_cache(maxsize=None)
def moment_control(delta):
 delta=tuple(delta)
 if delta not in tuple(tuple(sign*x for x in d) for _,d,_ in sector.SECTORS[:2] for sign in (1,-1)):
  raise ValueError('one of the two signed odd physical triples is required')
 w=sector.local_singlet(delta);v=w['state'];norm=w['norm_squared'];axes=[]
 for axis in range(3):
  plain=F(cubic_between_vectors(v,v,axis)/norm)
  gauged=F(cubic_between_vectors(v,v,axis,nontrivial_link())/norm)
  if plain!=F(-3,4) or gauged!=plain:raise ArithmeticError('full native composite cubic transport changed')
  perturbative=F(cubic_between_vectors(v,v,axis,electric_denominators=True)/norm)
  axes.append({'axis':axis,'identity_link_cubic_moment':plain,'nontrivial_native_link_cubic_moment':gauged,
               'third_order_decoupled_resolvent_transport_coefficient':perturbative})
 trial=sector.sector_trial(delta);packets=trial['hopping_packets'];g=trial['loop_norm_squared']
 variance=sum((r['norm_squared'] for r in packets),F())+5*g
 diagonal=sum((r['norm_squared']*r['relative_H0_energy'] for r in packets),F())+g*sum(map(sum,trial['loop_packet_matrix']))
 H0norm2=sum((r['norm_squared']*r['relative_H0_energy']**2 for r in packets),F())
 H0norm2+=g*sum((210*C)**2 for _,C,_ in sector.refinement.native_control())
 J=F(360)+F(6,7)
 second_norm_upper=sqrt_upper(H0norm2)+J*sqrt_upper(variance)
 remainder_factor=second_norm_upper**2
 tstar=F(3,2)/remainder_factor
 return {'native_occupation_increment':delta,'local_norm_squared':norm,'CAR_parity':1,
  'full_Wilson_axis_controls':axes,'first_nonzero_neighbor_moment_order':3,
  'third_order_coefficient_is_primary_endpoint_band':False,
  'neighbor_moments_through_order_three':[0,0,0,F(-3,4)],
  'bare_source_mean':trial['bare_trial_source_energy'],
  'centered_variance':variance,'centered_diagonal_third_moment':diagonal,
  'centered_second_action_norm_upper':second_norm_upper,
  'centered_H0_first_action_norm_squared':H0norm2,'full_perturbation_norm_upper':J,
  'amplitude_Taylor_remainder_upper_coefficient':remainder_factor/F(24),
  'positive_transfer_time_upper':tstar,
  'at_time_upper_transition_probability_lower':tstar**6/F(256),
  'positive_transfer_statement':'0<abs(t)<=time_upper implies Pr(target bare composite)>=t^6/256',
  'time_is_original_H_clock':True,'prepared_vector_is_original_selected_density':False,
  'physical_Gauss_support_and_native_occupations_preserved':True,
  'invariant_eight_state_compression_assumed':False}


@lru_cache(maxsize=1)
def cube_control():
 sites=tuple(product(range(2),repeat=3));A=s.zeros(8);W=s.zeros(8)
 for i,x in enumerate(sites):
  for j,y in enumerate(sites):
   A[i,j]=int(sum(a!=b for a,b in zip(x,y))==1)
   W[i,j]=(-1)**sum(a*b for a,b in zip(x,y))
 if W.T*W!=8*s.eye(8):raise ArithmeticError('complete cube character basis changed')
 D=W.T*A*W/8
 if D!=s.diag(*[3-2*sum(k) for k in sites]):raise ArithmeticError('cube adjacency spectrum changed')
 rows=[]
 for label,delta,_ in sector.SECTORS[:2]:
  control=moment_control(delta);mean=control['centered_diagonal_third_moment']
  rows.append({'sector_label':label,'Walsh_rows':[{'character':k,'graph_adjacency_eigenvalue':3-2*sum(k),
    'exact_centered_third_energy_moment':mean-F(3,4)*(3-2*sum(k))} for k in sites],
    'same_first_and_second_moments_on_all_rows':True,
    'largest_third_moment_difference':F(9,2)})
 return {'sites':sites,'nearest_neighbor_adjacency':[[int(x) for x in row] for row in A.tolist()],
  'orthogonal_Walsh_numerator':[[int(x) for x in row] for row in W.tolist()],
  'normalization_squared':8,'sector_rows':rows,'Walsh_label_is_conserved_physical_momentum':False}


@lru_cache(maxsize=1)
def multiplet_transport_control():
 """Exact complete local gauge-kernel transport, including both cube spins."""
 rows=[]
 for label,delta,_ in sector.SECTORS[:2]:
  for sign in (1,-1):
   d=tuple(sign*x for x in delta);local=spin.local_multiplet(d);B=local['kernel_basis']
   if B.cols!=6:raise ArithmeticError('the physical odd multiplet lost its six states')
   vectors=[{mask:B[j,i] for j,mask in enumerate(local['masks']) if B[j,i]} for i in range(B.cols)]
   matrices=[]
   for axis in range(3):
    M=s.Matrix(6,6,lambda j,i:cubic_between_vectors(vectors[i],vectors[j],axis))
    if M!=-s.eye(6)*s.Rational(3,4):raise ArithmeticError('full physical spin multiplet transport is not scalar')
    matrices.append([[F(v) for v in row] for row in M.tolist()])
   rows.append({'sector_label':label,'sign':sign,'native_occupation_increment':d,
    'local_gauge_kernel_dimension':B.cols,'three_exact_transport_matrices':matrices,
    'source_spin_Casimir_multiplicities':{str(k):int(v) for k,v in local['casimir_multiplicities'].items()},
    'both_G1_and_H_source_cube_subspaces_transported':True})
 return {'signed_multiplets':rows,'all_local_gauge_kernel_columns_retained':True,
  'spin_chart_is_Lorentz_little_group':False}


def bind_member(member):
 sector.branch._member_zero(member)
 bound=sector.bind_member(member);sector.verify_bound(member,bound)
 rows=[]
 for label,delta,_ in sector.SECTORS[:2]:
  for sign in (1,-1):
   rows.append({'sector_label':label,'sign':sign,**moment_control(tuple(sign*x for x in delta))})
 result=contract.serial({'member_id':member['member_id'],**member['ids'],
  'sector_threshold_receipt_sha256':contract.digest(bound),'transport_rows':rows,
  'native_mask_to_source_lexicographic_row_permutation':native_chart(),
  'exact_cube_moment_control':cube_control(),
  'full_local_spin_transport_control':multiplet_transport_control(),
  'preparation':'paired-memory physical-sector preparation with bare localized or normalized Walsh target',
  'complete_detector':'target physical rank-one projector and its full complementary outcome',
  'original_reference_retained_in_idle_memory':True,'compensating_opposite_sector_retained':True,
  'full_source_action_and_clock_unchanged':True,'selected_source_marginal_changes_by_explicit_preparation':True,
  'link_Haar_integration_exact_by_local_Gauss_invariance':True,
  'cutoff_cubic_identity_uses_first_action_admission':True,
  'isolated_dispersion_band_proved':False,'relativistic_spin_or_rest_mass_proved':False,
  'physical_mass_table_GeV':[]})
 return {'schema':'oph.primary_composite_transport.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
 if value!=bind_member(member):raise ValueError('composite transport source or result changed')


def certificate():
 paths=provenance.definition_inputs((*SOURCE_INPUTS,*sector.SOURCE_INPUTS))
 pins=provenance.definition_hashes(paths)
 result={'schema':'oph.primary_composite_transport_certificate.v1',
  'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
 if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during cubic transport replay')
 return result
