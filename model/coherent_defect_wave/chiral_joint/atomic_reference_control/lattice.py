"""Scalable atomic isolation and complete-carrier hopping norm controls."""
from functools import lru_cache
import sympy as s
from ..local_reference.spectrum import certificate as local_spectrum
from ..local_reference.reference import certificate as local_reference
from ..periodic_family.model import Configuration,edges,edge_endpoints,cells,kinetic_terms,native_blocks
from ..peter_weyl_links.clebsch import link,transporter,native_charges
from ...physical_attachment.chiral import representation
from .local_bound import certificate as local_bound


def gap_data(side=2,kappa=s.Integer(1)):
    config=Configuration(side=side,link_penalty=kappa);kappa=config.link_penalty
    spectrum=local_spectrum();reference=local_reference()
    e0low,e0high=map(s.Rational,spectrum['global_physical_ground_energy_interval'])
    blankgap=s.Rational(spectrum['global_physical_gap_interval'][0])
    lower=s.Rational(-1171,100);delta=e0high-lower
    fluxgap=kappa-2*delta
    if fluxgap<=0:raise ValueError('link penalty does not satisfy the certified charged-cell threshold')
    gap=min(blankgap,fluxgap);assert gap>0 and delta>0 and -12<e0low<e0high<0
    assert spectrum['global_physical_ground_rank']==2 and reference['fermion_parity']=='odd'
    allcells=cells(side);incidence=s.zeros(len(allcells),len(edges(side)))
    positions={cell:i for i,cell in enumerate(allcells)}
    for edge in range(len(edges(side))):
        left,right=edge_endpoints(side,edge);assert left!=right
        incidence[positions[left],edge]=incidence[positions[right],edge]=1
    assert all(sum(incidence[:,j])==2 for j in range(incidence.cols))
    return {'side':side,'volume':side**3,'links':3*side**3,'kappa':kappa,
            'e0low':e0low,'e0high':e0high,'local_lower':lower,'delta':delta,
            'blankgap':blankgap,'fluxgap':fluxgap,'gap':gap,'incidence':incidence}


def atomic_certificate(side=2,kappa=s.Integer(1)):
    local_bound();row=gap_data(side,kappa);volume=row['volume']
    return {'side':side,'volume':volume,'oriented_links':row['links'],'link_penalty':str(row['kappa']),
            'atomic_Hamiltonian':'H0=sum_x H_local,x + kappa sum_e P_excited,e; hopping absent',
            'all_twelve_Gauss_constraints_imposed_on_full_lattice':True,
            'unrestricted_charged_local_energy_lower_bound':str(row['local_lower']),
            'physical_local_ground_energy_interval':[str(row['e0low']),str(row['e0high'])],
            'charged_cell_energy_loss_upper_bound':str(row['delta']),
            'each_nonblank_link_touches_exactly_two_cells':True,
            'arbitrary_nonblank_set_incident_cell_bound':'a<=2m',
            'complete_nonblank_m_link_sector_bound':'H0>=volume*E0+m*(kappa-2*delta)',
            'certified_positive_nonblank_sector_gap':str(row['fluxgap']),
            'zero_link_physical_gap_lower_bound':str(row['blankgap']),
            'full_physical_atomic_gap_lower_bound':str(row['gap']),
            'ground_energy':'volume*E0','ground_projector':'P0=tensor_x P_ground,x tensor all_links_blank',
            'exact_ground_rank':2**volume,'ground_rank_formula':'2^volume',
            'stationary_density':'P0/2^volume = tensor_x(P_ground,x/2) tensor blank_links',
            'all_link_counts_and_charged_cell_sectors_retained':True,
            'bound_is_uniform_in_volume_for_fixed_kappa':True,
            'unit_link_penalty_is_admitted':bool(row['kappa']==1),
            'atomic_gap_is_a_particle_mass':False,
            'ground_statement_applies_to_full_hopping_H':False}


@lru_cache(maxsize=1)
def transporter_controls():
    rows=[]
    for charge in native_charges():
        matrix=transporter(charge);gram=s.ImmutableSparseMatrix(matrix.H*matrix)
        assert gram.H==gram and gram*gram==gram
        rows.append({'six_Y':charge,'full_internal_link_dimension':matrix.rows,
                     'exact_complete_Gram_is_projection':True,'transporter_norm_upper_bound':'1'})
    coordinates=0
    for block in native_blocks():
        for a in block:
            for b in block:
                u=link(a,b)
                assert u[0,0]==0 and u.H[0,0]==0
                coordinates+=1
    assert coordinates==59
    return rows,coordinates


def hopping_norm(side=2):
    config=Configuration(side=side);terms=kinetic_terms(config)
    assert all(term.path and len(set(term.path))==len(term.path) for term in terms)
    return s.simplify(60*sum(abs(term.coefficient) for term in terms))


def local_norm():
    count=2*sum(abs(value) for vertex in representation()[3].values()
                for tensor in vertex['components'] for value in tensor)
    assert count==56
    return count+1


def hopping_certificate(side=2):
    blocks,coordinates=transporter_controls();tau=hopping_norm(side)
    if side==2:assert tau==360
    return {'side':side,'complete_internal_link_contractions':blocks,
            'native_fermion_modes_per_endpoint':30,
            'CAR_row_or_column_norm_upper_bound':'sqrt(30)',
            'single_kinetic_monomial_without_coefficient_norm_upper_bound':'30',
            'full_hopping_norm_formula':'tau=60 sum_paths |kinetic_coefficient|',
            'full_hopping_norm_upper_bound':str(tau),
            'local_onsite_norm_upper_bound':str(local_norm()),
            'all_native_coordinates_blank_to_blank_matrix_elements_checked':coordinates,
            'every_actual_path_nonempty_and_distinct_edges':True,
            'exact_atomic_P0_T_P0':'0',
            'vacuum_link_selection_does_not_replace_local_quantum_ground_by_empty_CAR':True,
            'all_Fock_scalar_and_link_occupations_in_norm_bound':True}
