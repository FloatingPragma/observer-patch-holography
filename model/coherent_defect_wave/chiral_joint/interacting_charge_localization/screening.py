"""An exact screened remote-pair state in the actual native carrier."""
import sympy as s
from .geometry import graph,component
from ..periodic_family.model import Configuration,kinetic_terms,cells
from ..peter_weyl_links.clebsch import link
from ..interacting_charged_response.model import FAMILIES
from ..atomic_momentum_response.model import reference_data
from ...chiral_source.model import PAULI
from ...physical_attachment.chiral import representation


def certificate():
    edges,distance=graph();edge=next(e for e in edges if e[1:]==(3,7))
    term=next(t for t in kinetic_terms(Configuration()) if t.path==(edge[0],))
    assert edge==(9,3,7) and term.axis==0 and abs(term.coefficient)==s.Rational(1,4)
    assert 3 not in component((edge,)) and distance[3]==2
    reference,_=reference_data();records=[]
    for family in ('ec','L'):
        indices=FAMILIES[family];d=len(indices)
        assert reference['native_irrep_ground_occupations'][family]['exact_ground_mean_occupation']==str(d)
        # The spin-mixed, gauge-singlet actual local ground has one-body density
        # I/2 on each native irrep. This contracts the actual vacuum link Gram.
        norm=s.Integer(0)
        for i in range(2):
            for j in range(2):
                for a in indices:
                    for b in indices:
                        column=link(a,b).H[:,0]
                        norm+=s.conjugate(PAULI[term.axis][i,j])*PAULI[term.axis][i,j]*(column.H*column)[0]/64
        assert norm==s.Rational(d,32)
        q=int(6*representation()[2][-1][indices[0],indices[0]])
        records.append({'family':family,'native_internal_dimension':d,
            'remote_native_hop_squared_norm_on_actual_product_ground_density':str(norm),
            'combined_origin_probe_and_remote_hop_squared_norm':str(norm/2),
            'source_six_Y_by_changed_cell':{'0':-q,'3':-q,'7':q},
            'retained_origin_reference_six_Y':q,'total_six_Y':'0',
            'excited_link_count':1,'remote_hole_distance_from_origin':distance[3],
            'flux_path_from_origin_to_remote_hole_exists':False,
            'single_unscreened_charge_projectors_all_annihilate_state':True,
            'joint_Gauss_constraints_preserved':True,'full_H_hopping_expectation':'0',
            'atomic_energy_above_8E0_upper':'3(57-E0)+1',
            'total_energy_above_actual_ground_upper':'mu+3(57-E0)+1+v',
            'state_asserted_an_interacting_ground_or_asymptotic_scattering_state':False})
    return {'physical_input':'tensor of the eight actual local ground densities P_g/2, blank links, fixed reference vacuum; local spin purifications retained',
        'operator':'B_origin^r followed by the adjoint branch of actual edge 9, moving matter from cell 3 to cell 7',
        'changed_cells':{'origin':list(cells(2)[0]),'remote_hole':list(cells(2)[3]),'remote_particle':list(cells(2)[7])},
        'actual_edge':list(edge),'source_clock_and_energy_operator':'unchanged H_lambda+4 N_ref',
        'counterexample_target':'a charged matter cell at distance ell need not carry ell excited links from the static compensator when other charged matter is allowed',
        'records':records,'thermodynamic_string_breaking_claimed':False}
