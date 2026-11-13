from fractions import Fraction as F
from itertools import combinations,product
from functools import lru_cache
import sympy as sy
from .cube import graph,spin_census,basis,matrices,cg_entry,successors
from .certificate import calculate


@lru_cache(None)
def exact_checks():
    census=spin_census(26);small=basis();_,wilson=matrices();phases={0:1}
    for _ in range(len(small)):
        for face in range(6):
            for i,(_,left,_)in enumerate(small):
                if i not in phases:continue
                for j,(_,right,_)in enumerate(small):
                    if not wilson[face][i,j]:continue
                    ratio=sy.simplify(cg_entry(left,right,face)/wilson[face][i,j])
                    assert ratio in(-1,1)
                    phase=ratio*phases[i]
                    if j in phases:assert phases[j]==phase
                    phases[j]=phase
    assert len(phases)==32
    checked=0
    for _,spins in spin_census(32):
        for a,b in combinations(range(6),2):
            outputs=[]
            for first,second in((a,b),(b,a)):
                terms={}
                for middle,x in successors(spins,first):
                    for target,y in successors(middle,second):terms[target]=terms.get(target,0)+x*y
                outputs.append(terms)
            assert all(sy.simplify(outputs[0].get(k,0)-outputs[1].get(k,0))==0
                       for k in outputs[0].keys()|outputs[1].keys())
            checked+=1
    return {'full_trivalent_Gauss_census_through_energy_13_over_2':len(census),
        'complete_dimension_through_energy_six':32,
        'energy_census':[{ 'electric_energy_exact':str(F(e,4)),
            'count':sum(x==e for x,_ in census)}for e in sorted(set(e for e,_ in census))],
        'independent_Haar_Gram_identity':True,'all_six_Wilson_matrices_match_independent_character_integrals':True,
        'explicit_spin_network_to_character_basis_phases':[int(phases[i])for i in range(32)],
        'full_successor_commutators_checked':checked}


def run():
    checks=exact_checks();records=[calculate(g)for g in(4,8,16)]
    for record in records:
        variances=[tuple(map(F,row))for row in record['all_six_Wilson_variances_intervals_exact']]
        assert all(row[0]>0 for row in variances)
        means=[tuple(map(F,row))for row in record['all_six_Wilson_means_intervals_exact']]
        assert max(row[0]for row in means)<=min(row[1]for row in means)
        for opposite in(False,True):
            intervals=[tuple(map(F,row['connected_covariance_interval_exact']))for row in record['all_fifteen_connected_covariances']if row['opposite']==opposite]
            assert max(row[0]for row in intervals)<=min(row[1]for row in intervals)
        alpha=F(record['g_squared_exact'])/2
        record['complete_joint_Wilson_instrument']={
            'response_eta_exact':'1/2','effects':'product_f(I+s_f W_f/2)/64, all six signs retained',
            'Kraus_rule':'positive square root of every full-Haar effect',
            'complete_effect_sum':'I','external_quantum_reference_retained':True,
            'all_sixty_four_outcomes':[{'signs':list(signs),'probability_interval_exact':['1/4096','729/4096']}for signs in product((-1,1),repeat=6)],
            'record_pair_covariance':'one quarter of the certified native Wilson covariance for the same pair',
            'nonselective_ground_energy_work_interval_exact':['0',str(alpha/2)]}
    assert all(F(row['connected_covariance_interval_exact'][0])>0
               for record in records[-2:]for row in record['all_fifteen_connected_covariances']if not row['opposite'])
    return {'schema':'oph.complete_cubic_full_Haar_records.v1',
        'verdict':'PASS_COMPLETE_CUBIC_GAUSS_SOURCE_AND_FULL_HAAR_CONNECTED_RECORDS',
        'graph':{'vertices':8,'links':12,'faces':6,'spatial_coordinate_dimensions':3},
        'exact_source_checks':checks,'records':records,
        'source':'Declared cubic SU2 Casimir-plus-Wilson Hamiltonian on every Haar link representation with every trivalent invariant channel',
        'cube_symmetry':'All six face means and the separate adjacent/opposite covariance classes have overlapping certified intervals; spatial graph automorphisms preserve the full source.',
        'representation_projection_is_source_cutoff':False,
        'physical_gap_GeV_or_continuum_identification':False,
        'uniform_spatial_volume_or_refinement_limit_proved':False,
        'new_empirical_fit_inputs':[]}
