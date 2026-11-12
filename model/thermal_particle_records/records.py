"""Full number-projector law, retaining colour coherence and thermal complements."""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import logsumexp

from model.particle_thermal_history import thermodynamics as thermal
from model.particle_common_contract.provenance import import_closure

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
TEMPERATURES=(.1,1.,10.,100.)


def authenticate(path):
    value=json.loads(path.read_text());pins=value.get('source_bindings',value.get('input_sha256'))
    if not pins:raise ValueError('source bindings missing')
    for name,digest in pins.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('changed source: '+name)
    return value


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update({'model/particle_joint_mass_identity/results.json','model/particle_joint_mass_identity/PROOF.md',
        'model/particle_joint_mass_identity/joint.py','model/particle_thermal_history/results.json',
        'model/particle_thermal_history/PROOF.md','model/particle_thermal_history/thermodynamics.py',
        'model/conditional_particle_spectrum/results.json','model/conditional_particle_spectrum/PROOF.md'})
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(paths))}


@lru_cache(maxsize=1)
def inputs():
    identity=authenticate(ROOT/'model/particle_joint_mass_identity/results.json')
    gibbs=authenticate(ROOT/'model/particle_thermal_history/results.json')
    common='model/conditional_particle_spectrum/results.json';digest=hashlib.sha256((ROOT/common).read_bytes()).hexdigest()
    if identity['source_bindings'].get(common)!=digest or gibbs['source_bindings'].get(common)!=digest:
        raise ValueError('the group and Gibbs producers do not bind the same quadratic source')
    groups=identity['complete_calibrated_joint_instrument']['joint_groups']
    source=authenticate(ROOT/common)
    for sector,count in [('fermion',180),('boson',56)]:
        indices=[i for g in groups if g['source_sector']==sector for i in g['mode_indices']]
        if sorted(indices)!=list(range(count)):raise ArithmeticError('joint groups do not partition all modes')
    if len(groups)!=112 or any(len(g['candidate_species'])!=1 for g in groups):
        raise ArithmeticError('an exact joint group has no unique declared comparison label')
    return identity,gibbs,source,groups


def maximum_abs(matrix):
    return float(max(abs(matrix.data),default=0.)) if sparse.issparse(matrix) else float(np.max(abs(matrix)))


@lru_cache(maxsize=1)
def boson_rotation():
    """Actual circular-polarization isometry on all1653 capped Fock states."""
    modes=inputs()[2]['mode_inventory']['boson_modes'];blocks=defaultdict(list)
    for i,m in enumerate(modes):blocks[(m['species'],m['orientation'],tuple(m['momentum_GeV']),m['color_component'])].append(i)
    rotation=np.zeros((56,56),complex);helicity=[None]*56
    for indices in blocks.values():
        n=len(indices)
        if n==1:unitary=np.eye(1);spins=[0]
        elif n==2:unitary=np.array([[1,1],[-1j,1j]])/np.sqrt(2);spins=[-1,1]
        elif n==3:
            unitary=np.array([[1/np.sqrt(2),0,1/np.sqrt(2)],[-1j/np.sqrt(2),0,1j/np.sqrt(2)],[0,1,0]])
            spins=[-1,0,1]
        else:raise ArithmeticError('unknown physical polarization block')
        rotation[np.ix_(indices,indices)]=unitary
        if len({modes[i]['energy_GeV']['exact'] for i in indices})!=1:raise ArithmeticError('rotation mixes distinct exact energies')
        if len({modes[i]['electric_charge'] for i in indices})!=1:raise ArithmeticError('rotation changes electric charge')
        for i,spin in zip(indices,spins):helicity[i]=2*spin*modes[i]['momentum_GeV'][2]
    for group in inputs()[3]:
        if group['source_sector']=='boson' and any(helicity[i]!=group['helicity_twice'] for i in group['mode_indices']):
            raise ArithmeticError('source circular frame differs from joint helicity groups')
    words=[()]+[(i,) for i in range(56)]+list(combinations_with_replacement(range(56),2));index={word:i for i,word in enumerate(words)}
    entries={(0,0):1.+0j}
    for i in range(56):
        for a in np.flatnonzero(rotation[:,i]):entries[index[(int(a),)],index[(i,)]]=rotation[a,i]
    for i,j in combinations_with_replacement(range(56),2):
        targets=set()
        for a in np.flatnonzero(rotation[:,i]):
            for b in np.flatnonzero(rotation[:,j]):targets.add(tuple(sorted((int(a),int(b)))))
        for a,b in targets:
            value=(rotation[a,i]*rotation[b,j]+rotation[a,j]*rotation[b,i])/math.sqrt((1+int(i==j))*(1+int(a==b)))
            if abs(value)>0:entries[index[(a,b)],index[(i,j)]]=value
    rows,cols,values=zip(*[(a,b,v) for (a,b),v in entries.items()]);fock=sparse.csr_matrix((values,(rows,cols)),shape=(1653,1653))
    identity_error=maximum_abs(fock.getH()@fock-sparse.eye(1653))
    energy=np.array([sum(modes[i]['energy_GeV']['numeric'] for i in word) for word in words])
    charge=np.array([sum(float(F(modes[i]['electric_charge'])) for i in word) for word in words])
    energy_error=maximum_abs(sparse.diags(energy)@fock-fock@sparse.diags(energy))
    charge_error=maximum_abs(sparse.diags(charge)@fock-fock@sparse.diags(charge))
    if max(identity_error,energy_error,charge_error)>2e-12:raise ArithmeticError('complete boson Fock-frame intertwiner failed')
    return dict(matrix=fock,words=words,helicity=helicity,unitary_error=identity_error,
        energy_intertwiner_error=energy_error,charge_intertwiner_error=charge_error)


def gauge_block_control():
    """Independent colour blocks from the parent canonical mode inventory."""
    _,_,source,groups=inputs();checks=[]
    for sector in ('fermion','boson'):
        modes=source['mode_inventory'][sector+'_modes'];blocks=defaultdict(list)
        for i,m in enumerate(modes):
            if (sector=='fermion' and m['species'] in thermal.QUARKS) or (sector=='boson' and m['species']=='gluon'):
                spin=(m.get('spin_twice',m.get('helicity_twice')) if sector=='fermion' else boson_rotation()['helicity'][i])
                blocks[(m['species'],m['orientation'],tuple(m['momentum_GeV']),spin)].append(i)
        for key,indices in blocks.items():
            rank=3 if sector=='fermion' else 8
            if len(indices)!=rank or sorted(modes[i]['color_component'] for i in indices)!=list(range(rank)):
                raise ArithmeticError('canonical colour fiber is incomplete')
            for group in groups:
                if group['source_sector']!=sector:continue
                intersection=set(group['mode_indices'])&set(indices)
                if intersection and intersection!=set(indices):raise ArithmeticError('number observable measures a colour component')
            checks.append(dict(sector=sector,mode_indices=indices,color_rank=rank,whole_color_block_retained=True))
    if len(checks)!=52:raise ArithmeticError('full48 quark and4 gluon colour blocks required')
    return checks


def color_reference_control(rank):
    psi=np.zeros(rank*rank,complex)
    for i in range(rank):psi[i*rank+i]=1/math.sqrt(rank)
    rho=np.outer(psi,psi.conj());dephased=np.diag(np.diag(rho))
    distance=float(np.abs(np.linalg.eigvalsh(rho-dephased)).sum())
    exact=2*F(rank-1,rank)
    if abs(distance-float(exact))>1e-12:raise ArithmeticError('colour coherence negative control differs')
    return dict(rank=rank,group_count_one_Lueders_output_unchanged=True,
        arbitrary_reference_dimension_allowed=True,fine_color_dephasing_trace_distance_exact=str(exact),
        executed_trace_distance=distance)


def binomial_records(rank,energy,temperature):
    eps=energy/temperature;q=math.exp(-eps);normal=math.log1p(q)
    rows=[]
    for count in range(rank+1):
        log_probability=math.log(math.comb(rank,count))-count*eps-rank*normal
        value=math.exp(log_probability) if log_probability>-746 else 0.
        rows.append(dict(count=count,projector_rank=math.comb(rank,count),
            probability=value,log_probability=log_probability,positive_probability_underflowed=value==0.))
    return rows


def count_law(temperature):
    _,thermal_receipt,source,groups=inputs()
    control=thermal.finite_source_gibbs(temperature)
    stored=next(r for r in thermal_receipt['finite_source_Gibbs_controls'] if r['temperature_GeV']==temperature)
    if control!=stored:raise ArithmeticError('fresh finite Gibbs law differs from its frozen source')
    fgroups=[g for g in groups if g['source_sector']=='fermion'];bgroups=[g for g in groups if g['source_sector']=='boson']
    fmode=source['mode_inventory']['fermion_modes'];bmode=source['mode_inventory']['boson_modes']
    mode_to_group={i:j for j,g in enumerate(bgroups) for i in g['mode_indices']};classes=defaultdict(list)
    words=boson_rotation()['words']
    for i,word in enumerate(words):
        count=Counter(mode_to_group[j] for j in word);key=tuple(sorted(count.items()));classes[key].append(i)
    if len(classes)!=435:raise ArithmeticError('cap-two group count partition is incomplete')
    b_records=[];b_means=[0.]*28
    for key,indices in sorted(classes.items()):
        logs=[control['all_boson_occupation_records'][i]['log_probability'] for i in indices]
        log_probability=float(logsumexp(logs));probability=math.exp(log_probability) if log_probability>-746 else 0.
        rank=math.prod(math.comb(bgroups[g]['projector_rank']+n-1,n) for g,n in key)
        if rank!=len(indices):raise ArithmeticError('symmetric-Fock group rank is incorrect')
        b_records.append(dict(group_counts=[dict(joint_index=bgroups[g]['joint_index'],count=n) for g,n in key],
            original_occupation_indices=indices,projector_rank=rank,probability=probability,
            log_probability=log_probability,positive_probability_underflowed=probability==0.))
        for g,n in key:b_means[g]+=probability*n
    f_records=[];group_means=[]
    for group in fgroups:
        rank=group['projector_rank'];energy=fmode[group['mode_indices'][0]]['energy_GeV']['numeric']
        rows=binomial_records(rank,energy,temperature);mean=sum(r['count']*r['probability'] for r in rows)
        expected=sum(control['complete_factorized_fermion_occupations'][i] for i in group['mode_indices'])
        if abs(mean-expected)>1e-12:raise ArithmeticError('fermion binomial group law differs from complete mode Gibbs law')
        f_records.append(dict(joint_index=group['joint_index'],rank=rank,records=rows))
        group_means.append(dict(joint_index=group['joint_index'],source_sector='fermion',
            species=group['candidate_species'][0],mean_number=mean,
            energy_GeV=energy,charge=group['charge'],momentum_z_GeV=group['momentum_z_GeV']))
    for group,mean in zip(bgroups,b_means):
        group_means.append(dict(joint_index=group['joint_index'],source_sector='boson',
            species=group['candidate_species'][0],mean_number=mean,
            energy_GeV=bmode[group['mode_indices'][0]]['energy_GeV']['numeric'],
            charge=group['charge'],momentum_z_GeV=group['momentum_z_GeV']))
    mean_total=sum(g['mean_number'] for g in group_means);by_species=defaultdict(float)
    for row in group_means:
        row['fraction_of_expected_counts']=row['mean_number']/mean_total;by_species[row['species']]+=row['mean_number']
    energy=sum(g['mean_number']*g['energy_GeV'] for g in group_means)
    charge=sum(g['mean_number']*float(F(g['charge'])) for g in group_means)
    momentum=sum(g['mean_number']*g['momentum_z_GeV'] for g in group_means)
    raw_energy=control['fermion_mean_energy_GeV']+control['boson_mean_energy_GeV']
    if abs(energy-raw_energy)>2e-10 or max(abs(charge),abs(momentum))>2e-12:
        raise ArithmeticError('independent additive energy/charge/momentum control differs')
    probabilities=np.array([r['probability'] for r in control['all_boson_occupation_records']])
    rotation=boson_rotation()['matrix'];density=sparse.diags(probabilities)
    gibbs_rotation_error=maximum_abs(rotation.getH()@density@rotation-density)
    if gibbs_rotation_error>2e-14:raise ArithmeticError('Gibbs state not invariant under actual polarization rotation')
    qf=sum(math.exp(-m['energy_GeV']['numeric']/temperature) for m in fmode)
    qb=sum(math.exp(-m['energy_GeV']['numeric']/temperature) for m in bmode)
    log_vacuum=-control['fermion_log_partition']-math.log(control['boson_partition'])
    log_at_most_one=log_vacuum+math.log1p(qf+qb)
    vacuum=math.exp(log_vacuum);at_most_one=math.exp(log_at_most_one)
    initial_complement=1-at_most_one
    return dict(temperature_GeV=temperature,fermion_group_factor_records=f_records,
        full_fermion_joint_law='product of the84 listed binomial factor laws, with every count0..rank retained',
        full_fermion_joint_record_count_exact=str(math.prod(g['projector_rank']+1 for g in fgroups)),
        boson_joint_group_count_records=b_records,boson_probability_sum=sum(r['probability'] for r in b_records),
        original_boson_occupation_records=control['all_boson_occupation_records'],
        group_expected_counts=group_means,species_expected_counts=[dict(species=k,mean_number=v,
            fraction_of_expected_counts=v/mean_total) for k,v in sorted(by_species.items())],
        full_expected_particle_number=mean_total,
        fraction_convention='ratio of expected counts, measured by pooling all trials including empty trials; not expectation of the per-trial ratio',
        pre_probe_one_excitation_detector=dict(vacuum_probability=vacuum,exact_one_probability=at_most_one-vacuum,
            at_most_one_probability=at_most_one,many_excitation_complement_probability=initial_complement,
            log_vacuum_probability=log_vacuum,log_at_most_one_probability=log_at_most_one,
            identity='P(Ntotal<=1)=(1+sum_F exp(-E/T)+sum_B exp(-E/T))/(Z_F*Z_B_cap2)',
            evaluated_before_any_probe=True,claimed_to_be_the570_post_probe_record_law=False),
        independent_controls=dict(group_energy_equals_full_Gibbs_energy_error=abs(energy-raw_energy),
            total_charge_expectation=charge,total_momentum_expectation=momentum,
            full1653_polarization_rotation_Gibbs_invariance_error=gibbs_rotation_error,
            group_boson_mean_error=abs(sum(b_means)-control['boson_mean_number'])),
        measured_continuum_gas_or_cosmological_abundances=False)


def certificate():
    pins=bindings();identity,_,source,groups=inputs();rotation=boson_rotation();color=gauge_block_control()
    rows=[count_law(t) for t in TEMPERATURES]
    fgroups=[g for g in groups if g['source_sector']=='fermion'];fg_outcomes=math.prod(g['projector_rank']+1 for g in fgroups)
    result=dict(schema='oph.thermal_particle_records.v1',all_checks_pass=True,source_bindings=pins,
        source_Hilbert_dimension_exact=identity['complete_calibrated_joint_instrument']['source_Hilbert_dimension_exact'],
        joint_groups=groups,complete_number_projectors=dict(
            fermion_group_count=84,boson_group_count=28,fermion_count_outcomes_exact=str(fg_outcomes),
            boson_nonzero_rank_count_outcomes=435,total_nonzero_rank_joint_count_outcomes_exact=str(435*fg_outcomes),
            boson_zero_rank_count_vectors_symbolically_retained=3**28-435,
            zero_rank_boson_outcome_rule='every vector in{0,1,2}^28 with totalcount>2 has the exact zero effect and exact zero probability',
            total_factorized_joint_outcome_count_including_zero_effects_exact=str(3**28*fg_outcomes),
            total_projector_rank_exact=str((2**180)*1653),
            all_input_identity='product of the complete spectral partitions of commuting group number operators N_g; tensor identity on the474-dimensional original reference and any extra reference',
            number_operator='sum beta_j^dagger beta_j over the actual joint eigenspace; beta are the parent canonical quasiparticle fields',
            exact_projector_formula='P_g,n=product_(ell!=n)(N_g-ell I)/(n-ell), ell=0..r_g for fermions and0..2 for capped bosons',
            physical_or_fine_color_components_measured=False,
            many_body_total_color_Casimir_resolved=False,all_particle_occupations_retained=True),
        exact_color_block_partition=color,
        complete_boson_Fock_rotation=dict(dimension=1653,unitarity_error=rotation['unitary_error'],
            energy_intertwiner_error=rotation['energy_intertwiner_error'],charge_intertwiner_error=rotation['charge_intertwiner_error'],
            source_mode_indices_are_circular_polarization_columns=True,
            raw_occupation_words_reused_only_after_proving_Gibbs_rotation_invariance=True),
        arbitrary_reference_color_controls=[color_reference_control(r) for r in (3,8)],
        finite_temperature_count_laws=rows,
        scope=dict(same_calibrated_quadratic_source_and_original_global_boson_cap=True,
            full_group_number_readout_is_declared_service=True,native_pulse_compilation_or_hardware_accuracy=False,
            complete_all_input_reference_safe_Lueders_instrument=True,
            color_coherences_within_count_blocks_preserved=True,
            microscopic_colour_Gauss_law_or_confinement_derived=False,
            calibrated_parameters_are_source_predictions=False,
            native_temperature_or_Universe_abundances_identified=False,
            continuum_ideal_gas_substituted_for_finite_source=False,
            numerical_probabilities_have_certified_roundoff_intervals=False))
    if pins!=bindings():raise ArithmeticError('thermal count source changed during execution')
    return result


def main():
    p=argparse.ArgumentParser();m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=p.parse_args()
    result=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=result:raise ArithmeticError('finite thermal particle record receipt differs')
    print('PASS:complete finite thermal group-count law and pre-probe complement')


if __name__=='__main__':main()
