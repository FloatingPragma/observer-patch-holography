"""Flavour freedom of the actual finite conditional tensor construction."""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix

from model.common_source import calibrated_matter as matter
from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
TYPES=('Q','u^c','d^c','L','e^c')


def require(condition,message):
    if not condition:raise ValueError(message)


def zero(matrix):
    return all(s.expand(x)==0 or s.simplify(x)==0 for x in matrix)


def record(matrix):
    a=s.ImmutableSparseMatrix(matrix)
    return dict(shape=list(a.shape),entries=[[i,j,str(v)]for(i,j),v in sorted(a.todok().items())])


def antihermitian_basis():
    result=[]
    for i in range(3):
        a=s.zeros(3);a[i,i]=s.I;result.append(a)
    for i in range(3):
        for j in range(i+1,3):
            a=s.zeros(3);a[i,j]=1;a[j,i]=-1;result.append(a)
            a=s.zeros(3);a[i,j]=s.I;a[j,i]=s.I;result.append(a)
    require(len(result)==9 and all(zero(a.H+a)for a in result),'u3 basis failed')
    return result


def real_coordinates(matrices,symmetric_last=False):
    values=[]
    for k,a in enumerate(matrices):
        values.extend(a[i,j] for i in range(3) for j in range(3)
                      if not(symmetric_last and k==len(matrices)-1 and j<i))
    return s.Matrix([s.expand_complex(v).as_real_imag()[part] for part in (0,1) for v in values])


def tangent(yu,yd,ye,mnu=None):
    columns=[];basis=antihermitian_basis()
    for group in range(5):
        for generator in basis:
            q,u,d,l,e=[generator if k==group else s.zeros(3)for k in range(5)]
            values=[-q*yu+yu*u,-q*yd+yd*d,-l*ye+ye*e]
            if mnu is not None:values.append(l.T*mnu+mnu*l)
            columns.append(real_coordinates(values,mnu is not None))
    return s.Matrix.hstack(*columns)


def rotation(s12,s23,s13,phase=s.I):
    c12,c23,c13=[s.sqrt(1-x*x)for x in (s12,s23,s13)]
    r12=s.Matrix([[c12,s12,0],[-s12,c12,0],[0,0,1]])
    r23=s.Matrix([[1,0,0],[0,c23,s23],[0,-s23,c23]])
    r13=s.Matrix([[c13,0,s13*s.conjugate(phase)],[0,1,0],[-s13*phase,0,c13]])
    result=r23*r13*r12
    require(zero(result.H*result-s.eye(3)),'three-angle mixing is not unitary')
    return result


def example_matrices():
    v=rotation(s.Rational(3,5),s.Rational(5,13),s.Rational(8,17))
    yu=s.diag(1,2,4);yd=v*s.diag(2,3,5);ye=s.diag(1,3,7)
    mnu=s.Matrix([[2,1,s.I],[1,3,1],[s.I,1,5]])
    return yu,yd,ye,mnu,v


@lru_cache(maxsize=1)
def quotient():
    yu,yd,ye,mnu,_=example_matrices()
    charged=tangent(yu,yd,ye);weinberg=tangent(yu,yd,ye,mnu)
    ranks=[DomainMatrix.from_Matrix(a).rank()for a in (charged,weinberg)]
    require(ranks==[41,44],'generic orbit rank differs')
    u=rotation(s.Rational(5,13),s.Rational(8,17),s.Rational(7,25))
    rank_two=tangent(yu,yd,ye,u.conjugate()*s.diag(0,1,3)*u.H)
    rank_two_orbit=DomainMatrix.from_Matrix(rank_two).rank()
    require(rank_two_orbit==44,'generic rank-two orbit differs')
    baryon=s.zeros(45,1)
    for group in (0,1,2):
        for family in range(3):baryon[9*group+family]=1
    require(zero(charged*baryon)and zero(weinberg*baryon),'baryon stabilizer failed')
    leptons=[]
    for family in range(3):
        v=s.zeros(45,1);v[27+family]=v[36+family]=1;leptons.append(v)
        require(zero(charged*v),'charged-lepton stabilizer failed')
    kernel=s.Matrix.hstack(baryon,*leptons)
    require(kernel.rank()==4 and len(charged.nullspace())==4,'charged orbit nullspace differs')
    require(len(weinberg.nullspace())==1,'Majorana orbit nullspace differs')
    return dict(charged_raw_real_dimension=54,family_unitary_dimension=45,
        charged_exact_orbit_rank=ranks[0],charged_stabilizer_dimension=4,
        charged_invariant_parameters=13,
        charged_parameter_dictionary=dict(charged_singular_values=9,quark_mixing_angles=3,quark_CP_phase=1),
        Weinberg_additional_raw_real_dimension=12,Weinberg_exact_orbit_rank=ranks[1],
        Weinberg_stabilizer_dimension=1,Weinberg_generic_invariant_parameters=22,
        Weinberg_parameter_dictionary=dict(charged_singular_values=9,quark_mixing_angles=3,quark_CP_phase=1,
            neutrino_singular_values=3,lepton_mixing_angles=3,lepton_Dirac_CP_phase=1,lepton_Majorana_phases=2),
        rank_two_Weinberg_raw_real_dimension=64,rank_two_Weinberg_exact_orbit_rank=rank_two_orbit,
        rank_two_Weinberg_invariant_parameters=20,
        charged_tangent=record(charged),Weinberg_tangent=record(weinberg),
        charged_nullspace_basis=record(kernel),Weinberg_nullspace_basis=record(baryon),
        generic_rank_witness={k:record(v)for k,v in zip(['Yu','Yd','Ye','Mnu'],(yu,yd,ye,mnu))},
        domain='finite canonical family-basis equivalence with all preparations and detectors transported; no quantum chiral measure or theta-angle quotient asserted')


def components(yukawas):
    out={}
    for name,y in yukawas.items():
        rows=[]
        for internal in matter._oriented_components(name):
            directed=s.kronecker_product(internal,y.conjugate())
            rows.append(s.ImmutableSparseMatrix(directed+directed.T))
        out[name]=tuple(rows)
    return out


def kernel(yukawas):
    return sum((row[1]for row in components(yukawas).values()),s.zeros(45))


@lru_cache(maxsize=1)
def native_constraints():
    labels=action.ordering()['all_left_components']
    currents=action.generators();basis=antihermitian_basis();family_checks=0
    for kind in TYPES:
        for a in basis:
            f=s.diag(*[a if row['internal_type']==kind else s.zeros(3)for row in labels])
            for g in currents:
                require(zero(g*f-f*g),'family basis does not commute with an actual native current')
                family_checks+=1
    internal_checks=0
    for name in matter.CHANNELS:
        cs=matter._oriented_components(name)
        for primitive,g in zip(action._representation()[1],action.generators(families=False)):
            h=primitive[3:,3:]
            if action._representation()[4][name]['conjugate_Higgs']:h=-h.T
            for a in range(2):
                require(zero(g.T*cs[a]+cs[a]*g+sum((h[b,a]*cs[b]for b in range(2)),s.zeros(15))),
                        'oriented internal Higgs identity failed')
                internal_checks+=1
    # Distinct family coefficients have independent complete K45 images.
    columns=[];neutral_checks=0
    for name in matter.CHANNELS:
        for i in range(3):
            for j in range(3):
                e=s.zeros(3);e[i,j]=1
                y={k:e if k==name else s.zeros(3)for k in matter.CHANNELS}
                k=kernel(y)
                require(k==k.T and zero(action.electric_charge().T*k+k*action.electric_charge()),
                        'complete native family coefficient lost charge or transpose symmetry')
                require(zero(s.kronecker_product(action.SPIN_EPSILON,k).T+s.kronecker_product(action.SPIN_EPSILON,k)),
                        'full90-mode CAR creation coefficient lost antisymmetry')
                neutral_checks+=1;columns.append(s.Matrix(list(k)))
    gram=s.Matrix([[ (a.conjugate().T*b)[0] for b in columns]for a in columns])
    require(gram.is_diagonal() and all(x>0 for x in gram.diagonal()),'family coefficient map is not injective')
    yu,yd,ye,_,_=example_matrices();cs=components(dict(up=yu,down=yd,charged_lepton=ye));full_checks=0
    for name in matter.CHANNELS:
        for primitive,g in zip(action._representation()[1],currents):
            h=primitive[3:,3:]
            if action._representation()[4][name]['conjugate_Higgs']:h=-h.T
            for a in range(2):
                require(zero(g.T*cs[name][a]+cs[name][a]*g+sum((h[b,a]*cs[name][b]for b in range(2)),s.zeros(45))),
                        'full native complex-family Ward counterexample failed')
                full_checks+=1
    # The optional LH tensor has no family-coordinate restriction on a symmetric coefficient.
    h0,h1=s.symbols('H0 H1');lh=matter._LH_tensor();ell=lh*s.Matrix([h0,h1]);neutral_ward=0
    for primitive,g in zip(action._representation()[1],action.generators(families=False)):
        require(zero(g.T*ell+lh*primitive[3:,3:]*s.Matrix([h0,h1])),'native LH Ward identity failed')
        neutral_ward+=1
    return dict(family_generator_native_current_commutators=family_checks,
        oriented_internal_Higgs_Ward_identities=internal_checks,full45_complex_family_Ward_identities=full_checks,
        complete45_neutral_and90_CAR_coefficient_checks=neutral_checks,
        charged_family_embedding_complex_rank=27,charged_family_embedding_Gram=record(gram),
        extra_linear_family_constraints_from_native_gauge_Ward_equations=0,
        optional_LH_Ward_identities=neutral_ward,symmetric_Weinberg_complex_dimension=6,
        every_coefficient_in_the_declared_grammar_checked=True,
        statement='The actual native current/Higgs and canonical CAR identities admit every complex3by3 charged matrix and every complex symmetric Weinberg matrix; a selected source law may impose additional equations.')


def jarlskog(v):
    return s.simplify(s.im(v[0,0]*v[1,1]*s.conjugate(v[0,1]*v[1,0])))


def loop_readout(v):
    """Four complete probabilities measuring a rephasing-invariant current loop."""
    a=v[0,0]*v[1,1];b=v[0,1]*v[1,0]
    aa=s.simplify(s.conjugate(a)*a);bb=s.simplify(s.conjugate(b)*b)
    j=jarlskog(v)
    rows={'joint_Y_plus':s.simplify((aa+bb)/4-j/2),
          'joint_Y_minus':s.simplify((aa+bb)/4+j/2),
          'path0_failure':s.simplify((1-aa)/2),'path1_failure':s.simplify((1-bb)/2)}
    require(sum(rows.values())==1 and all(p>=0 for p in rows.values()),'complete current-loop probabilities failed')
    success=s.diag(a,b);fail0=s.diag(s.sqrt(1-aa),0);fail1=s.diag(0,s.sqrt(1-bb))
    require(zero(success.H*success+fail0.H*fail0+fail1.H*fail1-s.eye(2)),'loop instrument incomplete')
    require(s.simplify(rows['joint_Y_minus']-rows['joint_Y_plus']-j)==0,'loop inversion failed')
    # Preserve every output family pair of the literal V tensor V current
    # isometry, then both Y outcomes, rather than treating nontarget pairs
    # as an unrecorded filter failure.
    full=[];effect=s.zeros(2)
    for output0 in range(3):
        for output1 in range(3):
            a_rs=v[output0,0]*v[output1,1]
            b_rs=v[output0,1]*v[output1,0]
            for label,sign in [('plus',1),('minus',-1)]:
                kraus=s.Matrix([[a_rs,-sign*s.I*b_rs]])/s.sqrt(2)
                probability=s.simplify(abs((a_rs-sign*s.I*b_rs)/2)**2)
                effect+=kraus.H*kraus
                full.append(dict(output_family_pair=[output0,output1],Y=label,
                                 probability=str(probability),kraus=record(kraus)))
    require(zero(effect-s.eye(2)),'complete18-record current isometry failed')
    require(sum(s.sympify(row['probability'])for row in full)==1,'complete current records lost mass')
    target={row['Y']:s.sympify(row['probability'])for row in full if row['output_family_pair']==[0,1]}
    require(s.simplify(target['minus']-target['plus']-j)==0,'unconditional original-record contrast failed')
    return dict(probabilities={k:str(x)for k,x in rows.items()},J_exact=str(j),
        inverse='J=P(joint_Y_minus)-P(joint_Y_plus), using all attempted trials',
        all_failures_retained=True,normalization='no postselection',
        complete_native_current_records=full,complete_native_current_record_count=18,
        complete_path_effect_identity=True,
        four_record_filter_and_eighteen_record_current_are_separate_instruments=True,
        equality_or_coarse_graining_of_their_complete_reference_channels_claimed=False,
        native_current_protocol='path0 prepares families(0,1), path1 prepares(1,0) in two distinct colour-register copies; apply V tensor V, retain every output family pair and both pathY outcomes',
        native_current_inverse='J=P(output_pair01,Yminus)-P(output_pair01,Yplus), with all18 records retained',
        total_charge_transfer_same_on_both_paths=True,
        Kraus_operators={k:record(x)for k,x in [('success',success),('path0_failure',fail0),('path1_failure',fail1)]},
        terminal='PauliY projectors on the successful path qubit, both outcomes retained',
        apparatus='declared coherent two-path composition of four family-current amplitudes; not a physical laboratory weak process')


@lru_cache(maxsize=1)
def relation_and_readout():
    yu,yd,ye,_,v=example_matrices();hu=yu*yu.H;hd=yd*yd.H;comm=hu*hd-hd*hu
    j=jarlskog(v);require(j!=0 and not zero(comm),'generic alignment counterexample failed')
    loop=loop_readout(v);conjugate=loop_readout(v.conjugate())
    require(jarlskog(v.conjugate())==-j,'CP conjugation failed')
    require(all(s.simplify(abs(v[i,k])**2-abs(v.conjugate()[i,k])**2)==0 for i in range(3)for k in range(3)),
            'CP conjugate changed incoherent family probabilities')
    # Execute the actual native weak-current mass-frame block on all45 coordinates.
    labels=action.ordering()['all_left_components'];blocks=[s.eye(3)for _ in labels]
    for name in matter.CHANNELS:
        for(i,k),sign in matter._oriented_components(name)[1].todok().items():
            blocks[k]=sign*s.eye(3)
            if name=='down':blocks[i]=v
    frame=s.diag(*blocks);k=kernel(dict(up=yu,down=yd,charged_lepton=ye))
    require(zero(frame.H*frame-s.eye(45)),'full family mass frame not unitary')
    diagonal=frame.T*k*frame
    # Internal singlet signs are removed by frame; compare each explicit mass pair.
    for name,diag in [('up',yu),('down',s.diag(2,3,5)),('charged_lepton',ye)]:
        for(i,k),_ in matter._oriented_components(name)[1].todok().items():
            require(zero(diagonal[3*i:3*i+3,3*k:3*k+3]-diag),'native mass pair diagonalization failed')
    plus=matter.weak_raising_current();rotated=frame.H*plus*frame
    pairs=[(i,k)for(i,k),z in matter.weak_raising_current(families=False).todok().items()
           if labels[i]['internal_type']=='Q']
    require(len(pairs)==3 and all(zero(rotated[3*i:3*i+3,3*k:3*k+3]-v)for i,k in pairs),
            'actual full native weak current does not supply the loop amplitudes')
    native=action.yukawa_matrices();hnu=native['up']*native['up'].H;hnd=native['down']*native['down'].H
    require(zero(hnu*hnd-hnd*hnu),'selected common-feature trial lost alignment')
    return dict(tested_relation='[Yu Yu-dagger,Yd Yd-dagger]=0',
        native_common_feature_trial_relation_exact=True,
        relation_for_entire_native_gauge_admissible_family=False,
        generic_counterexample=dict(Yu=record(yu),Yd=record(yd),Ye=record(ye),CKM=record(v),
            commutator=record(comm),commutator_determinant=str(s.factor(comm.det())),
            J_exact=str(j),all_three_mixing_sines=['3/5','5/13','8/17'],CP_phase='pi/2',
            measured_masses_or_angles_used=False),
        source_trial='Actual action.yukawa_matrices gives a F3,b F3,J F3; aligned squared kernels imply zero quark J. This is a property of that declared trial, not of the native gauge identities.',
        full45_weak_current_color_blocks_checked=3,
        CP_conjugate_same_singular_values_and_modulus_squared=True,
        loop_readout=loop,CP_conjugate_loop_readout=conjugate,
        smallest_extra_source_equation='A source-predicted value of Im(V11 V22 conjugate(V12 V21)) fixes one remaining real CP direction once nonzero mixing angles are fixed.',
        smallest_extra_observable='One unconditional signed jointY contrast for output family pair01 of the complete18-record two-current instrument.',
        physical_family_or_particle_assignment_claimed=False)


def neutrino_phase_rank():
    u=rotation(s.Rational(5,13),s.Rational(8,17),s.Rational(7,25),s.I)
    d=s.diag(0,1,3);p=s.Rational(3,5)+s.I*s.Rational(4,5)
    mass=lambda v:s.simplify(v.conjugate()*d*v.H)
    m=mass(u);common=u*s.diag(1,p,p);changed=mass(common)
    require(zero((p*s.eye(3)).T*changed*(p*s.eye(3))-m),'rank-two common Majorana phase redundancy failed')
    charged=s.diag(1,3,7)
    require(zero((p*s.eye(3)).H*charged*(p*s.eye(3))-charged),'charged lepton phase compensation failed')
    relative=u*s.diag(1,1,p);other=mass(relative)
    first=s.simplify(m[0,0]*s.conjugate(m[0,0]));second=s.simplify(other[0,0]*s.conjugate(other[0,0]))
    require(first!=second,'remaining Majorana phase did not affect L-violating invariant')
    require(zero(u*d*d*u.H-relative*d*d*relative.H),'oscillation generator depends on Majorana phases')
    return dict(full_rank_Majorana_phases=2,rank_two_Majorana_phases=1,
        massless_column_phase_invisible_to_mass_kernel=True,
        common_phase_of_two_massive_columns_compensated_by_charged_lepton_basis=True,
        remaining_relative_phase_changes_Lepton_number_violating_invariant=dict(
            observable='abs(Mnu_ee)^2 in the charged-lepton mass basis',before=str(first),after=str(second)),
        oscillation_generator_U_m_squared_U_dagger_unchanged_by_Majorana_column_phases=True,
        independent_phase_sensitive_readout_needed_for_remaining_Majorana_direction=True,
        normal_ordering_is_not_source_derived=True)


def source_bindings():
    paths={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in {'.py','.md'}}
    paths.update(matter.source_bindings())
    paths.update(['model/conditional_standard_model/PROOF.md','model/particle_joint_reconstruction/CANDIDATE.md'])
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in sorted(import_closure(paths))}


def certificate():
    pins=source_bindings();charged=matter.parameters();neutral=matter.neutrino_parameters()
    result=dict(schema='oph.source_flavour_identifiability.v1',all_checks_pass=True,source_bindings=pins,
        exact_family_quotient=quotient(),native_source_constraints=native_constraints(),
        tested_relation_and_complete_readout=relation_and_readout(),neutrino_phase_redundancy=neutrino_phase_rank(),
        candidate_ancestry=dict(
            charged_yukawa_singular_values={k:dict(values=v,origin='empirical common-scale packet')for k,v in charged['matched_SM']['yukawa'].items()},
            quark_mixing=charged['CKM'],neutrino_oscillation_inputs=neutral['oscillation_inputs'],
            neutrino_architecture_choices=neutral['unmeasured_architecture_choices'],
            independent_empirical_flavour_coordinates_in_rank_two_packet=19,
            independent_unmeasured_Majorana_choice_in_rank_two_packet=1,
            zero_lightest_mass_and_ordering_are_architecture=True,
            two_printed_Majorana_phase_inputs_contain_one_redundant_common_phase_at_rank_two=True,
            additional_inputs_outside_this_flavour_count=['Higgs norm and units','gauge and scalar coefficients','theta/global-group choices','chiral measure and source/physical dictionary']),
        physical_premise_promotions=[],scope=dict(
            finite_family_basis_equivalence_proved=True,physical_gauge_redundancy_or_chiral_measure_proved=False,
            empirical_yukawas_or_neutrino_inputs_derived=False,source_selects_all_flavour_parameters=False,
            complete_declared_current_loop_records=True,physical_particle_spectrum_completed=False))
    require(source_bindings()==pins,'source changed during flavour calculation')
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=certificate();data=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:(HERE/'results.json').write_text(data)
    if args.check:require((HERE/'results.json').read_text()==data,'flavour receipt differs')
    print(json.dumps(dict(charged_parameters=13,Weinberg_parameters=22,rank_two_parameters=20,
                         J=result['tested_relation_and_complete_readout']['generic_counterexample']['J_exact'])))


if __name__=='__main__':main()
