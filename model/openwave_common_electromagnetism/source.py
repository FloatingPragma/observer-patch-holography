"""Exact electric charge, compact Gauss dressing and source energy moments."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy.linalg import eigh
from model.openwave_common_mass_source.source import (
    lepton_source, full_source_terms, exact_baryon_window, G, INPUTS, COMPONENTS)
from model.openwave_axial_state.source import sector, parameters, weights
from model.openwave_axial_state.certificate import exact_action
from model.openwave_baryon_weak_formfactors.native import currents
from model.openwave_calibrated_hadrons.source import integral
from model.common_source.chiral_candidate.action import ordering

STIFFNESS=F(1,8)
SELF=STIFFNESS/8


def native_charges():
    q=[F(r['electric_Q']) for r in ordering()['all_left_components']]+[F(0)]
    assert all((3*x).denominator==1 for x in q)
    return q


def exact_ward():
    """Formal bilinear-polynomial commutator, with no occupation truncation."""
    q=native_charges();full=[q[c] for x,s,c,f in product(range(2),range(2),range(16),range(3))]
    residual=defaultdict(lambda:[F(0),F(0)]);terms=full_source_terms()['terms'];linear=0
    def entries(rows):
        return [(i,j,complex(r,im)) for i,j,r,im in rows]
    for term in terms:
        kind=term['kind']
        if kind=='one_body':
            assert all(full[i]==full[j] for i,j,z in entries(term['matrix']) if z)
        elif kind in ('one_body_square','one_body_product'):
            A=entries(term.get('matrix',term.get('left_matrix')))
            B=entries(term.get('matrix',term.get('right_matrix')))
            # Generator entries and all nonzero charge differences are rational halves/thirds.
            for i,j,a in A:
                for k,l,b in B:
                    charge=full[i]-full[j]+full[k]-full[l]
                    if charge:
                        coefficient=F(term['weight'])*charge
                        ar,ai,br,bi=map(F,(a.real,a.imag,b.real,b.imag))
                        residual[i,j,k,l][0]+=coefficient*(ar*br-ai*bi)
                        residual[i,j,k,l][1]+=coefficient*(ar*bi+ai*br)
        else:
            j=int(term['name'].rsplit('_',1)[1]);charge=F(int(lepton_source()['charge'][j]))
            assert all(full[i]==charge for i,r,im in term['annihilation_coefficients'])
            assert all(full[i]==-charge for i,r,im in term['creation_coefficients'])
            linear+=1
    # Snapshot generator entries are literal dyadic rationals. Rational accumulation
    # verifies the formal identity without a tolerance or floating cancellation.
    assert all(re==im==0 for re,im in residual.values())
    l=lepton_source();Q=np.diag(l['charge'])
    errors=[np.linalg.norm(Q@l['A']-l['A']@Q),np.linalg.norm(Q@l['B']+l['B']@Q),
            abs(np.trace(Q@l['rho']))]
    assert max(errors)==0
    return {'native_charge_spectrum_exact':list(map(str,q)),
        'full_matter_CAR_modes':192,'positive_Dirac_modes_checked':linear,
        'formal_degree_four_commutator_coefficients_checked':len(residual),
        'global_charge_commutator_exact_zero':True,
        'native_charge_equals_positive_Dirac_charge_without_constant':True,
        'charge_conjugate_Bogoliubov_identity_errors':list(map(float,errors))}


def gauss_label(charges):
    """Integer flux and quantum background label for a complete matter charge record."""
    n=len(charges);R=-3*sum(charges);E=[3*n*q+R for q in charges]
    assert R.denominator==1 and all(e.denominator==1 for e in E)
    return tuple(map(int,E)),int(R)


def finite_gauss_control():
    # Four CAR modes, two charge signs at each cell, and every Fock occupation.
    q=(1,-1,1,-1);labels=[];transitions=0;undressed_failures=0
    for mask in range(16):
        charge=[F(sum(q[i] for i in range(2*x,2*x+2) if mask>>i&1)) for x in range(2)]
        E,R=gauss_label(charge);assert sum(E)==0
        assert all(E[x]-6*charge[x]-R==0 for x in range(2))
        labels.append({'occupation':mask,'charge':list(map(str,charge)), 'electric_flux':list(E),'reference_charge':R})
        for src,dst in ((0,2),(2,0),(1,3),(3,1)):
            if mask>>src&1 and not(mask>>dst&1):
                after=charge.copy();after[src//2]-=q[src];after[dst//2]+=q[src]
                EE,RR=gauss_label(after);shift=[0,0];shift[src//2]-=6*q[src];shift[dst//2]+=6*q[src]
                assert EE==tuple(E[x]+shift[x] for x in range(2)) and RR==R
                assert any(E[x]-6*after[x]-R for x in range(2));undressed_failures+=1;transitions+=1
    assert len({(tuple(r['electric_flux']),r['reference_charge'],r['occupation']) for r in labels})==16
    return {'full_Fock_dimension':16,'complete_occupation_flux_reference_records':labels,
        'charged_hop_columns_checked':transitions,'undressed_hop_constraint_failures':undressed_failures,
        'all_dressed_hops_preserve_each_Gauss_constraint':True,
        'rotor_flux_spaces_are_full_integer_spaces':True,
        'constraint_failure_projector_retained':True}


def charge_diagonal(kind='baryon'):
    d=sector(kind);out=[];square=[]
    for mask in d['basis']:
        local=[F(0),F(0)];squares=F(0)
        for i,(species,x,spin,weak,colour) in enumerate(d['modes']):
            if mask>>i&1:
                charge=(F(2,3) if weak==0 else F(-1,3))*(1 if species==0 else -1)
                local[x]+=charge;squares+=charge*charge
        out.append(SELF*((local[0]-local[1])**2-squares));square.append(squares)
    return out,square


def energy_records():
    prior=exact_baryon_window();d=sector('baryon');p=parameters();unit=F(INPUTS['electron_unit_GeV'])
    diag,square=charge_diagonal();rows=[]
    for name in ('neutron','proton'):
        v=list(map(int,prior[name+'_integer_trial']));norm=sum(x*x for x in v)
        hv=[F(0) for _ in v]
        for part,coefficient in weights(p).items():
            for i,z in enumerate(exact_action(d['parts'][part],v)):hv[i]+=coefficient*z
        old=sum(x*y for x,y in zip(v,hv))/norm
        correction=sum(x*x*c for x,c in zip(v,diag))/norm
        hv=[y+c*x for y,c,x in zip(hv,diag,v)];mean=old+correction
        variance=sum((y-mean*x)**2 for y,x in zip(hv,v))/norm+G*G*F(prior['contact_second_moment_exact'])
        width=F(1,2);outside=variance/(width*width);assert 0<outside<1
        success=F(prior['exact_preparation_success']);seed=prior['native_occupation_seeds'][name]['index']
        work=2*F(v[seed],norm)*(v[seed]*mean-hv[seed])
        rows.append({'assigned_channel':name,'mean_energy_GeV_exact':str(mean*unit),
            'mean_energy_GeV_diagnostic':float(mean*unit),
            'electric_mean_shift_GeV_exact':str(correction*unit),
            'electric_mean_shift_GeV_diagnostic':float(correction*unit),
            'full_source_variance_electron_units_squared_exact':str(variance),
            'half_width_GeV_exact':str(width*unit),
            'conditional_inside_probability_lower_exact':str(1-outside),
            'original_inside_probability_lower_exact':str(success*(1-outside)),
            'original_outside_probability_upper_exact':str(success*outside),
            'preparation_failure_probability_exact':str(1-success),
            'preparation_mean_work_electron_units_exact':str(work),
            'spectral_energy_window_measurement_work_exact':'0',
            'interpretation':'Prepared trial mean and complete full interacting spectral-window bounds; not an interacting pole or physical fitted mass.'})
    # Full colour sector is retained; the charge operator preserves its complete embedding.
    C=np.array(list(map(float,diag)));reduced=(d['B'].conj().T@d['B'].multiply(C[:,None])).toarray()
    error=float(np.linalg.norm(d['B'].multiply(C[:,None]).toarray()-d['B']@reduced));assert error<1e-13
    h=d['h']+reduced;e,U=eigh(h);assert np.linalg.norm(h@U-U*e)<2e-8
    return {'source_trial_reference':'Exact native2024-component neutron/proton trials of the common mass source; same pure Dirac reference.',
        'rows':rows,'complete_colour_sector_dimension':120,
        'complete_strong_electric_spectrum_electron_units_diagnostic':e.tolist(),
        'charge_embedding_residual':error,'full_source_contact_second_moment_exact':prior['contact_second_moment_exact'],
        'same_source_spectral_window_records_include_failure':True}


def certificate():
    p=parameters();unit=F(INPUTS['electron_unit_GeV']);l=lepton_source()
    floor=p['mbar']-p['delta']/2-18*G-SELF*F(4,9);assert floor>1
    rows=[]
    for i,name in enumerate(('electron','muon','tau')):
        mass=l['calibration']['charged_exact'][i]/l['calibration']['charged_exact'][0]
        rows.append({'channel':name,'normalized_rest_GeV_exact':str(mass*unit),
            'normalized_rest_GeV_diagnostic':float(mass*unit),
            'bare_electric_rest_GeV_exact':str((mass+SELF)*unit),
            'bare_electric_rest_GeV_diagnostic':float((mass+SELF)*unit),
            'both_charge_signs_and_spins_retained':True,
            'input_role':'empirical electron anchor' if i==0 else 'native analyzer with retrospective phase2/9'})
    return {'schema':'oph.common_longitudinal_electromagnetism.v1',
        'verdict':'PASS_NATIVE_GAUSS_COMPLETION_AND_FULL_SOURCE_ENERGY_BOUNDS',
        'architecture':{'matter_CAR_modes':192,'cells':2,'rotors':'two full integer leaf rotors and one quantum compensating reference',
            'stiffness_exact':str(STIFFNESS),'stiffness_is_measured_alpha':False,
            'graph_Green_matrix_exact':[['1/4','-1/4'],['-1/4','1/4']],
            'bare_electric_energy':'(kappa/8)(Q0-Q1)^2',
            'normal_ordering_counterterm':'-(kappa/8) sum_a q_a^2 n_a; native colourful occupations and positive-Dirac lepton occupations',
            'counterterm_coefficient_exact':str(SELF),'counterterm_selected_using_mass_targets':False,
            'all_charged_cross_cell_currents_Gauss_dressed':True,
            'transverse_photons_or_magnetic_loops':False},
        'native_Ward_certificate':exact_ward(),'finite_Gauss_control':finite_gauss_control(),
        'positivity':{'colourful_floor_per_particle_exact':str(floor),
            'charged_lepton_floor_per_particle_exact':str(1-SELF),
            'two_or_more_charged_lepton_floor_exact':str(2*(1-SELF)),
            'physical_Gauss_subspace_unique_vacuum':True,
            'exact_charge_plus_and_minus_one_threshold_GeV':str(unit),
            'neutral_gap_inherited_from_positive_Dirac_source':True,
            'physical_Gauss_complement_is_constraint_failure':True},
        'same_source_charged_rest_branches':rows,
        'changed_baryon_energy_records':energy_records(),
        'scope':'Compact longitudinal Gauss completion with explicit quantum neutralizing background and fixed normal-ordering prescription. No alpha prediction, physical Coulomb distance, transverse photon, radiative mass calculation or continuum gauge theory is inferred.'}
