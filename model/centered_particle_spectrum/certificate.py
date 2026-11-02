"""Exact rational checks of untrusted finite eigensystem proposals.

Every witness coefficient is a signed integer. All admission norms, interval
endpoints and positive-weight conclusions use Python integers/Fraction; LAPACK
is used only to propose the separately retained witness.
"""
from fractions import Fraction as F
from math import isqrt, comb
import numpy as np
from . import source

VECTOR_BITS=60
ENERGY_BITS=56
ROOT_BITS=100


def sqrt_bounds(value,bits=ROOT_BITS):
    if type(value) not in (int,F) or value<0 or type(bits) is not int or not 1<=bits<=512:
        raise ValueError('nonnegative exact radicand and bounded positive bit count required')
    value=F(value);scale=1<<bits
    root=isqrt((value.numerator*scale*scale)//value.denominator)
    lower=F(root,scale)
    return lower,lower if lower*lower==value else F(root+1,scale)


def decimal_bounds(lower,upper,places):
    if lower>upper or type(places) is not int or places<0:raise ValueError('ordered exact interval required')
    scale=10**places
    a=(lower*scale).numerator//(lower*scale).denominator
    b=-((-upper*scale).numerator//(-upper*scale).denominator)
    def render(n):
        sign='-' if n<0 else '';whole,frac=divmod(abs(n),scale)
        return sign+str(whole)+('.'+str(frac).zfill(places) if places else '')
    return [render(a),render(b)]


def _integer_array(value,shape,name):
    value=np.asarray(value)
    if value.shape!=shape or value.dtype!=np.dtype(np.int64):
        raise ValueError('exact int64 '+name+' with shape '+str(shape)+' required')
    return value.astype(object)


def eigensystem(kind,columns,values,vectors):
    if np.asarray(columns).ndim!=1 or np.asarray(columns).dtype!=np.dtype(np.int64):
        raise ValueError('exact projector column archive required')
    chosen=tuple(int(i) for i in columns)
    block=source.reduced(kind,chosen);n=block['dimension']
    ei=_integer_array(values,(n,),'energy numerators')
    vi=_integer_array(vectors,(n,n),'vector numerators')
    energies=tuple(F(int(e),1<<ENERGY_BITS) for e in ei)
    if any(a>b or (kind=='local' and a==b) for a,b in zip(energies,energies[1:])):
        raise ValueError('strictly ordered simple-eigenvalue proposal required')
    y=block['basis'].astype(object)@vi
    ty=np.zeros(y.shape,dtype=object)
    for i,row in enumerate(block['coefficient_source']):
        for j in np.flatnonzero(row):ty[i]+=int(row[j])*y[j]
    r=ty*(1<<ENERGY_BITS)-16*y*ei[None,:]
    gram=y.T@(block['sizes'][:,None]*y)
    norm_numerators=[int(gram[i,i]) for i in range(n)]
    norms=[F(v,1<<(2*VECTOR_BITS)) for v in norm_numerators]
    if any(v<=0 for v in norms):raise ArithmeticError('zero witness vector')
    defect=gram-(1<<(2*VECTOR_BITS))*np.eye(n,dtype=object)
    delta=F(max(sum(abs(int(v)) for v in row) for row in defect),1<<(2*VECTOR_BITS))
    if not delta<F(1,1000000):raise ArithmeticError('complete-basis Gram defect exceeds declared domain')
    squares=[F(sum(int(d)*int(v)**2 for d,v in zip(block['sizes'],col,strict=True)),
               1<<(2*(VECTOR_BITS+ENERGY_BITS+4))) for col in r.T]
    radius=sqrt_bounds(sum(squares,F(0))/(1-delta))[1]
    if not radius<F(1,10**9):raise ArithmeticError('complete-spectrum residual exceeds declared domain')
    residuals=[];errors=[];intervals=[]
    for i,e in enumerate(energies):
        separation=min(abs(e-other) for j,other in enumerate(energies) if j!=i)-radius
        residual=sqrt_bounds(squares[i]/norms[i])[1]
        isolated=separation>2*radius and residual<separation
        if not isolated and (kind=='local' or i==0):
            raise ArithmeticError('required spectral disc is not separately certified')
        residuals.append(residual)
        errors.append(2*residual/separation if isolated else None)
        bound=residual if isolated else radius
        intervals.append((e-bound,e+bound))
    return {'block':block,'columns':chosen,'energy_numerators':ei,'vector_numerators':vi,
            'orbit_vectors':y,'norms':norms,'gram_defect':delta,'complete_radius':radius,
            'residuals':residuals,'vector_errors':errors,'energy_intervals':intervals,
            'energies':energies}


def orbit_dimensions(local):
    """Irreducible orbit dimensions from exact finite-group tight-frame traces."""
    block=local['block'];q=block['q'].tocoo()
    y=np.zeros((q.shape[0],local['orbit_vectors'].shape[1]),dtype=object)
    y[q.row]=q.data[:,None]*local['orbit_vectors'][q.col]
    moments=[]
    for shift in range(8):
        rows,sign=source.fock_action(4,3,(0,1,2),shift)
        values=np.sum(y*(sign[:,None]*y[rows]),axis=0)
        moments.append([F(int(v),1<<(2*VECTOR_BITS))/norm
                        for v,norm in zip(values,local['norms'],strict=True)])
    dimensions=[]
    for i,error in enumerate(local['vector_errors']):
        approximate=sum((row[i]**2 for row in moments),F(0))
        tolerance=32*error
        candidates=[d for d in range(1,9) if approximate-tolerance<=F(8,d)<=approximate+tolerance]
        if len(candidates)!=1:raise ArithmeticError('finite-group orbit dimension not uniquely certified')
        dimensions.append(candidates[0])
    return dimensions,moments


def verify(witness):
    required={kind+'_'+name for kind in ('ground','local') for name in ('columns','energies','vectors')}
    if set(witness)!=required:raise ValueError('missing or extra eigensystem witness arrays')
    systems={kind:eigensystem(kind,witness[kind+'_columns'],witness[kind+'_energies'],
                             witness[kind+'_vectors']) for kind in ('ground','local')}
    ground=systems['ground'];local=systems['local']
    transition=source.transition(ground['columns'],local['columns']).astype(object)
    amplitudes=local['vector_numerators'].T@(transition@ground['vector_numerators'][:,0])
    rows=[];exact_weights=[];exact_gaps=[]
    dimensions,moments=orbit_dimensions(local)
    for i,num in enumerate(amplitudes):
        rootlo,roothi=sqrt_bounds(ground['norms'][0]*local['norms'][i])
        if rootlo<=0:raise ArithmeticError('normalization interval touches zero')
        raw=F(abs(int(num)),1<<(2*VECTOR_BITS))
        alo,ahi=raw/roothi,raw/rootlo
        error=ground['vector_errors'][0]+local['vector_errors'][i]
        if alo<=error:raise ArithmeticError('nonzero weight is not certified for line '+str(i))
        wlo,whi=(alo-error)**2,(ahi+error)**2
        elo=local['energy_intervals'][i][0]-ground['energy_intervals'][0][1]
        ehi=local['energy_intervals'][i][1]-ground['energy_intervals'][0][0]
        if elo<=0:raise ArithmeticError('positive observed excitation energy not certified')
        exact_weights.append((wlo,whi));exact_gaps.append((elo,ehi))
        overlaps=[]
        for moment in moments:
            lo=max(F(-1),moment[i]-2*local['vector_errors'][i])
            hi=min(F(1),moment[i]+2*local['vector_errors'][i])
            overlaps.append(decimal_bounds(lo,hi,places=18))
        spatial=[]
        signs=source.shift_gram_signs()
        for x in range(8):
            line=[]
            for z in range(8):
                mlo,mhi=(F(v) for v in overlaps[x^z])
                if signs[x,z]<0:mlo,mhi=-mhi,-mlo
                products=(wlo*mlo,wlo*mhi,whi*mlo,whi*mhi)
                line.append(decimal_bounds(min(products),max(products),places=32))
            spatial.append(line)
        rows.append({'index':i,'excitation_energy_interval':decimal_bounds(elo,ehi,14),
                     'local_branch_weight_interval':decimal_bounds(wlo,whi,32),
                     'energy_diagnostic':float((elo+ehi)/2),
                     'weight_diagnostic':float(raw*raw/(ground['norms'][0]*local['norms'][i])),
                     'strictly_positive_weight_proved':True,
                     'visible_spatial_multiplicity':dimensions[i],
                     'normalized_spatial_overlap_intervals':overlaps,
                     'removal_spatial_residue_intervals':spatial})
    normalization=(sum((w[0] for w in exact_weights),F(0)),sum((w[1] for w in exact_weights),F(0)))
    if not normalization[0]<=F(1,2)<=normalization[1]:raise ArithmeticError('half-occupation normalization not enclosed')
    second=(2*sum((w[0]*e[0]**2 for w,e in zip(exact_weights,exact_gaps,strict=True)),F(0)),
            2*sum((w[1]*e[1]**2 for w,e in zip(exact_weights,exact_gaps,strict=True)),F(0)))
    if not second[0]<=F(67,256)<=second[1]:raise ArithmeticError('exact signed second-moment sum rule not enclosed')
    mean=(2*sum((w[0]*e[0] for w,e in zip(exact_weights,exact_gaps,strict=True)),F(0)),
          2*sum((w[1]*e[1] for w,e in zip(exact_weights,exact_gaps,strict=True)),F(0)))
    basis_counts={'ground_spatial_invariant':ground['block']['orbit_dimension'],
                  'ground_spin_zero_eta_zero':ground['block']['dimension'],
                  'local_stabilizer_invariant':local['block']['orbit_dimension'],
                  'local_spin_half_eta_half':local['block']['dimension']}
    visible=sum(dimensions);full=comb(16,7)
    result={'reference':'unique centered N8, hidden-spin0, eta0 ground; K=H-N/2, c=1/8,U=1',
            'energy_unit':'declared dimensionless source-clock unit; no MeV/eV calibration',
            'basis_dimensions':basis_counts,
            'ground_energy_interval':decimal_bounds(*ground['energy_intervals'][0],places=14),
            'positive_visible_lines_per_local_spin_charge_branch':len(rows),
            'signed_local_retarded_atoms':2*len(rows),'all_visible_weights_strictly_positive':True,
            'unresolved_local_weight':'0','rows':rows,
            'spatial_residue_convention':'B_E[x,y]=<c_xdown G,P_E c_ydown G>; the hidden-spin factor is I_2',
            'addition_spatial_residue_rule':'A_E=eta B_E eta; eta_x=(-1)^(x1+x2+x3)',
            'local_addition_weight':'1/2','local_removal_weight':'1/2',
            'weight_interval_sum':decimal_bounds(*normalization,places=14),
            'signed_zeroth_moment':'1','signed_first_moment':'0','signed_second_moment':'67/256',
            'second_moment_interval':decimal_bounds(*second,places=14),
            'conditional_mean_excitation_energy_interval':decimal_bounds(*mean,places=14),
            'conditional_second_excitation_moment':'67/256',
            'visible_spatial_dimension_per_fixed_spin_weight':visible,
            'reachable_N7_dimension_all_linear_CAR_probes':2*visible,
            'unobserved_N7_dimension':full-2*visible,
            'reachable_N9_dimension_all_linear_CAR_probes':2*visible,
            'unobserved_N9_dimension':full-2*visible,
            'multiplicity_counts':{str(d):dimensions.count(d) for d in sorted(set(dimensions))},
            'proof_bounds':{kind:{'complete_eigenbasis_Gram_defect_upper':str(system['gram_defect']),
                                  'complete_spectrum_enclosure_radius_upper':str(system['complete_radius']),
                                  'maximum_unit_eigenvector_error_upper':str(max(x for x in system['vector_errors'] if x is not None))}
                            for kind,system in systems.items()},
            'global_neutral_gap_certified':False,'all_N7_N9_eigenvalues_resolved':False,
            'finite_energies_identified_as_relativistic_rest_masses':False,
            'physical_particle_species_identified':False,'native_observer_measurement_executed':False,
            'continuum_interaction_or_gravity_proved_by_this_packet':False}
    return result
