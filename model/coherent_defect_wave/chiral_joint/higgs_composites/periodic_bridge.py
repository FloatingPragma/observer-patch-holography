"""Exact composite/elementary/hopping residue identity on one periodic source."""
from collections import defaultdict
from functools import lru_cache
import numpy as np
import sympy as s
from .model import screened_eom,adjoint,clean
from ..model import sequence,fermion
from ..two_cell.symmetries import oriented_hopping
from ..periodic_spectrum.basis import coordinate_rows
from ..periodic_spectrum.spectrum import resolved_block,numerical_block
from ..periodic_spectrum.readout import GROUND,probe_matrix


@lru_cache(maxsize=None)
def composite_column(key,kind,site,spin):
    op=screened_eom('charged_lepton',9,spin)
    if kind=='remove':op=adjoint(op)  # [H,-c]=[H,c^dagger]^dagger
    mapping={9:0,10:1,11:2,24:3,25:4,26:5};out=defaultdict(lambda:s.Integer(0))
    for (word,r,c),value in op.items():
        if key[3+site]!=c:continue
        operators=tuple((6*site+mapping[mode],create) for mode,create in word)
        result=sequence(key[0],operators)
        if result is not None:
            dest=list(key);dest[0]=result[0];dest[3+site]=r
            out[tuple(dest)]+=value*result[1]
    return clean(out)


def hopping_column(key):
    mask,l01,l10,hl,hr=key;out=defaultdict(lambda:s.Integer(0))
    for (target,link),value in oriented_hopping(mask,l01,0).items():out[target,link,l10,hl,hr]+=value
    for (target,link),value in oriented_hopping(mask,l10,1).items():out[target,l01,link,hl,hr]+=value
    return clean(out)


def probe_column(key,kind,site,spin):
    result=fermion(key[0],6*site+3*spin,kind=='add')
    return {} if result is None else {(result[0],*key[1:]):result[1]*(-1 if kind=='remove' else 1)}


def apply(state,action):
    out=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        for dest,weight in action(key).items():out[dest]+=value*weight
    return clean(out)


def hopping_eom_column(key,kind,site,spin):
    left=apply(probe_column(key,kind,site,spin),hopping_column)
    right=apply(hopping_column(key),lambda key:probe_column(key,kind,site,spin))
    return clean({k:left.get(k,0)-right.get(k,0) for k in set(left)|set(right)})


def project_operator(action,target):
    sourcekeys,source,_,_,_=resolved_block(*GROUND)
    targetkeys,destination,_,_,_=resolved_block(*target)
    positions={key:i for i,key in enumerate(targetkeys)};entries={}
    for j,key in enumerate(sourcekeys):
        for dest,value in action(key).items():
            assert dest in positions,('operator leaves specified charged Gauss weight sector',dest)
            entries[positions[dest],j]=value
    raw=s.ImmutableSparseMatrix(len(targetkeys),len(sourcekeys),entries)
    image=raw*source;rows=coordinate_rows(*target[:3]);restricted=image.extract(rows,range(image.cols))
    result=s.ImmutableSparseMatrix(restricted.rows,restricted.cols,
        {(i,j):s.expand(s.conjugate(s.I**((targetkeys[rows[i]][0]&63).bit_count()))*value)
         for (i,j),value in restricted.todok().items()})
    assert destination*result==image
    return result


@lru_cache(maxsize=None)
def transfer_matrices(kind='add',site=1,spin=0):
    target,elementary=probe_matrix(kind,site,spin)
    composite=project_operator(lambda key:composite_column(key,kind,site,spin),target)
    hopping=project_operator(lambda key:hopping_eom_column(key,kind,site,spin),target)
    h0=resolved_block(*GROUND)[4];ht=resolved_block(*target)[4]
    assert ht*elementary-elementary*h0==hopping+composite
    assert composite and hopping
    return target,elementary,composite,hopping


def residue_certificate():
    e0,v0=numerical_block(*GROUND);w=v0[:,0];rows=[]
    for kind in ('add','remove'):
        target,t,j,k=transfer_matrices(kind)
        energies,vectors=numerical_block(*target);metric=resolved_block(*target)[2]
        matrix=lambda m:np.asarray(m.tolist(),dtype=complex)
        gt=matrix(metric)
        amplitudes=[vectors.T@gt@matrix(op)@w for op in (t,j,k)]
        elementary,composite,hopping=amplitudes
        differences=(energies-e0[0])*elementary-composite-hopping
        residual=float(np.max(np.abs(differences)));assert residual<1e-11
        weights=[np.abs(a)**2 for a in amplitudes]
        total=[float(np.sum(x)) for x in weights]
        direct=[float(np.real(np.vdot(matrix(op)@w,gt@(matrix(op)@w)))) for op in (t,j,k)]
        assert np.max(np.abs(np.asarray(total)-direct))<1e-11
        # The spectral residue is an unnormalized operator two-point weight.
        # It is distinct from a probability for an implemented composite probe.
        from ..periodic_spectrum.readout import clusters
        groups,intervals=clusters(target)
        records=[{'energy_indices':indices,'rank':len(indices),
                  'source_energy_interval':[str(lo),str(hi)],
                  'energy_transfer':float(np.mean(energies[indices])-e0[0]),
                  'elementary_residue':float(np.sum(weights[0][indices])),
                  'quantum_composite_residue':float(np.sum(weights[1][indices])),
                  'quantum_hopping_residue':float(np.sum(weights[2][indices])),
                  'elementary_composite_cross_residue':float(np.real(np.vdot(elementary[indices],composite[indices]))),
                  'elementary_hopping_cross_residue':float(np.real(np.vdot(elementary[indices],hopping[indices]))),
                  'composite_hopping_cross_residue':float(np.real(np.vdot(composite[indices],hopping[indices])))}
                 for indices,lo,hi in groups]
        assert sum(row['rank'] for row in records)==len(energies)
        assert records[0]['quantum_composite_residue']>0
        rows.append({'kind':kind,'target':list(target),'exact_complete_Gauss_transfer_identity':True,
                     'source_dimension':t.cols,'target_dimension':t.rows,
                     'maximum_binary64_amplitude_identity_error':residual,
                     'complete_residue_sums':dict(zip(('elementary','composite','hopping'),total,strict=True)),
                     'complete_energy_residue_records':records,
                     'first_residue_enclosure':first_residue_enclosure(kind)})
    return {'source':'actual two-oriented-link periodic lepton/Higgs source',
            'source_ground_contract':'periodic_spectrum.readout.GROUND; neutral L2 twoSz0 translation+1',
            'site':1,'spin':0,'same_H_exact_matrix_identity':'H_target T-T H_ground=K+J',
            'exact_eigenstate_identity':'(E_n-E_0)<n|T|0>=<n|K|0>+<n|J|0>',
            'quantum_Higgs_matrix_used':True,'residues_are_instrument_probabilities':False,
            'residue_precision':'binary64 diagnostic on complete isolated spectral clusters; exact operator identity separately proved',
            'per_vector_amplitudes_retained_as_physical_evidence':False,
            'diagnostic_replay_absolute_tolerance':1e-10,
            'isolated_first_eigenstate_composite_overlap_certified':True,
            'charged_boundary_compensation_inherited':True,'branches':rows,
            'physical_rest_mass_identified':False,'GeV_conversion':False}


def first_residue_enclosure(kind='add'):
    """Outward rational bound, using the parent's complete spectral enclosures."""
    from ..periodic_spectrum.readout import (ground_approximation,approximation,
                                            clusters,projector_error,SCALE)
    w,norm,ground_error=ground_approximation()
    target,_,j,_=transfer_matrices(kind)
    z,_,_,_=approximation(target);groups,intervals=clusters(target)
    indices,lo,hi=groups[0]
    assert indices==[0], 'the first spectral projector must be isolated and rank one'
    g=resolved_block(*target)[2];q=s.ImmutableMatrix(z[:,indices].tolist())/SCALE
    x=j*w;amplitude=q.T*g*x
    estimate=s.factor((amplitude.H*amplitude)[0]/norm)
    assert estimate.is_Rational
    # ||J||<=2 by the two exact norm-one scalar/CAR summands.  The ground
    # ground_error bounds pure-state trace distance. The positive effect
    # J^dagger P J<=4I contributes at most 4*ground_error; 8 is conservative.
    error=8*ground_error+4*projector_error(target,indices,intervals)
    precision=1<<27
    lower=s.Rational(s.floor(max(0,estimate-error)*precision),precision)
    upper=s.Rational(s.ceiling(min(4,estimate+error)*precision),precision)
    assert lower>0
    return {'isolated_energy_indices':indices,'source_energy_interval':[str(lo),str(hi)],
            'quantum_composite_residue_interval':[str(lower),str(upper)],
            'quantum_composite_residue_interval_decimal':[float(lower),float(upper)],
            'operator_norm_upper_bound':'2',
            'nonzero_overlap_certified':True,
            'residue_is_instrument_probability':False}
