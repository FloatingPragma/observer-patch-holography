"""Native infinite-chain bands with complete spectral-record bounds."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from . import exact_native
from .certified_rest import full
from .parameters import degree_coefficients
from .positivity import uniform_bound,joint_stability
from .source import probe,chain
from .analytic import certify_formula


def strings(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:strings(v)for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [strings(v)for v in value]
    if isinstance(value,float):return round(value,10)
    return value


@lru_cache(None)
def calculate():
    native=exact_native.certificate();formula=certify_formula()
    fixed=full();selected=full(*degree_coefficients())
    mn=F('.93956542194');ga=F('1.2753')
    assert fixed['full_source_GeV_interval'][1]<mn
    assert fixed['full_source_axial_interval'][0]>ga
    assert selected['full_source_GeV_interval'][0]<mn<selected['full_source_GeV_interval'][1]
    assert selected['full_source_axial_interval'][0]<ga<selected['full_source_axial_interval'][1]
    assert selected['full_source_axial_interval'][0]>ga-F('0.0013')
    assert selected['full_source_axial_interval'][1]<ga+F('0.0013')
    controls=[probe(L)for L in(2,3,4)]
    assert (chain(2,True)['H']-chain(2,False)['H']).nnz==0
    records=[]
    for L in(2,3,4):
        dimension=comb(4*L+2,3);rank=L+4*(L-1)
        records.append({'cells':L,'physical_N3_colour_singlet_dimension':dimension,
            'compact_fixed_projection_rank':rank,'filter_success_on_uniform_physical_input':F(rank,dimension),
            'filter_failure_on_uniform_physical_input':1-F(rank,dimension),
            'failed_source_state':'The complete orthogonal complement in the physical colour-singlet source; no occupation poststate outside Gauss is used.',
            'all_other_Fock_inputs':'Identity failure outcome; prepared weight zero.',
            'mean_filter_work':'Zero on I/d by the complete unital two-projector instrument.',
            'reference':'Opposite native charge correlated with each initial charge sector; unchanged by the charge-preserving filter.'})
    delta=chain(2)['parameters']['delta']
    unit=F('.00051099895069')
    proton=[x-delta*unit for x in selected['full_source_GeV_interval']]
    return strings({'schema':'oph.hadron_refinement_flow.v1',
        'verdict':'PASS_NATIVE_BULK_REST_BAND_AND_AXIAL_ENCLOSURES_WITH_CALIBRATED_DEGREE_FLOW',
        'exact_native_reduction':native,'exact_Bloch_coupling':formula,
        'degree_coefficients':{'colour':degree_coefficients()[0],'number':degree_coefficients()[1]},
        'fixed_coefficient_bulk_control':fixed,'selected_bulk_neutron':selected,
        'selected_bulk_proton_GeV_interval':proton,'native_vector_current_exact':'1',
        'finite_open_chain_controls':controls,'complete_preparation_records':records,
        'full_Fock_capacity_independent_stability':uniform_bound(True),
        'C_completion_joint_stability_without_nuclear_term':joint_stability(),
        'complete_infinite_fibre_spectral_instrument':{
            'carrier':'All relative three-quark positions in the infinite one-dimensional native Gauss chain at k=0,S=I=1/2.',
            'apparatus':'Five isolated low spectral outcomes and the complete high-spectrum complement; the ground outcome supplies the displayed mass interval.',
            'native_ground_state_compact_complement_probability_upper':selected['full_Q_probability_upper'],
            'native_ground_state_compact_probability_interval':[1-selected['full_Q_probability_upper'],1],
            'mean_spectral_readout_work':'Zero: the complete spectral projectors commute with the autonomous source.',
            'arbitrary_external_quantum_references':True},
        'calibration_ancestry':'The seven two-cell calibration domains are unchanged. Two degree-dependent coefficients use the inherited neutron-mass and axial targets. These agreements are calibrated enclosures, not independent predictions.',
        'scope':'One-dimensional infinite-chain native colour-singlet rest bands. Other bulk hadrons, three-dimensional Lorentz covariance, scattering poles and a physical continuum dictionary are not established by this transfer.'})
