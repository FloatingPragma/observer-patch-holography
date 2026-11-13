"""Complete pinned quark/diquark internal spectra, including all source terms."""
from fractions import Fraction as F
from math import comb
import numpy as np
import sympy as sy
from scipy.linalg import eigh
from model.openwave_nuclear_distance.spatial import carrier
from model.openwave_common_hyperfine_spectrum import source as native
from .source import shared_length_attachment


def constrained_potential():
    d,p=carrier(2),native.parameters()
    attach=shared_length_attachment();kappa=F(attach['kappa']);a=F(attach['a_fm']);unit=F(native.INPUTS['electron_unit_GeV'])
    ids=np.flatnonzero((d['n'][:,0]==1)&(d['n'][:,1]==2))
    assert len(ids)==40
    S=sum(d['angular']((x,),True) for x in range(2))[ids][:,ids].toarray()
    I=sum(d['angular']((x,),False) for x in range(2))[ids][:,ids].toarray()
    counts=d['counts'][ids]
    nu=counts[:,[j for j,o in enumerate(d['orbitals']) if o[2]==0]].sum(axis=1)
    D=np.diag(3-2*nu)
    from math import factorial
    gram=np.array([6//np.prod([factorial(int(v)) for v in row]) for row in counts],dtype=np.int64)
    def rational_action(value):
        integers=np.rint(4*value*np.sqrt(np.outer(gram,gram))).astype(np.int64)
        assert np.max(abs(4*value*np.sqrt(np.outer(gram,gram))-integers))<2e-12
        return sy.Matrix(40,40,lambda i,j:sy.Integer(int(integers[i,j]))/(4*sy.sqrt(int(gram[i]*gram[j]))))
    base=sy.Rational(p['delta']/2)*rational_action(D)+(sy.Rational(p['Js']+2*p['spin_gradient'])*rational_action(S)
          +sy.Rational(p['Jw']+2*p['eta'])*rational_action(I))
    constant=3*p['mbar']+F(5,9)*p['A']+F(8,3)*(p['U']-kappa/2)+6*p['t']
    matrices=[];rows=[]
    for r in (2,3,4):
        # Pinned sites x=2,y=2+r in an eight-cell chain: both degree2,
        # no adjacent occupied sites, so all finite-range terms coincide.
        scalar=constant+F(4,3)*kappa*r
        H=base+sy.Rational(scalar)*sy.eye(40)
        matrices.append(H)
        spectrum=eigh(np.asarray(H.evalf(18),float),eigvals_only=True)
        rows.append({'eight_cell_pinned_sites':[2,2+r], 'separation_cells':r,
                     'separation_fm':str(a*r),'complete_internal_dimension':40,
                     'all_constrained_source_levels_GeV':(spectrum*float(unit)).tolist(),
                     'electric_energy_GeV':str(F(4,3)*kappa*r*unit)})
    for H,J in zip(matrices,matrices[1:]):
        assert J-H==sy.Rational(F(4,3)*kappa)*sy.eye(40)
    slope=F(4,3)*kappa*unit/a
    assert slope==1
    return {'complete_internal_singlet_dimension':40,
            'native_pinned_occupation_dimension':12*comb(12,2),
            'full_eight_cell_three_quark_dimension':comb(96,3),
            'singlet_success_on_uniform_pinned_input':str(F(40,12*comb(12,2))),
            'coloured_failure_on_uniform_pinned_input':str(1-F(40,12*comb(12,2))),
            'outside_pinned_number_sector_probability':'0',
            'outside_pinned_number_sector_effect':'Identity minus the full pinned occupation projector',
            'complete_reference':'Charge-balancing reference correlated by physical charge sector, retained in every physical spectral outcome.',
            'physical_spectral_instrument':{'input':'I/40 on the full pinned Gauss-dressed singlet','complete_energy_projectors':'Every one of40 internal spectral ranks is retained','probability_per_resolved_rank':'1/40','outside_physical_pinned_carrier':'Complementary identity effect with zero probability on this input','mean_work':'0 for the nonselective unital spectral instrument'},
            'colour_rank_count_scope':'The792-dimensional occupation colour filter is a separate kinematic preparation count; physical spectral poststates lie in the40-dimensional Gauss kernel.',
            'complete_spectra':rows,
            'exact_total_constrained_H_difference_per_cell':str(F(4,3)*kappa),
            'exact_slope_GeV_per_fm':str(slope),
            'all_source_terms_included':{'number_and_flavour':'3mbar+delta(3-2Nu)/2',
                'local_baryon_density':'5A/9','local_colour':'(8/3)(U-kappa/2)',
                'graph_hopping_compression':'6t identity; every off-diagonal hop leaves the pinned sector',
                'local_and_edge_spin':'(Js+2lambda)(S_quark^2+S_diquark^2)',
                'local_and_edge_weak_isospin':'(Jw+2eta)(I_quark^2+I_diquark^2)',
                'global_colour':'0 on all40 singlets','pair_triplet':'0 on all pure-quark states',
                'Gauss_electric':'(4kappa/3)r identity'},
            'interpretation':'Total-source constrained quark/diquark potential: every internal eigenvalue shifts by the same linear amount. External position constraints are declared; this is not the freely evolving unpinned global ground energy or a continuum scattering potential.',
            'calibration':'The seven strong inputs and the quadrupole/illustrative tension attachment are explicit;1GeV/fm is an attachment target, not an independent prediction.'}
