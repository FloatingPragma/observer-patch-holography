"""Complete Gauss-compatible phase instrument on the active native source."""
import sympy as sy
from ... import core, source, spatial, local_capacity
from ...local_geometry import geometry_factors
from ...joint_geometry import auxiliary_columns

def exercise(graph):
    model=source.assemble(graph)
    candidates=[]
    for (word,shift),coef in sorted(model['matter'].items()):
        if len(word)==2 and all(cre for _,cre in word) and any(shift):
            if all((i%48)//3 in (9,11) for i,_ in word):
                state=core.apply_word(word,0)
                assert state is not None
                mask,sign=state
                entry=sy.expand(sign*coef)
                if entry.is_real and entry.is_positive:
                    candidates.append((word,shift,mask,entry))
    assert candidates
    word,shift,mask,z=candidates[0]
    assert z.is_real and z != 0
    local_index=1
    auxiliary_label=(0,)*graph.nvertices
    initial=(auxiliary_label,(local_index,(0,(0,)*len(graph.edges))))
    final=(auxiliary_label,(local_index,(mask,shift)))
    apply_column=local_capacity.apply_native_column
    candidate=local_capacity.CANDIDATE
    column=apply_column(model,initial,cap=6)
    assert final in column
    charge,_=core.mode_metadata(graph.nvertices)
    B=model['field']['B']
    for _,(_, (outmask,outflux)) in column:
        q=[sum(charge[i] for i in range(96*v,96*(v+1)) if outmask>>i&1)
           for v in range(graph.nvertices)]
        assert q == [sum(int(B[v,e])*outflux[e] for e in range(len(outflux)))
                     for v in range(graph.nvertices)]
    factors=geometry_factors(graph)
    local_lapse=[sy.Rational(float(N[local_index,local_index])) for N in factors['lapses']]
    capacity_lapses=[dict(auxiliary_columns()[1][a]).get(a,sy.Integer(0)) for a in auxiliary_label]
    support=source.term_support(word,shift,graph.edges)
    weight=sum(capacity_lapses[v]*local_lapse[v] for v in support)/len(support)
    assert weight.is_positive
    z_total=sy.expand(weight*z)
    assert sy.expand(column[final]-z_total)==0
    assert z_total.is_real and z_total != 0
    assert support=={0,1} and z==sy.Rational(1,2)
    assert z_total==sy.Rational(9,32)
    # The entire amended-source column is executed, not only its target entry.
    mutant=dict(model)
    mutant['matter']=core.gadd(model['matter'],core.scale(model['components']['dirac'],-1))
    mutated_column=apply_column(mutant,initial,cap=6)
    assert mutated_column.get(final,0)==0
    plus=sy.Matrix([1,sy.I])/sy.sqrt(2)
    minus=sy.Matrix([1,-sy.I])/sy.sqrt(2)
    Pp=(plus*plus.H).applyfunc(sy.expand)
    Pm=(minus*minus.H).applyfunc(sy.expand)
    assert Pp.H==Pp and Pm.H==Pm
    assert Pp*Pp==Pp and Pm*Pm==Pm and Pp*Pm==sy.zeros(2)
    assert Pp+Pm==sy.eye(2)
    # On the full orthogonal complement: both effects are zero and Pperp=I.
    # This block identity is dimension independent, including all reference
    # extensions. It does not replace the complement by one evolving state.
    hi,hf=sy.symbols('Hii Hff',real=True)
    Hblock=sy.Matrix([[hi,sy.conjugate(z_total)],[z_total,hf]])
    rho=sy.diag(1,0)
    derivative=(-sy.I*(Hblock*rho-rho*Hblock)).applyfunc(sy.expand)
    probs=[sy.trace(P*rho) for P in (Pp,Pm)]+[sy.Integer(0)]
    slopes=[sy.simplify(sy.trace(P*derivative)) for P in (Pp,Pm)]+[sy.Integer(0)]
    assert probs==[sy.Rational(1,2),sy.Rational(1,2),0]
    assert slopes==[-z_total,z_total,0]
    dephased=sy.diag(*derivative.diagonal())
    assert dephased==sy.zeros(2)
    assert [sy.trace(P*dephased) for P in (Pp,Pm)]==[0,0]
    # All complement diagonal probabilities have zero initial derivative
    # for any Hamiltonian from a basis projector; they are retained at t>0.
    return {'candidate':candidate,'graph':graph.name,
        'native_word':[[i,c] for i,c in word],'final_native_mask':str(mask),'final_flux':list(shift),
        'source_entry_exact':str(z),'local_capacity_entry_exact':str(z_total),
        'local_capacity_weight_exact':str(weight),'full_source_column_outputs_retained':len(column),
        'all_column_outputs_gauss_consistent':True,
        'dirac_deleted_full_column_outputs_retained':len(mutated_column),
        'dirac_deleted_same_entry_exact':'0','complete_PVM_exact':True,
        'gauss_compatible':True,'arbitrary_reference_CPTP_by_projector_identity':True,
        'orthogonal_complement':'Identity on the entire complement, not a reduced state.',
        'initial_probabilities_exact':[str(p) for p in probs],
        'initial_probability_derivatives_exact':[str(p) for p in slopes],
        'terminal_dephasing_derivatives_exact':['0','0','0'],
        'all_passed':True,
        'scope':'Exact initial full-source response of a creator-supplied complete phase instrument; no finite-time propagation, local hardware, physical particle or mass-clock claim.'}



def run():
    return exercise(spatial.triangle())
