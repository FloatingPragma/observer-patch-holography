from fractions import Fraction as F
import pytest
from model.openwave_native_process import source as native
from .source import (endpoint_matrix_unit_certificate, graph_certificate,
                     population_column,coherence_multiplier,apply_population,coarse)

def test_complete_endpoint_matrix_units_and_nonclassical_coherence():
    certificate=endpoint_matrix_unit_certificate()
    assert certificate['complete_endpoint_operator_columns']==48
    # C retains a coherence discarded by E; averaging only populations losesit.
    assert native.unread_column(1,4,(0,1),'C')=={(1,4):F(1,2)}
    assert native.unread_column(1,4,(0,1),'E')=={}
    assert native.unread_column(1,4,(0,1))=={(1,4):F(1,4)}

def test_global_branch_columns_and_exact_memory_counterexample():
    row=graph_certificate(0)
    assert row['complete_offdiagonal_columns']==240*239
    classes=row['coherence_multiplier_classes']
    for degrees in ([5,6],[6,5]):
        assert sum(x['ordered_column_count'] for x in classes if x['degrees']==degrees)==180*60
    assert F(row['cell_projection_one_step_l1_max_exact'])==F(11,7560)
    assert F(row['two_step_uniform_cell_memory_l1_exact'])>0
    assert sum(map(F,row['two_step_cell_memory_vector_exact']))==0

def test_reference_safe_mass_conservation_and_wrong_rate_rejected():
    g=native.graph(0)
    values=[F(i==0) for i in range(g.ports)]
    actual=apply_population(g,values)
    assert sum(actual)==1 and min(actual)>=0
    assert actual[g.neighbours[0][0]]==F(1,2*len(g.seams))
    assert actual[g.neighbours[0][0]]!=F(1,len(g.seams))
    assert all(0<coherence_multiplier(g,0,j)<1 for j in range(1,g.ports))

def test_projection_cannot_hide_zero_cell_difference():
    g=native.graph(0)
    boundary=next(i for i in range(12) if any(j//12!=0 for j in g.neighbours[i]))
    interior=next(i for i in range(12) if all(j//12==0 for j in g.neighbours[i]))
    v=[F(i==boundary)-F(i==interior) for i in range(g.ports)]
    assert not any(coarse(v))
    assert any(coarse(apply_population(g,v)))
    with pytest.raises(ValueError):population_column(g,g.ports)
