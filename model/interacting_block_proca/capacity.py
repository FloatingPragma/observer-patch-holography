"""Outward log-dimension bounds for the explicit finite source and readout."""
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import mpmath as mp


def endpoint(interval,upper=False):
    sign,mantissa,exponent,_=interval._mpi_[int(upper)]
    return F((-1)**sign*mantissa)*F(2)**exponent


def display(value,upper=False):
    value=F(value)
    with localcontext() as context:
        context.prec=13
        context.rounding=ROUND_CEILING if upper else ROUND_FLOOR
        return format(Decimal(value.numerator)/Decimal(value.denominator),'.12E')


def logical_capacity(N,V,kicks):
    assert isinstance(N,int) and isinstance(V,int) and isinstance(kicks,int)
    assert N>0 and V>0 and kicks>=0
    mp.iv.dps=160;iv=mp.iv
    blocks=N*V
    block=240*iv.log(2)+8*iv.log(5)+24*iv.log(119)
    source=blocks*block
    reference=iv.log(3072);pointer=2*kicks*iv.log(2)
    total=source+reference+pointer
    lo,hi=endpoint(total),endpoint(total,True)
    return {'native_eight_cell_blocks':blocks,'native_cells':8*blocks,
        'fermion_modes':240*blocks,'five_state_scalar_registers':8*blocks,
        '119_state_link_registers':24*blocks,
        'reference_dimension':3072,'retained_pointer_and_control_dimension':4**kicks,
        'Hilbert_dimension_expression':f'2^{240*blocks} * 5^{8*blocks} * 119^{24*blocks} * 3072 * 2^{2*kicks}',
        'logHilbertDimension_nats_formula':f'{blocks}*(240*log(2)+8*log(5)+24*log(119))+log(3072)+{2*kicks}*log(2)',
        'one_block_log_dimension_nats_interval_exact':[str(endpoint(block)),str(endpoint(block,True))],
        'source_log_dimension_nats_interval_exact':[str(endpoint(source)),str(endpoint(source,True))],
        'reference_log_dimension_nats_interval_exact':[str(endpoint(reference)),str(endpoint(reference,True))],
        'pointer_log_dimension_nats_interval_exact':[str(endpoint(pointer)),str(endpoint(pointer,True))],
        'logHilbertDimension_nats_interval_exact':[str(lo),str(hi)],
        'logHilbertDimension_nats_display_interval':[display(lo),display(hi,True)],
        'display_interval_is_also_outward_rounded':True,
        'carrier_scope':'exact full kinematic native source times retained reference and readout registers; upper bound for its Gauss-physical subspace',
        'replica_number_origin':'declared architectural replica count selected by sufficient record-error inequalities; not derived from a memory budget or scheduler law',
        'capacity_is_log_dimension_not_actual_state_entropy':True,
        'capacity_law_scope':'finite-dimensional entropy/record ceiling log(dim H); no replica count derived from a repair-law scheduler or cosmological budget',
        'additional_gate_synthesis_or_ground_preparation_workspace_included':False,
        'below_10_power_122_nats':hi<10**122,
        'comparison_to_10_power_122_scope':'numeric comparison only; no universe-capacity or observed-continuum identification'}
