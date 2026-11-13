"""Source-bound repeated preparation/reset/readout certificate."""
from fractions import Fraction as F
from .algebra import certificate as exact_certificate
from .controller import run as controller_run, SIDE, INTERNAL, SCALE, J, REGIONS, WEIGHTS
from model.openwave_passive_clock_records.clock import native_attachment

def run():
    masses=native_attachment()
    exact=exact_certificate()
    cases=[controller_run(F(row['source_rest_energy_electron_units_exact']),g) for row in masses for g in (2,3)]
    twirl=controller_run(F(1),2,True)
    return {'schema':'oph.repeated_native_mass_clock.v1','verdict':'PASS_REPEATED_REFERENCE_RESET_AND_COMPLETE_FINITE_CONTROLLER',
        'native_mass_attachments':masses,'exact_two_cycle_resource_theorem':exact,
        'static_controller':{'executed_reducing_dimension':SIDE*INTERNAL,'controller_positions':SIDE,'internal_occupations':INTERNAL,
            'scale_native_units_exact':str(SCALE),'hopping_over_scale_exact':str(J),'regions':{k:list(v) for k,v in REGIONS.items()},
            'initial_clock_integer_amplitudes':list(WEIGHTS),'initial_clock_support':[4,37],
            'initial_clock_phase':'i^x','normalization':'square root of the exact sum of squared integer amplitudes',
            'preparation_reference_reused':True,'phaseful_reference_count':3,'initially_empty_reset_battery_count':1,
            'gap_schedule':'Preparation gap m on x<64 and256<=x<272; readout1 gap m on x>=160; readout2 gap m on x>=400; battery/source gap m everywhere.',
            'reset_exchange_strength_over_scale_exact':'22/7',
            'semiboundedness':'H_total/L >= -max(abs(g),2,22/7); all source/reference energies and open-chain kinetic are positive.',
            'native_reduction':'The inherited Gauss-dressed charged pair has N=diag(0,0,1,1) and A=|3><0|. These are restrictions of the full native Fock Hamiltonian and bounded partial isometry. The reducing R carrier is executed, all 32768 occupation effects retained, and R-perp is an additional zero-probability effect on this input.',
            'complete_native_mass_and_changed_strength_controls':cases,'twirled_three_reference_control':twirl},
        'current_law':'For a clock corridor P without exchange, d<P A>/dt=i*m<P A>+i*L<[K,P]A>; the boundary contribution is bounded by2*L*J times the two adjacent-cell amplitude products.',
        'frequency_scope':'Both ideal isolated flights have exact frequency m independent of the collision strength. Actual finite-controller window frequencies have certified intervals and finite contact/boundary bias; exact strength independence of the complete repeated detector is not established.',
        'measurement_scope':'All terminal clock, native and four-resource occupations are recorded. Readout1 is a quantum record until the terminal instrument; its possible overwrite is bounded, not assumed absent. Clock-conditioned flight summaries retain every rejected clock outcome in the complete record.',
        'resource_scope':'Three initial phase references, one empty reset battery, one finite moving clock packet and its position-dependent gap architecture are supplied. Complete closed energy and terminal pinching work are retained. No reset of the battery, reference replenishment, memory erasure or self-starting sustained cycle is claimed.',
        'scope':{'two_complete_preparation_readout_cycles_executed':True,'same_preparation_reference_reused':True,
            'three_native_rest_energies_executed':True,'finite_time_propagation_error_enclosed_under_stated_arithmetic':True,
            'static_controller_exact_strength_independent_frequency':False,'sustained_replenished_clock':False,
            'native_twelve_port_process_attachment':False,'physical_seconds_identified':False,
            'additional_nature_validations':0,'literal_openwave_row3_validated':False}}
