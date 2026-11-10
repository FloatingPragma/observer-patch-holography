"""Full-complement graph and isolated spectral-projector residue errors."""
import sympy as s
from ..atomic_even_remainder.bounds import uniform_quartic_remainder
from ..family_cell.quotient import sqrt_interval


def residue_errors(gap, hopping_norm, coupling, ground_gap_coefficient, cluster_separation):
    gap, hopping_norm, coupling, ground_gap_coefficient, cluster_separation = map(s.sympify,
        (gap, hopping_norm, coupling, ground_gap_coefficient, cluster_separation))
    if any(v.is_Rational is not True or v <= 0 for v in
           (gap,hopping_norm,coupling,ground_gap_coefficient,cluster_separation)):
        raise ValueError('strictly positive exact rational residue controls required')
    v = coupling*hopping_norm
    if v > gap/5:
        raise ValueError('outside the full-source graph interval')
    graph = s.cancel(2*v/(gap-v))
    embedding = s.cancel(graph+graph**2/2)
    remainder = uniform_quartic_remainder(gap,hopping_norm,coupling)
    scaled_error = s.cancel(remainder/coupling**2)
    if scaled_error >= min(ground_gap_coefficient,cluster_separation)/4:
        raise ValueError('rescaled spectral clusters are not strictly isolated')
    ground_projector = s.cancel(scaled_error/(ground_gap_coefficient-scaled_error))
    cluster_projector = s.cancel(2*scaled_error/cluster_separation)
    full_ground_vector = s.cancel(embedding+2*ground_projector)
    full_cluster_projector = s.cancel(2*embedding+cluster_projector)
    amplitude = s.cancel(full_ground_vector+full_cluster_projector)
    return {'full_hopping_perturbation_norm':v,'invariant_graph_norm':graph,
        'graph_isometry_minus_atomic_embedding_norm':embedding,
        'full_Q_effective_H_remainder':remainder,'rescaled_effective_H_error':scaled_error,
        'effective_ground_projector_error':ground_projector,
        'effective_triplet_projector_error':cluster_projector,
        'full_ground_vector_error_after_phase_alignment':full_ground_vector,
        'full_triplet_projector_error':full_cluster_projector,
        'line_and_outside_record_amplitude_error':amplitude,
        'total_success_and_failure_amplitude_error':full_ground_vector}


def probability_interval(atomic_probability, amplitude_error):
    lo, hi = map(s.Rational, atomic_probability)
    amplitude_error=s.sympify(amplitude_error)
    if not 0 <= lo <= hi <= 1 or amplitude_error.is_Rational is not True or amplitude_error < 0:
        raise ValueError('a probability interval and nonnegative rational amplitude error are required')
    a,b=sqrt_interval((lo,hi),80)
    lower=max(s.Integer(0),a-amplitude_error)
    upper=min(s.Integer(1),b+amplitude_error)
    return s.cancel(lower**2),s.cancel(upper**2)
