# Omit exact identity recoveries

The existing supplied moving-code recovery decodes and re-encodes every
positive continuing flag, including flags whose native image is the
target code. This separate variant removes those redundant words and retains
the original moving-defect and biased-reset receipts as baseline evidence.

Run `python -m pytest -q model/recovery_elision` and
`python -m model.recovery_elision.variant` to check the certificate.

For a fixed known code `V_c`, every native flag is a positive weighted partial
injection with equal-size disjoint logical-column supports. When the two image
supports equal the target code supports, their uniform positive amplitudes
give the exact column identity

`K_h V_c = sqrt(p_h) V_c'`.

Therefore identity is a valid flagged recovery. This equality extends
linearly to arbitrary logical states entangled with an untouched reference.
The decision reads only the code label and complete flag, and no branch is
conditioned or normalized. It is an identity on the admitted flagged image;
the original and optimized full-register unitaries need not agree elsewhere.

The complete 43,200-slot census contains 11,580 positive flags and 3,120
eligible identity recoveries. All eligible flags keep the core label fixed.
Moving flags retain the existing recovery. Boundary reset services, zero and
abort slots, native repair words, record retention, and core transition
probabilities are unchanged.

For the existing biased reset process with uniform seams and equal C/E tags,
the exact stationary saving in supplied source pulses per native attempt is

`4737999067212019707204144087023030301441659146 / 97909756773858989776401485931436511086245`,

approximately 48,391.49, or 7.9 percent of the retained baseline cost. The
largest omitted recovery contains 196,460 pulses. These are counts of ideal
supplied service words. They establish no physical duration, energy, apparatus
error allowance, new particle interpretation, or reduction of required resets.

The optimized recovery is exposed in `variant.recover` and its matching
`variant.recovery_resources`. Existing callers continue using the original
benchmark until they explicitly adopt this variant and bind its complete
execution certificate.
