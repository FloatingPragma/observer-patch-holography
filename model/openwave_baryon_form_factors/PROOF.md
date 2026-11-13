# Charge profiles on a common strong source

Use the seven-input hyperfine action, its exact native three-quark colour
singlet and the full positive Gauss-law attachment in
`openwave_common_hyperfine_spectrum`. The two-cell operator identity preserves
the complete singlet Hamiltonian and its native currents. The quadrupole
and illustrative string-tension comparisons select one common length; no
proton or neutron charge-radius value selects that length or an action
coefficient.

## Exact isospin transport

Let v be the integer neutron trial from the complete even spin-half,
isospin-half reduction in `openwave_common_hyperfine_spectrum`. Its normalized
distance from the unique selected ground ray is bounded by the retained
residual divided by the next spectral distance, denoted s. Apply the
literal native isospin raising operator I+ before performing any reduction.
The resulting integer proton vector has the same squared norm. Every term
other than the flavour number splitting is an isospin scalar. Consequently

    H I+ v = I+ H v - delta I+ v.

The producer verifies this equality with rational arithmetic on all 2024
native occupation coordinates. The two trial variances agree exactly and
the same ground-ray bound applies. Their energy difference is delta times
the unchanged electron GeV unit, with no extra energy calibration. These
are strong-source calibrated bands, not weak-contact or continuum pole
masses.

## Intrinsic records and one physical length

For each native occupation let X=(x1+x2+x3)/3 be its number centre. Record
the actual native charges, 2/3 and -1/3, at distances |xi-X|. The possible
distances are 0,1/3,2/3 in graph units. Summing the exact integer-vector
probabilities gives signed radial weights w(r), with sum w equal to the
charge. On the neutron ray the inner nonzero radius has positive charge
and the outer radius has negative charge. The radius-zero contribution is
zero. These signs have certified margins under the ground-ray error.

For a common length a, the charge second moment is

    R2 = a^2 sum_r w(r) r^2.

Every three-quark occupation has sum of absolute charges at most two and
radius at most2/3. Thus the moment observable has norm at most8/9. The
pure-state trace distance is at most2s, giving second-moment error at
most16a^2s/9 and each radial-weight error at most4s. Signed neutron moments
are kept signed; no square root of a negative moment is called a radius.

For an explicitly uniform observation average over the orientation of
this graph, a plane-wave charge probe gives

    F(k) = sum_r w(r) sinc(k a r),
    -6 dF/d(k^2)|0 = R2.

The producer encloses each sinc by its rational Taylor polynomial and the
remainder |z|^(2M+2)/(2M+3)!, then adds the state error4s. The momentum
conversion uses the declared displayed value of hbar c. This orientation
average is an observation law. It does not construct an isotropic spatial
Hamiltonian or identify the electromagnetic Sachs form factor of a
relativistic nucleon.

## Complete physical preparation and measurement

Let B be the normalized native epsilon embedding of all120 colour singlets
into all2024 occupations. Every occupation row of B has support in at most
one orbital column. The positive effects

    E_i = B* |i><i| B

sum to the identity on the complete physical singlet carrier. Their
positive square roots define a complete instrument on that carrier and
transport to the Gauss kernel by the proved dressing isometry. This uses
compressed occupation effects, not coloured postmeasurement states.

Start with I/120 on this physical carrier and retain a charge-balancing
quantum reference correlated with each charge sector. The rank-one trial
filter succeeds with probability1/120. Conditional occupation probability
p_i is v_i^2/||v||^2. If d_i=(BB*)_ii, the complete original-input records are

    P(success,i)=p_i/120,  P(failure,i)=(d_i-p_i)/120.

The producer verifies every nonnegative weight, sum d_i=120 and total
probability one. Every unlisted native label has both probabilities zero;
the full external Fock/Gauss complement has its identity outcome. A
successful preparation loses no unspecified colour sector. Both filter and
compressed occupation instruments are unital and charge preserving. The
declared joint maximally mixed/reference input is unchanged nonselectively,
so its mean source and reference work is exactly zero. The conditional
success and failure states need not have zero work separately.

Replacing the trial by the exact ground ray changes the full unconditional
classical record distribution in total variation by at most2s/120.
The conservative half-channel bound4s controls arbitrary external quantum
references and subsequent common instruments.

## Independent radius comparison

The common length cancels from -R_n^2/R_p^2. Exact interval division gives
a source interval disjoint from the box formed by the supplied CODATA2022
proton radius and PDG2025 neutron mean-square radius. The comparison box
uses their separate quoted standard uncertainties; no joint statistical
likelihood is asserted. Consequently no choice of a common length matches
both quoted radius intervals on this selected source. Positive-core and
negative-shell records and the neutron/proton mass ordering remain exact
properties of the same declared member. They do not establish its full
physical nucleon identification.
