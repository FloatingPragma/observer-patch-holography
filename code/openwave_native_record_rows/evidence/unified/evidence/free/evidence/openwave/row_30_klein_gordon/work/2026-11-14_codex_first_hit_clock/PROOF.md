# The native first-hit clock realizes the Euclidean slice transfer

This theorem supplies an actual stopping-time readout of the slice transfer
already reconstructed from occupation correlations. It does not equate a
unit's raw event trajectory with a Lorentzian worldline. The Lorentzian free
field is the effective reconstruction of the full record correlations, as
explained in row29's interpretation audit.

## 1. Finite source and observer

Take the twelve-port graph at each point of a finite three-dimensional
periodic spatial array, and vertical coordinates z=-R,...,1, R>=0. Internal
seam weights are1 and the port-preserving spatial and vertical weights are
w=a^(-2). At the bottom there is no seam to z=-R-1. The top plane z=1 is a
stopping set. One native unit starts at z=0 with its internal port and
spatial position disclosed. Every event is the pinned endpoint rule, with
the chosen seam and returned coin/current retained.

The observer reconstructs the current location solely by continuity. It
stops at the first retained record that places the unit at z=1 and reports
the spatial/internal endpoint. This is an intervention on the preparation
and a record-defined stopping rule. It does not read hidden endpoint state.

Two disclosed schedulers give the same stopped-location law:

* Global IID: choose any seam with probability weight/W, W the total weight.
* Incident: from the record-reconstructed location choose an incident seam
  proportional to its weight, then apply the unchanged fair endpoint rule.

For a one-unit source, every nonincident attempt is a zero-current hold. In
the global scheduler, conditional on selecting an incident seam its law is
exactly the incident scheduler's law. Deleting those holds preserves every
successive location and first-hit endpoint. This is an exact random time
change, not equality of event counts or full coin tapes. An incident-scheduler
run is a registered adaptive N+ protocol; it cannot be relabelled a full
global-IID tape. Its shorter record is sufficient for the stopping observable.

## 2. Exact harmonic transfer

Let f be an eigenfunction of the internal-plus-spatial Laplacian with

    A=lambda+w q_spatial >=0,

where lambda is a native internal eigenvalue and q_spatial is the dimensionless
cubic symbol. For the hitting observable f(endpoint), the conditional expected
value from height z is f(start_spatial, start_port) h_z. The first-step equation
of either scheduler cancels its waiting probability and gives

    (A+2w)h_z-w h_(z-1)-w h_(z+1)=0,  -R<z<1,
    (A+w)h_(-R)-w h_(-R+1)=0,
    h_1=1.

These equations follow from current-event probabilities; a relativistic
Hamiltonian is not inserted. For A>0 define

    omega=2 asinh(sqrt(A/w)/2).

The unique solution is

    h_z=cosh(omega(z+R+1/2))/cosh(omega(R+3/2)).

The two interior exponential roots obey r+r^(-1)=2+A/w. The half-step
cosh satisfies the reflecting-bottom equation by direct substitution, and
the denominator fixes the top value. Uniqueness follows from the finite
maximum principle, or from almost-sure absorption in this connected finite
chain. At A=0 the solution is identically1, including the constant mode.

Thus the directly observed one-plane first-hit transfer multiplier is

    h_R(A)=cosh(omega(R+1/2))/cosh(omega(R+3/2)).              (1)

In particular it is strictly positive. For A>0,

    0 <= h_R(A)/exp(-omega)-1
      = [exp(-2omega(R+1/2))-exp(-2omega(R+3/2))]
         /[1+exp(-2omega(R+3/2))]
      <= exp(-2omega(R+1/2)).                               (2)

The infinite-depth limit is the exact slice transfer exp(-omega) of the
native occupation Green function. With physical slice length a,

    omega/a=(2/a)asinh(a sqrt(lambda+q_a(p))/2),

the same native dispersion and masses as the reviewed row30 theorem. This
proves a model-to-transfer correspondence, beyond merely naming -log of an
arbitrary matrix an energy.

## 3. Successive clock readings on the same history

On the infinite-depth cylinder, define tau_j as the first native hit of plane
j, starting below it. The vertical embedded walk is unbiased and recurrent,
so each tau_j is finite almost surely. This cylinder law can be defined by
the record-controlled incident scheduler, or as the consistent first-hit
limit of finite global-IID strips. A nonexistent uniform distribution over
infinitely many global seams is not assumed.

The strong Markov property and vertical translation symmetry give the same
endpoint transition H at each clock increment. Its internal/spatial
eigenvalues are exp(-omega), hence the n-clock endpoint transfer is H^n with
eigenvalues exp(-n omega). These clock readings are successive stopping
events of one history; they are not independent drawn copies or a cosine
supplied to a decoder. Every intermediate current is part of the history.

This is an operational Euclidean clock. It is entirely consistent that its
transfer is positive and contractive. The positive reconstruction uses this
same generator to describe Lorentzian waves. Requiring the primitive stop
sequence itself to have unitary probabilities would be a different demand.

## 4. Resource and censoring bounds

The ideal infinite-depth clock has no finite mean attempt cost. Finite strips
must therefore carry their resource accounting. In the global-IID schedule,
the vertical transition probabilities are w/(2W) upward and downward in the
interior, and only upward at the bottom. Solving the first-step equations
for the expected hit cost, starting at z=0, gives

    E_global tau=2W(R+1)/w.                                (3)

For the incident scheduler each interior vertex has weighted degree
d=5+8w, and each bottom vertex d-w. The same difference equations give

    E_incident tau=2d(R+1)/w-2.                            (4)

For w=1,R=2 this is76 actual endpoint attempts. For a finite cutoff T,
Markov's inequality gives P(tau>T)<=E tau/T. A sharper deterministic tail
can be computed from the finite substochastic pre-hit transfer Q:

    P(tau>T)=e_start^T Q^T 1.

For a bounded endpoint observable |f|<=B, reporting zero on a censored run
has unconditional bias at most B P(tau>T). Conditioning only on completed
runs changes the measure and needs its own normalization; it is not used
silently. All censored attempts and their probability are to be retained.

Equation(2) controls the finite-bottom discrepancy, and these tail bounds
control bounded runtime. Neither bound replaces the actual endpoint sampling
uncertainty. The architecture, slice choice and scheduler are disclosed, while
the transfer law and its spectrum follow from native record probabilities.
