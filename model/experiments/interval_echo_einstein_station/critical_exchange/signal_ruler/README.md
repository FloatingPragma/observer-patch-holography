# A signal ruler on the cooled field source

An emitter at port 2532 swaps its first ring station with the adjacent independent rail at 2499. One arm loads the rail, the other leaves it empty. The displaced sender is discarded. A receiver counts its own occupied outcome on every launch of either arm. Subtracting the empty-arm rate from the loaded-arm rate gives a signal with the same travel profile on the cold sea and on two excited preparations. One excitation increases the independently read bond energies at stations 51 and 55 by \(0.000265185\) and \(0.000244539\). For the other, an individual receiver-arm probability changes by as much as 0.00492929. [Regional energy measurements](../README.md) · [Signal and source controls](results.json)

The cold preparation replays the native work-comparison history starting with all 123 loads in the upper orbitals. It takes 11,828 proposals to reach the selected lower band. One held source replaces the registered coherent pair of orbitals. The other turns the occupied orbital of the signed energy-redistribution experiment through \(\pi/4\) toward the enrolled empty orbital. The signal experiment uses the same guide and one fixed propagation program on every launch. [Energy selector](../selection/README.md)

Write \(C\) for the one-particle correlation matrix of the incoming material and \(P=|s\rangle\langle s|\) for the sender mode. Placing the sender first in the occupation order makes the physical reset a tensor product of the chosen sharp load with the incoming marginal on every other station. In particular \(C_b=(I-P)C(I-P)+bP\) for a rail carrying \(b=0,1\), including entangled incoming material. For the number-conserving quadratic propagation with one-particle matrix \(U(t)\), the receiver's rate is \([U(t)C_bU(t)^\dagger]_{rr}\). Hence the exact signed signal is
\[
p_r(1,t)-p_r(0,t)=|\langle r|U(t)|s\rangle|^2 .
\]
This compares recorded occupancies in two cohorts and includes exclusion of a second load at an occupied binary station. It applies equally to arbitrary coherent, mixed or non-Gaussian inputs of this binary guide.

The common program is the native Fourier rotation, mode phases for \(\epsilon_k=-\cos k\), and inverse rotation. All two-mode rotations in this compilation use adjacent stations in an open ordering. It acts as the same fermionic one-particle unitary in every filling sector. On the enrolled 123-load source it equals the ring hopping action. A literal spin hopping gate on the bond closing a ring introduces the fermion parity factor \((-1)^{M-1}\), and a reset can change \(M\). Direct word evolution distinguishes this physically different implementation in the four-station check below.

The emitter and receivers use integer action-clock appointments from 0 through 80. The calibration assigns each receiver the peak of the signed signal on that clock grid. The three preparations give the following predictions. In the frozen count banks the maximum-rate appointment agrees on all three sources at every receiver.

| Separation in calibrated guide stations | Peak appointment | Peak occupied-rate excess |
| ---: | ---: | ---: |
| 1 | 2 | 0.332612 |
| 4 | 5 | 0.153063 |
| 8 | 10 | 0.101031 |
| 16 | 18 | 0.068164 |
| 51 | 54 | 0.032354 |
| 55 | 58 | 0.030774 |

All 2,916 individual Born rates were evaluated on the prepared covariances; subtracting each pair agrees with the independently summed Fourier propagation probability to \(1.4\times10^{-15}\). The additional synthetic scan contains 262,144 copies for each single-arm appointment. A 99 percent simultaneous Hoeffding bound over all elementary rates is 0.005033 per rate. An appointment within four such radii of the observed maximum remains compatible with being the population peak. The first two peaks are resolved on that test. At separations 16, 51 and 55 multiple appointments remain compatible on all three sources. Every original count and population probability is archived. Replay checks the frozen counts and validates every rate against the freshly calculated native output. [Count banks](counts.npz) · [Replay](measure.py)

For an interaction test enroll four consecutive binary stations, with only the three neighboring contacts active, the sender at 0, receiver at 2, and action time 3. Compare an empty background with a spectator at station 1. Exact evolution of all 16 binary words gives the following loaded-minus-empty readings.

| Four-station operation | Empty background | One spectator |
| --- | ---: | ---: |
| Quadratic open hopping | +0.367524 | +0.367524 |
| Hopping plus \(\sum_j n_j n_{j+1}\) | +0.367524 | +0.222847 |
| Three active native blank seam terms at \(\kappa=8\) | +0.003191 | −0.987056 |
| Literal spin hopping on a closed four-site ring | +0.990018 | +0.358187 |

The active blank term adds \(8.5P_{\ne}\) on every bond to the hopping; external federation contacts and reference services are parked in this four-station control. The large negative signed response means a loaded sender reduces the receiver's occupancy relative to the empty-sender arm in that background. It supplies a direct instance of matter-sensitive native propagation. The longer guided ruler measures one spatial ordering and a scheduled action clock. Identifying a common proper area, moving clocks or gravitational curvature requires other observer measurements. [Fixed-filling and independent 16-word checks](test_signal.py)
