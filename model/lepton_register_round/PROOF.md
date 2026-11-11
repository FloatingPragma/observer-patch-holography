# The charged-lepton register, its weak round and its read counts

All arithmetic below is exact unless a tolerance is given: integers and
rationals, the field `Q(sqrt5)` for the carrier chart, integer pairs
`(A, B)` for numbers `A + sqrt2 B` in the register, and the field `Q(w)`,
`w = exp(i pi/3)`, `w^2 = w - 1`, for fluxes. Symbolic identities use sympy;
evaluated phases use mpmath at 40 digits. Each statement is a named check in
[results.json](results.json); the prefix of the name is the section number.

Every check on the public program calls it: `Observer.repair`, `reflect`,
`step`, `inverse_step`, `quarter_turn`, `embedded_port` and
`embedded_return_measurement` of `model/coherent_seam_repair/program.py`, the
local completion `local_repair_unitary` of `model/coherent_seam_repair/common.py`,
the edge gate `exchange` and `dressed_swap` of
`model/coherent_defect_wave/spatial/family.py`, and the native current rule
`native` of `model/observer/native.py`. The construction items that the
program does not contain are the premises at the end.

## 1. The register

The half-turns `g_C` and `g_E` that exchange the ends of the seams
`C = (2, 3)` and `E = (2, 10)` compose to `g = g_E g_C` of order three, which
fixes exactly the faces `f = (4, 5, 11)` and `f' = (6, 7, 8)`. The third-turn
`rho_f`, counterclockwise seen from outside, equals `g^-1`. The six ports on
neither face form the hexagon `2 -> 3 -> 9 -> 1 -> 0 -> 10` in the screen
sense of `f`, the sense of `n_f . (x_a x x_b) > 0` with `n_f` the face
centre; `rho_f` advances the hexagon by two steps. The C seam runs `2 -> 3`
and the E seam `10 -> 2`. (`face_is_fixed_by_the_two_seam_half_turns`,
`hexagon_and_screen_orientations`)

A seam object has a load in `{tail, head, idle}`, a record `r` in `{0, 1}`
and an address `o` in `{+, -}`, with `+` the screen orientation
`tail -> head`. Its event `Q` is `Observer.repair` at the two addresses of
its seam, read basis state by basis state: address `(tail, head)` for `o = +`
and `(head, tail)` for `o = -`. A load on any other port is unchanged at both
addresses, and the C and E objects have the same event matrix
(`object_event_is_observer_repair_on_the_seam`). `Q^2 = I`
(`object_event_is_an_involution`). With `R_2` the exchange of the two
orientations, `W = Q R_2`, and the oriented record `V = Z_r` at `o = -`, the
integer matrices satisfy `(2W)^T (2W) = 4I`, `R_2^2 = V^2 = I`
(`object_factors_unitary`). The coupling `T` applies `R` to F when C repairs
at `+` and `R^-1` at `-`. The tick `U_0 = T V (W_C (x) W_E)` on
`12 x 12 x 3 = 432` states satisfies `(4 U_0)^T (4 U_0) = 16 I` in integer
arithmetic (`full_432_state_tick_unitary`). The clock
`K_theta = exp(i theta J)` is diagonal with unimodular entries, so
`U(theta) = K_theta U_0` is unitary.

Each object is prepared as `sqrt(1/3)|tail> + sqrt(2/3)|idle>` with a blank
record and the uniform address. After `V W` both objects have the weights
idle `2/3`, agree (record 0 on the seam) `1/6`, repair at `+` `1/12`,
repair at `-` `1/12` (`object_weights_2_3_1_6_1_12_1_12`). The analyzer is
`Phi = (phi_H + phi_+ + phi_-)/sqrt3` with `phi_H` = C idle and E in the
agreeing state `(|tail> + |head>)/sqrt2`, record 0, uniform address;
`phi_+-` = C in the agreeing state with record 1 at address `+-` and E idle.
The three flags are orthonormal (`analyzer_flags_orthonormal`).

The compression `M(theta)[x', x] = <Phi (x) x'| U(theta) |in (x) x>` is
computed from the integer matrix `4 U_0` and the integer forms
`6 |in_C (x) in_E> = r_0 + sqrt2 r_1` and `4 sqrt3 Phi = sqrt2 h + 2 (e_+ + e_-)`.
Writing `M = M_0 + z M_+ + z^-1 M_-`, `z = exp(i theta)`, each part is
`(A + sqrt2 B)/(96 sqrt3)` with integer `A, B`. The program finds
`A_0 = 32 I`, `B_+ = 16 R`, `B_- = 16 R^-1` and all other parts zero, that is

```
M(theta) = (sqrt2/3) C(theta),   C(theta) = sqrt(1/6) I + sqrt(1/12) (z R + z^-1 R^-1)
```

(`compression_is_sqrt2_over_3_C_theta`). Every part commutes with `R`
(`compression_commutes_with_R`), and `(b/a)^2 = 2 * 16^2/32^2 = 1/2`
(`balance_b_over_a_is_one_over_sqrt2`). On the characters
`v_k = 3^(-1/2) sum_j w3^(-jk) |x_j>`, `w3 = exp(2 pi i/3)`, one has
`R v_k = w3^k v_k` and `M v_k = lambda_k v_k` with
`lambda_k = (sqrt2/3)(1/sqrt6 + (1/sqrt3) cos(theta + 2 pi k/3))`
(`eigenvalues_on_the_characters`). For `C(theta)` the eigenvalues satisfy
`sum lambda_k^2 = 1` and `sum lambda_k = sqrt6/2` for every theta, so
`Q = sum m_k/(sum sqrt m_k)^2 = 2/3` with `m_k` proportional to `lambda_k^2`
whenever all `lambda_k >= 0`; the eigenvalue of `k = 1` vanishes at
`theta = pi/12`, the edge of the positive chamber
(`koide_two_thirds_in_the_positive_chamber`).

Control: with `V = I` the same computation gives `B_- = -16 R^-1`, the
quadrature form `(sqrt2/3)(sqrt(1/6) I + sqrt(1/12)(z R - z^-1 R^-1))`
(`plain_record_gives_quadrature_form`). The oriented record reverses the
sign of the repair at `-`, so both oriented repairs enter `M` with the same
sign.

Read counts of the full 432-state vector `U_0 |in, x_0>`:
`P(C repairs +) = 1/12`, `P(E agrees) = 1/6`, `P(both) = 1/72`
(`repair_plus_1_12_hold_1_6_joint_1_72`), so the exclusive read has
`P(R+ xor H) = 1/12 + 1/6 - 2/72 = 2/9` (`exclusive_read_count_is_2_9`).
The three parts of `M` have disjoint supports and `A_c B_c = 0` entrywise,
so for every real theta the read effect `M^dag M = (2/9) C(theta)^2` has the
trace `sum_c (A_c^2 + 2 B_c^2)/(96^2 * 3) = 2/9`, the read count
(`trace_of_the_read_effect_equals_the_count`).

Register outputs: at `theta = 2/9`, inside the chamber, with electron, muon
and tau the characters `w3^1`, `w3^2`, `w3^0` and masses proportional to
`lambda_k^2`, the ratios are `m_mu/m_e = 206.770316` and
`m_tau/m_e = 3477.4728`, within half a unit of the last stated digit
(`register_outputs_at_theta_2_9`). No measured value enters this package.

## 2. Double loads on the hexagon

The hexagon object carries two loads on the six hexagon ports, the
four-valued record of the public completion and its own address register over
the twelve ordered hexagon seams, reflected about its uniform state. In the
completion, record 2 marks one of two loads on the first port of the ordered
seam moving to the second, record 3 the reverse. The tick is `W = Q R` and the
inverse tick `R Q`. The plain law applies the completion in the order of the
address, as `Observer.repair` does. The oriented law applies it at both
addresses of a seam in the one order whose first port lies in the third-turn
orbit `{2, 9, 0}`; on the C seam this is the screen orientation, the order
that the oriented record of section 1 selects.

The completion is a real symmetric involution
(`completion_is_a_real_symmetric_involution`). Let `D_p` be the blank double
at port `p` with the uniform address and `S = (D_2 + D_10)/sqrt2`. In both
laws one tick from `S` gives no emission `2/3`, emission across the E seam
(in-dimer, the hold H) `1/6`, emission from port 2 across `(2, 3)`
(forward, R+) `1/12` and from port 10 across `(0, 10)` (backward, R-)
`1/12`; for two independent objects `P(R+ xor H) = 2/9`
(`symmetric_double_tick_weights_2_3_1_6_1_12_1_12_and_count_2_9`).
`D_2` emits `1/6` forward and `1/6` in-dimer, `D_10` `1/6` backward and
`1/6` in-dimer (`one_port_doubles_emit_1_6_across_each_of_their_two_seams`).
The weights `2/3, 1/6, 1/12, 1/12` equal the object weights of section 1 term
by term (`double_load_weights_equal_the_register_object_weights`), so the
analyzer of section 1 on the normalized flag branches gives `(sqrt2/3) C(theta)`
on these objects with no activity input.

## 3. The weak round with the inverse tick

One round on two independent objects A and B in `S (x) S`: the tick on both,
the kick `exp(-i h J (x) Qx)` with `Qx` the projector on (A forward) xor
(B in-dimer), the inverse tick on both, and a return test that projects on
`S (x) S`. With `z = exp(-i h j)` for the eigenvalue `j` in `{-1, 0, 1}` of
`J`, the round maps `|SS>` to `|SS> + (z - 1)|X>`,
`X = (W^-1 (x) W^-1) Qx (W (x) W) |SS>`. In both laws:

- `W^-1 W S = S` (`oriented_inverse_tick_restores_the_pair`,
  `plain_inverse_tick_restores_the_pair`).
- `<S|W^-1 Pi_F W|S> = 1/12`, `<S|W^-1 Pi_I W|S> = 1/6`, `P = 2/9`. A zero
  return leaves `S (x) S`, so every round starts as the first and reads the
  tick-1 weights; the six-tick cycle of the next paragraph is never entered
  (`every_round_reads_the_tick_one_weights`).
- `<SS|X> = <X|X> = 2/9`: the zero-return amplitude is `1 - P + P z`
  (`zero_return_amplitude_is_1_minus_P_plus_P_z`) and the non-return weight is
  `|1 - z|^2 (<X|X> - <SS|X>^2) = 4 sin^2(h j/2) * 14/81`
  (`non_return_weight_is_14_81_per_unit_abs_1_minus_z_squared`).
- The non-return vectors `W^-1 Pi W S - s S` have the squared norms `11/144`
  (A) and `5/36` (B) (`non_return_norms_11_144_and_5_36`). Both lie on blank
  doubles. A has `10/144` at port 2 with a non-uniform address and `1/144`
  with the uniform address, which is the antisymmetric double
  `(D_2 - D_10)/sqrt2`; B has `5/36` with a non-uniform address
  (`non_return_vectors_are_blank_doubles`). A record read detects no
  non-return; the return test reads the address register and the symmetric
  location.
- The two laws give the same round numbers
  (`plain_and_oriented_rounds_agree`).

With N rounds per unit clock time (`h = 1/N`, `j = 1`) the amplitude is
`(1 - P + P exp(-i/N))^N`. The phase per unit clock time is
`0.205334030970278` at `N = 1`, `0.222062099670503` at `N = 10` and
`0.222222206218563` at `N = 1000`, within `1e-12`
(`phase_per_unit_clock_time_at_N_1_10_1000`). With no return read between
rounds, `W W^-1 = I` joins the N kicks into one kick of angle `t`, which gives
the `N = 1` value at `t = 1`. The exact series is

```
N atan(P sin(1/N) / (1 - P + P cos(1/N))) = 2/9 - 35/(2187 N^2) + O(N^-4)
N log|1 - P + P exp(-i/N)|               = -7/(81 N) + O(N^-3)
```

(`phase_series_2_9_minus_35_over_2187_N2`, `modulus_series_minus_7_over_81_N`),
and `|(1 - P + P exp(-i j/N))^N - exp(-i P j)| <= 1/N` at
`N = 1, 2, 10, 100, 1000` (`distance_to_the_limit_at_most_1_over_N`).

Control without the inverse tick, oriented law: the counts at ticks 1 to 6
are `2/9, 35/81, 38/81, 35/81, 2/9, 0`, with mean `8/27`, and `W^6 S = S`
(`forward_only_counts_2_9_35_81_38_81_35_81_2_9_0`,
`forward_only_mean_8_27_and_period_six`).

The strong round on the public one-load program, with its own methods only:
`embedded_port(0)`, `step`, `quarter_turn` on each of the 20 coordinates with
record 1 and nonzero amplitude (the kick with `h j = -pi/2` on the flag
"record 1"), `inverse_step`, `embedded_return_measurement(0)`. It gives
`inverse_step(step(x)) = x`, `P(record 1) = 1/12`, the return amplitude
`11/12 + i/12 = 1 - p + p i`, yes `61/72` and no `11/72 = 2 p (1 - p)`
(`public_one_load_round_returns_11_12_plus_i_12`).

## 4. The clock unit

The declared clock is
`H = sum_x (I - Re U_x) + sum_edges (pi/4) dGamma(I - S_xy)` with evolution
`exp(-i t H)`, and the tick's edge gate is `E = a I + b S`, `a = (1 - i)/2`,
`b = (1 + i)/2`. Since `S` is a Hermitian involution,
`exp(-i tau (pi/4)(I - S)) = c_0 I + c_1 S` with
`c_0 = (1 + exp(-i pi tau/2))/2` and `c_1 = (1 - exp(-i pi tau/2))/2`. This
equals `E` exactly when `exp(-i pi tau/2) = -i`, that is `tau` in `1 + 4Z`.

On the public `dressed_swap` and `exchange`, applied to a two-cell state on
actual link words with nontrivial frame signs: `S` is a nontrivial involution
and `E^2 = S` (`swap_is_a_nontrivial_involution_and_gate_squares_to_it`);
`(c_0, c_1)` at `tau = 1` is `(a, b)`
(`clock_at_unit_time_has_the_gate_coefficients`); among
`tau = 1, 5, -3, 9, -7, -1, 3, 2, 0, 4, 1/2, 4/5, pi/2, 2/pi` the clock state
equals the gate state exactly for `1, 5, -3, 9, -7`
(`clock_equals_the_edge_gate_exactly_on_1_plus_4Z`), the inverse gate for
`-1, 3` and `S` for `2`
(`clock_gives_inverse_gate_and_swap_on_minus_1_and_2_plus_4Z`);
`exp(-i pi (1 + 4n)/2) = -i` for integer `n`
(`exp_minus_i_pi_1_plus_4n_over_2_is_minus_i`). The one-body spectrum of
`tau (pi/4)(I - S)` is `{0, tau pi/2}`, and the principal logarithm of `E`
has its phases in `(-pi, pi]`, so `tau = 1` is the one member of `1 + 4Z`
with `|tau| pi/2 < pi` (`principal_branch_member_is_tau_1`).

Attempt count on `Observer.step`, with `A` the embedding of a port with a
blank record and the uniform address (`embedded_port`): for every port
`1 - <x|W x>/<x|x> = P(record 1 in W x) = 1/12`
(`one_minus_overlap_is_record_one_weight_1_12_on_every_port`); the identity
holds on six random rational consensus inputs
(`identity_on_random_consensus_inputs`) and fails for a blank record with a
non-uniform address (`identity_needs_the_consensus_address`). Its domain is
the consensus sector. `A^dag W A = K_12 = I - L/60` with `L` the port
Laplacian of the 30 seams (`consensus_block_is_I_minus_L_over_60`), and
`x - (W x + W^-1 x)/2 = A (L/60) e_p` on all twelve ports
(`onsite_term_keeps_the_consensus_sector`). The onsite term of the clock
therefore keeps the consensus sector and acts there as the generator
`I - K_12` of the attempt semigroup
`exp(-t (I - K_12)) = sum_n Poisson(n; t) K_12^n`, one attempt per unit time.

Reading: the edge term reproduces the tick's root-SWAP at clock time one, and
the onsite term runs one attempt, one tick, per unit time. With the rounds of
section 3 counted on this clock, the oriented current carries the phase
`theta = P = 2/9` per tick.

## 5. The kick and the orientation sign

On every seam `{a, b}` of the carrier, `Observer.repair` on the basis
`(a, 0), (b, 0), (a, 1), (b, 1)` is `1 - K^2 + K` at the address `(a, b)` and
`1 - K^2 - K` at `(b, a)`, with `K = |p,1><m,0| + |m,0><p,1|`,
`m = (|a> - |b>)/sqrt2`, `p = (|a> + |b>)/sqrt2` and `K^3 = K`; the images of
these four basis states stay in their span
(`repair_is_1_minus_K2_plus_sigma_z_K_on_every_seam`). The public event at
the two addresses of a seam is `Q = 1 - K^2 + sigma_z(o) K`, and so is the
object event of section 1 (`object_event_is_1_minus_K2_plus_sigma_z_K`). The
oriented record removes the sign: `V Q = 1 - K^2 + K` on every column with a
blank record, and the two differ on written records
(`oriented_record_gives_1_minus_K2_plus_K_on_blank_records`).

The oriented current is `J = |1><1|_r sigma_z(o) = (V_- - V_+)/2`, with
`V_+-` the record reflection `Z_r` at orientation `+-` and `Z_r = V_+ V_-`
(`current_is_half_difference_of_oriented_records`). `J` takes the values
`0, +1, -1` on the flags C idle, C repairs `+`, C repairs `-`
(`flags_carry_current_0_plus_1_minus_1`).

The named interaction is the oriented toggle
`T(eta) = exp(i eta J (x) Z_F Z_I)`, `Z_F = 1 - 2 Pi_F` on A and
`Z_I = 1 - 2 Pi_I` on B. In both laws
`(W^-1 (x) W^-1)(Z_F (x) Z_I)(W (x) W)|SS> = |SS> - 2|X>`
(`pair_toggle_identity_SS_minus_2X`). The standing pair carries no flag, so
`T(eta)^-1` multiplies it by `u = exp(-i eta j)`; after the tick `T(eta)`
multiplies a branch with flags `(f, i)` by `u^(-(1 - 2f)(1 - 2i))`. The
program computes the round `W^-1 T(eta) W T(eta)^-1 |SS>` as a polynomial in
`u` and finds, after the inverse tick, the coefficient `|SS> - |X>` at `u^0`,
the coefficient `|X>` at `u^2` and no other power. Hence
`W^-1 T(eta) W T(eta)^-1 = W^-1 exp(-2 i eta J (x) Qx) W` on the standing
pair: the kick of section 3 with `h = 2 eta`
(`kick_is_the_toggle_commutator_with_h_2_eta`). Its zero-return amplitude is
`1 - P + P exp(-2 i eta j)` with `P = <SS|X> = 2/9`
(`toggle_round_amplitude_1_minus_P_plus_P_exp_minus_2i_eta_j`).

Crossing: a seam with its own two-orientation address register, whose
reflection exchanges the orientations. In the plain law two ticks carry the
blank double from one end to the other with amplitude one, and four ticks
return it, on all 30 seams from both ends and both starting orientations, 120
starts (`double_crosses_its_seam_in_two_plain_ticks`). With one fixed order
at both addresses, two ticks return the double to its port
(`fixed_order_law_returns_the_double`).

With the crossing of A between rounds, plain law: B stays a double at port 2
read for the in-dimer flag, and A alternates between ports 2 and 10. The
round acts as `1 - c + c z` with `c = 5/18` for A at port 2 and `c = 1/6`
for A at port 10, mean `2/9`
(`one_port_round_coefficients_5_18_and_1_6_mean_2_9`); the no weights over
`|1 - z|^2` are `65/324` and `5/36` (`one_port_no_weights_65_324_and_5_36`).
Every non-return component keeps the blank double at its port with a zero
uniform-address part, so the return test is the address consensus at that
port, the form of `Observer.embedded_return_measurement`
(`one_port_non_returns_stay_at_the_port`). This round needs no superposition
of the two locations.

## 6. The sixty-address carrier

The carrier walk is the plain law with the sixty ordered seams of
`Observer.addresses`, the uniform reflection and the completion at the
addressed seam. On one load it equals `Observer.repair` and `Observer.reflect`
on all 1440 basis states
(`walk_equals_observer_repair_and_reflect_on_1440_basis_states`).

A double at any port with the uniform address gives after one tick: blank
`5/6`, record 2 `1/12`, record 3 `1/12`
(`double_load_blank_5_6_record_2_and_3_1_12_at_every_port`). Each port is
the first port of five addresses and the second of five. Every emission
amplitude is `+1/sqrt60` (`emission_amplitudes_all_positive`). Record 2 is
the native current `+1` (a load moves from the first to the second port of
the address) and record 3 is `-1`; the configuration after the event equals
the public native rule's for both coins on all `12 x 60` cases
(`records_2_and_3_are_the_native_current_signs`).

Two cells with independent registers: X holds the defect's register and
reads `R+` = record 2 and `R-` = record 3; a second cell Y reads `H` = a
written record. By an explicit sum over the `60 x 60` address pairs,
`P(R+) = 1/12`, `P(H) = 1/6` and `P(R+ and H) = 1/72` at every port
(`two_cells_read_1_12_1_6_joint_1_72_at_every_port`). The count is
`2/9 = 3/n - 4/n^2` with `n = 12` ports, from `P(R+) = 1/n` and
`P(H) = 2/n` (`two_cell_count_is_2_9_equal_to_3_over_n_minus_4_over_n2`), and
`H : R+ : R- = 2 : 1 : 1`.

A cell with one load on each of two non-adjacent ports (all 36 such pairs):
the address touches no load `2/3`, the agreeing half `1/6`, a repair with the
load at the first port `1/12` and at the second port `1/12`
(`native_pairs_give_2_3_1_6_1_12_1_12`). The activity `1/3` of section 1 is
the address count `20/60`. The repair amplitude is `+1/2` with the load at
the first port and `-1/2` at the second
(`second_endpoint_repair_has_the_minus_sign`); this is the sign the oriented
record removes.

One register: the entries of the completion are `0`, `1` and `+-1/2`
(`completion_entries_are_0_1_and_plus_minus_one_half`), so from a definite
configuration with a blank record and the uniform address of one sixty-address
register every branch weight is a multiple of `1/240`. The program checks
this for all 78 two-load configurations, and `240 * 2/9 = 160/3`
(`one_register_weights_are_multiples_of_1_240`). Two doubles on non-adjacent
ports of one register give `P(R+) = 1/12`, `P(H) = 1/6` and
`P(R+ and H) = 0` for all 36 placements, count `1/4`
(`one_register_two_doubles_count_1_4`). The count `2/9` needs two independent
address registers.

The weak round on two cells that each hold a double at port 2, and again at
port 10: the inverse tick is exact, the overlaps are `1/12` and `1/6`, and
`<start|X> = <X|X> = 2/9` (`weak_round_on_sixty_addresses_has_P_2_9`). The
non-return vectors have the squared norms `11/144` and `5/36`, a blank record,
the initial configuration and a zero uniform-address part, so a reset of the
address consensus restores the standing state and the return test is the
address consensus alone
(`sixty_address_non_returns_11_144_and_5_36_need_only_the_address_consensus`).

## 7. The class of the defect selects the pair

A class-k flux through the hexagon is the link phase `w^k` on the hexagon
seam `(2, 3)`: at both addresses of that seam the completion is conjugated by
`diag(w^(k n_3))`, `n_3` the load on port 3. The separated sector holds two
loads on distinct ports with record 0 or 1 (30 states); the pair sector holds
the blank doubles and the pairs on the two ends of a hexagon seam with
record 2 or 3 (18 states). The discriminant `D = A^dag Q A` is the address
average of the phased completion. For every class both sectors are closed
under `D` and `D` is Hermitian
(`sectors_closed_and_discriminant_hermitian_for_every_class`). Its fixed
vectors are the consensus states, the uniform-address states that the tick
leaves unchanged. Exact Gaussian elimination in `Q(w)` gives

| Class | Separated loads | Pair on one port |
| --- | --- | --- |
| 0 | 2 | 1 |
| 1, 2, 4, 5 | 0 | 0 |
| 3 | 0 | 1 |

(`consensus_states_by_class_0_2_1_class_3_0_1_others_none`). A single load
sees the holonomy `w^k`, a double load its square, and `w^(2k) = 1` only for
classes 0 and 3.

Around a class-3 flux the consensus state of two loads is unique. All 18 of
its amplitudes have equal modulus
(`class_3_consensus_state_is_uniform_over_18_states`); its blank-record part
is the uniform double over the six hexagon ports, with weight `1/3`
(`its_blank_record_part_is_the_uniform_double_with_weight_1_3`); the tick
leaves it unchanged (`class_3_consensus_state_is_fixed_by_the_tick`). One
tick from the uniform double, with the class-3 sign, gives no emission `2/3`,
emission across an E-type seam `1/6`, across a C-type seam from the orbit
`{2, 9, 0}` (forward) `1/12` and from the other orbit (backward) `1/12`; two
independent objects count `2/9`, and the inverse tick restores the double
(`uniform_double_tick_weights_2_3_1_6_1_12_1_12_count_2_9`).

## 8. Many rounds on the sixty-address carrier

Two cells X and Y hold a blank double at port 2, and again at port 10, with independent registers over the sixty
addresses. The flags are the records of section 6: `R+` is record 2 of X and `H` is a written record of Y. The
event maps each state `|D, 0, k>` of the blank double to one basis state and back
(`event_maps_each_blank_double_address_to_one_basis_state`). The checks below are in `rounds.py`.

### The complete record of the return tests

The round of section 3 maps the standing pair to `|start> + (z - 1)|X>`, `z = exp(-i h j)`. Each cell is tested by
the projector on its uniform address state. `|X>` lies on blank doubles with `<start|X> = <X|X> = 2/9`
(`yes_amplitude_overlap_and_norm_are_2_9`). Its four test components have the weights

| X | Y | Weight in `X` |
| --- | --- | --- |
| yes | yes | `4/81` |
| no | yes | `11/324` |
| yes | no | `125/1296` |
| no | no | `55/1296` |

(`test_outcome_weights_4_81_11_324_125_1296_55_1296`), and the three "no" weights sum to `14/81 = P (1 - P)`
(`no_weights_sum_to_14_81`). The vector `|X>` does not depend on `j`. Each outcome therefore acts on the clocked
system as a function of `J`: yes on both cells as `1 - P + P z_J`, and an outcome with a "no" as
`(z_J - 1) sqrt(weight)`. A "no" removes the branch `j = 0`. With every outcome kept and the address register of a
cell prepared again after its "no", the coherence between the branches `j` and `j'` is multiplied per round by

```
(1 - P + P z)(1 - P + P/w) + P (1 - P)(z - 1)(1/w - 1) = 1 - P + P z/w,   z = exp(-i h j), w = exp(-i h j')
```

(`complete_record_multiplier_is_1_minus_P_plus_P_z_over_w`). With `N` rounds per unit clock time the multiplier is
`(1 - P + P exp(-i (j - j')/N))^N`, which tends to `exp(-i P (j - j'))`. The limit phase `P = 2/9` per unit clock
time holds on the complete record, with no selection of the returned branch. The expected number of rounds with a
reset per unit clock time is `N 4 sin^2(1/(2N)) 14/81` on a branch with `abs(j) = 1`.

### The round without the reset

A cell whose test said no keeps the blank double at its port with a non-uniform address. Its next tick reads the
flag of X with weight `11/12` and the flag of Y with weight `5/6`, the complements of `1/12` and `1/6`
(`not_returned_cells_read_11_12_and_5_6`). The count of the pair in the four test states (both returned, X not, Y
not, neither) is `2/9, 7/9, 7/9, 2/9` (`counts_without_reset_are_2_9_7_9_7_9_2_9`). Each single non-return changes
the count to `1 - 2/9`. The reset of the address consensus is part of the round.

### Partners that share the kick

With `n` partner cells and the generator `G = (1/n) sum_l Qx_l`, `Qx_l` the exclusive read of X and the partner `l`,
the yes amplitude of one round is

```
(11/12)(5/6 + u/6)^n + (1/12)(1/6 + 5u/6)^n,   u = exp(-i h j/n)
```

(`shared_kick_polynomial_is_the_product_form`; for two partners by a sum over all `60^3` address triples,
`two_partner_polynomial_by_enumeration_of_216000_triples`). The mean of `G` is `2/9` for `n = 1, 2, 3, 6`
(`shared_kick_generator_has_mean_2_9_for_1_2_3_6_partners`), so the phase per unit clock time tends to `2/9` for
every number of partners. Reading each partner in turn inside one round of X is the same operator, because the
kicks are diagonal in the records and the test of a partner commutes with every operation on the other cells.
Partners read at once, with `H` the event that at least one partner is written, give the count
`1/12 + (5/6)(1 - (5/6)^n)`: `2/9, 73/216, 563/1296` for `n = 1, 2, 3`
(`partners_read_at_once_count_73_216_and_563_1296`). The phases per unit clock time at `N = 1` are
`0.205334030970` (`n = 1`), `0.215732955808` (`n = 2`) and `0.217602641007` (`n = 3`)
(`shared_kick_phases_at_N_1_and_10`).

### Forward only on sixty addresses

Without the inverse tick the count `p + q - 2 p q`, `p` the weight of record 2 and `q` the weight of a written
record after `t` ticks, is `2/9`, `3544109/7290000`, `21968948741/41006250000` and
`1396558119030851/2989355625000000` at ticks 1 to 4 (`forward_only_counts_on_sixty_addresses`), that is `0.2222`,
`0.4862`, `0.5357`, `0.4672`.

### The register tick with record flags

On the two cells the analyzer of section 1 reads `phi_H` = Y written and X blank, `phi_+-` = record 2 (3) of X and
Y blank. Inside each branch all amplitudes are equal and positive. The compression is
`M = a I + b (z R + z^-1 R^-1)` with `3 a^2 = P(X blank) P(Y written)` and `3 b^2 = P(X record 2) P(Y blank)`:
`a^2 = 5/108`, `b^2 = 5/216`, `M = sqrt(5/18) C(theta)` and `tr M^dag M = 5/18`
(`double_load_register_has_a2_5_108_b2_5_216_trace_5_18`, `compression_on_sixty_addresses_is_sqrt_5_18_C` within
`1e-14`). For a native pair on two non-adjacent ports in each cell with the oriented record, `a^2 = 1/27`,
`b^2 = 1/54`, `M = (sqrt2/3) C(theta)` and `tr M^dag M = 2/9`
(`native_pair_register_has_a2_1_27_b2_1_54_trace_2_9`). Both have `(b/a)^2 = 1/2`
(`both_registers_are_balanced`). The round "tick, gate `exp(-i h |Phi><Phi|)`, inverse tick, return test" acts on
the register F in the returned branch as `1 + (exp(-i h) - 1) M^dag M`, which is diagonal on the three characters
(`mass_round_return_operator_is_1_plus_e_minus_1_M_dag_M`, within `1e-13`).

### Ten thousand rounds

The round is executed on complex arrays over the `60 x 60` addresses of the pair, with the projected state carried
from round to round. The first five rounds agree with the same rounds through `Walk.step` and `Walk.inverse_step`
on dictionaries within `1e-12`. After `10^4` rounds at `N = 100` the phase per unit clock time equals
`100 atan(P sin(1/100)/(1 - P + P cos(1/100))) = 0.222220621848` within `1e-11`, and the count read by every
checked round is `2/9` within `1e-14` (`ten_thousand_rounds_keep_the_phase_and_the_count`).

## Premises

The exact checks above use these declared constructions; the program does not
derive them.

- P1. The register: the C and E seam objects with their own loads, records
  and two-orientation addresses; the neutral register F carrying the
  third-turn `R` of `f`; the coupling `T`; the clock
  `K_theta = exp(i theta J)` on the oriented current; the oriented record
  `V`.
- P2. The preparation with activity `1/3` and the load at the screen tail,
  and the analyzer `Phi`, the Fourier zero over three flags.
- P3. The hexagon objects with two loads and their own twelve-address
  registers, and the standing state `S`.
- P4. The weak round: the kick on the exclusive read of fresh records, the
  inverse tick, a return test projecting on the standing state, N rounds per
  unit clock time.
- P5. The declared clock of `model/coherent_defect_wave/energy` with
  coefficient one on the onsite term and `pi/4` on the edge term. The
  coefficient one of the kick in that generator is part of this premise.
- P6. The oriented toggle `T(eta)` across the clocked register and the pair.
- P7. Two cells with independent sixty-address registers, each with two
  loads moved by the completion at the addressed seam.
- P8. The class-k flux as the link phase `w^k` on one hexagon seam.
- P9. The return test of a cell is the projector on its uniform address
  state, and a cell is prepared in that state again after a "no". With
  partners, the kick is the generator `(1/n) sum_l Qx_l`.

The public program executes one load per carrier; the two-load walks apply
its local completion at the addressed seam. A source operation that puts two
loads on one port is open, the coherent law keeps the sectors of section 7,
and the angle `eta` of the toggle per application, which sets the lepton
scale, is open.
