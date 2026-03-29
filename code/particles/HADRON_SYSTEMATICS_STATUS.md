# Hadron Pipeline Status

This file tracks the execution and publication contract for the active hadron
prediction path separately from the main particle table.

## Execution Receipt
- artifact: `runtime_schedule_receipt_N_therm_and_N_sep`
- receipt_kind: `non_theorem_external_runtime_receipt`
- status: `receipt_filled_waiting_backend_dump`
- `N_therm`: `2048`
- `N_sep`: `512`
- kernel_id: `fixed_schedule_rhmc_hmc`
- initial_configuration: `U_cold = 1`
- next_single_residual_object_after_execution: `oph_hadron_stable_channel_sequence_evaluator`

## Stable-Channel Execution State
| ensemble_id | cfg_ids | n_cfg | n_src_per_cfg | t_extent | arrays_written |
| --- | --- | ---: | ---: | ---: | --- |
| `qcd_2p1_seed_n0` | `qcd_2p1_seed_n0__cfg0,qcd_2p1_seed_n0__cfg1` | `2` | `2` | `291` | `no` |
| `qcd_2p1_seed_n1` | `qcd_2p1_seed_n1__cfg0,qcd_2p1_seed_n1__cfg1` | `2` | `2` | `581` | `no` |
| `qcd_2p1_seed_n2` | `qcd_2p1_seed_n2__cfg0,qcd_2p1_seed_n2__cfg1` | `2` | `2` | `1162` | `no` |

## Stable-Channel Numerical State
| ensemble_id | channel | forward_window_limit_exists | am_ground_candidate | stat_err | sys_err | mass_gev_candidate |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `qcd_2p1_seed_n0` | `pi_iso` | `False` | `None` | `None` | `None` | `None` |
| `qcd_2p1_seed_n0` | `N_iso` | `False` | `None` | `None` | `None` | `None` |
| `qcd_2p1_seed_n1` | `pi_iso` | `False` | `None` | `None` | `None` | `None` |
| `qcd_2p1_seed_n1` | `N_iso` | `False` | `None` | `None` | `None` | `None` |
| `qcd_2p1_seed_n2` | `pi_iso` | `False` | `None` | `None` | `None` | `None` |
| `qcd_2p1_seed_n2` | `N_iso` | `False` | `None` | `None` | `None` | `None` |

## Published Systematics Budget
| ensemble_id | channel | statistics | continuum | volume | chiral | publishable |
| --- | --- | --- | --- | --- | --- | --- |
| `qcd_2p1_seed_n0` | `pi_iso` | `pending` | `pending` | `pending` | `pending` | `no` |
| `qcd_2p1_seed_n0` | `N_iso` | `pending` | `pending` | `pending` | `pending` | `no` |
| `qcd_2p1_seed_n1` | `pi_iso` | `pending` | `pending` | `pending` | `pending` | `no` |
| `qcd_2p1_seed_n1` | `N_iso` | `pending` | `pending` | `pending` | `pending` | `no` |
| `qcd_2p1_seed_n2` | `pi_iso` | `pending` | `pending` | `pending` | `pending` | `no` |
| `qcd_2p1_seed_n2` | `N_iso` | `pending` | `pending` | `pending` | `pending` | `no` |

## Residual Open Object
- current: `backend_correlator_dump.production.json from real production RHMC/HMC execution on the theorem-emitted seeded family`
- after-execution: `StableChannelForwardWindowConvergence`
- after-convergence: `oph_hadron_rho_scattering_readout`

## Working Rule
Keep hadron rows out of the public particle table until the active pipeline
has a real production backend dump plus published statistical and
continuum/volume/chiral systematics on the seeded `2+1`, QED-off stable-channel
branch.
