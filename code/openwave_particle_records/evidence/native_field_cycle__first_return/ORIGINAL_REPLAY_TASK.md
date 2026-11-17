# Archived failed-gate replay

Read the exact archived probe bytes and execute them with the original
logical file location, so unchanged provider paths resolve correctly.
Override only check_pins to verify the archived original freeze and all
229 original source pins. No scientific operation, q0 gate or preparation
is changed. Catch the original AssertionError and compare its raw
traceback-local diagnostic exactly with the archived FAILURE_GATE.log.
Require that original gate failure to reproduce; success of the original
gate would fail this replay. This script never writes evidence.

Freeze the reader/task and archived code/task/design/pins/prereg/log before
execution. The amended headline is a separate read-only replay and does
not call or change this historical diagnostic.
