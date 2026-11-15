#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
charge_packet="$(cd "$(dirname "$0")" && pwd)"
python3 -u "$charge_packet/embedding.py" > "$charge_packet/EMBEDDING_RUN.log" 2>&1
python3 -u "$charge_packet/topology.py" > "$charge_packet/TOPOLOGY_RUN.log" 2>&1
python3 -u "$charge_packet/dynamics.py" > "$charge_packet/DYNAMICS_RUN.log" 2>&1
python3 -u "$charge_packet/native_flux.py" > "$charge_packet/NATIVE_RUN.log" 2>&1
python3 "$charge_packet/freeze.py"
python3 "$charge_packet/freeze.py" --verify
