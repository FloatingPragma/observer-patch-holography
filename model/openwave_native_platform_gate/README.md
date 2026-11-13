# Headless native-process review

```
python -m model.openwave_native_platform_gate.headless
python -m model.openwave_native_platform_gate --write --openwave-root /path/to/pinned/openwave
python -m model.openwave_native_platform_gate --check --openwave-root /path/to/pinned/openwave
python -m pytest -q model/openwave_native_platform_gate
```

Use the exact revision in `model/openwave_m10_review_gate/contract.json`.
This author-side reproduction has no official score credit and sends no
upstream submission. The source-to-physical-observable join remains explicit.
