# Native class freezeout and boundary response

The unchanged native population law gives exact probabilities for all 46
terminal classes, their conditional lifetimes and the full six-class
limiting distribution. Uniform zero-flux input relaxes in an average
`5557709/1125000` native attempts; the largest mean over all input states
is 15. Complete finite populations quantify the response to declared
asymmetric preparation and reset controls.

The [proof](PROOF.md) derives the equations and uniform error bounds. The
[receipt](results.json) keeps every state and terminal class. The word
freezeout means cessation of native class changes; identifying a thermal
cosmic epoch requires separate physical assumptions.

```sh
python3 -m pytest -q model/native_freezeout_response
python3 -m model.native_freezeout_response --write
python3 -m model.native_freezeout_response --check
```
