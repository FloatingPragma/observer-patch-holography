# Exact native response

The complete one-load C/E process has population response I-L/(2M) and
explicit real coherence decay factors on 240, 960 and 3840 source ports.
All endpoint operator columns and all global one-load columns are covered.
The natural cell reduction has exact l1 residual 11/(12M); retaining its
port-memory complement gives a closed finite convolution identity.

This establishes the response of the actual repair process. It does not
identify that response with the separately declared M10.1 cosine operator.

Run `python -m model.openwave_native_response --write`, then
`python -m model.openwave_native_response --check`.
