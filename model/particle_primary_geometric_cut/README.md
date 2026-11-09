# Regional entropy and boundary charges

An observer can restrict this finite source to the matter in a region and
its internal links, retaining all crossing links outside. This cut rotates
with the source cube, so its entropy, boundary-charge records and modular
response rotate together. The exact boundary representation term obeys
`S_edge <= ell * crossing_links * log(6)`. It is part of the regional
entropy, without an identification as proper area or gravitational entropy.

The [proof](PROOF.md) specifies the cut, actual Wilson-loop charge witnesses,
complete center records and the separate domains of regional energy and
entropy. Run `python -m model.particle_primary_geometric_cut --check` for
the source-bound receipt and `pytest model/particle_primary_geometric_cut`
for the finite algebra and covariance controls. These controls do not
allocate the full spatial density.
