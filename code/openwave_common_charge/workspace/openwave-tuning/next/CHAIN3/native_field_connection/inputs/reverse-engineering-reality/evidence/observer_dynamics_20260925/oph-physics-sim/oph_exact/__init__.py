"""Exact-carrier lanes of the OPH simulator.

Every lane in this package runs the declared twelve-port icosahedral carrier
with the canonical seam-mean repair law, ``T = I - L/60``, and reports exact
receipts.  The package sits beside ``oph_fpe`` so that the CR-0 python-tree
snapshot of ``oph_fpe`` is untouched by lane additions.  Lanes import from
``oph_fpe`` where an exact object is already committed there.
"""
