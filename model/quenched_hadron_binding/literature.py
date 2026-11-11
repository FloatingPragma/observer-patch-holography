"""Imported numbers: measured masses, one convention and one published lattice value of the same theory. None enters
a fit, except the measured ratio m_pi / m_rho, which fixes the light quark mass at the physical point."""

MEASURED = {
    'source': 'Particle Data Group, Phys. Rev. D 110, 030001 (2024), https://pdg.lbl.gov/',
    'm_nucleon_MeV': 938.919,          # mean of proton 938.272 and neutron 939.565
    'm_rho_MeV': 775.26,
    'm_pion_MeV': 138.039,             # isospin average (2 x 139.570 + 134.977) / 3
    'mN_over_mrho': 1.2111,
    'mpi_over_mrho': 0.1780,
}

CONVENTION = {
    'source': 'conventional string tension from Regge slopes and quarkonia, about 5 per cent wide',
    'sqrt_sigma_MeV': 440.0,
    'mrho_over_sqrtsigma': 1.762,
    'mN_over_sqrtsigma': 2.134,
}

LITERATURE = {
    'plaquette_beta_5.7': {
        'value': 0.5492, 'source': 'G. P. Lepage and P. B. Mackenzie, Phys. Rev. D 48, 2250 (1993), arXiv:hep-lat/9209022'},
}
