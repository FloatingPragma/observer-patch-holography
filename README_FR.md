# Observer Patch Holography (OPH)

> L'OPH part d'une idée simple : aucun observateur ne voit le monde entier d'un seul coup. Chaque observateur n'accède qu'à un patch local, et les patchs voisins doivent s'accorder sur leur recouvrement. L'OPH demande quelle part de la physique peut être reconstruite à partir de cette seule contrainte.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

L'OPH est un programme de reconstruction. Espace-temps, structure de jauge, particules, enregistrements et synchronisation des observateurs y apparaissent comme des conséquences de la cohérence de recouvrement sur un écran holographique fini.

## Par où commencer

Pour le noyau technique compact, commencez par **Paper 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**. **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** porte les dérivations particules. **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** développe l'image consensus-réparation. **Paper 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)** couvre l'architecture d'écran finie, les enregistrements et la machinerie observateur. Ce README, le Paper 1 et le livre sont les grandes vues d'ensemble.

## Ce que l'OPH apporte

- Un paquet à cutoff fixe pour les patches d'observateurs, les collerettes, la réparation de recouvrement, la jauge supérieure, les enregistrements et le checkpoint/restauration.
- Une reconstruction de la géométrie lorentzienne, du temps modulaire, de la dynamique d'Einstein de type Jacobson et de la cosmologie de Sitter en patch statique.
- Une reconstruction de jauge compacte vers le quotient réalisé du Modèle Standard `SU(3) x SU(2) x U(1) / Z_6`, avec le réseau exact des hypercharges, le triplet de couleur réalisé `N_c = 3` et le comptage des générations `N_g = 3`.
- Un programme particules avec porteurs structurels exactement sans masse, secteur électrofaible, surface Higgs/top, masses de quarks avec Yukawas, structure neutrino et constantes locales.
- Une architecture microphysique d'écran concrète qui met mesure, enregistrements et observateurs à l'intérieur de la physique.

L'OPH utilise une seule entrée quantitative externe, la capacité totale de l'écran `N_scr = log dim H_tot`, lue depuis l'horizon de Sitter, ainsi qu'un ratio local de pixel `P = a_cell / l_P^2`. Pour la constante cosmologique observée, le ratio nu d'aire d'horizon vaut `N_patch = (R_dS / l_P)^2 ≈ 1.05e122`, tandis que la capacité entropique utilisée par l'OPH vaut `N_scr = pi N_patch ≈ 3.31e122`. La même cellule d'écran est décrite deux fois : vue de l'extérieur comme un pixel situé légèrement au-dessus de l'équilibre auto-similaire exact `φ = (1 + sqrt(5)) / 2`, et vue de l'intérieur comme la plus petite échelle d'observation électromagnétique disponible pour les observateurs de l'univers simulé. `P` est la valeur pour laquelle ces deux descriptions coïncident.
Le papier de synthèse écrit cette autoréférence comme un problème de point fixe. L'OPH trouve la constante de structure fine en demandant quel décalage non nul d'une cellule d'écran holographique fait coïncider le déplacement géométrique extérieur depuis l'équilibre auto-similaire parfait avec l'échelle d'observation électromagnétique émise par l'univers vivant sur ce même écran. Le décalage extérieur vaut `α_ext(P) = (P - φ) / sqrt(pi)`, ou encore `P = φ + α_ext(P) sqrt(pi)`, tandis que le côté intérieur est le couplage électromagnétique émis par cette même cellule. Pour la valeur centrale externe 2022 `α⁻¹(0) = 137.035999177`, cette formule extérieure donne `P = 1.630968209403959...`. La même géométrie de point fixe est aussi testée par une note matérielle en cavité optique. Sur cette petite surface quantitative, l'OPH émet des prédictions concrètes pour les couplages, les masses et les grandeurs gravitationnelles.

## Surface locale d'unification

L'OPH place une surface locale d'unification autour de l'entrée UV locale calibrée. La même échelle pilotée par `P` porte la voie bosonique électrofaible et Higgs ainsi que la voie entropique gravitationnelle, tandis que la géométrie lorentzienne fournit la vitesse causale invariante et que le paquet local de lecture fournit l'affichage SI. La présentation produit relevée donne `ellbar_shared = ellbar_SU(2) + ellbar_SU(3)` ; la même loi locale fixe `ellbar_shared = P/4`, et la lecture locale en unités SI est `G_SI = c^3 a_cell / (hbar P)` relativement au datum microscopique déclaré `a_cell`.
Sur la surface publique des constantes, `hbar` et `k_B` appartiennent à cette couche aval de lecture en unités familières plutôt que d'apparaître comme des constantes OPH émises de manière autonome.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260415" alt="Schéma d'unification OPH" width="92%">
  </a>
</p>

Les sorties particules détaillées vivent dans [code/particles/RESULTS_STATUS.md](code/particles/RESULTS_STATUS.md) et [code/particles/EXACT_NONHADRON_MASSES.md](code/particles/EXACT_NONHADRON_MASSES.md).

**Pile générale des théorèmes et dérivations**

<p align="center">
  <a href="assets/prediction-chain.svg?v=20260412" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg?v=20260412" alt="Pile théorématique et de dérivation OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La pile OPH complète, des axiomes jusqu'à la relativité, la structure de jauge, les particules et les observateurs. Cliquez pour ouvrir le SVG complet.</sub></p>

## Lignes quantitatives sélectionnées

Ce tableau condensé garde les lignes OPH les plus faciles à comparer directement avec les valeurs
PDG/NIST. Les résultats structurels comme la géométrie lorentzienne `3+1D`, le quotient de jauge
du Modèle Standard `SU(3) x SU(2) x U(1) / Z_6`, le réseau exact des hypercharges, le triplet
de couleur réalisé `N_c = 3` et le comptage des générations `N_g = 3` sont énoncés dans les papiers et ne sont pas répétés ici.

| Quantité | Symbole | OPH | PDG/NIST | Δ |
| --- | --- | --- | --- | --- |
| Constante gravitationnelle | G | 6.6742999959e-11 | 6.67430(15)e-11 | 0.00003σ |
| Vitesse de la lumière | c | 299792458 | 299792458 (exact) | match |
| Structure fine (inv.) | α⁻¹(0) | fermeture `P` donnant 137.035999177 | 137.035999177(21) | match |
| Masse du photon | m_γ | 0 eV | <1e-18 eV | sous la borne |
| Masse du gluon | m_g | 0 GeV | 0 GeV | match |
| Masse du graviton | m_grav | 0 eV | <1.76e-23 eV | sous la borne |

**Secteur des quarks**

| Quark | Symbole | OPH | PDG | Δ |
| --- | --- | --- | --- | --- |
| Bottom | m_b(m_b) | 4.183 GeV | 4.183 ± 0.007 | match |
| Charm | m_c(m_c) | 1.273 GeV | 1.2730 ± 0.0046 | match |
| Strange | m_s(2 GeV) | 93.5 MeV | 93.5 ± 0.8 | match |
| Down | m_d(2 GeV) | 4.70 MeV | 4.70 ± 0.07 | match |
| Up | m_u(2 GeV) | 2.16 MeV | 2.16 ± 0.07 | match |

`Δ` donne l'écart en sigma lorsque le PDG ou le NIST publie une incertitude à un sigma. Sinon, il
indique `match` ou `sous la borne`.

Pour les quarks, le PDG utilise ses conventions standard : `u`, `d` et `s` à `2 GeV`, et `c` et
`b` dans le schéma `MS` à leur propre échelle de masse.
Les papiers contiennent aussi les dérivations structurelles du Modèle Standard listées plus haut
ainsi qu'une famille neutrino, qui n'apparaît pas dans ce tableau faute de ligne de comparaison
PDG/NIST directe à un seul nombre.

La surface électrofaible porte aussi une valeur Higgs `m_H = 125.1995304097179 GeV`
et une valeur top compagnon `m_t = 172.3523553288312 GeV`.
Les masses hadroniques demandent la dynamique QCD non perturbative des états liés; elles ne sont
pas de simples entrées de quarks dans un tableau de comparaison rapide.

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse large de la pile OPH.
- **Papier 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique compact pour la relativité, la gravité et la structure du Modèle Standard.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivations particules, masses et couplages.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : formulation point fixe, réparation, et consensus.
- **Papier 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)** : architecture d'écran finie, enregistrements, et machinerie observateur.

## Plus

- **Site officiel :** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Page theory of everything :** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Page simulation theory :** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Livre :** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Application d'étude guidée :** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions et explications détaillées :** OPH Sage sur [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage) ou [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab :** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Objections courantes :** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **Note IBM Quantum :** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Guide du dépôt

- **[`paper/`](paper)** : PDF, sources LaTeX et métadonnées de release.
- **[`book/`](book)** : source du livre OPH.
- **[`code/`](code)** : sorties calculatoires, surface particules et expériences.
- **[`assets/`](assets)** : schémas et figures publics.
- **[`extra/`](extra)** : notes publiques maintenues, objections, comptes rendus expérimentaux et quelques essais de support.
