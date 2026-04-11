# Observer Patch Holography (OPH)

> L'OPH part d'une idée simple : aucun observateur ne voit le monde entier d'un seul coup. Chaque observateur n'accède qu'à un patch local, et les patchs voisins doivent s'accorder sur leur recouvrement. L'OPH demande quelle part de la physique peut être reconstruite à partir de ce point de départ une fois le ledger complet des axiomes et des branches rendu explicite.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

L'OPH est un programme de reconstruction. Espace-temps, structure de jauge, particules, enregistrements et synchronisation des observateurs y apparaissent comme des conséquences du paquet OPH enraciné dans la cohérence de recouvrement sur un écran holographique fini, avec les prémisses de branche explicites énoncées dans les papiers.

## Ce que l'OPH apporte

- Un paquet théorématique à cutoff fixe pour les patches d'observateurs, les collerettes, la réparation de recouvrement, la jauge supérieure, les enregistrements et le checkpoint/restauration.
- Une voie conditionnelle vers la géométrie lorentzienne, le temps modulaire, la dynamique d'Einstein de type Jacobson et la cosmologie de Sitter en patch statique sur le sous-réseau géométrique premier extrait ; la branche d'Einstein utilise la stationnarité à cap fixe, le pont modulaire sur les surfaces nulles et la branche projective séparée sur intervalles bornés, tandis que le scaffold UV/BW restant est la réalisation de la paire de caps géométriques sur ce sous-réseau puis la rigidité des paires de coupures ordonnées, avec le plancher commun éventuel de transport modulaire sur collerette locale fixe comme plus petit bloqueur inférieur.
- Une voie conditionnelle de jauge compacte dans la branche bosonique vers le quotient réalisé du Modèle Standard `SU(3) x SU(2) x U(1) / Z_6`, sous les prémisses de reconstruction par secteurs transportables et sous MAR, avec le réseau exact des hypercharges et la chaîne de comptage réalisée `N_g = 3`, `N_c = 3`.
- Un programme particules avec porteurs structurels exactement sans masse, une branche de calibration électrofaible de Phase II émise vers l'avant avec une surface théorématique publique `W/Z` target-free fermée plus une paire gelée exacte utilisée seulement comme validation compare-only, un étage quantitatif Higgs/top, des sidecars exacts non hadroniques et des voies de continuation explicites pour la saveur et les hadrons.
- Une architecture microphysique d'écran concrète qui met mesure, enregistrements et observateurs à l'intérieur de la physique.

## Surface locale d'unification

L'OPH place une surface locale d'unification autour de l'entrée UV locale calibrée. La même échelle pilotée par `P` porte la voie bosonique électrofaible et Higgs ainsi que la voie entropique gravitationnelle, tandis que la branche lorentzienne fournit la vitesse causale invariante et que le paquet local de lecture fournit l'affichage SI.
Sur la surface publique actuelle des constantes, `hbar` et `k_B` restent dans cette couche aval de lecture en unités familières plutôt que d'apparaître comme des constantes OPH émises de manière autonome.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260407" alt="Schéma d'unification OPH" width="92%">
  </a>
</p>

Les constantes, chaînes de théorèmes et fronts de preuve ouverts pour cette surface sont suivis dans [extra/OPH_PHYSICS_CONSTANTS.md](extra/OPH_PHYSICS_CONSTANTS.md).

**Pile générale des théorèmes et dérivations**

<p align="center">
  <a href="assets/prediction-chain.svg?v=20260412" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg?v=20260412" alt="Pile théorématique et de dérivation OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La pile OPH complète, des axiomes jusqu'à la relativité, la structure de jauge, les particules, les observateurs et les fronts encore ouverts. Cliquez pour ouvrir le SVG complet.</sub></p>

## Points forts côté particules

### Résultats théorématiques et structurels

- Zéros structurels exacts pour le photon, les gluons et le graviton.
- Sortie électrofaible sur la branche de calibration D10, avec lignes publiques `W/Z` target-free fermées et paire gelée exacte utilisée seulement comme validation compare-only
  `W = 80.377 GeV`, `Z = 91.18797809193725 GeV`.
- Étage quantitatif Higgs/top en aval du coeur électrofaible, avec une graine forward scalaire unique fermée qui porte les lignes publiques
  `H = 125.218922 GeV`, `t = 172.388646 GeV`.

### Surface exacte non hadronique

| Voie | Sortie exacte | Note de statut |
| --- | --- | --- |
| Porteurs structurels | `m_photon = m_gluon = m_graviton = 0` | exactitude structurelle de rang théorème |
| Sidecar électrofaible | `W = 80.377 GeV`, `Z = 91.18797809193725 GeV` | surface de réparation gelée exacte |
| Sidecar exact Higgs/top | `(H, t) = (125.1995304097179, 172.3523553288311) GeV` | tranche inverse exacte mais compare-only sur le même Jacobien Higgs/top |
| Témoin chargé | `(e, mu, tau) = (0.00051099895, 0.1056583755, 1.7769324651340912) GeV` | témoin exact sur une chaîne de lecture quadratique fermée à trois points |
| Témoin quark | `(u, d, s, c, b, t) = (0.00216, 0.00470, 0.0935, 1.273, 4.183, 172.3523553288311) GeV` | témoin exact sur une chaîne de lecture quadratique fermée à trois points |
| Branche théorème neutrino | `(m1, m2, m3) = (0.017454720257976796, 0.019481987935919015, 0.05307522145074924) eV` avec `Δm21²`, `Δm31²`, `Δm32²` émis sur la branche à cycle pondéré | famille absolue à cycle pondéré de rang théorème |

Les lignes publiques Higgs/top sont portées par la graine forward scalaire unique fermée. La paire inverse exacte ci-dessus reste un sidecar compare-only sur le même Jacobien et ne remplace pas la branche forward publique.

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

### Succès de continuation

- La voie de continuation quark émet des lignes publiques pour `u`, `d`, `s`, `c` et `b` sur la branche de continuation sélectionnée.
- La branche neutrino à cycle pondéré atteint le régime PMNS et hiérarchie observé avec
  `theta12 = 34.2259°`, `theta23 = 49.7228°`, `theta13 = 8.68636°`, `delta = 305.581°`,
  et `Δm21² / Δm32² = 0.03072111`.
- La surface exacte non hadronique est regroupée dans
  [code/particles/EXACT_NONHADRON_MASSES.md](code/particles/EXACT_NONHADRON_MASSES.md).

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : papier de synthèse de l'ensemble OPH.
- **Papier 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : papier de dérivation relativité/structure du Modèle Standard.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivation des particules, surface exacte, et carte des continuations.
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
