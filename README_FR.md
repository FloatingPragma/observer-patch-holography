# Observer Patch Holography (OPH)

> L'OPH part d'une idée simple : aucun observateur ne voit le monde entier d'un seul coup. Chaque observateur n'accède qu'à un patch local, et les patchs voisins doivent s'accorder sur leur recouvrement. L'OPH demande quelle part de la physique peut être reconstruite à partir de ce point de départ une fois le ledger complet des axiomes et des branches rendu explicite.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

L'OPH est un programme de reconstruction. Espace-temps, structure de jauge, particules, enregistrements et synchronisation des observateurs y apparaissent comme des conséquences du paquet OPH enraciné dans la cohérence de recouvrement sur un écran holographique fini, avec les prémisses de branche explicites énoncées dans les papiers.

## Règle d'autorité et de lecture

Pour le statut théorématique du noyau reconstruit et le niveau de preuve des claims, consultez d'abord **Paper 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**. Le niveau de preuve propre à chaque voie reste dans les papiers compagnons correspondants, notamment **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**, **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** et **Paper 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)**. Ce README, le Paper 1 et le livre sont des surfaces de synthèse synchronisées : ils résument et organisent les résultats sans en rehausser le niveau de preuve.

## Ce que l'OPH apporte

- Un paquet théorématique à cutoff fixe pour les patches d'observateurs, les collerettes, la réparation de recouvrement, la jauge supérieure, les enregistrements et le checkpoint/restauration.
- Une voie conditionnelle vers la géométrie lorentzienne, le temps modulaire, la dynamique d'Einstein de type Jacobson et la cosmologie de Sitter en patch statique sur le sous-réseau géométrique premier extrait ; la branche d'Einstein utilise la stationnarité à cap fixe, le pont modulaire sur les surfaces nulles et la branche projective séparée sur intervalles bornés, tandis que le scaffold UV/BW restant est la réalisation de la paire de caps géométriques sur ce sous-réseau puis la rigidité des paires de coupures ordonnées, avec le plancher commun éventuel de transport modulaire sur collerette locale fixe comme plus petit bloqueur inférieur.
- Une voie conditionnelle de jauge compacte dans la branche bosonique vers le quotient réalisé du Modèle Standard `SU(3) x SU(2) x U(1) / Z_6`, sous les prémisses de reconstruction par secteurs transportables et sous MAR, avec le réseau exact des hypercharges et la chaîne de comptage réalisée `N_g = 3`, `N_c = 3`.
- Un programme particules avec porteurs structurels exactement sans masse, une branche électrofaible de fermeture quantitative de Phase II émise vers l'avant avec une surface théorématique publique `W/Z` target-free fermée plus une paire gelée exacte utilisée seulement comme validation compare-only, un théorème exact de séparation Higgs/top à source seule sur la surface électrofaible déclarée de running, matching et seuils avec une tranche inverse exacte Higgs/top gardée comme validation compare-only, une fermeture quark exacte sur classe publique sélectionnée avec Yukawas forward exactes explicites, des surfaces exactes non hadroniques et des voies de continuation explicites aux frontières théorématiques ouvertes.
- Une architecture microphysique d'écran concrète qui met mesure, enregistrements et observateurs à l'intérieur de la physique.

L'OPH utilise une seule entrée quantitative externe, la capacité totale de l'écran `N_scr = log dim H_tot`, déduite de la constante cosmologique, ainsi qu'un ratio local de pixel `P = a_cell / l_P^2`. La même cellule d'écran est décrite deux fois : vue de l'extérieur comme un pixel situé légèrement au-dessus de l'équilibre auto-similaire exact `φ = (1 + sqrt(5)) / 2`, et vue de l'intérieur comme la plus petite échelle d'observation électromagnétique disponible pour les observateurs de l'univers simulé. `P` est la valeur pour laquelle ces deux descriptions coïncident.
Le papier de synthèse écrit cette autoréférence comme un problème de point fixe. Le décalage extérieur vaut `α_ext(P) = (P - φ) / sqrt(pi)`, ou encore `P = φ + α_ext(P) sqrt(pi)`, tandis que le côté intérieur est le couplage électromagnétique émis par cette même cellule sur la branche quantitative déclarée. L'équation de fermeture affirme que ces deux lectures d'un même objet coïncident. La preuve par point fixe montre que cette fermeture sélectionne un `P` admissible sur l'intervalle physique. Sur cette petite surface quantitative, l'OPH émet des prédictions concrètes pour les couplages, les masses et les grandeurs gravitationnelles.

## Surface locale d'unification

L'OPH place une surface locale d'unification autour de l'entrée UV locale calibrée. La même échelle pilotée par `P` porte la voie bosonique électrofaible et Higgs ainsi que la voie entropique gravitationnelle, tandis que la branche lorentzienne fournit la vitesse causale invariante et que le paquet local de lecture fournit l'affichage SI. Sur la surface locale d'extension déclarée, la présentation produit relevée de la branche quotient réalisée donne `ellbar_shared = ellbar_SU(2) + ellbar_SU(3)` ; la même loi D10 sur cette surface fixe `ellbar_shared = P/4`, et la lecture locale en unités SI est `G_SI = c^3 a_cell / (hbar P)` relativement au datum microscopique déclaré `a_cell`.
Sur la surface publique des constantes, `hbar` et `k_B` restent dans cette couche aval de lecture en unités familières plutôt que d'apparaître comme des constantes OPH émises de manière autonome.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260415" alt="Schéma d'unification OPH" width="92%">
  </a>
</p>

Les surfaces de statut pour les particules dans ce dépôt sont [code/particles/RESULTS_STATUS.md](code/particles/RESULTS_STATUS.md) et [code/particles/EXACT_NONHADRON_MASSES.md](code/particles/EXACT_NONHADRON_MASSES.md).

**Pile générale des théorèmes et dérivations**

<p align="center">
  <a href="assets/prediction-chain.svg?v=20260412" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg?v=20260412" alt="Pile théorématique et de dérivation OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La pile OPH complète, des axiomes jusqu'à la relativité, la structure de jauge, les particules, les observateurs et les fronts encore ouverts. Cliquez pour ouvrir le SVG complet.</sub></p>

## Lignes quantitatives sélectionnées

Ce tableau condensé ne garde que les lignes OPH avec égalité exacte, accord en sigma publié,
respect clair d'une borne supérieure par rapport aux valeurs de référence PDG/NIST utilisées dans
les papiers, ou une ligne explicitement marquée comme candidate de fermeture de Phase II. Les
résultats structurels comme la branche lorentzienne `3+1D`, le quotient de jauge
du Modèle Standard `SU(3) x SU(2) x U(1) / Z_6`, le réseau exact des hypercharges et la chaîne de
comptage `N_g = 3`, `N_c = 3` sont énoncés dans les papiers et ne sont pas répétés ici. La voie
bosonique `W/Z/H` appartient à la branche de fermeture quantitative de phase II; elle est donc discutée dans
les papiers mais omise de ce tableau rapide. La ligne de structure fine reste une candidate de
fermeture `P` de Phase II tant que la voie de transport à source seule n'a pas publié son
certificat de fermeture lisible par machine.

| Quantité | Symbole | OPH | PDG/NIST | Δ |
| --- | --- | --- | --- | --- |
| Constante gravitationnelle | G | 6.6742999959e-11 | 6.67430(15)e-11 | 0.00003σ |
| Vitesse de la lumière | c | 299792458 | 299792458 (exact) | match |
| Structure fine (inv.) | α⁻¹(0) | candidate de fermeture `P` de Phase II | 137.035999177(21) | certificat en attente |
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
indique `match`, `sous la borne`, ou le statut d'une ligne candidate déclarée.

Pour les quarks, le PDG utilise ses conventions standard : `u`, `d` et `s` à `2 GeV`, et `c` et
`b` dans le schéma `MS` à leur propre échelle de masse.
Les papiers contiennent aussi les dérivations structurelles du Modèle Standard listées plus haut
ainsi qu'une famille neutrino de rang théorème, qui n'apparaissent pas dans ce tableau faute de
ligne de comparaison PDG/NIST directe à un seul nombre.
La surface neutrino publique inclut aussi une paire physique de phases de Majorana de
rang théorème sur la branche de transport weighted-cycle en base partagée ; voir
`code/particles/RESULTS_STATUS.md`.

La surface électrofaible déclarée de fermeture quantitative porte aussi un théorème exact à source seule pour
le Higgs, avec `m_H = 125.1995304097179 GeV`, ainsi qu'une coordonnée compagnon pour le top
`m_t = 172.3523553288312 GeV` sur la même surface jacobienne.
À la précision effectivement publiée par le PDG, la ligne Higgs tombe sur la moyenne 2025.
La ligne top publique exacte sur la surface quark sélectionnée utilise l'entrée PDG 2025 en masse
extraite de section efficace `Q007TP4`.
Le pont vers la moyenne auxiliaire en mesures directes
`Q007TP = 172.56 ± 0.31 GeV` reste ouvert et est suivi dans
[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207).

Le secteur des leptons chargés suit un niveau de preuve distinct. Le dépôt contient un témoin exact sur une même famille, la lecture des mêmes labels `q_e`, un caractère déterminant côté source défini pour un vecteur de multiplicité fixé, un relèvement conditionnel de la ligne déterminant sur les données chargées physiques, puis une lecture algébrique des masses à partir d'une échelle absolue chargée de rang théorème. La voie théorématique publique n'émet ni vecteur exposant déterminant sectorialisé de rang théorème, ni identification du caractère déterminant côté source avec la ligne déterminant physique chargée. Les masses publiques de l'électron, du muon et du tau ne sont donc pas émises depuis `P`.

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : papier de synthèse de l'ensemble OPH ; il organise la suite sur une seule surface et hérite du niveau de preuve du papier compact sur le noyau reconstruit ainsi que des ledgers des papiers compagnons correspondants, au lieu de les rehausser.
- **Papier 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : surface faisant autorité pour le noyau reconstruit et le niveau de preuve sur la chaîne Lorentz/gravité et la branche structurelle réalisée du Modèle Standard.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivation des particules, surface exacte, et carte des frontières théorématiques.
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
