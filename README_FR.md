# Observer Patch Holography (OPH)

> L'OPH part d'une idée simple : aucun observateur ne voit le monde entier d'un seul coup. Chaque observateur n'accède qu'à un patch local, et les patchs voisins doivent s'accorder sur leur recouvrement. L'OPH demande quelle part de la physique peut être reconstruite à partir de cette seule contrainte.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

L'OPH est un programme de reconstruction pour la physique fondamentale. Il part d'observateurs finis sur un écran holographique fini et progresse vers l'extérieur. Sa base de travail est algébrique-quantique : algèbres de patchs, états, probabilités de type trace/Born sur les surfaces d'enregistrement déclarées et entropie généralisée font partie du point de départ formel. Le programme ne cherche pas à dériver chaque ingrédient mathématique à partir de premiers principes. Son objectif est de construire une théorie du tout cohérente et complète en utilisant cette base d'information algébrique pour reconstruire l'univers effectif observé : espace-temps, structure de jauge, particules, enregistrements et synchronisation des observateurs y apparaissent comme des conséquences de la cohérence de recouvrement, et non comme des primitives.

La thèse opérationnelle est plus précise que "l'information est fondamentale".
L'OPH modélise la réalité comme un processus de consensus à point fixe fondé
sur les observateurs. Des patches d'observateurs finis portent des
enregistrements locaux, ne comparent que les données visibles dans leurs
recouvrements, réparent les désaccords par des mouvements de récupération
déclarés et convergent vers des points fixes stables qui survivent au
raffinement. Le monde public est la sortie stable par recouvrement de ce
processus. En ce sens, l'OPH traite la réalité comme un processus
computationnel, et non comme une scène statique sur laquelle le calcul se
produirait simplement.

## Par où commencer

Pour le noyau technique compact, commencez par **Paper 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**. **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** porte les dérivations particules. **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** développe l'image consensus-réparation. **Paper 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)** couvre l'architecture d'écran finie, les enregistrements et la machinerie observateur. Ce README, le Paper 1 et le livre sont les grandes vues d'ensemble.

## Ce que l'OPH apporte

- Un cadre à cutoff fixe pour les patches d'observateurs, les collerettes, la réparation de recouvrement, la jauge supérieure, les enregistrements et le checkpoint/restauration.
- Une reconstruction de la géométrie lorentzienne, du temps modulaire, de la dynamique d'Einstein de type Jacobson et de la cosmologie de Sitter en patch statique.
- Une reconstruction de jauge compacte vers le quotient réalisé du Modèle Standard `SU(3) x SU(2) x U(1) / Z_6`, avec le réseau exact des hypercharges, le triplet de couleur réalisé `N_c = 3` et le comptage des générations `N_g = 3`.
- Un programme particules à niveaux de revendication explicites : porteurs structurels exactement sans masse, ligne `W/Z` en comparaison seule, surface Higgs/top déclarée, masses de quarks sélectionnées avec Yukawas, branche neutrino pondérée, constantes locales candidates, et frontières no-go courantes pour les leptons chargés, la classification globale des cadres de quarks et le codomaine top direct auxiliaire.
- Une architecture microphysique d'écran concrète qui met mesure, enregistrements et observateurs à l'intérieur de la physique.

La mécanique quantique est traitée comme le langage algébrique d'information porté par l'architecture OPH. Le test de reconstruction est de savoir si cette base retrouve de manière cohérente l'univers effectif, et non si chaque ingrédient mathématique a été dérivé depuis un point de départ vide. Une reconstruction des espaces de Hilbert, des `C*`-algèbres ou des algèbres de von Neumann, des probabilités de Born, de la trace et de l'entropie depuis de simples enregistrements opérationnels serait une autre question, pas une condition préalable au programme.

Le mécanisme est la boucle de consensus à point fixe. Les observateurs locaux
n'accèdent pas à un état global depuis l'extérieur. Ils portent des états de
patch finis, échangent les données visibles dans les recouvrements, rejettent
les prolongements incohérents et conservent les motifs stables qui peuvent être
synchronisés. Géométrie, particules, lois et enregistrements sont les points
fixes à grande échelle de ce calcul en réseau d'observateurs.

L'OPH utilise une seule entrée quantitative externe, la capacité totale de l'écran `N_scr = log dim H_tot`, lue depuis l'horizon de Sitter, ainsi qu'un ratio local de pixel `P = a_cell / l_P^2`. Pour la constante cosmologique observée, le ratio nu d'aire d'horizon vaut `N_patch = (R_dS / l_P)^2 ≈ 1.05e122`, tandis que la capacité entropique utilisée par l'OPH vaut `N_scr = pi N_patch ≈ 3.31e122`. La même cellule d'écran est décrite deux fois : vue de l'extérieur comme un pixel situé légèrement au-dessus de l'équilibre auto-similaire exact `φ = (1 + sqrt(5)) / 2`, et vue de l'intérieur comme la plus petite échelle d'observation électromagnétique disponible pour les observateurs de l'univers simulé. `P` est la valeur pour laquelle ces deux descriptions coïncident.
Le papier de synthèse écrit cette autoréférence comme un problème de point fixe. L'OPH trouve la constante de structure fine en demandant quel décalage non nul d'une cellule d'écran holographique fait coïncider le déplacement géométrique extérieur depuis l'équilibre auto-similaire parfait avec l'échelle d'observation électromagnétique émise par l'univers vivant sur ce même écran. Le décalage extérieur vaut `α_ext(P) = (P - φ) / sqrt(pi)`, ou encore `P = φ + α_ext(P) sqrt(pi)`, tandis que le côté intérieur est le couplage électromagnétique émis par cette même cellule. Le tronc courant émet `P = 1.63097210492078846050203640439` et `α_cand⁻¹ = 136.994020662724205139718642793`, mais cette racine reste candidate tant que le point final Thomson projeté par Ward et le certificat RG/seuils ne sont pas fermés. La valeur CODATA/NIST `α⁻¹(0) = 137.035999177` reste donc une comparaison, non une ligne promue.

## Surface locale d'unification

L'OPH place une surface locale d'unification autour de l'entrée UV locale calibrée. La même échelle pilotée par le `P` candidat porte la voie électrofaible en comparaison, la voie Higgs ainsi que la voie entropique gravitationnelle, tandis que la géométrie lorentzienne fournit la vitesse causale invariante et que la couche locale de lecture fournit l'affichage SI. La présentation produit relevée donne `ellbar_shared = ellbar_SU(2) + ellbar_SU(3)` ; la même loi locale fixe `ellbar_shared = P/4`, et la lecture locale en unités SI est `G_SI = c^3 a_cell / (hbar P)` relativement au datum microscopique déclaré `a_cell`.
Sur la surface publique des constantes, `hbar` et `k_B` appartiennent à cette couche aval de lecture en unités familières plutôt que d'apparaître comme des constantes OPH émises de manière autonome.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260415" alt="Schéma d'unification OPH" width="92%">
  </a>
</p>

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
| Structure fine (inv.) | α⁻¹(0) | tronc candidat donnant 136.994020662724 | 137.035999177(21) | point final/certificat ouverts |
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
| Top | m_t, ligne section efficace | 172.35235532883115 GeV | 172.3523553288312 | match sélectionné |

`Δ` donne l'écart en sigma lorsque le PDG ou le NIST publie une incertitude à un sigma. Sinon, il
indique `match` ou `sous la borne`.

Pour les quarks, le PDG utilise ses conventions standard : `u`, `d` et `s` à `2 GeV`, et `c` et
`b` dans le schéma `MS` à leur propre échelle de masse.
Les papiers contiennent aussi les dérivations structurelles du Modèle Standard listées plus haut
ainsi qu'une famille neutrino, qui n'apparaît pas dans ce tableau faute de ligne de comparaison
PDG/NIST directe à un seul nombre.

La surface particules courante porte aussi des valeurs `W/Z` en comparaison seule
`80.377 GeV` et `91.18797809193725 GeV`, une valeur Higgs
`m_H = 125.1995304097179 GeV` et une valeur top sélectionnée
`m_t = 172.35235532883115 GeV` sur le codomaine section efficace `Q007TP4`.
La ligne top directe auxiliaire `Q007TP` reste en comparaison seule. Les masses
des leptons chargés sont des témoins ancrés sur cible, pas des prédictions
promues; l'issue #201 enregistre le no-go courant pour l'attachement manquant
à la ligne de déterminant. La branche neutrino pondérée émet
`(0.017454720257976796, 0.019481987935919015, 0.05307522145074924) eV` sur sa
branche déclarée. Les masses hadroniques demandent un backend hadronique OPH;
elles ne sont pas de simples entrées de quarks dans un tableau de comparaison
rapide. Les issues #153 et #157 sont fermées comme hors périmètre / bloquées
computationnellement, non comme résolues; il faut les rouvrir seulement
lorsqu'un backend hadronique OPH comme GLORB/Echosahedron produit des données
spectrales et des systématiques de production.

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse large de la pile OPH.
- **Papier 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique compact pour la relativité, la gravité et la structure du Modèle Standard.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivations particules, masses et couplages.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : formulation point fixe, réparation, et consensus.
- **Papier 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)** : architecture d'écran finie, enregistrements, et machinerie observateur.

## Articles supplémentaires

- **[Observer-Patch Holography and the Dark Matter Phenomenon](extra/oph_dark_matter_paper.pdf)** : note de recherche supplémentaire sur la phénoménologie de la matière noire, la limite galactique de type MOND, les diagnostics SPARC, le transport par réparation et les tests cosmologiques compressés. Le code de reproductibilité vit dans [code/dark_matter](code/dark_matter/).

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
