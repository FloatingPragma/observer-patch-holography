# Observer Patch Holography (OPH)

> L'OPH part d'une idée simple : aucun observateur ne voit le monde entier d'un seul coup. Chaque observateur n'accède qu'à un patch local, et les patchs voisins doivent s'accorder sur leur recouvrement. L'OPH demande quelle part de la physique peut être reconstruite à partir de cette seule contrainte.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [OPH Textbooks](https://learn.floatingpragma.io/) | [Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [site](https://floatingpragma.io/oph/) | [OPH Lab](https://oph-lab.floatingpragma.io)

L'OPH est un programme de reconstruction pour la physique fondamentale. Il part d'observateurs finis sur une géométrie finie d'écran holographique. Sa base de travail est algébrique-quantique : algèbres de patchs, états, probabilités de type trace/Born sur les surfaces d'enregistrement déclarées et entropie généralisée font partie du point de départ formel. Le programme ne cherche pas à dériver chaque ingrédient mathématique à partir de premiers principes. Son objectif est de construire une théorie du tout cohérente et complète en utilisant cette base d'information algébrique pour reconstruire l'univers effectif observé : espace-temps, structure de jauge, particules, enregistrements et synchronisation des observateurs y apparaissent comme des conséquences de la cohérence de recouvrement, et non comme des primitives.

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

Pour le noyau technique compact, commencez par **Paper 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**. **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** porte les dérivations particules. **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** développe l'image consensus-réparation. **Paper 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** couvre l'architecture fédérée de patches, les enregistrements, les règles de preuve matérielle publique et la machinerie observateur. Ce README, le Paper 1 et le livre sont les grandes vues d'ensemble.

## Ce que l'OPH apporte

- Un cadre à cutoff fixe pour les patches d'observateurs, les collerettes, la réparation de recouvrement, la jauge supérieure, les enregistrements et le checkpoint/restauration.
- Une reconstruction de la géométrie lorentzienne, du temps modulaire, de la dynamique d'Einstein de type Jacobson et de la cosmologie de Sitter en patch statique.
- Une reconstruction de jauge compacte vers le quotient réalisé du Modèle Standard \(SU(3)\times SU(2)\times U(1)/\mathbb Z_6\), avec le réseau exact des hypercharges, le triplet de couleur réalisé \(N_c=3\) et le comptage des générations \(N_g=3\).
- Sur la branche compacte de jauge visible par support déclarée, avec échelle
  euclidienne en quatre dimensions, positivité par réflexion, complétude de
  réparation et extraction continue visible par support, une forme euclidienne
  Yang-Mills et une identification du gap de Yang-Mills avec le gap de
  réparation. Le statut face au problème de Clay dépend de l'acceptation de
  cette branche OPH comme construction axiomatique Yang-Mills en quatre
  dimensions.
- Un programme particules à niveaux de revendication explicites : porteurs structurels exactement sans masse, ligne \(W/Z\) en comparaison seule, surface Higgs/top déclarée, masses de quarks sélectionnées avec Yukawas, branche neutrino pondérée, constantes locales candidates, et frontières no-go pour les leptons chargés, la classification globale des cadres de quarks et le codomaine top direct auxiliaire.
- Une architecture microphysique d'écran concrète qui met mesure, enregistrements et observateurs à l'intérieur de la physique.

La mécanique quantique est traitée comme le langage algébrique d'information porté par l'architecture OPH. Le test de reconstruction porte sur la cohérence avec laquelle cette base retrouve l'univers effectif. Une reconstruction des espaces de Hilbert, des \(C^*\)-algèbres ou des algèbres de von Neumann, des probabilités de Born, de la trace et de l'entropie depuis de simples enregistrements opérationnels constitue une question distincte du programme.

Le mécanisme est la boucle de consensus à point fixe. Les observateurs locaux
n'accèdent pas à un état global depuis l'extérieur. Ils portent des états de
patch finis, échangent les données visibles dans les recouvrements, rejettent
les prolongements incohérents et conservent les motifs stables qui peuvent être
synchronisés. Géométrie, particules, lois et enregistrements sont les points
fixes à grande échelle de ce calcul en réseau d'observateurs.

## Géométrie, symétrie et simulateurs

Le langage de sphère dans l'OPH est un langage de géométrie. Dans les cartes
régulatrices symétriques, une coupe accessible à un observateur peut être
représentée par la deux-sphère \(S^2\). Ces cartes décrivent une géométrie
angulaire visible par support. Le simulateur n'est pas une sphère réelle. Le
substrat n'est pas un ordinateur quantique sphérique littéral. Le simulateur
est un modèle fini de patches et de recouvrements pour les contraintes
algébriques exposées par cette géométrie.

L'ancre de symétrie finie est \(A_5\), le groupe des rotations de
l'icosaèdre. Elle fournit le squelette icosaédrique derrière le langage des
porteurs de patches échosaédriques : une façon finie et très symétrique
d'organiser ports, recouvrements et données locales de comparaison sans traiter
le porteur comme une boule lisse.

L'ancre exceptionnelle est le groupe de Lie \(E_8\) et sa structure de réseau
de racines. \(E_8\) donne le langage de fermeture exceptionnelle utilisé dans
le côté symétries supérieures et représentations de la pile OPH. Le groupe
icosaédrique binaire et le \(E_8\) affine se rencontrent par la correspondance
de McKay. C'est pourquoi le langage \(A_5\)-icosaédrique et le langage de type
\(E_8\) appartiennent à une même histoire de symétrie. Ces noms désignent des
contraintes de symétrie et une structure de régulateur.

L'OPH utilise une seule entrée quantitative externe, la capacité totale de l'écran lue depuis l'horizon de Sitter, ainsi qu'un ratio local de pixel \(P\), l'aire d'une cellule d'écran en unités de Planck. Pour la constante cosmologique observée, le ratio nu d'aire d'horizon vaut environ \(1.05\times10^{122}\), tandis que la capacité entropique utilisée par l'OPH vaut environ \(3.31\times10^{122}\). La même cellule d'écran est décrite deux fois : vue de l'extérieur comme un pixel situé légèrement au-dessus de l'équilibre auto-similaire au nombre d'or, et vue de l'intérieur comme la plus petite échelle d'observation électromagnétique disponible pour les observateurs de l'univers encodé.
Le papier de synthèse écrit cette autoréférence comme un problème de point fixe. Le calcul de structure fine part de l'équilibre entropique au nombre d'or et de la normalisation gaussienne de bord, puis transporte l'ancre électrofaible jusqu'au point de Thomson. La surface publique utilise \(P\simeq1.6309682094\). La lecture donne \(\alpha^{-1}(0)=137.035999177(21)\) et \(\alpha(0)\simeq0.00729735256433\). Le même pixel local porte la gravité, la structure fine, la jauge réalisée, les lignes de masses particules, les enregistrements et la synchronisation des observateurs.

## Surface locale d'unification

L'OPH place une surface locale d'unification autour du pixel public \(P\simeq1.6309682094\). La même échelle porte la voie électrofaible en comparaison, la voie Higgs ainsi que la voie entropique gravitationnelle, tandis que la géométrie lorentzienne fournit la vitesse causale invariante et que la couche locale de lecture fournit l'affichage SI. Les formules détaillées vivent dans les papiers.
Sur la surface publique des constantes, \(\hbar\) et \(k_B\) appartiennent à cette couche aval de lecture en unités familières plutôt que d'apparaître comme des constantes OPH émises de manière autonome.

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
PDG/NIST. Les résultats structurels comme la géométrie lorentzienne en \(3+1\) dimensions, le quotient de jauge
du Modèle Standard \(SU(3)\times SU(2)\times U(1)/\mathbb Z_6\), le réseau exact des hypercharges, le triplet
de couleur réalisé \(N_c=3\) et le comptage des générations \(N_g=3\) sont énoncés dans les papiers et ne sont pas répétés ici.

| Quantité | Symbole | OPH | PDG/NIST | Δ |
| --- | --- | --- | --- | --- |
| Constante gravitationnelle | G | 6.6742999959e-11 | 6.67430(15)e-11 | 0.00003σ |
| Vitesse de la lumière | c | 299792458 | 299792458 (exact) | match |
| Structure fine (inv.) | α⁻¹(0) | 137.035999177 | 137.035999177(21) | match |
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

La surface particules porte aussi des valeurs \(W/Z\) en comparaison seule
\(80.377\,\mathrm{GeV}\) et \(91.18797809193725\,\mathrm{GeV}\), une valeur Higgs
\(m_H=125.1995304097179\,\mathrm{GeV}\) et une valeur top sélectionnée
\(m_t=172.35235532883115\,\mathrm{GeV}\) selon la convention PDG de masse top par section
efficace. La ligne top directe auxiliaire reste en comparaison seule. Les
masses des leptons chargés sont des témoins ancrés sur cible, hors du rang des
prédictions promues. La branche neutrino pondérée émet
`(0.017454720257976796, 0.019481987935919015, 0.05307522145074924) eV` sur sa
branche déclarée. Les masses hadroniques demandent un backend hadronique OPH;
elles ne sont pas de simples entrées de quarks dans un tableau de comparaison
rapide. Les lignes hadroniques gardent deux surfaces séparées : les masses
source demandent un backend hadronique OPH comme GLORB/Echosahedron avec
données spectrales et systématiques de production; la fermeture empirique
utilise une classe séparée de charge utile \(e^+e^-\to\mathrm{hadrons}\).

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse du programme OPH, des observateurs finis jusqu'à l'univers effectif reconstruit.
- **Papier 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique pour la relativité, la gravité, la structure réalisée du Modèle Standard et la forme/gap Yang-Mills sur la branche compacte de jauge visible par support sous ses hypothèses déclarées.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivations particules, lignes de masses, structure des couplages et surface quantitative.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : dynamique de réparation à point fixe, stabilité des enregistrements et consensus public.
- **Papier 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** : architecture fédérée de patches, cadrage de symétrie \(A_5\)-icosaédrique et de type \(E_8\), règles de preuve matérielle publique, enregistrements, mouvements de récupération et synchronisation observateur.
- **Note de calibration. [Digital Calibration Note for OPH Screen Microphysics](paper/screen_microphysics_digital_calibration_note.pdf)** : suite de simulation octaédrique \(\mathbb Z_2/S_3\), limitée à une calibration numérique des contraintes de patches. Ce n'est pas un modèle de substrat sphérique littéral.

## Articles supplémentaires

- **[Breaking SHA-256 with Physics](extra/breaking_sha256_with_physics.pdf)** ([TeX](extra/breaking_sha256_with_physics.tex)) : recherche physique par contraintes pour la preuve de travail SHA-256d.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)** ([TeX](extra/fine_structure_constant_derivation.tex)) : dérivation à point fixe de la ligne de structure fine.
- **[A Narrow Note on the OPH Gravity Chain](extra/jacobson_gravity_note.pdf)** ([TeX](extra/jacobson_gravity_note.tex)) : isolation de la route de Jacobson pour la branche gravitationnelle OPH.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)** ([TeX](extra/yang_mills_gap_clay_problem.tex)) : route OPH visible par support vers le problème de Yang-Mills de Clay, limitée à la branche compacte de jauge déclarée, avec identification du gap au gap de réparation.
- **[Observer-Patch Holography and the Dark Matter Phenomenon](extra/oph_dark_matter_paper.pdf)** ([TeX](extra/oph_dark_matter_paper.tex)) : phénoménologie de la matière noire et limite galactique de type MOND.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)** ([TeX](extra/thinking_as_patch_net_fixed_point_search.tex)) : cognition et qualia comme consensus récurrent de patches.

## Plus

- **Site officiel :** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Page theory of everything :** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Page simulation theory :** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Démo OPH du problème à trois corps :** [3body.floatingpragma.io](https://3body.floatingpragma.io) : un simulateur et parcours de preuve supplémentaire pour la formulation OPH en réseau fini de patches du problème à trois corps, présentée comme un exemple de recollement par holonomie de boucle, sans revendication de solution élémentaire fermée.
- **Livre :** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Application d'étude guidée :** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions et explications détaillées :** OPH Sage sur [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage) ou [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab :** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Objections courantes :** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **Note IBM Quantum :** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Table de statut

La ligne publique de structure fine utilise \(\alpha^{-1}(0)=137.035999177(21)\) et
\(P\simeq1.6309682094\). Les enregistrements d'audit source et de résidu au point
final vivent dans le papier particules, séparés de la ligne publique de point fixe.

La paire faible est une ligne de validation. Les masses absolues des leptons
chargés sont des témoins ancrés sur cible. La moyenne directe auxiliaire du top
est une ligne de validation. Les lignes hadroniques source demandent le backend
OPH de liaison forte avec données spectrales de production et systématiques; la
fermeture empirique utilise une classe séparée de charge utile
\(e^+e^-\to\mathrm{hadrons}\).

La CP forte est en cours de travail dans le théorème des quarks sur classe
sélectionnée : le corpus disponible ne dérive pas l'angle theta de QCD,
n'émet pas l'angle CP fort physique et ne prouve pas que la phase CP forte
physique s'annule. Le pont requis est la descente de phase, d'anomalie et
d'angle topologique sur la branche réalisée.

## Guide du dépôt

- **[`paper/`](paper)** : PDF, sources LaTeX et métadonnées de release.
- **[`book/`](book)** : source du livre OPH et PDF téléchargeable généré.
- **[`code/`](code)** : sorties calculatoires, surface particules et expériences.
- **[`assets/`](assets)** : schémas et figures publics.
- **[`extra/`](extra)** : notes publiques maintenues, objections, comptes rendus expérimentaux et quelques essais de support.
