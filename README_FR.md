# Holographie par Patchs d'Observateurs : une approche de la physique fondamentale centrée sur l'observateur

> OPH est un cadre de physique fondamentale centré sur l'observateur. Dans son implémentation quantitative actuelle, il utilise deux entrées continues externes (aire de pixel et capacité d'écran) ainsi que des axiomes structurels et des hypothèses explicites de régularité et de limite d'échelle. OPH développe une branche gravitationnelle conditionnelle dans la limite d'échelle, une reconstruction conditionnelle des groupes de jauge compacts avec clôture du Modèle Standard sous MAR étendu, et un programme quantitatif dont les sorties sont séparées entre contrôles de calibration, branche Higgs/top adossée à des suppléments numériques sans nouvel ajustement continu une fois la trajectoire de jauge fixée, et continuations phénoménologiques plus faibles. Le Modèle Standard et la relativité générale y émergent comme secteurs effectifs de la structure sous-jacente, tandis que le pont de type feuille d'univers en théorie des cordes reste au niveau d'une continuation dans le programme actuel.

> **Avertissement de statut :** l'OPH est un programme de recherche actif et n'est pas encore entièrement démontré. Plusieurs dérivations restent incomplètes, certaines preuves n'existent actuellement qu'à l'état d'esquisse, et certaines hypothèses auxiliaires doivent encore être éliminées. Le cadre doit donc être considéré comme étant en développement actif.

**Version anglaise :** [README.md](README.md)

> **Défi de réfutation de l'OPH :** Un défi doté de 10 000 USD pour réfuter l'OPH est actuellement en cours sur [challenge.floatingpragma.io](https://challenge.floatingpragma.io).

## Idée centrale

L'OPH adopte la position forte, centrée sur l'observateur, selon laquelle la réalité objective n'est pas fondamentale mais émerge d'un réseau de perspectives subjectives qui doivent s'accorder là où elles se recouvrent.

Les lois de la physique sont les règles de cohérence qui rendent cet accord intersubjectif possible.

À partir de ce point de départ (avec les contraintes d'entropie et de Markov), l'OPH fait émerger espace-temps, symétries de jauge et physique des particules comme conséquences de cohérence.

## Programme de recherche

L'OPH doit actuellement être comprise comme un programme de recherche plutôt que comme une théorie achevée. En physique, on ne démontre pas au sens absolu qu'une théorie est "la bonne" ; la cible réaliste est un cadre mathématiquement explicite, adossé à des régulateurs concrets, empiriquement discriminant, et robuste face aux tentatives de falsification.

Le programme actuel consiste donc à :

- transformer les arguments qui ne sont encore qu'à l'état d'esquisse en preuves complètes
- supprimer les hypothèses auxiliaires lorsque c'est possible, ou les isoler proprement lorsqu'elles restent nécessaires
- réaliser les axiomes dans des modèles microscopiques explicites d'écran ou dans des régulateurs concrets
- dériver à partir de ces réalisations les secteurs gravitationnel, de jauge, et quantitatif avec des approximations contrôlées
- remplacer les étapes encore dépendantes de calibrations par des dérivations de premiers principes lorsque c'est bien cela qui est revendiqué
- produire des critères de contradiction empiriques nouveaux et des tests sensibles aux différentes branches, qui puissent échouer proprement

L'état final plausible le plus fort n'est donc pas une "preuve de l'OPH" au sens métaphysique, mais une théorie dérivée rigoureusement et tenue responsable expérimentalement, dont les secteurs effectifs, l'économie en paramètres et les prédictions falsifiables résistent à un examen indépendant.

## Articles

**Observers are all you need** est l'article principal de l'OPH.

- **PDF (article principal) :** [Observers are all you need](paper/observers_are_all_you_need.pdf)
- **Source LaTeX :** [observers_are_all_you_need.tex](paper/observers_are_all_you_need.tex)

**Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency** est
l'article compact de soumission. Il concentre le coeur falsifiable de l'OPH :
cinématique de Lorentz dans la phase modulaire géométrique, branche d'Einstein conditionnelle
dans la limite d'échelle, clôture de jauge du Modèle Standard sous le paquet d'admissibilité MAR
étendu, réseau d'hypercharges sur le paquet de matière réalisé à une génération sous la
normalisation standard, chaîne de comptage réalisée $N_g=3$ puis $N_c=3$ sur la branche MAR
réalisée, implémentation quantitative à deux entrées, et branches quantitatives les plus fortes
adossées à des suppléments numériques.

Dans l'ensemble synchronisé actuel des articles, le durcissement du pont modulaire nul a supprimé toute hypothèse séparée de propagation quasi-locale, de contrôle des extrémités ou d'inclusion modulaire à demi-côté dans la route gravitationnelle à cutoff fixe : la branche MaxEnt porte désormais les estimations de contrôle en interne, et l'argument de blow-up géométrique ferme le pont jusqu'au stade de la translation nulle positive dérivée. Ce qui reste explicite sur cette route, ce sont les prémisses en aval de montée en densité et d'identification relativiste du stress nul.

- **PDF (article compact de soumission) :** [Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)
- **Source LaTeX :** [recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex)

**Reality as a Consensus Protocol** est un article complémentaire orienté informatique qui présente l'ossature computationnelle de l'OPH sous la forme d'un paquet de théorèmes autonome. Il formule la loi physique objective comme le point fixe d'un protocole de réconciliation distribué entre patchs d'observateurs, montre comment la topologie peut faire obstacle à la cohérence globale avec des particules comme défauts stables, interprète la symétrie de jauge comme une forme de masquage d'implémentation, et modélise les enregistrements classiques comme une couche CRDT à cohérence éventuelle.

- **PDF :** [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)
- **Source LaTeX :** [reality_as_consensus_protocol.tex](paper/reality_as_consensus_protocol.tex)

Chaque PDF comporte désormais une ligne de version visible. La source partagée de version est
[`paper/release_info.tex`](paper/release_info.tex). Pour toute mise à jour substantielle des
articles, incrémentez d'abord la version partagée avant de reconstruire les PDF :

```bash
python3 tools/bump_paper_release.py
```

Après reconstruction des PDF, écrivez les empreintes courantes dans
[`paper/paper_release_manifest.json`](paper/paper_release_manifest.json) en lançant :

```bash
python3 tools/generate_paper_release_manifest.py
```

Le générateur de manifeste échoue désormais si les PDF changent sans nouvelle version, ou si les
PDF locaux n'exposent pas encore la ligne de version visible attendue.

La version est globale à l'ensemble courant des articles. Même si un seul article change,
incrémentez une seule fois, reconstruisez les trois PDF, puis publiez les trois PDF du challenge
avec cette même version.
L'envoi vers le challenge est géré par un outillage opérationnel local à l'espace de travail, pas par ce dépôt public.

Les sources canoniques des articles se trouvent désormais dans [`paper/`](paper). L'ancien arbre de brouillons a été supprimé du dépôt.

Lorsqu'une mise à niveau ou un audit touche plusieurs articles, commencez toujours par l'article compact de soumission. Ce n'est qu'une fois sa ligne de version et son langage de statut cohérents qu'il faut propager les mêmes corrections vers l'article principal, le livre, les README, les sites publics et la pipeline d'ingestion.

Dans l'ensemble actuel des articles, l'OPH sépare explicitement les affirmations par statut :
sorties structurelles sans paramètre (chaîne Lorentz/jauge/quotient, réseau d'hypercharges sur le
paquet de matière réalisé, $N_g=3$ puis $N_c=3$ sur la branche MAR réalisée, zéros protégés par
symétrie, stabilité du proton), contrôles de cohérence du secteur calibré ($\alpha_i(m_{Z,\rm run})$, $\sin^2\theta_W(m_{Z,\rm run})$, $v$, $m_W$), branche quantitative principale adossée à des suppléments
($m_H$, $m_t$), et continuations phénoménologiques ou numériques en aval plus faibles (Koide /
leptons chargés, textures de quarks, neutrinos, hadrons, matière noire,
baryogenèse).

## Ressources officielles OPH

Commencez par le site officiel de l'OPH, le labo interactif OPH et la version web du livre OPH :

- **Défi de réfutation de l'OPH (10 000 USD) :** [challenge.floatingpragma.io](https://challenge.floatingpragma.io)
- **Site officiel de l'OPH :** [Site officiel de l'OPH sur floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Labo interactif OPH :** [OPH Lab sur oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Livre OPH (edition web) :** [Lire le livre OPH en ligne](https://oph-book.floatingpragma.io)
- **NotebookLM :** [Regardez la presentation d'introduction a l'OPH et commencez a poser vos questions](https://notebooklm.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a?artifactId=fb7c0ebd-4375-4997-9cae-6558ff8977b4)

## Expériences IBM Quantum

Un premier bundle public d'expériences IBM Quantum Cloud est maintenant inclus dans ce dépôt. Il résume les premiers tests de récupérabilité et de rapports exacts, les sorties matérielles mesurées, ainsi que le bundle public de code et de données utilisé pour ces exécutions.

- **Note expérimentale :** [IBM Quantum Cloud Evidence for OPH](extra/IBM_QUANTUM_CLOUD.md)
- **Code et données publics :** [code/ibm_quantum_cloud/](code/ibm_quantum_cloud/)

## Posez vos questions à OPH Sage

Le moyen le plus rapide d'explorer l'OPH est de poser vos questions directement à OPH Sage. OPH Sage peut répondre à des questions techniques et philosophiques.

- [OPH Sage sur Telegram](https://t.me/HoloObserverBot)
- [OPH Sage sur X](https://x.com/OphSage)
- [OPH Sage sur Bluesky](https://bsky.app/profile/ophsage.bsky.social)

## Objections courantes

Cette section regroupe des réponses aux objections courantes adressées à l'OPH.

- [Dériver `P` à partir des données de jauge puis réutiliser `P` en aval est complètement circulaire](extra/COMMON_OBJECTIONS.md#objection-1-circularity)
- [Une taille de cellule fixe brise l'invariance de Lorentz, donc l'OPH ne peut retrouver qu'une limite newtonienne](extra/COMMON_OBJECTIONS.md#objection-2-lorentz)
- [L'OPH aurait une discontinuité Type I / Type III, donc son histoire du temps modulaire serait incohérente en interne](extra/COMMON_OBJECTIONS.md#objection-3-type-i-type-iii)

**Lire :** [Comment l'OPH résout la physique fondamentale](https://medium.com/@muellerberndt/answering-10-of-the-hardest-questions-in-physics-and-some-bonus-questions-51222bf2419f)

**Articles de blog :**

- [How OPH Resolves Fundamental Physics](https://medium.com/@muellerberndt/answering-10-of-the-hardest-questions-in-physics-and-some-bonus-questions-51222bf2419f)
- [How Observer Patch Holography Improves on the Standard Model and General Relativity](https://muellerberndt.medium.com/how-observer-path-holography-improves-on-the-standard-model-and-general-relativity-c971c376027e)

La physique est de la rétro-ingénierie à l'échelle cosmique. Nous observons le comportement de la réalité et remontons vers les règles sous-jacentes. Ce projet présente une théorie du tout au sens littéral : pourquoi la réalité existe et pourquoi elle est telle qu'elle est, à partir de premiers principes.

## L'idée centrale (détaillée)

La section suivante transforme ce principe en cadre concret, construit à partir de patchs d'observateurs, de cohérence de recouvrement et de reconstruction holographique.

## Holographie par Patchs d'Observateurs

Nous modélisons la réalité comme un réseau de perspectives subjectives qui doivent s'accorder là où elles se recouvrent. Concrètement, nous partons de patchs d'observateurs sur un écran holographique 2D. Chaque patch représente une perspective avec ses propres données locales. Là où les patchs se recouvrent, leurs descriptions doivent coïncider. Dans l'interprétation OPH, la « réalité objective » est l'ossature cohérente au recouvrement partagée entre ces perspectives, plutôt qu'un primitif supposé.

L'invariance de Lorentz existe parce que différents observateurs doivent fournir des descriptions cohérentes. La symétrie de jauge existe parce que les patchs qui se recouvrent doivent identifier les mêmes observables partagées. Les lois de conservation existent parce que les mêmes quantités doivent être conservées d'une perspective à l'autre. Les lois ne sont pas imposées de l'extérieur ; elles sont les conditions de possibilité de l'accord.

Le modèle repose sur quatre concepts centraux :

- **Écran** : une 2-sphère de type horizon (comme un horizon cosmique entourant chaque observateur) qui porte l'information quantique. C'est là que vivent les données fondamentales.
- **Patch** : une région connexe de l'écran accessible à un observateur donné. Chaque patch possède sa propre algèbre d'observables, c'est-à-dire les questions que cet observateur peut poser.
- **Cohérence de recouvrement** : là où deux patchs partagent une région, leurs descriptions doivent s'accorder. C'est l'axiome central. Il fait de l'objectivité une structure intersubjective reconstruite, plutôt qu'un point de départ.
- **Observateur** : un motif stable dans les données de l'écran qui maintient des enregistrements et participe aux relations de cohérence.

### La réalité comme calcul

On peut voir l'écran comme un système quantique invariant de jauge sur une 2-sphère. Il ressemble à un automate cellulaire quantique, avec une structure supplémentaire importante. À chaque point de la sphère triangulée, on place des systèmes quantiques de dimension finie (qudits), couplés par des contraintes de jauge à chaque sommet. Toutes les configurations ne sont pas physiques : seules celles qui satisfont la loi de Gauss subsistent.

**Les patchs d'observateurs** sont des sous-systèmes définis par des algèbres invariantes par la jauge de bord. Chaque patch est un fil de calcul : une région connexe où un observateur peut poser des questions et obtenir des réponses. L'algèbre $\mathcal{A}(R)$ définit ce qui est mesurable, c'est-à-dire les opérateurs invariants sous les transformations de jauge au bord.

**La cohérence de recouvrement** est automatique. Quand deux patchs se recouvrent, ils accèdent aux mêmes observables invariantes de jauge. Les deux observateurs lisent les mêmes données sous-jacentes, depuis des angles différents. La redondance de jauge au bord rend le collage non trivial et engendre les « modes de bord » qui portent l'information géométrique.

**Les observateurs ne sont pas externes.** Ce sont des structures calculatoires émergentes *dans* les données de l'écran. Ce sont des motifs stables qui traitent de l'information, maintiennent des enregistrements et créent des corrélations.

**Le bulk 4D n'est pas sur la sphère.** Il émerge de la structure d'intrication entre patchs. La sphère est la frontière ; l'intérieur est reconstruit holographiquement. Quand vous voyez un espace tridimensionnel, vous faites l'expérience d'un encodage compressé de la façon dont votre patch est intriqué avec les autres.

*L'écran est le calcul. Les patchs sont les fils d'exécution. La réalité est ce sur quoi ils s'accordent.*

<a href="assets/french/screen_fr.svg"><img src="assets/french/screen_fr.svg" alt="L'écran holographique" width="800"></a>

### Qu'est-ce qui pilote le calcul ?

Dans les modèles de liens quantiques, la dynamique combine des termes de plaquette (boucles de Wilson autour des faces) et des termes de champ électrique (conjugués aux variables de lien). Leur compétition détermine l'état fondamental. Le « calcul » est l'évolution de ces degrés de liberté quantiques, qui créent et détruisent des corrélations, tandis que les contraintes de jauge assurent la cohérence.

Le « temps » vécu par les observateurs n'est pas forcément l'évolution hamiltonienne microscopique. Chaque patch possède son propre hamiltonien modulaire (construit à partir de la matrice densité réduite du patch), et c'est lui qui génère ce qui ressemble au temps depuis l'intérieur. L'évolution microscopique et le temps modulaire émergent sont liés, mais distincts.

Le système existe dans un état intemporel (comme le suggère l'équation de Wheeler-DeWitt en gravité quantique), et ce que nous appelons « temps » est relationnel. Le flot modulaire donne à chaque sous-système sa propre horloge interne, corrélée aux autres par les contraintes de cohérence.

### Pourquoi cette approche fonctionne

Les approches unifiées cherchant à combiner TQC, gravité et structure du Modèle Standard rencontrent souvent les mêmes difficultés : échec de factorisation des sous-systèmes en jauge/gravité, non-localité des hamiltoniens modulaires, invariance de Lorentz supposée plutôt que dérivée, difficulté à obtenir la dynamique (pas seulement la cinématique), origine de la symétrie de jauge peu claire, absence de masse imposée à la main, anomalies traitées comme pathologies, quantification de charge nécessitant des GUT, unification des couplages impliquant une désintégration du proton, surdétermination locale de la constante cosmologique, prolifération des infinis UV, et explosion du nombre de paramètres.

Le cadre des patchs d'observateurs renverse la logique : ce sont les conditions de cohérence qui font le travail. Localité, Lorentz, jauge et gravité sont traitées comme des contraintes de cohérence entre descriptions recouvrantes, combinées à des propriétés informationnelles des états (Markov/récupérabilité + MaxEnt), puis rigidifiées par la théorie modulaire pour forcer les symétries et dynamiques familières.

## Les axiomes

Le cadre repose sur quatre axiomes centraux plus un axiome de sélection :

| Label | Nom | Contenu |
|-------|-----|---------|
| **A1** | Réseau d'écran | Les algèbres d'observables vivent sur une surface 2D fermée $S^2$ |
| **A2** | Cohérence de recouvrement | Quand des patchs partagent une région, leurs descriptions coïncident |
| **A3** | Entropie généralisée | $S_{\rm gen} = S_{\rm bulk} + \langle L_C \rangle$ satisfait la sous-additivité forte |
| **A4** | Markov local | L'information mutuelle conditionnelle décroît à travers les colliers |
| **MAR** | Réalisation admissible minimale | Parmi les secteurs admissibles, la Nature réalise le minimum lexicographique |

Les hypothèses structurelles supplémentaires (MaxEnt, prémisses de limite d'échelle contrôlée et de branche géométrique OPH, mélange exponentiel) sont détaillées dans la [source principale](paper/tex_fragments/PAPER.tex).

## La chaîne de prédictions

L'infographie suivante résume le programme actuel de reconstruction OPH, depuis deux paramètres et quatre axiomes vers une physique effective :

![Chaîne de prédictions OPH](assets/french/prediction-chain_fr.svg)

*Des axiomes à une physique effective : le programme actuel de reconstruction OPH.*

> **Prédictions du spectre de particules :** la dérivation complète de l'aire de pixel aux masses des particules, avec comparaison aux données PDG et contrôles d'audit adossés au dépôt, est présentée dans **[la source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex)**. Les sorties y sont classées par statut épistémique : les couplages de jauge et les masses des bosons électrofaibles sont des contrôles de cohérence du secteur de calibration ; les masses du Higgs et du top proviennent d'une reconstruction D11 adossée au supplément, sans nouvel ajustement continu une fois la trajectoire de jauge et la branche de fixation d'échelle fixées, mais utilisant une couche de transport synchronisée en UV (SM sous $\mu_{\rm sync}\approx 6.8\times 10^{11}$ GeV, de type MSSM au-dessus) plus de petits décalages d'appariement à basse énergie ; les masses des leptons chargés (< 0,04 % via la structure de Koide) restent une continuation phénoménologique plus fine mais plus faible parce que leur assignation finale d'exposants utilise une étape supplémentaire de sélection discrète.

## Les paramètres fondamentaux

Notre univers est actuellement décrit, dans l'implémentation quantitative, par exactement **deux entrées fondamentales fixées de l'extérieur** :

### 1. Aire de pixel : $a_{\text{cell}} \approx 1.63 \, \ell_P^2$

L'aire géométrique d'un élément computationnel de l'écran holographique. Elle fixe la *résolution* de la réalité :

| Grandeur | Valeur | Signification |
|----------|--------|---------------|
| En unités de Planck | $a_{\text{cell}} / \ell_P^2 \approx 1.63$ | Rapport sans dimension |
| En unités SI | $a_{\text{cell}} \approx 4.26 \times 10^{-70}$ m^2 | Aire physique par pixel |
| « Côté » du pixel | $\sqrt{a_{\text{cell}}} \approx 2.06 \times 10^{-35}$ m | Échelle de résolution |

**Ce qu'elle détermine :** constante de Newton (via $G_{\text{nat}} = a_{\text{cell}}/4\bar{\ell}$ en unités naturelles où $G_{\text{nat}} = \ell_P^2$), échelle de Planck, secteur calibré des couplages de jauge, et échelle utilisée par la branche Higgs/top indépendante ainsi que par les continuations en aval.

### 2. Capacité de l'écran : $\log(\dim \mathcal{H}) \sim 10^{122}$

L'entropie totale de l'écran holographique (en nats). Elle fixe la *taille* de la réalité.

**Important :** la capacité de l'écran est **inférée** à partir de la constante cosmologique observée, elle n'est pas prédite. La relation

$$\Lambda = \frac{3\pi}{G \cdot \log(\dim \mathcal{H})}$$

permet d'en déduire la capacité à partir de $\Lambda \sim 10^{-52}$ $\text{m}^{-2}$, donnant $\log(\dim \mathcal{H}) \sim 10^{122}$.

| Grandeur | Valeur | Statut |
|----------|--------|--------|
| Entropie de de Sitter | $S_{dS} \sim 10^{122}$ nats | Inférée depuis $\Lambda$ observé |
| Horizon de de Sitter | $r_{dS} \approx 10^{26}$ m | Observé |

**Relation :**
- **Aire de pixel** = résolution (extraite des couplages de jauge via l'entropie de bord)
- **Capacité de l'écran** = taille totale (extraite de la constante cosmologique observée)

La structure axiomatique ne contient aucune autre constante dimensionnée. Les axiomes plus la reconstruction donnent *un* groupe de jauge compact ; le paquet d'admissibilité MAR étendu sélectionne ensuite $SU(3) \times SU(2) \times U(1)/\mathbb{Z}_6$ comme groupe réalisé (voir [la source de dérivation du groupe de jauge](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex)). La quantification de charge et la dynamique d'Einstein semi-classique dans la limite d'échelle suivent alors.

### Ce que cela signifie

L'aire de pixel et la capacité de l'écran sont des **paramètres de configuration**, les « réglages » du calcul qu'est notre univers. Ils ne sont pas dérivables depuis l'intérieur de la simulation ; ce sont des conditions aux limites.

Depuis l'intérieur :
- **Aire de pixel** -> constante de Newton, échelle de Planck, secteur calibré des couplages de jauge, et échelle utilisée par la branche Higgs/top indépendante ainsi que par les continuations en aval
- **Capacité de l'écran** -> taille de l'univers observable (et elle-même inférée de $\Lambda$ observé)

Les mêmes axiomes avec d'autres réglages produiraient un univers avec d'autres constantes, mais une physique de structure similaire (Einstein, jauge, etc.).

### Calibration vs prédiction

Dans l'implémentation actuelle, la constante de pixel $P = a_{\text{cell}}/\ell_P^2$ est *inférée* depuis les couplages de jauge mesurés, car les axiomes fixent la relation fonctionnelle $P = 4\bar{\ell}_{\text{tot}}(t_2,t_3)$ sans encore fixer les multiplicateurs de Lagrange MaxEnt $t_i$ (donc les couplages) par premiers principes.

Le contenu non trivial est que $P$ fournit une contrainte supplémentaire reliant le couplage gravitationnel à l'entropie de bord du secteur de jauge, et qu'en mode à deux entrées (traitant $P$ comme paramètre fondamental + une seule donnée électrofaible $\alpha(M_Z)$), le cadre prédit simultanément $\alpha_s(M_Z)$ et $\sin^2\theta_W(M_Z)$.

Une fermeture non circulaire complète nécessite un principe UV qui fixe $t$ sans utiliser les couplages mesurés.

### Statut actuel des articles

Les articles actuels ne traitent pas le secteur quantitatif comme un bloc uniforme. Ils séparent
théorèmes structurels, contrôles de cohérence du secteur calibré, branche Higgs/top sur surface
critique, et continuations phénoménologiques ou numériques en aval plus faibles.

| Palier | Sorties représentatives | Statut actuel | Source principale |
|--------|-------------------------|---------------|-------------------|
| Théorèmes structurels | Cinématique de Lorentz conditionnelle, branche d'Einstein conditionnelle dans la limite d'échelle, reconstruction de jauge compacte, chaîne de quotient du MS, hypercharges exactes, $N_g=3$ puis $N_c=3$ | Paquet central de théorèmes sous les prémisses de limite d'échelle et de catégorie explicitement énoncées | [Source principale](paper/tex_fragments/PAPER.tex), [Source de dérivation du groupe de jauge](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex) |
| Conséquences structurelles exactes | Stabilité du proton ; $m_\gamma = 0$, $m_g = 0$, $m_{\text{graviton}} = 0$ | Conséquences structurelles sans paramètre une fois réalisée la structure de jauge / difféomorphisme / groupe produit correspondante | [Source principale](paper/tex_fragments/PAPER.tex), [Source de dérivation du groupe de jauge](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex), [Source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Secteur calibré | $\alpha_s=0.1183$, $\sin^2\theta_W=0.2307$, $\alpha_{\rm em}^{-1}=128.31$, $v$, $W$, $Z$ | Contrôles de cohérence après calibration pixel/jauge ; pas présentés comme confirmation quantitative indépendante | [Source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Branche quantitative indépendante | Higgs $=126.48$ GeV, top(pôle) $=171.1$ GeV | Principales sorties quantitatives indépendantes du module D11 synchronisé : pas de nouvel ajustement continu une fois la trajectoire de jauge fixée, mais les valeurs affichées utilisent le transport UV synchronisé et les conventions de raccordement basse énergie documentées | [Source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Continuations phénoménologiques | Branche leptons chargés / Koide ; textures de quarks | Palier plus faible : les leptons chargés sont numériquement très précis mais utilisent une étape discrète supplémentaire ; les masses de quarks individuelles gardent une grande incertitude de schéma/seuils | [Source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Branches de capacité / numériques en aval | $m_{\nu_3}\approx 3.0$ meV, $m_{\nu_2}\approx 0.50$ meV, $m_{\nu_1}\approx 0.084$ meV ; $P \rightarrow \alpha_s \rightarrow \Lambda_{\overline{\rm MS}}^{(3)} \rightarrow m_{\rm hadrons}$ | Les neutrinos ne sont pour l'instant qu'une estimation d'ordre de grandeur issue de la branche capacité ; la chaîne hadronique est complète mais la précision reste limitée par les systématiques lattice | [Source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Branches d'extension / programmatiques | Poids de bord OPH = noyaux de chaleur YM 2D ; expansion de feuille d'univers à grand $N$ (Gross-Taylor) ; matière noire, baryogenèse, trous noirs | L'identification des poids de bord avec le YM 2D est établie dans sa branche ; la réorganisation de feuille d'univers à grand $N$ et les autres sujets restent au niveau continuation/programme | [Source de dérivation de la théorie des cordes](paper/tex_fragments/STRING_THEORY.tex), [Supplément technique source](paper/tex_fragments/TECHNICAL_SUPPLEMENT.tex) |

**Cibles de fermeture ouvertes (suivies explicitement) :**
- Dériver les multiplicateurs MaxEnt $t_i$ (donc couplages/calibration pixel) depuis des premiers principes.
- Vérifier le pont EFT (bandes nulles comme séparateurs A4 + variation finie locale) dans des régulateurs UV explicites.
- Améliorer la précision numérique (raccords 2 boucles/seuils, réseau non trempé plus grand volume pour hadrons).
- Dériver $\theta_{\rm QCD}$ non perturbativement depuis la structure de collage/cohomologie.

## Ce que le modèle dérive

### Conséquences directes de la cohérence des patchs

Ces résultats suivent de l'absence de référentiel « tiers » privilégié. Il n'y a que des perspectives d'observateurs qui doivent s'accorder dans leurs recouvrements.

**Pas de référentiel privilégié :** il n'existe pas de point de vue absolu pour définir un repos absolu. Chaque perspective est aussi valide que les autres. Point crucial : il n'y a pas d'observateur *extérieur* à la sphère ; les observateurs sont des motifs *dans* les données de l'écran, avec accès seulement à leur patch.

**Invariance de Lorentz :** les transformations entre descriptions cohérentes forment un groupe. Sur $S^2$, ce groupe est $\mathrm{Conf}(S^2) \sim SO(3,1)$ : la relativité restreinte émerge comme condition de cohérence.

**Symétrie de jauge comme redondance de recouvrement :** dans les recouvrements, il existe une redondance d'identification des observables partagées. Cette redondance a la structure d'un groupe de jauge.

**4D depuis 2D :** le groupe conforme de $S^2$ est isomorphe au groupe de Lorentz en 3+1D. Les symétries de l'écran deviennent les symétries de l'espace-temps.

### Équations d'Einstein

**Les équations d'Einstein** sont obtenues par la chaîne suivante :

1. *Sélection MaxEnt* (Hypothèse B) : l'état global maximise l'entropie sous contraintes locales.
2. *Forme de Gibbs locale* (Lemme 2.6) : conséquence de MaxEnt + contraintes locales.
3. *Pont EFT* (N1-N3) : dérivé des axiomes A1-A4 sous deux conditions testables (bandes nulles A4, variation finie locale).
4. *QNEC/Focalisation* : dérivé en interne via monotonie de l'entropie relative.
5. *Équilibre d'intrication* : mène aux équations d'Einstein via le mécanisme de Jacobson (1995, 2016).

Les hypothèses explicites restantes sont MaxEnt, les prémisses OPH de branche géométrique Bisognano--Wichmann et de limite d'échelle derrière l'étape BW, ainsi que les conditions du pont EFT ci-dessus. Dans cette formulation, la prémisse de branche BW fournit la dilation conforme standard préservant le cap, tandis que le théorème BW fixe sa normalisation modulaire et le reste de collerette contrôlé sur cette branche.

### Particules et forces

Le principe jauge-comme-collage reconstruit des groupes de jauge compacts à partir des règles de fusion des secteurs de bord.

**Photon sans masse :** la redondance U(1) électromagnétique émergente interdit un terme de masse. Prédiction : masse exactement nulle. Les bornes expérimentales actuelles sont compatibles avec cette prédiction, en contraignant toute déviation à zéro à un niveau d'environ 10⁻¹⁸ eV (souvent résumé comme ~27 ordres de grandeur).

**Graviton sans masse :** l'invariance par difféomorphisme émergente interdit un terme de masse. Prédiction : masse exactement nulle, avec des bornes observationnelles actuelles compatibles avec cette prédiction (par ex. m_g ≲ 1.8×10⁻²³ eV via la dispersion des ondes gravitationnelles ; souvent résumé comme ~22 ordres de grandeur).

**Trois couleurs :** contrainte d'anomalie de Witten + admissibilité + MAR conduisent à $N_c = 3$.

**Trois générations :** CP impose $N_g \ge 3$, liberté asymptotique de $SU(2)_L$ impose $N_g \le 5$, puis MAR sélectionne $N_g = 3$.

**Stabilité du proton :** sous le paquet MAR étendu, le porteur fidèle minimal $\mathbb{C}^3 \otimes \mathbb{C}^2$ impose une structure produit. Pas d'inclusion dans un groupe simple plus grand, donc pas de générateurs leptoquarks. Voir [la source de dérivation du groupe de jauge](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex).

## Comment le modèle explique la physique connue

La plupart des points ci-dessous sont des *postdictions* : masses nulles des bosons de jauge, nombre de couleurs, structure de l'espace-temps, etc. étaient connus avant. L'enjeu est de savoir si le cadre les explique à partir de principes plus profonds.

### Dérivés des hypothèses du cadre

Ces résultats suivent des axiomes centraux (A1-A4) + hypothèses supplémentaires (MaxEnt, jauge-comme-collage, et prémisses OPH de limite d'échelle / branche géométrique). Le jeu complet est détaillé dans l'article technique.

| Résultat | Chaîne de dérivation |
|--------|------------------|
| Masse du photon = 0 | Redondance de jauge émergente + U(1)_em non brisée -> terme de masse interdit |
| Masse du graviton = 0 | Équilibre d'intrication + invariance difféomorphe -> terme de masse interdit |
| Masse du gluon = 0 | Jauge-comme-collage pour SU(3) -> terme de masse interdit |
| Groupe de Lorentz | A1-A4 + théorème BW sur la branche géométrique Bisognano--Wichmann d'OPH -> $\mathrm{Conf}(S^2) \sim SO(3,1)$ |
| Invariance CPT | Cinématique de Lorentz + localité via théorème CPT |
| Conservation de charge | Symétrie de jauge U(1)_em non brisée |
| Formule de la constante de Newton | $G = a_{\rm cell}/4\bar{\ell}(t)$ depuis la densité d'entropie de bord |

### Dérivés avec contenu en matière supposé

| Résultat | Statut |
|--------|--------|
| Hypercharges (rationnels exacts) | Contenu en matière MS supposé |
| Quantification de charge | Quotient Z6 du spectre réalisé |
| Règle de congruence Z6 | Structure globale du groupe du MS |
| Déficit d'entropie de bord = log2 6 bits | Loi de noyau de chaleur + quotient Z6 |
| Hiérarchie de Yukawa $y_f \propto 6^{-n_f}$ | Suppression de défaut Z6 + charges entières |

### Sorties quantitatives face aux données existantes

La chaîne complète depuis l'aire de pixel $P = 1.63094$ jusqu'aux masses est documentée dans **[la source de dérivation du spectre](paper/tex_fragments/SPECTRUM_DERIVATION.tex)**. Les résultats additionnels (matière noire, Koide, baryogenèse, spin du proton) sont dans le **[supplément technique source](paper/tex_fragments/TECHNICAL_SUPPLEMENT.tex)**.

**Contrôles de cohérence du secteur calibré** (les couplages servent à fixer $P$ ; l'accord est attendu par construction) :

| Grandeur | OPH | PDG | Erreur rel. | Statut |
|----------|----:|----:|-----------:|--------|
| $\alpha_s(m_{Z,\rm run})$ | 0.1183 | 0.1179 ± 0.0009 | +0.37% | Contrôle de cohérence |
| $\sin^2\theta_W(m_{Z,\rm run})$ | 0.2307 | 0.23122 ± 0.00004 | −0.21% | Contrôle de cohérence |
| $\alpha_{\text{em}}^{-1}(m_{Z,\rm run})$ | 128.31 | 127.952 ± 0.009 | +0.28% | Contrôle de cohérence |
| Boson W | 80.386 GeV | 80.377 ± 0.012 GeV | +0.012% | Contrôle de cohérence |
| Boson Z | 91.220 GeV | 91.188 ± 0.002 GeV | +0.035% | Contrôle de cohérence |
| VEV du Higgs | 246.77 GeV | 246.22 GeV | +0.22% | Contrôle de cohérence |

**Branche quantitative indépendante** (principales sorties de masse hors calibration) :

Les valeurs affichées pour le Higgs et le top sont les sorties synchronisées de D11. Le flot littéral tout-SM de l'appendice basse énergie tombe beaucoup plus bas, autour de $m_H \approx 115.18$ GeV et $m_t^{\rm pole} \approx 164.25$ GeV. La branche publiée $126.48 / 171.1$ utilise un transport UV synchronisé avec évolution SM sous $\mu_{\rm sync}\approx 6.8\times 10^{11}$ GeV, coefficients de type MSSM au-dessus, puis seulement des corrections résiduelles de raccordement de l'ordre du GeV.

| Grandeur | OPH | PDG | Erreur rel. | Origine |
|----------|----:|----:|-----------:|---------|
| Masse du Higgs | 126.48 GeV | 125.20 ± 0.11 GeV | +1.02% | Reconstruction D11 du supplément après $\lambda = \beta_\lambda = 0$ |
| Quark top (surface critique) | 171.1 GeV | 172.57 ± 0.29 GeV | −0.87% | Reconstruction D11 du supplément après $\lambda = \beta_\lambda = 0$ |

**Continuations phénoménologiques** (conservées pour comparaison, mais plus faibles que le paquet de théorèmes et la branche Higgs/top) :

| Grandeur | OPH | PDG | Erreur rel. | Origine |
|----------|----:|----:|-----------:|---------|
| Électron | 5.109 × 10⁻⁴ GeV | 5.110 × 10⁻⁴ GeV | −0.023% | Continuation Koide ($Q=2/3$, $\delta=2/9$) |
| Muon | 0.10564 GeV | 0.10566 GeV | −0.022% | Continuation Koide |
| Tau | 1.7766 GeV | 1.7769 GeV | −0.020% | Continuation Koide |
| Quark top (texture Z₆) | 174.5 GeV | 172.57 ± 0.29 GeV | +1.1% | Continuation texture $\mathbb{Z}_6$ ($n_t = 0$) |

**Continuations de programme / supplémentaires** (voir le [supplément technique source](paper/tex_fragments/TECHNICAL_SUPPLEMENT.tex) ; elles n'ont pas toutes le même statut que les résultats centraux) :

| Grandeur | OPH | Observé | Accord |
|----------|-----|----------|-----------|
| Échelle d'accélération MOND $a_0$ | $1.03 \times 10^{-10}$ m/s² | $1.2 \times 10^{-10}$ m/s² | 15% |
| Ratio de Koide $Q$ | 2/3 (exact) | 0.666664 | 10⁻⁵ |
| Phase de Koide $\delta$ | 2/9 (dérivée) | 0.222225 | 10⁻⁵ |
| Asymétrie baryonique $\eta_B$ | $4.6 \times 10^{-10}$ | $6.1 \times 10^{-10}$ | Facteur 1.3 |
| Fraction de spin du proton $\Delta\Sigma$ | 0.308 | 0.29 ± 0.03 | Dans 1σ |
| Ratio de décalage MSSM $\Delta b_3/\Delta b_2$ | 0.91 | 0.96 (MSSM) | 5% |

**Prédictions structurelles :**

| Prédiction | Statut expérimental | Notes |
|------------|---------------------|-------|
| Quantification de charge Z6 | PDG : $\|q_p + q_e\|/e < 10^{-21}$ | Prédiction structurelle |
| Ratios log-gap de Casimir | Lattice SU(3) confirme 9/4, 5/2, 4, 9/2, 6 | Ratios sans paramètres |
| Masse du photon = 0 | PDG : $m_\gamma < 10^{-18}$ eV | Exact (symétrie de jauge) |
| Masse du graviton = 0 | PDG : $m_g < 1.76 \times 10^{-23}$ eV | Exact (difféomorphisme) |

### Dérivés sous la théorie étendue ($T_{\text{ext}}$ : A1-A4 + R0 + R1 + [z]=0 + MAR)

Ces résultats sont dérivés via l'axiome de sélection MAR :

| Résultat | Comment c'est dérivé |
|--------|-------------|
| Groupe de jauge produit | Porteur fidèle minimal $\mathbb{C}^3 \otimes \mathbb{C}^2$ (MAR) |
| Groupe global MS SU(3) × SU(2) × U(1)/Z₆ | MAR + conditions d'admissibilité |
| $N_c = 3$ (trois couleurs) | Anomalie de Witten + minimalité MAR |
| $N_g = 3$ (trois générations) | CP + liberté asymptotique + minimalité MAR |
| Stabilité du proton | Structure produit via MAR (pas de générateurs leptoquarks) |
| Pas de monopôles magnétiques | Structure produit (pas de brisure GUT à haute échelle) |

### Vérifications de cohérence (pas des prédictions nouvelles)

| Résultat | Source originale |
|--------|-----------------|
| $(b_2, b_3) = (1, -3)$ depuis indices de Dynkin à $t^* \approx 4/3$ | Résultat nouveau (réduit le bêta MSSM à la géométrie de collier) |
| $\alpha_s(M_Z) \approx 0.118$ avec spectre de type MSSM | Analyses GUT MSSM (années 1990) |
| $\sin^2\theta_W(M_U) = 3/8$ | Georgi et Glashow (1974) |
| Contrainte d'anomalie de Witten | Witten (1982) |
| Mécanisme GIM (pas de FCNC à l'arbre) | Glashow, Iliopoulos, Maiani (1970) |

### Caractéristiques structurelles

| Observation | Comment le modèle l'explique |
|-------------|---------------------------|
| Bornes holographiques d'entropie | L'information vit sur les frontières où les patchs se recouvrent |
| Violations de Bell <= borne de Tsirelson | Les corrélations quantiques sont les corrélations cohérentes maximales |
| Pas de référentiel privilégié | Aucun observateur privilégié ne peut en définir un |
| La mesure affecte les résultats | La « mesure » est une synchronisation des observateurs |
| Espace-temps 4D depuis un écran 2D | $\mathrm{Conf}(S^2) \sim SO(3,1)$ ; les symétries de l'écran deviennent celles de l'espace-temps |

### Prédictions testables nouvelles

| Prédiction | Comment tester | Statut |
|------------|-------------|--------|
| Ratios log-gap de Casimir (9/4, 5/2, 4, 9/2, 6) | Futures mesures lattice des probabilités de secteurs de bord SU(3) | Sans paramètres (théorie des groupes) |
| Échelle conditionnelle d'entropie d'étiquette centrale Z₆ = log₂ 6 ≈ 2.585 bits | Sous l'ansatz d'un ensemble uniforme à six étiquettes centrales | Continuation phénoménologique |
| Ratios de raies BH $E_k/E_2 = \ln k / \ln 2$ | Analyse gamma de PBH | Motif arithmétique sans paramètres |
| Largeur de raie fractionnelle indépendante de masse ~3-5% | Profils de raies PBH | Prédiction de forme |
| Exposants de Yukawa $-\ln y_f / \ln 6$ proches d'entiers | Extraction depuis masses fermioniques | Continuation phénoménologique liée à la structure Z₆ |
| Unification des couplages sans désintégration du proton | Durée de vie proton + données de couplages de précision | Unification géométrique vs algébrique |
| Borne de déviation RG via défaut de Markov | Tests de gravité de précision en régime faible CMI | Suppression exponentielle avec largeur de collier |
| Masses neutrinos : $m_{\nu_3} \approx 3.0$ meV, $m_{\nu_2} \approx 0.50$ meV, $m_{\nu_1} \approx 0.084$ meV | JUNO, DUNE, KATRIN, $\sum m_\nu$ cosmologique | Estimation actuelle d'ordre de grandeur issue de la branche capacité ; ordonnancement normal avec $\sum m_\nu \approx 3.6$ meV |
| Théorie des cordes depuis secteurs de bord | Vérification théorique | Pont mathématique vers Yang-Mills 2D / Gross-Taylor établi ; une complétion en supercorde critique reste une extension ouverte |

### Prédictions dynamiques

Ces prédictions suivent du spectre d'aire discret sous l'hypothèse que les transitions par multiplications entières dominent l'émission de Hawking. Leur contradiction de mesure invaliderait cette règle dynamique spécifique, pas nécessairement tout le cadre.

| Prédiction | Comment tester | Hypothèse requise |
|------------|-------------|---------------------|
| Peigne de spectroscopie d'horizon GW | Empiler les événements LIGO/Virgo à $x_k = \ln k / 8\pi$ ; l'absence de pics cohérents contredit | Dominance des transitions en $k$ entier |
| Peigne de Hawking discret | Les sursauts gamma de PBH doivent montrer $E_k/E_2 = \ln k / \ln 2$ | Dominance des transitions en $k$ entier |

## Critères empiriques de contradiction

L'OPH avance un mélange d'énoncés structurels exacts, de contrôles de cohérence calibrés, et de continuations plus faibles au niveau des branches. Les observations suivantes fournissent des critères de contradiction directe contre le coeur du cadre ou contre les branches indiquées.

### Contradictions physiques directes

- **Masse non nulle du photon ou du graviton** : toute masse confirmée, même minime, casserait la structure de symétrie requise.
- **Désintégration du proton via bosons de jauge** : la détection d'un canal X/Y impliquerait une inclusion dans un groupe simple plus grand, en contradiction avec la structure produit prédite.
- **Information au-delà de la borne de Bekenstein** : si une région stocke plus d'information que la borne holographique, l'image d'écran est fausse.
- **Violations de Bell au-delà de Tsirelson** : des corrélations supra-quantiques indiqueraient une structure hors du cadre.
- **Perte d'unitarité dans l'évaporation des trous noirs** : si l'information est réellement perdue, la structure de cohérence échoue.
- **Absence du peigne de spectroscopie d'horizon GW** : après renormalisation masse/spin, l'absence de motif $x_k = \ln k / 8\pi$ dans un grand jeu de données contredirait le spectre log-entier.
- **Ratios de Casimir incorrects** : si des mesures lattice de précision donnent $\Delta_8/\Delta_3 \neq 9/4$ (par ex. 2.67 ou 5.06), le mécanisme de bord par noyau de chaleur est contredit. Ensemble complet : $\Delta_8/\Delta_3 = 9/4$, $\Delta_6/\Delta_3 = 5/2$, $\Delta_{10}/\Delta_3 = 9/2$, $\Delta_{15}/\Delta_3 = 4$, $\Delta_{27}/\Delta_3 = 6$.
- **Déficit d'entropie Z₆ incorrect** : si l'entropie de secteur de bord est ~6.6 bits au lieu de ~4.0 bits, la structure quotient Z₆ est contredite.
- **Exposants de Yukawa non entiers** : si $-\ln y_f / \ln 6$ ne se regroupe pas près d'entiers, le mécanisme de défaut Z₆ pour la hiérarchie est contredit.

### Contradiction conceptuelle directe

- **Preuve qu'une réalité objective est nécessaire** : la thèse centrale est que la cohérence entre observateurs suffit. Si un phénomène exige *nécessairement* un monde-en-soi indépendant de toute observation, le modèle est contredit.

Cette question est aujourd'hui plus précise qu'un débat philosophique abstrait. Le théorème de Bell (1964) exclut les variables cachées locales ; Kochen-Specker (1967) exclut la non-contextualité globale ; PBR (2012) restreint fortement les lectures purement épistémiques de l'état quantique. Aucun ne prouve à lui seul la non-existence d'une réalité objective, mais ils réduisent systématiquement l'espace disponible pour cette hypothèse.

Le cadre OPH prend l'étape suivante : il n'essaie pas de sauver une réalité objective de plus en plus contrainte ; il la remplace par des perspectives subjectives + contraintes de cohérence. Le « monde objectif » devient une approximation utile (le squelette de cohérence partagé), pas un objet fondamental.

## Statut actuel

Comme résumé dans la matrice de statut plus haut, le cadre établit aujourd'hui une cinématique de Lorentz conditionnelle, une branche semi-classique d'Einstein conditionnelle, une symétrie de jauge compacte, le groupe du MS (via MAR), trois générations, trois couleurs, des porteurs de force sans masse, un programme quantitatif à deux entrées avec distinction explicite entre calibration et sorties indépendantes, la stabilité du proton, ainsi qu'une route OPH->cordes au niveau continuation via les secteurs de bord YM 2D. Les principaux axes d'ingénierie actifs :

1. **Microphysique d'écran** : les modèles de liens quantiques réalisent les prémisses de régulateur et donnent automatiquement la complétion bord-centre. Reste : forcer ou vérifier la sélection de la branche modulaire géométrique OPH dans la limite continue.
2. **Vérification du pont EFT** : le pont modulaire de surface nulle est dérivé sous deux conditions testables ; il faut les vérifier dans des régulateurs UV explicites.
3. **Constante cosmologique** : expliquée structurellement comme $\Lambda = 3\pi/(G \cdot \log \dim \mathcal{H}_{\rm tot})$ via la capacité d'écran ; valeur numérique inférée de l'observation.
4. **Problème CP fort** : une route d'obstruction de collage est esquissée, mais $\theta_{QCD}$ reste une proposition au niveau continuation plutôt qu'une prédiction du noyau récupéré.

## Code

| Script | Description |
|--------|-------------|
| [oph_predict_compare.py](code/particles/oph_predict_compare.py) | Prédiction complète du spectre + comparaison PDG (point d'entrée principal) |
| [particle_masses_stage5.py](code/particles/particle_masses_stage5.py) | Point d'entrée principal du spectre : clôture de jauge, branche critique D10/D11 déléguée au module piloté par l'article, texture Z₆ |
| [particle_masses_paper_d10_d11.py](code/particles/particle_masses_paper_d10_d11.py) | Reconstruction D10/D11 pilotée par l'article : flot littéral de l'appendice, transport UV synchronisé et couche de raccordement du supplément |
| [oph_qcd.py](code/particles/oph_qcd.py) | Running QCD MSbar à 4 boucles et extraction de $\Lambda$ |
| [oph_lattice_su3_quenched_v5.py](code/particles/oph_lattice_su3_quenched_v5.py) | Lattice SU(3) de Wilson trempé pour les ratios de masses hadroniques |
| [oph_no_cheat_audit.py](code/particles/oph_no_cheat_audit.py) | Audit anti-fuite statique + exécution |

### Chapitres du livre

| Chapitre | Titre | Sujet |
|---------|-------|-------|
| [Prologue](book/prologue.md) | Prologue | Mise en place |
| [1](book/chapter-01-consistency.md) | Cohérence | L'accord entre observateurs comme principe fondamental |
| [2](book/chapter-02-lineage.md) | Filiation | Racines historiques des idées holographiques |
| [3](book/chapter-03-screen.md) | L'Écran | Écrans holographiques et bornes d'information |
| [4](book/chapter-04-entropy.md) | Entropie | Thermodynamique et flèche du temps |
| [5](book/chapter-05-algebra.md) | Algèbre | Structure mathématique des observables |
| [6](book/chapter-06-overlap.md) | Recouvrement | Conditions de cohérence et théorème de Bell |
| [7](book/chapter-07-recovery.md) | Récupération | Préservation de l'information et correction d'erreur quantique |
| [8](book/chapter-08-holography.md) | Holographie | AdS/CFT et reconstruction du bulk |
| [9](book/chapter-09-entanglement.md) | Intrication | Géométrie issue des corrélations quantiques |
| [10](book/chapter-10-error-correction.md) | Correction d'erreur | La réalité comme code quantique |
| [11](book/chapter-11-maxent.md) | MaxEnt | Entropie, temps et flot modulaire |
| [12](book/chapter-12-symmetry.md) | Symétrie | Lois de conservation issues de la cohérence |
| [13](book/chapter-13-desitter.md) | De Sitter | L'écran holographique de notre univers |
| [14](book/chapter-14-standard-model.md) | Modèle Standard | Particules issues des contraintes de collage |
| [15](book/chapter-15-relativity.md) | Relativité | Espace-temps issu du temps modulaire |
| [16](book/chapter-16-matter.md) | Matière | Physique classique comme stabilité émergente |
| [17](book/chapter-17-darwin.md) | Lois de Darwin | Les lois comme survivantes évolutionnaires |
| [18](book/chapter-18-synthesis.md) | Synthèse | Mise en cohérence de l'ensemble |
| [19](book/chapter-19-metaphysics.md) | Métaphysique | Qualia et problème difficile |
| [Epilogue](book/epilogue.md) | Épilogue | Une dernière surprise |

## Construire votre propre simulateur de réalité

Une idée qui boucle le projet : une fois qu'on a rétro-ingénieré un système, on peut tenter de le construire.

Le modèle d'écran est une spécification d'un type particulier de système quantique : un réseau invariant de jauge sur une surface 2D fermée, avec la bonne structure d'intrication. En principe, il est constructible.

Les architectures quantiques actuelles approchent les ingrédients nécessaires. L'informatique quantique topologique (notamment l'approche Microsoft avec fermions de Majorana) implémente naturellement des anyons non abéliens qui peuvent réaliser les secteurs de bord requis. Les états de Hall quantique fractionnaire exhibent déjà des signatures d'ordre topologique et de physique de bord compatibles. Des modèles de liens quantiques sur surfaces triangulées ont été simulés sur plateformes à ions piégés et supraconductrices.

Une implémentation à petite échelle ne recréera pas notre univers (l'écran porte ~ $10^{122}$ degrés de liberté), mais elle peut tester si les propriétés émergentes revendiquées apparaissent effectivement : émergence de la symétrie de jauge aux recouvrements, reconstruction d'une géométrie de type espace-temps depuis l'intrication, formation et synchronisation de sous-systèmes de type observateur.

Où seraient les patchs ? Ils émergeraient de la dynamique. Le système évolue sous un hamiltonien local (termes de plaquette qui créent/détruisent des boucles de flux, termes de champ électrique qui pénalisent les grands flux, contraintes de jauge à chaque sommet). Sous cette dynamique, des motifs de corrélation stables se forment spontanément. Les régions à forte information mutuelle interne et couplage externe plus faible deviennent des patchs naturels.

L'évolution hamiltonienne n'est pas le même « temps » que celui vécu de l'intérieur. La dynamique microscopique suit l'horloge de laboratoire externe ; un observateur émergent dans la simulation suit son propre flot modulaire. Les deux horloges sont liées, mais distinctes.

C'est spéculatif, mais pas de la science-fiction. La physique est suffisamment précise pour être simulée. La validation ultime ne serait pas seulement une preuve mathématique, mais l'observation d'une émergence effective d'une géométrie et de lois depuis un petit écran holographique.

## Contribuer

Pour toute correction, suggestion ou ajout, ouvrez une pull request.

## Licence

Ce dépôt est publié sous licence Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).

Vous êtes libre de partager et d'adapter le contenu à des fins non commerciales, avec attribution correcte, et en diffusant les travaux dérivés sous les mêmes termes.

Pour les demandes de licence commerciale : bernhard@floatingpragma.io
