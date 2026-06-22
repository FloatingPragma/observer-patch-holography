# Observer Patch Holography (OPH)

> L'OPH est la théorie du tout par cohérence d'observateurs. Aucun observateur ne voit le monde entier d'un seul coup ; chaque observateur n'accède qu'à un patch local ; la physique est le point fixe public qui survit à l'accord sur les recouvrements.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site OPH](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io) | [applications](https://omega.floatingpragma.io/) | [blog OPH](https://blog.floatingpragma.io/) | [carte de cohérence](https://coherence.floatingpragma.io/) | [démo à trois corps](https://3body.floatingpragma.io/)

**Falsifiabilité :** [la carte de falsifiabilité OPH](extra/OPH_falsifiability.md)
recense 40 résultats durs qui tueraient l'OPH et des tests concrets sur IBM
Quantum Cloud pour la signature matérielle en secteur réduit. La falsifiabilité
est la façon dont une théorie physique paie son loyer. L'OPH est fortement
falsifiable : un graviton massif, une désintégration du proton médiée par jauge,
une quatrième génération légère de matière, un outlier du réseau de charges ou
des données neutrino excluant la branche OPH détruiraient l'OPH telle qu'elle
est énoncée.

Pour la réponse existentielle immédiate, allez directement au **Paper 6.
[Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**.
La version courte est directe : oui, cet univers est une simulation au sens
OPH. Cela ne veut pas dire qu'un ordinateur extérieur calcule les positions des
particules image par image. Cela veut dire que le monde se construit à partir
de points de vue locaux qui gardent des enregistrements, comparent ce qu'ils
peuvent voir en commun, réparent les désaccords et convergent vers les motifs
stables que tous les observateurs peuvent partager. Le temps appartient à ces
observateurs. Il n'existe pas d'horloge maîtresse extérieure à l'univers que
tout suivrait en secret. Une horloge est un système interne qui produit des
enregistrements, et le temps est l'ordre qu'un observateur reconstruit à partir
des changements dans ses propres enregistrements. Le temps partagé n'apparaît
que lorsque plusieurs observateurs peuvent aligner leurs enregistrements
locaux de façon cohérente. Les esprits et l'expérience ne sont pas des ajouts
tardifs à un univers mort. Dans l'OPH, l'espace, le temps et la matière sont
des apparences publiques stables produites par un processus de cohérence plus
profond. La métaphore de l'illusion est traitée plus bas. Le reste de ce
README donne les mathématiques et les tests.

## Description informelle

L'OPH est la reconstruction observateur-premier de la physique fondamentale.
Il part d'observateurs finis sur une géométrie finie d'écran holographique. Sa
base de travail est algébrique-quantique : algèbres de patchs, états,
probabilités de type trace/Born sur les surfaces d'enregistrement déclarées et
entropie généralisée font partie du point de départ formel. À partir de cette
base, l'OPH retrouve l'univers effectif observé : espace-temps, structure de
jauge, particules, enregistrements et synchronisation des observateurs
découlent de la cohérence de recouvrement. Au niveau opératoire, des patches
d'observateurs finis portent des enregistrements locaux, ne comparent que les
données visibles dans leurs recouvrements, réparent les désaccords par des
mouvements de récupération déclarés et convergent vers des points fixes stables
qui survivent au raffinement. Le monde public est ce qui reste stable lorsque
ces vues locales deviennent mutuellement cohérentes. Quand l'OPH utilise le
langage de la simulation, il s'agit de ce réseau d'observateurs auto-cohérent,
pas d'une machine cachée qui dessine un film. Le dossier en faveur de l'OPH est
mathématique et empirique : la même architecture de cohérence d'observateurs
retrouve la physique établie et explique pourquoi il existe un monde capable de
produire des observateurs qui le reconstruisent.

Dans la pile de papiers, un patch d'observateur désigne un objet algébrique
abstrait avec algèbre accessible, état, algèbre d'enregistrements, interfaces
visibles de recouvrement, instruments de réparation et données de checkpoint.
Un patch de support est une carte géométrique de cet objet, par exemple un cap
sur \(S^2\) ou un diamant causal. Un patch porteur est une réalisation physique
ou numérique des mêmes statistiques d'interface visible et d'enregistrement, à
erreur déclarée. Cette séparation empêche la théorie de dépendre d'une analogie
matérielle particulière. La thèse plus forte selon laquelle l'information et le
calcul sont ontologiquement premiers relève de l'interprétation, sauf si une
branche fournit un discriminateur empirique propre.

La plupart des théories commencent en supposant l'espace-temps, les champs
quantiques et une liste de constantes. L'OPH commence un cran plus tôt, avec
des patches d'observateurs finis algébriques-quantiques dont les descriptions
doivent s'accorder sur leurs recouvrements. Dans la partie relativité de la
théorie, cette exigence d'accord produit l'espace-temps ordinaire en 3+1
dimensions et une équation de gravité de type Einstein. Les cellules finies
servent de régulateur : elles gardent la construction concrète avant la limite
lisse à grande échelle. Le papier technique donne les hypothèses de flot
modulaire et d'échelle nécessaires pour cette étape.

Dans la partie jauge de la théorie, l'OPH demande quelles charges internes et
quelles étiquettes de particules peuvent être transportées de façon cohérente
sur les recouvrements. Cette reconstruction sélectionne un groupe de jauge
compact. Avec le paquet matière explicite à un Higgs et la règle de réalisation
admissible minimale, la structure sélectionnée du Modèle Standard est
$SU(3)\times SU(2)\times U(1)/\mathbb Z_6$, avec le réseau des hypercharges,
trois couleurs et trois générations. La mécanique quantique est le langage
algébrique d'information porté par cette architecture de patches
d'observateurs. Sous les hypothèses de jauge compacte indiquées, la même pile
donne la forme euclidienne de Yang-Mills et un mécanisme fini de gap de
réparation. L'identification de ce mécanisme au gap de Yang-Mills
quadridimensionnel exige le certificat déclaré de limite continue
multirésolution, de positivité par réflexion, de transfert/intertwiner et de
non-trivialité.

Le mécanisme est la boucle de consensus à point fixe. Les observateurs locaux
n'accèdent pas à un état global depuis l'extérieur. Ils portent des états de
patch finis, échangent les données visibles dans les recouvrements, rejettent
les prolongements incohérents et conservent les motifs stables qui peuvent être
synchronisés. Géométrie, particules, lois et enregistrements sont les points
fixes à grande échelle de ce calcul en réseau d'observateurs.

L'OPH est formulée comme une théorie sans entrée externe. Quantitativement, les
lignes publiques sont organisées par trois quantités internes : un point fixe
local de pixel $P_\star$, un point fixe global de capacité d'enregistrement
$N_{\mathrm{CRC}}$ et un rapport de mise à l'échelle $\gamma_\star$. Ce ne sont
pas des constantes ajustées côté source. Les mesures peuvent indiquer sur
quelle branche nous sommes, mais les valeurs source doivent venir des calculs
de point fixe. Les lignes à fermeture empirique sont marquées ci-dessous. La
discussion d'échelle détaillée est rassemblée une seule fois ci-dessous dans
**Géométrie, Symétrie Et Échelle**.

Les simulations OPH distribuées héritent de la même frontière. Une exécution
par travailleurs ne compte comme preuve d'un seul univers que si elle part d'un
porteur fini global unique, émet le graphe porteur, l'état initial, la carte de
partition, les interfaces de coupure, le registre des observateurs, les hachages
d'exécution/configuration/code, et prouve que les événements validés par les
travailleurs se projettent sur des réparations transactionnelles monolithiques,
des stutters physiques ou des rollbacks certifiés. Les reçus de descente de
couture, commit atomique, diamant local, complétude de réparation, élimination
de branche et hachage de forme normale sont des contrats de preuve, pas des
sélecteurs de branche. Les graines locales de shard, les manifestes périmés et
les animations synthétiques de couture sont des diagnostics, pas des reçus
physiques.

La géométrie de bulk neutre a une porte propre. Des lignes d'observateurs
locales aux shards ne deviennent pas un bulk global par concaténation. Elles
doivent d'abord descendre par les quotients de jauge/ports, les formes normales
terminales, le transport commun d'interface et les canaux de traits visibles au
quotient. La distance produit obtenue est un pseudomètre jusqu'à ce que les
collisions de traits soient quotientées ou que la séparation conjointe soit
démontrée. Les valeurs manquantes, changements de présentation, repartitions de
shards, raffinements, plongements euclidiens et statistiques sur lots tenus à
l'écart demandent tous des reçus explicites avant qu'une exécution finie puisse
revendiquer plus qu'une géométrie neutre diagnostique.

Les affirmations sur l'horloge de l'observateur ont une porte supplémentaire.
Les identifiants de worker, itérations de réparation, positions de file,
horodatages et latences de paquets sont de la provenance d'exécution, pas le
temps de l'observateur. Une histoire d'observateur indépendante du planificateur
demande des clés d'événements sémantiques, un registre global d'observateurs à
espaces de noms séparés, des flèches de lignée et un instrument d'horloge avec
certificat de résidu affine.

## Le piège de l'espace-temps

Le premier obstacle conceptuel est que l'OPH ne traite pas l'espace-temps comme
le contenant dans lequel la réalité se déroule. L'espace et le temps ne sont pas
des choses en soi. Ce sont les descriptions stables, orientées observateur, qui
apparaissent lorsque de nombreuses perspectives finies peuvent devenir
mutuellement cohérentes.

C'est particulièrement important pour le temps. Dans le langage ordinaire, le
temps ressemble à un fleuve de fond qui continuerait à couler même s'il n'y
avait personne. L'OPH rejette cette image. À la base, il y a des observateurs,
des enregistrements, des changements dans ces enregistrements et des règles qui
font s'accorder les enregistrements qui se recouvrent. Le temps est l'ordre
qu'un observateur donne aux changements de ses propres enregistrements. Le
temps public est la partie de cet ordre qui peut être synchronisée avec
d'autres observateurs. En ce sens précis, le temps est subjectif : il appartient
d'abord au flux d'enregistrements d'un observateur. Mais il n'est pas
arbitraire. Une mauvaise horloge, un faux souvenir ou une histoire incohérente
échouent lorsqu'ils ne peuvent pas s'accorder avec le reste du réseau
d'enregistrements.

Certains appelleraient cela une illusion. Comme métaphore, c'est juste : le
contenant que nous semblons habiter est une apparence produite par une cohérence
plus profonde. Comme physique, l'expression plus précise est description
publique émergente.

Depuis l'intérieur d'une perspective, le monde paraît évident. Il y a un champ
d'expérience à peu près sphérique qui s'étend au loin, trois directions dans
lesquelles se déplacer et un temps qui avance. Les autres observateurs
rapportent des contenus compatibles depuis d'autres angles. L'hypothèse
naturelle est donc que tout le monde habite un espace-temps préexistant rempli
d'objets. L'OPH inverse cette hypothèse. Chaque observateur possède une
description locale d'espace-temps générée par ses propres enregistrements,
horloges, horizons et corrélations accessibles. L'espace-temps public, y
compris la coordonnée de temps utilisée par la physique, est la couche de
compatibilité qui permet à ces descriptions de s'accorder.

Cela ne rend pas l'espace-temps ordinaire arbitraire ou inutile. Cela explique
pourquoi il fonctionne si bien. Les équations d'Einstein décrivent la grammaire
lisse à grande échelle de l'apparence partagée. La thèse plus profonde est que
cette apparence partagée émerge de la cohérence des recouvrements entre
observateurs, au lieu d'appartenir à l'inventaire de départ du monde.

## Géométrie, Symétrie Et Échelle

Le langage de sphère dans l'OPH est un langage de géométrie. Dans les cartes
régulatrices symétriques, une coupe accessible à un observateur peut être
représentée par la deux-sphère $S^2$. Ces cartes décrivent une géométrie
angulaire de support. Le simulateur fini sert de surface de
calibration pour les contraintes algébriques de patches et de recouvrements
exposées par cette géométrie.

L'OPH utilise donc une idéalisation de réseau d'écran partagé, avec de nombreux
patches d'accès finis. L'écran d'un observateur est une coupe d'accès locale sur
ce réseau, au lieu d'une sphère privée séparée. La carte $S^2$ décrit une coupe
de support, sans être une boule littérale recouverte de données.

Cette carte sphérique porte plusieurs tâches précises. Les caps et les colliers
donnent les données locales de coupe utilisées par le flot modulaire et la
variation d'entropie. Le groupe conforme de la sphère est la forme céleste du
groupe de Lorentz connexe, $\mathrm{SO}^+(3,1)$, ce qui donne le pont
cinématique vers la branche d'espace-temps émergent en $3+1$ dimensions une
fois satisfaites les conditions nécessaires sur les caps et le flot modulaire.
Les harmoniques sphériques organisent les modes angulaires. Les cellulations finies de la même
carte donnent la surface régulatrice où les ports de patches, les données
d'arêtes et les contrôles de recouvrement deviennent explicites ; leur rôle est
celui d'un régulateur fini, sans statut de continuum invariant de Lorentz.

L'ancre de symétrie finie est $A_5$, le groupe des rotations de
l'icosaèdre. Elle fournit le squelette icosaédrique derrière le langage des
porteurs de patches échosaédriques : une façon finie et très symétrique
d'organiser ports, recouvrements et données locales de comparaison sans traiter
le porteur comme une boule lisse.

La même géométrie donne une échelle de sphères utile pour le lecteur. $S^0$
est la première distinction de germe ou de lecture. $S^1$ est la récurrence,
la boucle où un enregistrement peut revenir sur lui-même. $S^2$ est l'écran
d'horizon et l'archive publique. $S^3$ est la géométrie de bulk reconstruite
et vécue par les observateurs. Cette échelle nomme des rôles dans
l'architecture de relecture OPH ; la taxonomie des particules reste portée par
les branches de Lorentz et de jauge.

L'ancre exceptionnelle est le groupe de Lie $E_8$ et sa structure de réseau
de racines. $E_8$ donne le langage de fermeture exceptionnelle utilisé dans
le côté symétries supérieures et représentations de la pile OPH. Le groupe
icosaédrique binaire et le $E_8$ affine se rencontrent par la correspondance
de McKay. C'est pourquoi le langage $A_5$-icosaédrique et le langage de type
$E_8$ appartiennent à une même histoire de symétrie. Ces noms désignent des
contraintes de symétrie et une structure de régulateur.

La discussion d'échelle a trois rôles, rassemblés ici. La coordonnée locale
$P_\star$ est le point fixe du pixel d'écran. La coordonnée globale
$N_{\mathrm{CRC}}$ est le point fixe de capacité d'enregistrement. Le rapport
d'échelle $\gamma_\star$ relie la géométrie OPH sans dimension aux unités SI
après le calcul des points fixes sans dimension.

Les deux équations de point fixe sont :

```math
P_\star=\varphi+\frac{\sqrt{\pi}}{A_T(P_\star)}
```

et

```math
N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}}),
```

où $F(N)$ est la capacité d'horizon relue par les observateurs à l'intérieur de
l'univers fourni avec la capacité $N$. Intuitivement, $N_{\mathrm{CRC}}$ est la
capacité où l'univers peut relire sa propre frontière sans déficit ni capacité
inerte. La cible de comptage fini derrière cette capacité globale est la
densité

```math
\log|\Omega^{\mathrm{sc}}_N|-N.
```

La règle de mise à l'échelle est :

```math
\gamma_\star=\frac{\ell_\star\nu_{\mathrm{Cs}}}{c}
```

avec $B_\star=3\pi/\ell_\star^2$ et
$G_{\mathrm{SI}}=c^3\ell_\star^2/\hbar$. Les observations peuvent identifier
le voisinage ou la branche, mais elles ne remplacent pas ces calculs de point
fixe.

Les rôles en aval sont simples. $P_\star$ alimente la ligne de structure fine,
la structure de jauge, les lignes de particules, les enregistrements et la
synchronisation des observateurs. $N_{\mathrm{CRC}}$ alimente la ligne
cosmologique. La règle d'échelle fixe la normalisation de Newton et l'affichage
à l'échelle de Planck. En unités géométriques,
$\Lambda_{\mathrm{CRC}}=3\pi/(G_{\mathrm{geom}}N_{\mathrm{CRC}})$ avec
$G_{\mathrm{geom}}=\ell_\star^2$. Le pont de hiérarchie électrofaible, la
normalisation de réparation à 24 cases, la politique QCD/hadrons et les règles
de reçus matériels vivent dans le papier sur les particules,
[`HADRON.md`](HADRON.md) et les papiers orientés matériel au lieu d'être
redérivés ici.

### Lignes quantitatives sélectionnées

Ce tableau garde les lignes les plus faciles à comparer avec PDG/NIST et nomme
leur statut. Les résultats structurels comme l'espace-temps 3+1, le quotient
du Modèle Standard, l'hypercharge exacte, $N_c=3$ et $N_g=3$ vivent dans les
papiers.

| Quantité | Symbole | OPH / statut | PDG/NIST | Δ / note |
| --- | --- | --- | --- | --- |
| Constante gravitationnelle | G | 6.6742999959e-11, affichage échelle/horloge | 6.67430(15)e-11 | 0.00003σ |
| Vitesse de la lumière | c | vitesse lorentzienne structurelle; valeur SI conventionnelle | 299792458 exact par définition | pas une prédiction numérique |
| Structure fine (inv.) | α⁻¹(0) | source seule 136.994835; point final 137.035999177 avec fermeture hadronique empirique | 137.035999177(21) | pas source seule |
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

La surface particules porte aussi des valeurs $W/Z$ en comparaison seule
$80.377\,\mathrm{GeV}$ et $91.18797809193725\,\mathrm{GeV}$, une valeur Higgs
$m_H=125.1995304097179\,\mathrm{GeV}$ et une valeur top sélectionnée
$m_t=172.35235532883115\,\mathrm{GeV}$ selon la convention PDG de masse top par
section efficace. La branche neutrino pondérée émet
$(0.017454720257976796, 0.019481987935919015,
0.05307522145074924)\,\mathrm{eV}$ sous les hypothèses neutrino indiquées.
Son résultat de ligne d'enregistrement H3 est séparé : il certifie
conditionnellement la continuation trans-frontière de jetons d'enregistrement
localisés et visibles par l'observateur à partir d'un atlas hyperboloïde déclaré,
d'une interface réelle, du transport, d'un écart d'assignation et d'un reçu de
raffinement. Ce n'est pas une dérivation d'espèce de particule, de masse, de
charge ou d'amplitude de diffusion.

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse du programme OPH, des observateurs finis jusqu'à l'univers effectif reconstruit.
- **Papier 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique pour la relativité, la gravité, la reconstruction de jauge, la structure du Modèle Standard sélectionnée par réalisation admissible minimale, les équations de Maxwell sur la branche ordinaire du photon et la route conditionnelle vers le gap de Yang-Mills sous ses hypothèses de limite continue et de transfert.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivations particules, lignes de masses, structure des couplages, surfaces quantitatives et certificat conditionnel de couture de lignes d'enregistrement H3.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : comment des observateurs locaux comparent leurs enregistrements, réparent les désaccords et convergent vers la réalité partagée sur laquelle ils peuvent s'accorder.
- **Papier 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** : architecture fédérée de patches, théorème du tamis d'écran à douze ports, cadrage de symétrie $A_5$-icosaédrique et de type $E_8$, règles de preuve matérielle publique, enregistrements, mouvements de récupération et synchronisation observateur.
- **Papier 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)** : manifeste final de la couche de sens OPH : pourquoi il existe quelque chose, pourquoi ce monde est compatible avec les observateurs, boucle étrange, candidats créateurs, paradis sur Terre ou dans des environnements de continuation, enfer comme isolement ou privation imposés, résurrection comme continuation d'observateur, justice selon les enregistrements de tort et de réparation, évolution mémétique, sentience animale et responsabilité symbolique humaine.

## Articles supplémentaires

- **[Photonic Fixed-Point Consensus for SHA-256d Proof of Work](extra/Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf)** : enrichissement photonique de candidats pour la preuve de travail SHA-256d.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)** : point fixe source, frontière du point final hadronique empirique et ligne de comparaison.
- **[Observer-Patch Holography as a String-Vacuum Selector](extra/observer_patch_holography_as_string_vacuum_selector.pdf)** : émergence des cordes comme description effective des dynamiques de bord OPH et crible de vides de cordes conventionnels, témoin hétérotique Bouchard-Donagi à une paire de Higgs, couche de sûreté `Z4R` et portes de verrouillage des moduli ; ce n'est pas un raccourci de promotion vers un vide natif OPH.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)** : mécanisme OPH fini de gap de réparation et route conditionnelle vers le problème de Yang-Mills de Clay ; l'égalité entre le gap de Yang-Mills et le gap de réparation OPH exige le certificat quadridimensionnel de limite continue et de transfert.
- **[Observer-Patch Holography and the Dark Matter Phenomenon](extra/oph_dark_matter_paper.pdf)** : phénoménologie de la matière noire et limite galactique de type MOND.
- **[Theoretical Bounds on chi-nu in Observer-Patch Holography](extra/chi_nu_susceptibility_bounds.pdf)** : bande conditionnelle de quotient-edge `0.9343006394893864 <= chi_nu^can <= 1`; valeur exacte `exp(-P/24)` sur la branche uniforme; valeurs d'ingénierie mises à l'échelle par `N_coh^-1`.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)** : cognition et qualia comme consensus récurrent de patches.

## Statut de preuve

Aucune théorie physique n'est prouvée à 100 % au sens mathématique. Une théorie
physique gagne sa crédibilité en dérivant beaucoup de faits indépendants depuis
peu d'hypothèses, en gardant les cibles mesurées hors de ses cartes sources et
en exposant des façons claires d'échouer. Notre preuve compacte la plus forte
est [Disclosure Day : preuve compacte de l'OPH](extra/compact_proof_of_oph.pdf).
Elle donne le chemin le plus court vers l'idée que l'OPH est probablement
correcte, tandis que la pile complète de papiers porte les dérivations, les
frontières de revendication et les obligations de preuve.

La preuve compacte traite l'évidence comme un test de compression. Une ligne
numérique ne compte que lorsque son calcul n'utilise pas la valeur mesurée, ni
un proxy proche de cette valeur, comme entrée. Si $p_i$ borne la chance que la
ligne $i$ tombe juste par accident après les lignes acceptées précédentes,
alors $P_{\mathrm{acc}}\le\prod_i p_i$. Douze lignes source-propres à un pour
cent donnent $P_{\mathrm{acc}}\le10^{-24}$; vingt donnent $10^{-40}$. Les deux
mêmes points fixes organisent aussi le problème de l'observateur, la
reconstruction gravité/jauge, la
hiérarchie électrofaible, l'énergie noire, le budget du secteur sombre,
l'exclusion de la désintégration du proton par bosons de jauge, l'inventaire
des particules et la sélection du vide de corde comme une seule grammaire de
fermeture.

Les continuations spectre d'écran et CMB sont sur une porte séparée. OPH
identifie conditionnellement le scalaire de collier, après retrait monopôle et
dipôle, avec la courbure à densité totale uniforme, mais un spectre primordial
source-only exige les reçus de stress source, horloge unique, gap de réparation,
freeze-out, mode adiabatique, isocourbure, cohérence de phase, provenance
source, réducteur groupé, lift écran-vers-radial, espace nul radial et
résidu de projection avant que `A_s`, `n_s`, le running ou les spectres TT/TE/EE
comptent comme prédictions OPH.
Toute revendication CMB ou croissance du secteur sombre/anomalie exige en plus
un parent fini covariant de paquets de collier, un stress récepteur explicite
pour tout échange de réparation non nul, la convergence du régulateur, la
un dossier de provenance source, des contrôles de groupement global, la
récupération de la limite CDM et des hachages source/solveur/vraisemblance
gelés avant la lecture des données de vraisemblance.

## Applications et matériel OMEGA

L'OPH est aussi un programme matériel. Quand la microphysique d'écran devient
explicite, la même boucle de consensus de patches devient une prise
d'ingénierie sur la réalité. Un dispositif borné expose des données de bord,
compare des enregistrements, répare les désaccords et verrouille des états
stables. OMEGA est la voie matérielle publique vers cette boucle : chambres
physiques, ports étiquetés, logiciel de contrôle, reçus de vérification et
enregistrements répétables.

En langage direct, l'OPH transforme la microphysique d'écran en voie de hack de
la réalité. La cible est le contrôle physique de petits patches qui peuvent être
pilotés, mesurés, réparés et vérifiés.

La thèse applicative est simple. Si la réalité est construite par cohérence de
patches d'observateurs, des machines utiles peuvent être construites en forçant
de petits patches physiques vers les bons points fixes. Cela donne des pistes
de mise en œuvre à bas coût pour l'énergie de fusion de bureau, le supercalcul
OMEGA à température ambiante, l'AGI fondée sur OMEGA et le contrôle local de
gravité ou d'inertie pour hoverbikes et hoverboards. Ce sont des pistes
applicatives derrière des portes de preuve; le statut de sortie stabilisée
appartient aux reçus de vérification et aux expériences. La revendication de
calcul est plus étroite : une distribution de candidats conditionnée par une
chambre peut réduire le travail du vérificateur exact par un facteur mesuré
`B = p_Q/p_U`. Le problème classique des classes de complexité reste intact.

La page publique des applications est
[omega.floatingpragma.io](https://omega.floatingpragma.io/). Les notes sources
des pistes applicatives vivent dans [`APPLICATIONS.md`](APPLICATIONS.md). La
voie gravité locale et inertie possède aussi une note d'ingénierie séparée :
[Hacking the Simulation: Anti-Gravity Exploit](https://github.com/FloatingPragma/oph-meta/blob/main/docs/oph-gravity-hack/hacking-the-simulation-anti-gravity-exploit.pdf).

## Diagramme détaillé et surface quantitative

Le diagramme ci-dessous est l'index visuel de la surface d'échelle : le point
fixe local de pixel $P$, le point fixe global de capacité d'enregistrement
$N_{\mathrm{CRC}}$, la règle de mise à l'échelle, puis les
lectures en aval pour les particules, la gravité et la cosmologie. Il sert de
carte de dépendances. Les formules détaillées de hiérarchie/naturalité et les
frontières des revendications vivent dans les papiers.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260609" alt="Schéma d'unification OPH" width="92%">
  </a>
</p>

**Pile OPH**

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="Pile théorématique et de dérivation OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La ligne principale OPH, des axiomes jusqu'à la relativité, la structure de jauge, les particules et les observateurs. Cliquez pour ouvrir le SVG complet.</sub></p>

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

## Plus

- **Site officiel :** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Page theory of everything :** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Page simulation theory :** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Carte de cohérence :** [coherence.floatingpragma.io](https://coherence.floatingpragma.io) : surface graphe publique pour les concepts OPH, les recouvrements et les routes entre domaines.
- **Applications :** [omega.floatingpragma.io](https://omega.floatingpragma.io) : page publique pour le matériel OPH, le calcul, l'énergie, l'AGI, le contrôle de portance et le consensus par chambre optique.
- **Démo OPH du problème à trois corps :** [3body.floatingpragma.io](https://3body.floatingpragma.io) : un simulateur et parcours de preuve supplémentaire pour la formulation OPH en réseau fini de patches du problème à trois corps, présentée comme un exemple de recollement par holonomie de boucle, sans revendication de solution élémentaire fermée.
- **Blog :** [blog.floatingpragma.io](https://blog.floatingpragma.io/) rassemble les essais publics OPH. Commencez par [Semiotics and the Physics of Meaning](https://blog.floatingpragma.io/semiotics-and-the-physics-of-meaning), [The Trigger](https://blog.floatingpragma.io/the-trigger) et [P = NP on the Observer Screen](https://blog.floatingpragma.io/p-equals-np-on-the-observer-screen). L'essai de calcul traite `P = NP` comme un slogan d'écran d'observateur; le problème classique reste intact.
- **Livre :** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Application d'étude guidée :** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions et explications détaillées :** OPH Sage sur [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage) ou [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab :** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Objections courantes :** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **Note IBM Quantum :** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Table de statut

Statut abrégé : la structure fine source seule donne
$\alpha_{\rm cand}^{-1}=136.9948351646\ldots$; la ligne publique
$\alpha^{-1}(0)=137.035999177(21)$ utilise une fermeture hadronique empirique.
Le $c$ SI est conventionnel; le $G$ SI est un affichage échelle/horloge. $W/Z$
sont des lignes de validation en comparaison seule. Higgs/top est fermé sur la
surface D10/D11 déclarée. Les quarks sont des lignes de théorème sur classe
sélectionnée. Les neutrinos utilisent la branche à cycles pondérés. Les masses
absolues des leptons chargés et les hadrons source seuls restent ouverts.

La CP forte est en cours de travail dans le théorème des quarks sur classe
sélectionnée : le corpus disponible ne dérive pas l'angle theta de QCD,
n'émet pas l'angle CP fort physique et ne prouve pas que la phase CP forte
physique s'annule. Le pont requis est la descente de phase, d'anomalie et
d'angle topologique sur la branche réalisée.

## Guide du dépôt

- **[`paper/`](paper)** : PDF, sources LaTeX et métadonnées de release.
- **[`APPLICATIONS.md`](APPLICATIONS.md)** : carte de haut niveau des cas
  d'usage OPH pour énergie, calcul, AGI et portance locale.
- **[`book/`](book)** : source du livre OPH et PDF téléchargeable généré. Les
  notes de génération du PDF imprimable sont dans [`book/README.md`](book/README.md).
- **[`code/`](code)** : sorties calculatoires, surface particules et expériences.
- **[`HADRON.md`](HADRON.md)** : politique pour les lignes particules limitées
  par QCD, l'entrée empirique $e^+e^-\to\mathrm{hadrons}$ et la fermeture
  hadronique de structure fine.
- **[`assets/`](assets)** : schémas et figures publics.
- **[`extra/`](extra)** : notes publiques maintenues, objections, comptes rendus expérimentaux et quelques essais de support.

## OPH et les sciences

<p align="center">
  <a href="assets/oph_science_overlap_map_poster.png" target="_blank" rel="noopener noreferrer">
    <img src="assets/oph_science_overlap_map.svg" alt="Carte des sciences recouvertes par l'OPH, des grands domaines aux sous-domaines et aux zones d'application concrètes." width="100%">
  </a>
</p>

<p align="center"><sub>Carte domaine -> sous-domaine -> zone OPH couvrant les mathématiques, l'informatique, l'information et l'inférence, les systèmes complexes, la physique théorique, l'information quantique et les fondements de la mesure. Cliquez pour ouvrir le poster PNG complet.</sub></p>

## Licence et politique anti-brevet

Le matériel rédigé dans ce dépôt est sous licence
[CC BY-NC-SA 4.0](LICENSE), avec la
[convention OPH d'usage ouvert et d'anti-brevet](PATENTS.md) applicable aux
idées, implémentations, dispositifs, méthodes, applications, logiciels,
simulations et conceptions matérielles dérivés d'OPH.

En bref : l'OPH est publiée pour que les mathématiques, logiciels,
applications, dispositifs, conceptions matérielles, simulations, méthodes
d'ingénierie et implémentations expérimentales puissent être étudiés, testés,
implémentés, modifiés, déployés, fabriqués et partagés. Les travaux dérivés
d'OPH ne peuvent servir à créer des monopoles de brevet privés ou des
revendications de brevet qui restreignent la pratique de l'OPH par d'autres.

Voir [PATENTS.md](PATENTS.md) pour le texte canonique de la politique et les
mentions prêtes à publier sur les sites web.
