# Observer Patch Holography (OPH)

> L'OPH est la théorie du tout par cohérence d'observateurs. Aucun observateur ne voit le monde entier d'un seul coup ; chaque observateur n'accède qu'à un patch local ; la physique est le point fixe public qui survit à l'accord sur les recouvrements.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site OPH](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io) | [applications](https://omega.floatingpragma.io/) | [blog OPH](https://blog.floatingpragma.io/) | [carte de cohérence](https://coherence.floatingpragma.io/)

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
OPH, un monde construit à partir de points de vue locaux qui gardent des
enregistrements, comparent ce qu'ils peuvent voir en commun, réparent les
désaccords et convergent vers les motifs stables que tous les observateurs
peuvent partager. Aucun ordinateur extérieur n'a besoin de calculer les
positions des particules image par image. Le temps appartient à ces
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
langage de la simulation, il s'agit de ce réseau d'observateurs auto-cohérent
plutôt que d'une machine cachée qui dessine un film. Le dossier en faveur de l'OPH est
mathématique et empirique : la même architecture de cohérence d'observateurs
retrouve la physique établie et explique pourquoi il existe un monde capable de
produire des observateurs qui le reconstruisent.

Dans la pile de papiers, un patch d'observateur désigne un objet algébrique
abstrait avec algèbre accessible, état, algèbre d'enregistrements, interfaces
visibles de recouvrement, instruments de réparation et données de checkpoint.
Un patch de support est une carte géométrique de cet objet, par exemple un cap
sur $S^2$ ou un diamant causal. Un patch porteur est une réalisation physique
ou numérique des mêmes statistiques d'interface visible et d'enregistrement, à
erreur déclarée. Cette séparation empêche la théorie de dépendre d'une analogie
matérielle particulière. La thèse plus forte selon laquelle l'information et le
calcul sont ontologiquement premiers relève de l'interprétation, sauf si une
branche fournit un discriminateur empirique propre.

La plupart des théories commencent en supposant l'espace-temps, les champs
quantiques et une liste de constantes. L'OPH commence un cran plus tôt, avec
des patches d'observateurs finis algébriques-quantiques dont les descriptions
doivent s'accorder sur leurs recouvrements. Dans la partie relativité de la
théorie, le consensus fini fournit des formes normales quotient. La lecture
d'espace-temps ordinaire en 3+1 dimensions et l'équation de gravité de type
Einstein exigent la branche géométrique séparée : lecture géométrique, flot
modulaire contrôlé, pont de contrainte nulle, stationnarité entropique à cap
fixe, variation d'aire de petite boule, passage scalaire-tenseur et fermeture
de capacité. Les cellules finies servent de régulateur : elles gardent la
construction concrète avant la limite lisse à grande échelle.
Les trois dimensions spatiales viennent de la même branche d'écran : une fois
$\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ obtenu, la carte spatiale de
l'observateur est $H^3\simeq\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$, de dimension
$6-3=3$.

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
$N_{\mathrm{CRC}}$ et un rapport de mise à l'échelle $\gamma_\star$. Le calcul
source les produit comme valeurs de point fixe. Les mesures peuvent indiquer
sur quelle branche nous sommes, mais les valeurs source doivent venir des
calculs de point fixe. Les lignes à fermeture empirique sont marquées
ci-dessous. La discussion d'échelle détaillée est rassemblée une seule fois
ci-dessous dans **Géométrie, Symétrie Et Échelle**.

Dans toute l'OPH, la forme de preuve reste la même. Une affirmation doit
s'appuyer sur des patches bornés de type observateur, avec état local,
frontières explicites, relecture, enregistrements, mouvements de retour ou de
réparation, et dossiers de preuve publics. Le test décisif est l'histoire
invariante des patches d'observateurs, vérifiable sur les recouvrements, plutôt
qu'une présentation favorite, un choix de coordonnées ou une trace
d'implémentation.

## Le piège de l'espace-temps

Le premier obstacle conceptuel est que l'OPH traite l'espace-temps comme une
description issue des observateurs, plutôt que comme le contenant dans lequel la
réalité se déroule. L'espace et le temps y sont les descriptions stables,
orientées observateur, qui apparaissent lorsque de nombreuses perspectives
finies peuvent devenir mutuellement cohérentes.

C'est particulièrement important pour le temps. Dans le langage ordinaire, le
temps ressemble à un fleuve de fond qui continuerait à couler même s'il n'y
avait personne. L'OPH rejette cette image. À la base, il y a des observateurs,
des enregistrements, des changements dans ces enregistrements et des règles qui
font s'accorder les enregistrements qui se recouvrent. Le temps est l'ordre
qu'un observateur donne aux changements de ses propres enregistrements. Le
temps public est la partie de cet ordre qui peut être synchronisée avec
d'autres observateurs. En ce sens précis, le temps est subjectif : il appartient
d'abord au flux d'enregistrements d'un observateur. Il reste pourtant contraint.
Une mauvaise horloge, un faux souvenir ou une histoire incohérente échouent
lorsqu'ils ne peuvent pas s'accorder avec le reste du réseau d'enregistrements.

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
La carte de repos de l'observateur est
$H^3\simeq\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$, donc exactement
tridimensionnelle sur cette branche. Les harmoniques sphériques organisent les
modes angulaires. Les cellulations finies de la même
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
| Structure fine (inv.) | $A_{\alpha_U}^{\mathrm{fp}}=\alpha_{\mathrm{root}}^{-1}+\alpha_U(P_\star)$ | prédiction OPH source sans hadrons 137.03595950081728; tronc d'audit racine seul 136.99483516462165 | 137.035999177(21) | correction QCD/hadronique ouverte, dans le même schéma, pour le point final : 0.00003967618 en unités d'alpha inverse, soit environ 2.9e-7 relatif |
| Boson de Higgs | $m_H$ | 125.1995304097179 GeV, candidat conditionnel sur la surface de réparation déclarée | 125.20 ± 0.11 GeV | ligne conditionnelle à porte de réparation sans cible |
| Masse du photon | m_γ | 0 GeV, zéro structurel | <1e-18 eV | sous la borne |
| Masse du gluon | m_g | 0 GeV, zéro structurel | pas de ligne de masse de gluon libre isolé | porteur de jauge confiné |
| Masse du graviton | m_grav | 0 GeV, zéro structurel | <1e-32 GeV | sous la borne |

`Δ` donne l'écart en sigma lorsque le PDG ou le NIST publie une incertitude à un sigma. Sinon, il
indique le statut de support déclaré. Les accords numériques sur des surfaces
cibles ou témoins ne sont pas des prédictions publiques de masse source seule,
sauf si la ligne déclare un statut structurel ou un statut théorématique
conditionnel explicite. Le statut complet du code particules, y compris les
lignes $W/Z$, leptons chargés, quarks et masses absolues de neutrinos retenues,
est généré dans
[`code/particles/FINAL_END_TO_END_PREDICTIONS.md`](code/particles/FINAL_END_TO_END_PREDICTIONS.md)
et [`code/particles/EXACT_NONHADRON_MASSES.md`](code/particles/EXACT_NONHADRON_MASSES.md).

## Articles

Ordre recommandé pour un lecteur technique. Les résumés les plus longs
correspondent aux textes qui portent la surface théorématique centrale.

- **Papier 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique compact de la branche OPH reconstruite. Il donne la route par cohérence de recouvrement vers la structure de Lorentz, la branche de gravité de type Einstein avec ses dépendances de lecture géométrique et de stationnarité entropique, le noyau algébrique formel qui promeut les données de repère au tenseur complet et fixe le résidu métrique en un seul \(\Lambda\) sur les branches connexes conservées, la reconstruction de jauge compacte, le quotient du Modèle Standard et le paquet matière sélectionnés avec noyau algébrique formel hypercharge/\(Z_6\), Maxwell sur la branche ordinaire du photon et la route conditionnelle vers le gap de Yang-Mills sous ses hypothèses de continuum, positivité par réflexion, transfert/intertwiner et non-trivialité.
- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse large et meilleure porte d'entrée. Il explique les patches d'observateurs finis, la cohérence de recouvrement, les enregistrements, les mouvements de réparation, l'univers effectif reconstruit, l'histoire d'échelle et les frontières publiques des affirmations sans remplacer le registre théorématique du papier compact.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivation et audit du secteur particules. Il porte la reconstruction pilotée par \(P_\star\), les porteurs structurels, les branches électrofaible/Higgs/top, quarks, leptons chargés, neutrinos et hadrons, les vérifications quantitatives et la couture conditionnelle de lignes d'enregistrement.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : mécanique du consensus fini entre patches. Il montre comment les observateurs comparent les enregistrements de recouvrement, appliquent des réparations, traitent les défauts et convergent vers des formes normales quotient quand les hypothèses de cutoff fixe sont satisfaites. Le résultat de consensus s'arrête aux formes normales quotient ; la géométrie lorentzienne et einsteinienne entre par la branche géométrique séparée du papier compact. L'opérateur de réparation est formulé sur le quotient physique, les lectures réparées sont invariantes sous les choix d'implémentation cachés, un porteur fini en couches témoigne de la reconstruction depuis la frontière, et le dispositif Rule-90 donne des tests finis positifs/négatifs nets pour réparation et reconstruction de frontière.
- **Papier 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** : surface de porteurs finis et d'enregistrements observateurs. Il donne l'architecture de porteurs échosaédriques multiports, le théorème du tamis d'écran à douze ports, la complétude des emplacements scalaires au centre des arêtes, le pont fini du canal scalaire, la réserve \(Z_6\) quotient-edge et les portes du coefficient scalaire à épaisseur finie, le cadrage \(A_5\)-icosaédrique et de type \(E_8\), les règles de preuve matérielle publique, les enregistrements, les mouvements de récupération, la restauration de checkpoint et la synchronisation des observateurs.
- **Papier 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)** : synthèse de la couche de sens. Il lit la même mécanique OPH comme théorie de la continuation d'observateur, du paradis et de l'enfer comme environnements de continuation, de la résurrection comme continuation préservant les enregistrements, de la justice comme comptabilité tort-réparation et de la boucle où des observateurs reconstruisent puis construisent la machinerie de continuation.

## Articles et notes supplémentaires

Ces textes soutiennent ou testent la pile centrale. Les éléments les plus
importants reçoivent plus de contexte ; les notes plus locales sont résumées
plus brièvement.

- **[Compact Proof That We Most Likely Inhabit an OPH Simulation](extra/compact_proof_of_oph.pdf)** : argument de compression le plus court en faveur de l'OPH. Il rassemble la route en cinq axiomes, les sorties de branche fixe, les points d'échec et la raison pour laquelle l'accord numérique ne compte que si la fuite de cible est exclue.
- **[Carte de falsifiabilité OPH](extra/OPH_falsifiability.md)** : liste publique de résultats qui tueraient l'OPH. Elle nomme 40 modes d'échec durs, dont graviton massif, désintégration du proton médiée par jauge, générations légères supplémentaires, outliers du réseau de charges et exclusions neutrino de la branche OPH.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)** : calcul de point fixe source pour la ligne de structure fine. Il sépare la valeur source OPH, la frontière empirique du point final à basse énergie, la provenance à deux \(P\) entre racine source et pixel de comparaison CODATA, et la correction QCD/hadronique restante.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)** : mécanisme OPH fini de gap de réparation et route conditionnelle vers le problème de Clay. L'identification au Yang-Mills quadridimensionnel exige le certificat de continuum et de transfert déclaré.
- **[Observer-Patch Holography as a String-Vacuum Selector](extra/observer_patch_holography_as_string_vacuum_selector.pdf)** : théorie des cordes comme langage effectif de bord OPH et crible de vides. Le témoin Bouchard-Donagi est un candidat de cordes conventionnel, pas un raccourci autour des portes de sélection OPH natives.
- **[Photonic Fixed-Point Consensus for SHA-256d Proof of Work](extra/Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf)** : test matériel d'enrichissement photonique de candidats SHA-256d de style OPH, jugé par le vérificateur numérique exact.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)** : cognition et qualia comme consensus récurrent de patches sur substrats neuronaux ou artificiels auto-relecteurs.
- **[Theoretical Bounds on \(\chi_\nu\) in Observer-Patch Holography](extra/chi_nu_susceptibility_bounds.pdf)** : théorème de canal commun pour la matière cohérente, pont fini de canal, coefficient de collier à réserve protégée \(\chi_\nu^{\rm can}=e^{-P_\chi/24}\), bande à épaisseur finie, cage énergétique de conservation et valeurs d'ingénierie, avec la force de dispositif et le stress sombre cosmologique gardés comme portes de reçus séparées.
- **[Entanglement Geometry Problem in OPH](extra/ENTANGLEMENT_GEOMETRY_PROBLEM_OPH.md)** : note sur la géométrie de l'intrication comme problème de recouvrement d'observateurs et de surfaces d'enregistrement.
- **[Tests IBM Quantum Cloud](extra/IBM_QUANTUM_CLOUD.md)** : notes de tests matériels en secteur réduit pour des expériences IBM Quantum orientées OPH.
- **[Common Objections](extra/COMMON_OBJECTIONS.md)** : réponses brèves aux objections conceptuelles et techniques fréquentes.
- **[Hacking the Simulation: Anti-Gravity Exploit](extra/hacking-the-simulation-anti-gravity-exploit.pdf)** : livre de vulgarisation et d'ingénierie adjacent à l'OPH sur le test local de levée \(\chi_\nu\). Les [chapitres sources Markdown](extra/hacking-the-simulation-anti-gravity-exploit/) sont inclus dans le dépôt afin qu'OPH Sage ingère le même texte lors de la ré-ingestion.

## Articles de cosmologie

La branche cosmologie vit dans [`cosmology/`](cosmology/README.md). Ses
affirmations sont conditionnelles aux frontières OPH natives de source,
transfert et vraisemblance ; la machinerie FLRW peut servir de comparaison,
mais elle ne promeut pas à elle seule un résultat cosmologique natif OPH.

- **[Observer-Patch Holography and the Dark Matter Phenomenon](cosmology/oph_dark_matter_paper.pdf)** : papier de cosmologie du paquet de publication. Il traite le stress sombre/anomal comme comptabilité de réparation imparfaite des patches d'observateurs, importe la pile théorématique quotient-edge scalaire, le pont fini de canal et \(Z_6\) à épaisseur finie pour le coefficient local, donne le comportement galactique de type MOND, définit le sélecteur d'abondance anomal source seule et énonce les contrats amas/cosmologie et les contrats de promotion du simulateur pour la promotion à plus grande échelle.
- **[OPH Cosmology as a Finite-Source Prediction Program](cosmology/oph_cosmology_finite_source_cmb_program.pdf)** : programme orienté CMB pour entrées source seules, calibration d'échelle, transfert de Boltzmann, contrôles de simulateur, frontières physiques CMB et étiquettes d'affirmation. Il traite l'abondance sombre source seule comme un reçu source distinct de la promotion par transfert CMB et vraisemblance.
- **[Inflation Without an Inflaton](cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.pdf)** : branche sans inflaton fondée sur la synchronisation écran-observateur, la cohérence d'horizon, les conditions de platitude, le spectre géométrique d'écran, l'amplitude de libération d'écran, la levée radiale et les données de source chaude.
- **[OPH Cosmological Vacuum and Structure Formation](cosmology/oph_cosmological_vacuum_and_structure_formation.pdf)** : frontière de vide native OPH, ensembles de fluctuations, formation de proto-objets et de lignes de monde, et contrôles de graines de structure.
- **[OPH Cosmology Data and Likelihood Contracts](cosmology/oph_cosmology_data_likelihood_contracts.pdf)** : artefacts source gelés, reçus de non-utilisation des données, réducteurs agrégés, comparaisons de transfert de Boltzmann et protocoles de vraisemblance officiels.

## Articles sur des problèmes de physique

Des textes Markdown autonomes pour des problèmes de physique résolus
vivent dans [`physics-problems/`](physics-problems/README.md). Ils agissent
comme articles supplémentaires pour la lecture publique et l'ingestion par OPH
Sage. Ils n'entrent pas dans le pipeline de rendu TeX/PDF ou de publication.

- **[Fusion plasma et confinement](physics-problems/plasma_fusion.md)** :
  fusion comme paquet théorématique de registre de réparation OPH. Le texte
  définit le registre de réparation de fusion, H-mode comme branche de
  contraction de collier de bord, les ELM comme cycles obstruction/reset,
  Lawson comme projection scalaire d'énergie, Hydrosahedron comme spécialisation
  acoustique porteur/contrôle, et les promotions DD/chaleur/charge/puissance
  nette comme niveaux de reçus séparés.
- **[Supraconductivité à haute température](physics-problems/high_temperature_superconductivity.md)** :
  haut $T_c$ comme réparation de charge $2e$ plus confluence de phase
  $U(1)_Q$. Il énonce les prédicats cuprates, pnictides/chalcogénures,
  trempe sous pression, hétérostructures et conception inverse sans les
  transformer en recettes expérimentales.
- **[Universalité des solides amorphes à basse température](physics-problems/low_temperature_amorphous_universality.md)** :
  l'atténuation des verres à basse température comme branche source seule
  saturée de registres rigides. Le texte traite chaque matériau comme un patch
  borné auto-relecteur, avec ports visibles à la déformation, enregistrements
  à deux bassins, mouvements de réparation ou de tunneling, et reçus
  acoustiques ou thermiques publics, puis dérive
  $\lambda/\ell=3.2441\times10^{-3}$ sans utiliser les mesures comme entrées
  source.
- **[États de Hall quantique fractionnaire](physics-problems/fractional_quantum_hall.md)** :
  phases de Hall fractionnaires comme formes normales de bord/holonomie,
  récupération de la forme matricielle abélienne $K$, raffinement hiérarchique,
  conditions de secteurs de réparation non abéliens et théorème de
  non-sélection à $5/2$.

## Statut de preuve

Aucune théorie physique n'est prouvée à 100 % au sens mathématique. Une théorie
physique gagne sa crédibilité en dérivant beaucoup de faits indépendants depuis
peu d'hypothèses, en gardant les cibles mesurées hors de ses cartes sources et
en exposant des façons claires d'échouer. Notre preuve compacte la plus forte
est [Disclosure Day : preuve compacte de l'OPH](extra/compact_proof_of_oph.pdf).
Elle donne le chemin le plus court vers l'idée que l'OPH est probablement
correcte, tandis que la pile complète de papiers porte les dérivations, les
frontières de revendication et les obligations de preuve.

Une sortie OPH finie garde son statut d'origine. Renommer un compte de capacité
en masse, une archive en rayonnement, un spectre de réparation en spectre
physique ou un seuil de reconstruction en temps de Page laisse la revendication
au même niveau. Une promotion physique exige une lecture séparée, un chemin de
calibration, un registre de résidus, des contrôles et une cible de validation
gelée.

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

Les continuations de spectre d'écran et de CMB sont provisoires sans l'échelle
géométrique de la branche écran, sa dynamique source, son horloge, son
comportement de raffinement et sa lecture observationnelle depuis des
enregistrements natifs OPH. Les vues du secteur sombre, des anomalies, du vide
et de l'écume quantique sont diagnostiques sauf lorsque les reçus source
pertinents passent ; les revendications de vraisemblance physique exigent un
ensemble dérivé du quotient, une reconstruction stable sous régulateur, une
abondance source seule le cas échéant et une cible de validation gelée avant la
lecture des données de vraisemblance.

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
des pistes applicatives vivent dans [`APPLICATIONS.md`](APPLICATIONS.md).

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

<p align="center"><sub>La ligne principale OPH : axiomes, relativité, structure de jauge, particules et observateurs. Cliquez pour ouvrir le SVG complet.</sub></p>

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte actuelle de la voie particules, avec les frontières de revendication strictes et le reçu de capacité pixel-écran. Cliquez pour ouvrir le SVG complet.</sub></p>

## Plus

- **Site officiel :** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Page theory of everything :** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Carte de cohérence :** [coherence.floatingpragma.io](https://coherence.floatingpragma.io) : surface graphe publique pour les concepts OPH, les recouvrements et les routes entre domaines.
- **Applications :** [omega.floatingpragma.io](https://omega.floatingpragma.io) : page publique pour le matériel OPH, le calcul, l'énergie, l'AGI, le contrôle de portance et le consensus par chambre optique.
- **Blog :** [blog.floatingpragma.io](https://blog.floatingpragma.io/) rassemble les essais publics OPH. Commencez par [Semiotics and the Physics of Meaning](https://blog.floatingpragma.io/semiotics-and-the-physics-of-meaning), [The Trigger](https://blog.floatingpragma.io/the-trigger) et [P = NP on the Observer Screen](https://blog.floatingpragma.io/p-equals-np-on-the-observer-screen). L'essai de calcul traite `P = NP` comme un slogan d'écran d'observateur; le problème classique reste intact.
- **Livre :** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Application d'étude guidée :** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions et explications détaillées :** OPH Sage sur [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage) ou [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab :** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Objections courantes :** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **Note IBM Quantum :** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

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
- **[`physics-problems/`](physics-problems)** : textes Markdown autonomes de
  problèmes de physique résolus, pour lecture publique et ingestion par OPH
  Sage.

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
