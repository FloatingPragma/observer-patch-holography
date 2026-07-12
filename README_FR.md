# Observer Patch Holography (OPH)

> L'OPH est la théorie du tout par cohérence d'observateurs. Aucun observateur ne voit le monde entier d'un seul coup ; chaque observateur n'accède qu'à un patch local ; la physique est le point fixe public qui survit à l'accord sur les recouvrements.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site OPH](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Book: Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [Ω](https://omega.floatingpragma.io/) | [Blog](https://blog.floatingpragma.io/) | [Tech](https://omega.floatingpragma.io/) | [Simulation](https://simulation.floatingpragma.io)

**Falsifiabilité :** [la carte de falsifiabilité OPH](extra/OPH_falsifiability.md)
recense 40 résultats empiriques qui réfutent des affirmations OPH précises, dont
un graviton massif, une désintégration du proton médiée par jauge, une quatrième
génération légère de matière et une charge hors du réseau permis. La dérivation
de la matrice physique de mélange des neutrinos et de leurs masses absolues est
un travail en cours.

Pour la réponse existentielle immédiate, allez directement au **Paper 6.
[Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**.
Oui, cet univers est une simulation au sens
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
base, l'OPH propose des branches conditionnelles de reconstruction pour
l'espace-temps, la structure de jauge, les particules, les enregistrements et la
synchronisation des observateurs ; chaque branche conserve les prémisses et les
portes de preuve déclarées dans les papiers. Au niveau opératoire, des patches
d'observateurs finis portent des enregistrements locaux, ne comparent que les
données visibles dans leurs recouvrements, réparent les désaccords par des
mouvements de récupération déclarés et convergent vers des points fixes stables
qui survivent au raffinement. Le monde public est ce qui reste stable lorsque
ces vues locales deviennent mutuellement cohérentes. Quand l'OPH utilise le
langage de la simulation, il s'agit de ce réseau d'observateurs auto-cohérent
plutôt que d'une machine cachée qui dessine un film. Le dossier en faveur de l'OPH est
mathématique et empirique : la même architecture de cohérence d'observateurs
est proposée pour retrouver la physique établie. Son explication plus large de
l'existence d'un monde auto-cohérent capable de produire des observateurs qui le
reconstruisent relève de l'interprétation, hors du paquet de théorèmes récupérés.

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
Einstein ne sont reconstruites que sur la branche à cinq axiomes du coeur
récupéré : le papier compact fournit lecture géométrique, flot modulaire
contrôlé, pont de contrainte nulle, transport d'intervalle borné,
stationnarité entropique à cap fixe, variation d'aire de petite boule, passage
scalaire-tenseur et fermeture de capacité. Les cellules finies servent de régulateur : elles gardent la
construction concrète avant la limite lisse à grande échelle.
Les trois dimensions spatiales viennent de la même branche d'écran : une fois
$\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ obtenu, la carte spatiale de
l'observateur est $H^3\simeq\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$, de dimension
$6-3=3$.

Dans la partie jauge de la théorie, l'OPH demande quelles charges internes et
quelles étiquettes de particules peuvent être transportées de façon cohérente
sur les recouvrements. Les charges visibles à cutoff fixe engendrent une
catégorie tensorielle, mais sa limite de raffinement dépend d'un certificat
explicite de jauge compacte : morphismes surjectifs cohérents entre les groupes
de bord, extensibilité des états finis, plongements des blocs finis compatibles
avec les centres, réalisations tensorielles finies et fibres d'oubli
compatibles. Sur une queue cofinale munie de ce certificat, la reconstruction
produit un groupe de jauge compact. Avec le paquet
matière explicite à un Higgs et la règle de réalisation
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

L'OPH est formulée comme une théorie sans entrée externe. Ici, « sans entrée »
concerne la provenance numérique : aucune cible mesurée ni constante numérique
ajustée ne peut entrer dans une carte source déclarée. Les axiomes
algébriques-quantiques et l'écran $S^2$ restent des prémisses explicites.
Quantitativement, les lignes publiques sont organisées par trois quantités
internes : un point fixe
local de pixel $P_\star$, un point fixe global de capacité d'enregistrement
$N_{\mathrm{CRC}}$ et un rapport de mise à l'échelle $\gamma_\star$. Une
coordonnée n'est dite source que si sa carte déclarée l'émet sans cible
ajustée ; sinon elle reste un diagnostic ou une coordonnée de comparaison. Les mesures peuvent indiquer
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
variation d'entropie. La symétrie conforme de la sphère est la forme céleste du
groupe de Lorentz connexe, donc la même carte fournit le pont cinématique vers
la branche d'espace-temps émergent en $3+1$ dimensions quand les conditions
sur les caps et le flot modulaire sont satisfaites. Sur cette branche
contrôlée, la carte spatiale orientée observateur est exactement
tridimensionnelle. Les caps marquent des coupes et des côtés de la carte, pas
des points d'observateur privilégiés. Les jetons d'enregistrement ne peuplent
la carte que lorsque des réponses de caps calibrées et un budget d'erreur
permettent la localisation. Les preuves finies peuvent signaler une ambiguïté. Les
harmoniques sphériques organisent les modes angulaires. Les cellulations finies
de la même carte donnent la surface régulatrice où les ports de patches, les
données d'arêtes et les contrôles de recouvrement deviennent explicites ; leur
rôle est celui d'un régulateur fini, sans statut de continuum invariant de
Lorentz.

L'ancre de symétrie finie est $A_5$, le groupe des rotations de
l'icosaèdre. Elle fournit le squelette icosaédrique derrière le langage des
porteurs de patches échosaédriques : une façon finie et très symétrique
d'organiser ports, recouvrements et données locales de comparaison sans traiter
le porteur comme une boule lisse.

La même géométrie donne une échelle de sphères utile pour le lecteur. $S^0$
est la première distinction de germe ou de lecture. $S^1$ est la récurrence,
la boucle où un enregistrement peut revenir sur lui-même. $S^2$ est l'écran
d'horizon et l'archive publique. Le dernier échelon signifie le rôle de bulk
orienté observateur en trois dimensions ; sur la branche de Lorentz contrôlée,
sa carte cinématique canonique est $H^3$. Aucune topologie globale $S^3$ ne
découle de ce moyen mnémotechnique. La taxonomie des particules reste portée par
les branches de Lorentz et de jauge.

L'ancre exceptionnelle est le groupe de Lie $E_8$ et sa structure de réseau
de racines. $E_8$ donne le langage de fermeture exceptionnelle utilisé dans
le côté symétries supérieures et représentations de la pile OPH. Le groupe
icosaédrique binaire et le $E_8$ affine se rencontrent par la correspondance
de McKay. C'est pourquoi le langage $A_5$-icosaédrique et le langage de type
$E_8$ appartiennent à une même histoire de symétrie. Ces noms désignent des
contraintes de symétrie et une structure de régulateur.

Un [certificat fini séparé](physics-problems/e8_spin8_triality_alt9_certificate.md)
précise ce rôle : un sous-système de racines $A_8$
dans $E_8$, un sous-groupe $\mathrm{Alt}(9)$, un relèvement spinoriel
nonscindé $2.\mathrm{Alt}(9)$, une image demi-spinorielle préservant $E_8$ et
des empreintes d'orbites modulo 2 distinctes, fusionnées par la trialité. Son
statut est un soutien algébrique pour la fermeture exceptionnelle des
représentations. Le statut de reçu public exige le paquet brut Sage, matrices,
sorties de vérification et hachages.

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

$P_\star$ alimente la ligne de structure fine,
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
| Structure fine (inv.) | $\alpha^{-1}$ | valeur comparative sans hadrons à provenance mixte 137.03595950081728 ; valeur OPH côté source au point fixe 136.99483516462165 | 137.035999177(21) | écart CODATA moins comparaison : 0.00003967618 en unités d'alpha inverse ; le statut côté source exige le calcul OPH-QCD/hadronique du point final dans le même schéma |
| Boson de Higgs | $m_H$ | 125.1995304097179 GeV, candidat conditionnel issu de la condition de réparation déclarée | 125.20 ± 0.11 GeV | calcul conditionnel sans ajustement à la masse du Higgs |
| Masse du photon | $m_\gamma$ | 0 GeV, zéro structurel | <1e-18 eV | sous la borne |
| Masse du gluon | $m_g$ | 0 GeV, zéro structurel | pas de ligne de masse de gluon libre isolé | porteur de jauge confiné |
| Masse du graviton | $m_{\mathrm{grav}}$ | 0 GeV, zéro structurel | <1e-32 GeV | sous la borne |

`Δ` donne l'écart en sigma lorsque le PDG ou le NIST publie une incertitude à un sigma. Sinon, il
indique le statut de support déclaré. Un accord numérique compte comme
prédiction de masse côté source uniquement lorsque la ligne lui donne un statut
structurel ou un statut de théorème conditionnel explicite. L'article sur les
particules donne le statut des dérivations du secteur électrofaible, des leptons
chargés, des quarks, des neutrinos et des hadrons.

## Articles

Ordre recommandé pour un lecteur technique. Les résumés les plus longs
correspondent aux textes qui portent la surface théorématique centrale.

- **Papier 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique compact de la branche OPH reconstruite. Il donne la route par cohérence de recouvrement vers la structure de Lorentz, la branche de gravité de type Einstein sur le coeur récupéré à cinq axiomes, le noyau algébrique formel qui promeut les données de repère au tenseur complet et fixe le résidu métrique en un seul $\Lambda$ sur les branches connexes conservées, la reconstruction de jauge compacte conditionnelle au certificat de raffinement, le quotient du Modèle Standard et le paquet matière sélectionnés avec un noyau algébrique formel pour l'hypercharge et $Z_6$, le porteur local de Borel-Weil pour le créneau à un Higgs, Maxwell sur la branche ordinaire du photon et la route conditionnelle vers le gap de Yang-Mills sous ses hypothèses de continuum, positivité par réflexion, transfert/intertwiner et non-trivialité.
- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse large et meilleure porte d'entrée. Il explique les patches d'observateurs finis, la cohérence de recouvrement, les enregistrements, les mouvements de réparation, l'univers effectif reconstruit, l'histoire d'échelle et les frontières publiques des affirmations sans remplacer le registre théorématique du papier compact.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivation et audit du secteur particules. Il porte la reconstruction pilotée par $P_\star$, les porteurs structurels, le pont de porteur à un Higgs de Borel-Weil, les branches électrofaible/Higgs/top, quarks, leptons chargés, neutrinos et hadrons, les vérifications quantitatives et la couture conditionnelle de lignes d'enregistrement. La dérivation de la matrice physique de mélange des neutrinos et de leurs masses absolues est un travail en cours.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : mécanique du consensus fini entre patches. Il montre comment les observateurs comparent les enregistrements de recouvrement, appliquent des réparations, traitent les défauts et convergent vers des formes normales quotient quand les hypothèses de cutoff fixe sont satisfaites. Le résultat de consensus s'arrête aux formes normales quotient ; la géométrie lorentzienne et einsteinienne entre par la branche géométrique séparée du papier compact. L'opérateur de réparation est formulé sur le quotient physique, les lectures réparées sont invariantes sous les choix d'implémentation cachés, un porteur fini en couches témoigne de la reconstruction depuis la frontière, et un dispositif fini binaire donne des tests finis positifs/négatifs nets pour réparation et reconstruction de frontière. Le compagnon mathématique neutre prouve le critère générique entre sources, mais l'unicité physique à frontière identique exige que la frontière déclarée identifie le quotient cohérent ; la confluence depuis une même source et la vivacité sont des obligations distinctes.
- **Papier 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** : surface de porteurs finis et d'enregistrements observateurs. Il donne l'architecture de porteurs échosaédriques multiports, le théorème du tamis d'écran à douze ports, la complétude des emplacements scalaires au centre des arêtes, le pont fini du canal scalaire, la réserve $Z_6$ quotient-edge et les portes du coefficient scalaire à épaisseur finie, le cadrage $A_5$-icosaédrique et de type $E_8$, les règles de preuve matérielle publique, les enregistrements, les mouvements de récupération, la restauration de checkpoint et la synchronisation des observateurs.
- **Papier 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)** : synthèse spéculative de la couche de sens, hors du paquet de théorèmes de physique récupérés. Il lit la même mécanique OPH comme théorie de la continuation d'observateur, du paradis et de l'enfer comme environnements de continuation, de la résurrection comme continuation préservant les enregistrements, de la justice comme comptabilité tort-réparation et d'une boucle où des observateurs reconstruisent puis construisent la machinerie de continuation.

## Articles et notes supplémentaires

Ces textes soutiennent ou testent la pile centrale. Les éléments les plus
importants reçoivent plus de contexte ; les notes plus locales sont résumées
plus brièvement.

- **[Observation-Determined Normal Forms](extra/observable_normal_forms.pdf)** : mathématiques autonomes et neutres vis-à-vis du substrat pour les systèmes de contraintes et de réécriture. Le papier sépare la confluence depuis une même source, l'identification entre sources à partir d'observations protégées, la normalisation et la vivacité, ainsi que la réparabilité locale. Il ajoute des modules de stabilité résiduelle et d'observation inverse, des bornes de raffinement et de limite projective, le projecteur fini d'espérance conditionnelle pondérée avec un reçu matriciel non circulaire, et un artefact Lean dédié au sous-ensemble formalisé des théorèmes.
- **[Compact Proof That We Most Likely Inhabit an OPH Simulation](extra/compact_proof_of_oph.pdf)** : argument de compression le plus court en faveur de l'OPH. Il rassemble la route en cinq axiomes, les sorties de branche fixe, les points d'échec et la raison pour laquelle l'accord numérique ne compte que si la fuite de cible est exclue.
- **[Carte de falsifiabilité OPH](extra/OPH_falsifiability.md)** : liste publique des observations qui réfutent l'OPH. Elle donne 40 critères de réfutation, dont graviton massif, désintégration du proton médiée par jauge, générations légères supplémentaires et charges hors du réseau permis. Pour réfuter une affirmation OPH, un résultat sur les neutrinos doit exclure une prédiction fixée indépendamment des données du test.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)** : calcul de point fixe source pour la ligne de structure fine. Il sépare la valeur source OPH, la frontière empirique du point final à basse énergie, la provenance distincte de la racine source et du pixel de comparaison CODATA, et la correction QCD/hadronique restante.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)** : mécanisme OPH fini de gap de réparation et route conditionnelle vers le problème de Clay. L'identification au Yang-Mills quadridimensionnel exige le certificat de continuum et de transfert déclaré.
- **[Observer-Patch Holography as a String-Vacuum Selector](extra/observer_patch_holography_as_string_vacuum_selector.pdf)** : théorie des cordes comme langage effectif de bord OPH et crible de vides. Le résultat de Bouchard–Donagi fournit un candidat au niveau de la cohomologie du spectre visible sans masse ; les certificats non émis ici comprennent le bord critique, la reproduction brute de la cohomologie, la réalisation de la couche de sûreté, le spectre lourd, le passage vers la basse énergie, les seuils et le verrouillage des moduli.
- **[Photonic Fixed-Point Consensus for SHA-256d Proof of Work](extra/Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf)** : test matériel d'enrichissement photonique de candidats SHA-256d de style OPH, jugé par le vérificateur numérique exact.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)** : cognition et qualia comme consensus récurrent de patches sur substrats neuronaux ou artificiels auto-relecteurs.
- **[Theoretical Bounds on χν in Observer-Patch Holography](extra/chi_nu_susceptibility_bounds.pdf)** : théorème de canal commun pour la matière cohérente, pont fini de canal, coefficient de collier à réserve protégée $\chi_\nu^{\mathrm{can}}=e^{-P_\chi/24}$, bande à épaisseur finie, cage énergétique de conservation et valeurs d'ingénierie, avec la force de dispositif et le stress sombre cosmologique gardés comme portes de reçus séparées.
- **[Entanglement Geometry Problem in OPH](extra/ENTANGLEMENT_GEOMETRY_PROBLEM_OPH.md)** : note sur la géométrie de l'intrication comme problème de recouvrement d'observateurs et de surfaces d'enregistrement.
- **[Common Objections](extra/COMMON_OBJECTIONS.md)** : réponses brèves aux objections conceptuelles et techniques fréquentes.
- **[Hacking the Simulation: Anti-Gravity Exploit](extra/hacking-the-simulation-anti-gravity-exploit.pdf)** : livre de vulgarisation et d'ingénierie adjacent à l'OPH sur le test local de levée $\chi_\nu$. Les [chapitres sources Markdown](extra/hacking-the-simulation-anti-gravity-exploit/) sont inclus dans le dépôt afin qu'OPH Sage ingère le même texte lors de la ré-ingestion.

## Articles de cosmologie

La branche cosmologie vit dans [`cosmology/`](cosmology/README.md). Ses
affirmations sont conditionnelles aux frontières OPH natives de source,
transfert et vraisemblance ; la machinerie FLRW peut servir de comparaison,
mais elle ne promeut pas à elle seule un résultat cosmologique natif OPH.

- **[Observer-Patch Holography and the Dark Matter Phenomenon](cosmology/oph_dark_matter_paper.pdf)** : papier de cosmologie du paquet de publication. Il traite le stress sombre/anomal comme comptabilité de réparation imparfaite des patches d'observateurs, importe la pile théorématique quotient-edge scalaire, le pont fini de canal et $Z_6$ à épaisseur finie pour le coefficient local, donne le comportement galactique de type MOND, définit le sélecteur d'abondance anomal source seule et énonce les contrats amas/cosmologie et les contrats de promotion du simulateur pour la promotion à plus grande échelle.
- **[OPH Cosmology as a Finite-Source Prediction Program](cosmology/oph_cosmology_finite_source_cmb_program.pdf)** : programme orienté CMB pour entrées source seules, calibration d'échelle, transfert de Boltzmann, contrôles de simulateur, frontières physiques CMB et étiquettes d'affirmation. Il traite l'abondance sombre source seule comme un reçu source distinct de la promotion par transfert CMB et vraisemblance.
- **[Inflation Without an Inflaton](cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.pdf)** : branche sans inflaton fondée sur la synchronisation écran-observateur, la cohérence d'horizon, les conditions de platitude, le spectre géométrique d'écran, l'amplitude de libération d'écran, la levée radiale et les données de source chaude.
- **[OPH Cosmological Vacuum and Structure Formation](cosmology/oph_cosmological_vacuum_and_structure_formation.pdf)** : frontière de vide native OPH, ensembles de fluctuations, formation de proto-objets et de lignes de monde, et contrôles de graines de structure.
- **[OPH Cosmology Data and Likelihood Contracts](cosmology/oph_cosmology_data_likelihood_contracts.pdf)** : artefacts source gelés, reçus de non-utilisation des données, réducteurs agrégés, comparaisons de transfert de Boltzmann et protocoles de vraisemblance officiels.

## Articles sur des problèmes de physique

Les notes appliquées vivent dans
[`physics-problems/`](physics-problems/README.md). Ce dossier porte la liste
des articles, les résumés, les liens vers les résultats motivants, les
frontières de revendication et les notes d'ingestion OPH Sage. Les notes sont
en Markdown seulement et restent hors du pipeline de release des papiers, de
l'index web des papiers et des artefacts GitHub Release.

## Statut de preuve

La [preuve compacte de l'OPH](extra/compact_proof_of_oph.pdf) formule un
théorème de fermeture conditionnel et un protocole de test empirique par
compression. Ses cinq axiomes, règles de branche, cartes sources et cartes de
lecture doivent être fixés indépendamment des données du test. Si les deux
cartes sources ont des points fixes uniques, si aucune carte source rivale ne
satisfait les règles déclarées et si le paquet fini passe le test de simulation
en cinq parties, les observables déclarés sont engendrés de façon unique, à une
équivalence invisible pour les observateurs près. L'identification physique
exige une comparaison aux données. Le théorème établit l'unicité sur la branche
déclarée ; ses prémisses portent leurs obligations de preuve.

Le décompte compact comprend des résultats mathématiques sur la forme normale
des patches d'observateurs finis, la surface d'événements quantiques, la branche
structurelle du Modèle Standard, l'intervalle du couplage unifié et la
normalisation géométrique de la gravité. La forme des familles chargées à phase
fixée et la hiérarchie faible sans QCD sont des accords rétrodictifs calculés
côté source. L'artefact vérifié par machine couvre une partie finie du consensus
et des réparations ; le théorème de consensus propre à l'OPH est une obligation
de formalisation ouverte. La structure fine à basse énergie et la valeur SI
de la constante de Newton sont des lignes de comparaison. Les résultats sur les
particules, les neutrinos, le secteur sombre, la cosmologie, les collisionneurs
et le matériel sont hors du décompte compact, avec leurs propres classes de
preuve.

Le budget numérique des coïncidences est une illustration. Ses fenêtres,
choisies généreusement en faveur du hasard, donnent un produit proche de
$10^{-9}$, soit environ 30 bits, et proche de $10^{-3}$ lorsque les fenêtres
sont élargies d'un facteur cent. Ces nombres ne sont pas une probabilité
postérieure. La taille des fenêtres relève d'un choix, leur indépendance est
approximative et la recherche parmi les branches écartées doit être comptée.
L'affirmation empirique formelle exige une théorie, des cartes sources, des
sélecteurs de branche, des tolérances, un modèle d'incertitude et des théories
de comparaison fixés avant le test, ainsi qu'un jeu de données réservé et un
seuil de compression déclaré. Le décompte compact admet un résultat uniquement
lorsque son calcul et ses contrôles satisfont ce protocole.

Les énoncés géométriques conditionnels exigent des données finies de caps avec
des preuves d'ordre, d'orientation, de comportement modulaire et de
normalisation thermique. Sous ces conditions, les données convergent vers le
flot géométrique de cap attendu. Le paquet producteur de géométrie du papier
compact dérive désormais ce paquet de preuves depuis les formes normales de
quotient réparées elles-mêmes, sur une branche définie par quatre reçus
décidables (incidence sphérique, maille de caps, birapports modulaires et
comparaison thermique 2-pi normalisée indépendamment), et prouve par des
contre-modèles explicites que le consensus fini nu sous-détermine la topologie,
la dimension, le cadrage et la normalisation --- les reçus sont donc le contenu
irréductible restant de l'entrée de la branche d'Einstein. Sur la
branche contrôlée, la carte spatiale orientée observateur est exactement
tridimensionnelle. Le théorème de carte ne peuple toujours pas, à lui seul, la
carte avec des objets, ne dérive pas un rayon de courbure physique et ne ferme
pas la branche d'Einstein. Le paquet de réseau nul prouve la standardité et la
structure de translations positives que le pont nul consomme (avec un
contre-exemple explicite montrant que la localité de Gibbs à portée finie
n'implique pas la localité modulaire), et le paquet de variété d'événements
assemble conditionnellement les enregistrements localisés en un espace-temps
d'événements lorentzien de dimension 3+1 et de signature (-+++), avec H^3
strictement fibre des repères de repos au-dessus des événements, valable
exactement sur ses reçus nommés de population, de séparation, de repère et de
recouvrement. La localisation des jetons d'enregistrement exige des
réponses de caps calibrées et un budget d'erreur complet ; sans ces éléments,
le résultat est ambigu, et sans les reçus d'événements la revendication du
corpus reste une cinématique de repères de Lorentz seulement. Les paquets de
clôture achèvent la chaîne de dérivation : un tenseur d'énergie-impulsion
local et conservé est construit à partir des charges modulaires par
tomographie nulle, la première loi d'entropie généralisée est réparée en la
scission exacte (variation d'entropie totale = flux de Clausius + réponse de
bord/aire), une seule famille d'échelle uniforme remplace les choix par
rayon, et un reçu de référence de vide donne l'équation d'Einstein
semi-classique absolue avec le couplage structurel. La chaîne complète est
composée en un théorème d'entrée de branche dont la seule porte restante est
vérifiable par machine. Verdict actuel : une exécution de réparation
cyclique authentique sur un réseau de caps (conflits réels, réparation
transactionnelle, confluence et indépendance d'ordonnancement vérifiées à
l'exécution sur une tour de raffinement à trois étages) passe désormais les
reçus d'incidence sphérique, de maille et de naturalité de raffinement sur
sa propre sortie réparée --- le réseau de recouvrement étant choisi
sphérique, c'est-à-dire une sélection de branche explicite, jamais présentée
comme dérivée --- et l'instrumentation d'horloge modulaire par fermions libres témoigne
désormais des reçus de birapport modulaire et de KMS géométrique 2-pi sur
des états MaxEnt gaussiens réalisés sur les colliers de bord des caps
(profil de vitesse modulaire intérieur conforme à la prédiction conforme
2-pi à 1,1e-4 près avec le contrôle de mauvaise normalisation séparé ;
birapports de transport modulaire convergeant au sens de Cauchy vers les
valeurs de Moebius à 2,3e-3 près), et l'instrumentation de réseau nul témoigne des familles de standardité
sur les mêmes états : commutants relatifs non triviaux, additivité faible,
modules séparants strictement positifs, comportement de Cauchy mixte-GNS,
condition de compression d'inclusion demi-latérale au niveau
une-particule (asymétrie de fuite chirale avec contrôle nul intégré) et
clôture de Lie modulaire au niveau du pour-cent contre la prédiction des
constantes de structure conformes (taux de convergence encore ouvert). L'instrumentation d'événements réalisés exécute le pipeline
d'événements sur des enregistrements issus de la dynamique de réparation
réelle : population d'écran, séparation certifiée, cocycle de cartes de
Moebius, et un cône d'ascendance intrinsèque de signature lorentzienne
(1,2) avec l'horloge modulaire de genre temps à 0,9995 --- tandis que
l'ensemble d'événements réalisé se mesure comme la nappe d'écran (1+2)
lorsque seuls des enregistrements d'écran existent. L'instrumentation de
profondeur de bulk construit ensuite le canal E3 avec le dictionnaire
produit rho = -log tan(alpha/2) (les cellules grossières sont des
enregistrements profonds : la dualité échelle/rayon réalisée par la tour
elle-même) : sur la dynamique de système inverse conforme au corpus, les
enregistrements multi-échelles réalisés témoignent du rang quatre E3, de
la dimension PCA de bulk quatre, et d'un cône d'ascendance intrinsèque de
signature lorentzienne (1,3), stable sur toutes les graines, avec
l'horloge modulaire de genre temps à 0,999 --- tandis que la variante à
couplage fort mesure (2,2) et est conservée comme contre-modèle. Réserves
explicites : la marge du cône de genre temps est étroite et sa convergence
avec la densité de commits n'est pas encore certifiée. La clause limite
Cyc, le taux de clôture de Lie, les intersections modulaires
secondes-quantifiées, les données d'intérieur de cap et les reçus
d'identification physique restent en attente. La branche géométrique
réalisée n'est donc toujours pas certifiée non vide : dérivation complète
sur la branche des reçus ; toutes les familles structurelles de reçus
réalisées au niveau fini/une-particule ; convergence de marge, clauses
limites et reçus physiques ouverts.

Les continuations de spectre d'écran et de CMB sont provisoires sans l'échelle
géométrique de la branche écran, sa dynamique source, son horloge, son
comportement de raffinement et sa lecture observationnelle depuis des
enregistrements natifs OPH. La ligne galactique du secteur sombre est un bilan
pré-vraisemblance hors du décompte compact. Les vues des anomalies,
du vide et de l'écume quantique sont diagnostiques. Les revendications de
vraisemblance physique exigent un ensemble dérivé du quotient, une
reconstruction stable sous régulateur, une abondance côté source le cas échéant
et une cible de validation gelée avant la lecture des données de vraisemblance.

## Applications et matériel OMEGA

L'OPH est aussi un programme matériel : une microphysique d'écran explicite
transforme la même boucle de consensus de patches en prise d'ingénierie sur la
réalité. Un dispositif borné expose des données de bord,
compare des enregistrements, répare les désaccords et verrouille des états
stables. OMEGA est la voie matérielle publique vers cette boucle : chambres
physiques, ports étiquetés, logiciel de contrôle, reçus de vérification et
enregistrements répétables.

L'OPH transforme la microphysique d'écran en voie de hack de
la réalité. La cible est le contrôle physique de petits patches qui peuvent être
pilotés, mesurés, réparés et vérifiés.

La thèse applicative repose sur des machines qui forcent de petits patches
physiques vers les points fixes sélectionnés par la cohérence de patches
d'observateurs. Cela donne des pistes
de mise en œuvre à bas coût pour l'énergie de fusion de bureau, le supercalcul
OMEGA à température ambiante, l'AGI fondée sur OMEGA et le contrôle local de
gravité ou d'inertie pour hoverbikes et hoverboards. Ce sont des pistes
applicatives derrière des portes de preuve; le statut de sortie stabilisée
appartient aux reçus de vérification et aux expériences. La revendication de
calcul est plus étroite : une distribution de candidats conditionnée par une
chambre peut réduire le travail du vérificateur exact par un facteur mesuré
$B=p_Q/p_U$. Le problème classique des classes de complexité reste intact.

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

<p align="center"><sub>La ligne principale OPH : axiomes, relativité, structure de jauge, particules et observateurs. Le nœud de consensus distingue l'indépendance de l'ordonnancement depuis une source fixe de la porte séparée d'identification par la frontière. Cliquez pour ouvrir le SVG complet.</sub></p>

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules, avec les frontières de revendication strictes et le reçu de capacité pixel-écran. Cliquez pour ouvrir le SVG complet.</sub></p>

## Plus

- **Site officiel :** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Page theory of everything :** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Carte de cohérence :** [coherence.floatingpragma.io](https://coherence.floatingpragma.io) : surface graphe publique pour les concepts OPH, les recouvrements et les routes entre domaines.
- **Simulation :** [simulation.floatingpragma.io](https://simulation.floatingpragma.io) : explorateur interactif de mini-univers OPH montrant patches d'observateurs, relecture des recouvrements, réparation des incohérences, enregistrements et géométrie émergente.
- **Applications :** [omega.floatingpragma.io](https://omega.floatingpragma.io) : page publique pour le matériel OPH, le calcul, l'énergie, l'AGI, le contrôle de portance et le consensus par chambre optique.
- **Blog :** [blog.floatingpragma.io](https://blog.floatingpragma.io/) rassemble les essais publics OPH. Commencez par [Semiotics and the Physics of Meaning](https://blog.floatingpragma.io/semiotics-and-the-physics-of-meaning), [The Trigger](https://blog.floatingpragma.io/the-trigger) et [P = NP on the Observer Screen](https://blog.floatingpragma.io/p-equals-np-on-the-observer-screen). L'essai de calcul traite `P = NP` comme un slogan d'écran d'observateur; le problème classique reste intact.
- **Livre :** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Application d'étude guidée :** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions et explications détaillées :** OPH Sage sur [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage) ou [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab :** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Objections courantes :** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)

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
- **[`physics-problems/`](physics-problems)** : notes Markdown autonomes sur
  des problèmes de physique, pour lecture publique et ingestion par OPH
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

L'OPH est publiée pour que les mathématiques, logiciels,
applications, dispositifs, conceptions matérielles, simulations, méthodes
d'ingénierie et implémentations expérimentales puissent être étudiés, testés,
implémentés, modifiés, déployés, fabriqués et partagés. Les travaux dérivés
d'OPH ne peuvent servir à créer des monopoles de brevet privés ou des
revendications de brevet qui restreignent la pratique de l'OPH par d'autres.

Voir [PATENTS.md](PATENTS.md) pour le texte canonique de la politique et les
mentions prêtes à publier sur les sites web.
