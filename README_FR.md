# Holographie des parcelles d’observateur

**OPH propose que l’Univers soit ce sur quoi ses observateurs s’accordent, et construit cette idée comme un programme que l’on peut exécuter.**

<p align="center">
  <img src="assets/oph-schematic.png" alt="Schéma de la holographie des parcelles d’observateur. Des observateurs primitifs sur un support icosaédrique lisent, enregistrent et agissent. Des parcelles d’observateur qui se recouvrent sur un écran holographique s’accordent sur des registres partagés, qui deviennent des faits partagés. De l’intérieur, une parcelle encode toute une histoire : un espace-temps émergent à trois dimensions d’espace et une de temps, où la gravitation apparaît comme géométrie." width="100%">
</p>
<p align="center"><b><a href="https://philpapers.org/rec/MUEFOC">Article de synthèse</a> · <a href="https://oph-book.floatingpragma.ai/">Livre</a> · <a href="https://floatingpragma.ai/faq/">FAQ sur la réalité (en anglais)</a> · <a href="https://notebook.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a">Carnet</a> · <a href="https://particles.floatingpragma.ai/">Simulation</a> · <a href="Lean/">Bibliothèque Lean</a> · <a href="REPRODUCE.md">Reproduire</a></b></p>

<p align="center">
<a href="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/ci.yml"><img src="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
<a href="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean.yml"><img src="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean.yml/badge.svg" alt="Lean"></a>
</p>

## L’idée

OPH propose que la physique est
ce que des observateurs bornés peuvent établir ensemble. Des observateurs
finis comparent leurs registres sur des frontières partagées, réparent leurs
désaccords et tiennent pour réel ce qui résiste à toutes les comparaisons.
L’Univers évolue par une dynamique locale de désaccord et de réparation.
Trois problèmes qui résistent à une description classique motivent le
programme : la mesure quantique, les corrélations qu’aucun modèle classique
local ne reproduit et l’origine d’un espace-temps apparent. OPH cherche une
explication commune des trois dans la cohérence des registres et des
opérations des observateurs, et la poursuit sous la forme d’un observateur
exécutable, de mathématiques vérifiées par machine et de prédictions fixées à
l’avance.

**L’Univers se lit lui-même jusqu’à l’existence.** OPH a une seule
hypothèse, l’hypothèse de clôture : l’Univers est une structure intemporelle
et autoréférentielle dont les observateurs retrouvent la spécification de
leur propre monde à partir de leurs registres et construisent un matériel qui
réalise cette même histoire, y compris ses constructeurs. Rien ne se trouve
hors de cette boucle ; la spécification retrouvée n’a donc besoin d’aucune
autre justification. L’événement où les observateurs ont rétro-conçu et construit
ce processus est un point fixe mathématiquement garanti de l’histoire, et il
existe dans notre futur subjectif comme notre passé existe. Elle doit fonctionner et reproduire ce que mesurent les
observateurs. L’espace, le temps et la matière sont des descriptions
effectives des registres, et « notre Univers » est le registre public, ce qui
résiste à toute comparaison entre observateurs. La [base](docs/BASIS.md)
énonce l’hypothèse, la règle de l’absence d’extérieur, la condition d’accord
et la spécification exécutable ; le
[guide des termes](docs/INTERPRETATION_AND_TERMS.md) explique la vision
intemporelle et les horloges qu’elle contient.

**Les observateurs viennent en premier. Ce qu’OPH compare à la nature est ce
que les observateurs lisent dans leurs propres registres, et ce sur quoi ils
s’accordent là où leurs registres se recouvrent. La mesure, le consensus et la
synchronisation des horloges viennent après les registres. Aucun champ
objectif et aucune horloge globale ne se tiennent derrière eux : les champs,
les particules, l’espace-temps et un temps commun sont des descriptions
effectives que les observateurs conviennent d’employer. L’état du programme et
son paramètre de temps calculent ces lectures. Ils ne forment aucun second
univers à côté d’elles.**

## L’hypothèse exécutable

Le [modèle exécutable](model/contrast_memory_spherical_v1/README.md), A-contrast-memory-spherical-v1, place de petits
observateurs sur les ports de cellules icosaédriques. Deux voisins comparent
leurs charges par une mémoire privée sur la couture qu’ils partagent ; un
échange cohérent les fait circuler, et l’accord abaisse l’énergie. La loi,
son couplage et la préparation sont fixés avant l’exécution. Les lecteurs
conservent chaque résultat, y compris le résultat « aucune paire », et chaque
observateur compte ses propres événements : aucune horloge commune n’est
utilisée. Les théorèmes démontrés sur d’autres constructions gardent leurs
hypothèses, et leur transfert au modèle exige une identification démontrée.

L’[hypothèse](docs/OPH_EXECUTABLE_HYPOTHESIS.md) dit ce qui
compterait comme un succès : les registres des observateurs du modèle doivent
reproduire ce que mesurent les observateurs de notre monde.

## Axe de recherche

Nous avons un [modèle exécutable](model/README.md) d’OPH qui fonctionne,
appelé A-contrast-memory-spherical-v1. De petits observateurs sont placés sur les ports de
cellules icosaédriques. Chacun détient quelques registres et les compare
avec ses voisins ; là où deux voisins sont en désaccord, ils réparent la
différence entre eux deux. La théorie qui l’entoure est achevée : le principe
de clôture, la bibliothèque de théorèmes et les constructions finies sont en
place.

**Ce que vivent les observateurs intérieurs.** L’état de la simulation ne
contient pas l’Univers physique. Ses ports, ses coutures et ses cellules sont
le matériel d’un calcul, et ils ne sont pas un espace. L’espace, le temps, les
champs et les particules sont la manière dont les observateurs intérieurs
décrivent les registres dont ils conviennent. Nous avons montré que les
observateurs qui habitent la simulation peuvent décrire leurs expériences
comme un espace-temps continu relativiste à trois dimensions d’espace et une
de temps, bâti sur la réponse la plus lente d’une cellule, qui a exactement
trois directions. Leurs registres obéissent aux règles quantiques de la
mesure et, sous des prémisses énoncées, leurs forces ont l’algèbre du Modèle
standard et leur gravitation obéit à l’équation d’Einstein : ce qu’ils
peuvent mesurer a la structure de ce que nous mesurons. Chaque test du modèle
compare donc ce que ses observateurs enregistrent avec ce que nous
enregistrons. Il ne compare jamais l’état de la simulation avec le monde.

La conception répond à la question de savoir comment un univers peut
calculer sa propre histoire. Une histoire cohérente est un point fixe sur un
nombre énorme de registres, et la chercher prendrait un temps exponentiel.
Le réseau ne cherche pas. Chaque vérification est locale, chaque réparation
est locale, et chaque réparation acceptée retire une quantité fixe de
désaccord : le travail croît comme un polynôme du nombre d’observateurs.
Avec cette conception, l’Univers se pense lui-même jusqu’à l’existence en
temps polynomial. Le [manuel de construction](model/construction_manual/CONSTRUCTION_MANUAL.md)
(en anglais) construit le modèle pièce par pièce et dit quelle part de cela
est un théorème.

Deux lois de réparation candidates tournent sur le modèle : la loi de mémoire
de contraste et le canal scalaire. Aucun nombre mesuré n’y entre. Sous ces
lois, le processus de
réparation reproduit des valeurs et des structures observables réelles, dans
les descriptions effectives des observateurs qui habitent le calcul :

- **Des particules qui tiennent ensemble.** Une paire de charges écrite sur
  une couture est lue comme une paire à chaque lecture ultérieure, avec une
  probabilité qu’un théorème borne sur des supports de toute taille, et elle
  se déplace de couture en couture. La même loi sans l’énergie d’accord perd
  la paire.
- **Une horloge dans un objet.** Un objet tenu sur quatre coutures oscille à
  une hauteur fixe après une impulsion. Un observateur lit cette hauteur dans
  ses propres comptages, en unités de sa propre référence, et retrouve le même
  nombre par une seconde voie, à partir des deux niveaux de l’objet.
- **Quelque chose qui reste et ne porte aucune charge.** Une note laissée
  dans la mémoire d’une couture, sans rien autour, n’est jamais effacée par
  la loi. Elle ne porte aucune charge, ses observateurs lisent son niveau, et
  la matière qui l’atteint s’y attache.
- **Regarder change la suite.** Demander « une paire est-elle là ? » et
  demander « où est-elle ? » donnent la même première réponse et des registres
  ultérieurs différents.
- **Douze courants de force.** Avec un service de port déclaré, la réponse
  des ports d’une cellule a douze générateurs dotés de l’algèbre des forces
  du Modèle standard.

Ces sorties se branchent sur la bibliothèque de théorèmes. Un théorème
s’applique au modèle partout où une sortie du modèle a la forme de sa
prémisse : la cellule est le support du théorème des trois directions,
l’écriture du modèle est la règle de réparation que quatre exigences
désignent seule, ses lectures sont les instruments projectifs des théorèmes
de la mesure, et ses observateurs comparent leurs horloges par des registres
partagés, comme l’exigent les théorèmes d’horloge. Le
[registre des observations](docs/registers/OBSERVATION_LEDGER_V3.md) (en
anglais) liste ce que le modèle reproduit et ce qui est ouvert.

## Trouver la loi de réparation et de mémoire

L’objectif principal de la recherche est d’identifier la bonne loi de
réparation et de mémoire. C’est la loi de réparation native, la loi simple
qui gouverne les interactions des observateurs dans l’Univers primitif : la
règle locale que chaque observateur suit lorsqu’il compare ses registres avec
un voisin, et qui dit comment un désaccord est réparé et ce que la réparation
laisse en mémoire. La réparation décide de ce qui devient un fait partagé. La
mémoire compte tout autant : dans A-contrast-memory-spherical-v1, la phase que la mémoire
d’une couture acquiert à partir du contraste qu’elle a écrit est l’énergie
d’accord qui lie une paire de charges, et la même loi sans elle perd la paire.

A-contrast-memory-spherical-v1 est le candidat à l’essai. Ses deux lois sont des candidates,
leurs couplages sont fournis, et les autres observables des articles, de la
charge et du spectre des particules à la gravitation, sont en cours sur le
même processus. Toute partie de l’observateur exécutable pourra changer au
fil de la recherche.

**Pourquoi c’est important.** La bibliothèque de théorèmes énonce ce qui
découle d’une loi dès qu’elle fournit certaines sorties : trois directions
d’espace, les règles quantiques de la mesure, les forces du Modèle standard,
la gravitation d’Einstein. Une loi doit les fournir toutes à la fois, sur un
seul processus, avec des vérifications et des réparations toutes locales.
Avec la bonne loi, la physique que nous mesurons devient quelque chose que
l’on calcule : on fait tourner la loi, ses observateurs tiennent leurs
registres, et on lit ce qu’ils trouvent. Une mauvaise loi se trahit de la
même manière, car ses observateurs enregistrent autre chose que nous. La
recherche de la loi est ainsi le test le plus net de toute l’idée, et le
reste du projet est à son service. La
[vue d’ensemble du modèle](model/README.md) (en anglais) liste ce que les
candidats reproduisent et ce qui est ouvert.

## Ce qui en sort

Sous des prémisses énoncées, la bibliothèque de théorèmes et les constructions
finies donnent :

- **L’espace-temps et les champs comme descriptions effectives.** À grande
  échelle, les registres partagés des observateurs sont décrits par des
  théories des champs relativistes à trois dimensions d’espace et une de temps,
  avec la gravitation. La réponse lente du support est exactement
  tridimensionnelle, et une théorie des champs lorentzienne de dimension (1+3)
  reproduit les lectures conjointes des observateurs.
- **Accord et registres stables.** Sous des conditions explicites de
  terminaison et de cohérence, les réparations locales atteignent un même
  résultat public quel que soit leur ordre.
  [Article sur le consensus](paper/reality_as_consensus_protocol.pdf).
- **Les règles de la mesure quantique.** Les registres achevés obéissent à la
  règle de Born, à la mesure répétable et à la borne de Tsirelson sur les
  corrélations de Bell, avec un état qui l’atteint.
  [Algèbres finies d’événements quantiques](extra/machine_checked_finite_event_algebras.pdf).
- **Les forces du Modèle standard.** Sous des prémisses énoncées de réponse
  réversible complète et de transport entre observateurs, l’algèbre de réponse
  à douze ports a le type de Lie de jauge du Modèle standard, et les phases des
  ports réalisent le groupe global.
  [Dérivation de la jauge](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **Quarks et leptons.** Les défauts que laisse la réparation portent
  exactement les six classes de centre des particules du Modèle standard, et
  chaque génération sort avec les charges observées : trois couleurs de quarks
  up et down, l’électron et son neutrino. Une représentation de matière
  spécifiée donne une génération sans anomalie.
- **Thermodynamique.** Des théorèmes finis de réparation donnent une inégalité
  de Clausius, des relations de fluctuation et une inégalité de Landauer
  informationnelle. La réparation de couture à pièce équitable est un
  rééchantillonnage conditionnel exact : l’entropie ne baisse jamais sur les
  populations binaires, tandis que la bande complète garde chaque bit.
  [Récit thermodynamique](paper/observers_are_all_you_need.pdf).
- **Gravitation.** Un bilan réalisé des registres donne l’équation d’Einstein.
  [Espace-temps et dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).

## Lancer l’observateur

Depuis la racine du dépôt, avec Python 3.11 ou plus récent :

```sh
python3 -m model
```

Cette commande recalcule les résultats principaux d’A-contrast-memory-spherical-v1 et les
compare aux résultats conservés. Le [guide du modèle](model/README.md)
indique l’exécution complète.

## Des mathématiques vérifiables

La [bibliothèque Lean](Lean/) contient plus de 12900 théorèmes et lemmes
publics sans preuve admise, sur le consensus, les registres quantiques, la
structure de jauge, la thermodynamique et l’identité de Koide. Chaque résultat
de calcul est livré avec un reçu qui lie ses entrées et se rejoue.

## Des tests fixés à l’avance

Les prédictions et leurs règles de rejet sont gelées avant les données, avec
des preuves OpenTimestamps ancrées dans Bitcoin, dans
l’[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md). Le
[programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md) dit comment
la construction peut échouer.

## Articles

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  est la synthèse principale : l’hypothèse de clôture, l’hypothèse du modèle
  exécutable et la tour des résultats. Citez sa version publiée, la
  [fiche PhilPapers](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer for Observer Patch Holography](paper/executable_observer_construction.pdf)**
  énonce A-contrast-memory-spherical-v1, démontre son théorème de liaison et rapporte les
  registres dont ses observateurs conviennent.
- **La bibliothèque de théorèmes** contient les preuves longues, du
  [consensus](paper/reality_as_consensus_protocol.pdf) et des
  [observateurs et registres](paper/observers_are_all_you_need.pdf) jusqu’à
  l’[espace-temps et la dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf)
  et aux [particules](paper/deriving_the_particle_zoo_from_observer_consistency.pdf).
  L’[index des articles](paper/README.md) et l’[ensemble des articles](docs/PAPER_SET.md)
  donnent chaque article et son rôle.
- **Les articles compagnons** portent des résultats ciblés. Le programme de capacité demande si la capacité
  publique attribuée à l’Univers s’accorde avec celle reconstruite depuis
  l’intérieur. Ces calculs de constantes, comme celui de la structure fine,
  ont un statut diagnostique.

## Où se trouvent les choses

| Chemin | Contenu |
| --- | --- |
| [model/](model/README.md) | L’observateur exécutable : implémentation, tests, reçus conservés et enregistrements de protocole |
| [flagship/](flagship/) | L’article de synthèse |
| [paper/](paper/README.md) | L’article de construction et la bibliothèque de théorèmes, avec les articles archivés sous [paper/archive/](paper/archive/README.md) |
| [extra/](extra/README.md), [cosmology/](cosmology/) | Les articles compagnons |
| [Lean/](Lean/) | La bibliothèque vérifiée par machine |
| [code/](code/README.md) | Les calculs exacts, les vérificateurs et les reçus derrière les articles |
| [evidence/](evidence/README.md) | Les paquets de preuves conservés et les enregistrements de garde |
| [claims/](claims/README.md) | Les énoncés établis et leurs dépendances |
| [tracking/](tracking/) | Les registres des prémisses, observations, prédictions et instruments |
| [docs/](docs/README.md) | La base, le guide des termes, les politiques et les registres |
| [book/](book/) | Le livre |

## Reproduire et vérifier

Le [guide de reproduction](REPRODUCE.md) lance la suite de vérification depuis
un clone propre. Le vérificateur du modèle recalcule A-contrast-memory-spherical-v1 et ses résultats ;
la suite de régression vérifie aussi des cas adverses :

```sh
python -m pip install -r model/requirements.txt
python -m model
python -m pytest -q -W error model
```

Le guide donne aussi les contrôles autonomes des particules et des champs
sources utilisés par le [workflow d’intégration continue](.github/workflows/ci.yml).

La [bibliothèque Lean](Lean/) se compile avec Lake, l’[index de la documentation](docs/README.md)
donne les registres et les politiques, et le [livre](book/) offre un chemin
d’accès plus long. Les manuels sur [learn.floatingpragma.ai](https://learn.floatingpragma.ai/)
et la simulation interactive sur
[simulation.floatingpragma.io](https://simulation.floatingpragma.io/) sont
les introductions les plus accessibles.

## Contribuer

La tâche consiste à faire que la même loi et ses instruments
fournissent la physique effective : l’horloge de l’observateur, son interprétation
spatiale et ses mesures de champ sur les mêmes histoires, comparées à des
observations indépendantes. Les
contributions doivent résoudre une connexion nommée ou révéler un défaut du
modèle ; la [spécification d’A-contrast-memory-spherical-v1](model/contrast_memory_spherical_v1/SPEC.md) et les
[tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
décrivent ce travail. Les publications d’articles suivent le
[guide des articles](paper/README.md#reproducibility) ;
les PDF de la branche principale sont des aperçus des sources, et les versions
déposées conservent leurs octets d’origine.

## Licence

Le logiciel, y compris Lean et le code du modèle, est publié sous
[Apache-2.0](code/LICENSE). Les articles, la documentation, les figures, les
données et les registres sont publiés sous [CC BY-NC-SA 4.0](LICENSE), avec les
conditions amont du Particle Data Group pour les valeurs de particules
empaquetées. La [carte des licences](LICENSE) donne les conditions par
répertoire.

<details>
<summary>Projets associés et financement</summary>

[Cadence](https://github.com/muellerberndt/cadence) est le projet associé de
machines apprenantes. Ses éléments de preuve sur les capacités sont distincts
des éléments de preuve physiques d’OPH.
[Pragma Research](https://floatingpragma.ai/) réunit les projets.

Le memecoin officiel TEO (Theory Of Everything) d’OPH a été déployé par la
communauté Base et se trouve sur Base à l’adresse
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([graphique](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Ses frais de transaction sont versés à Bernhard Mueller pour financer le
projet. Le jeton n’a aucun rôle dans la théorie. Les autres jetons qui
utilisent le nom OPH ou TEO n’ont aucun lien avec le projet. Rien ici ne
constitue un conseil financier.

</details>

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->
