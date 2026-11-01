# Holographie des parcelles d’observateur

**La physique à partir des registres que les observateurs peuvent partager.**

[Article principal](https://philpapers.org/rec/MUEFOC) · [Livre](https://oph-book.floatingpragma.io/) · [Carnet](https://notebook.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a) · [Bibliothèque Lean](Lean/) · [Reproduire](REPRODUCE.md) · [Pragma Research](https://floatingpragma.io/) · [English](README.md)

[![Mandatory Suite](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml)
[![Lean CI](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml)
[![Runnable observer model](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml)

## Mission

La mission de l’Holographie des parcelles d’observateur (OPH) est de
construire une théorie de la réalité cohérente et complète, formelle et
rigoureuse. Cette théorie n’a pas de point de départ : c’est une boucle
fermée, et elle doit remplir quatre conditions.

1. **Cohérente.** La boucle se referme. Son point fixe existe, construit
   explicitement et vérifié par machine, et une théorie qui possède un modèle
   ne contient aucune contradiction.
2. **Rigide.** Toute grandeur sans dimension est imposée. Une seule valeur
   mesurée fixe les unités, et rien d’autre.
3. **Complète**, au sens
   qu’[Einstein, Podolsky et Rosen](https://doi.org/10.1103/PhysRev.47.777)
   ont donné à ce mot en 1935. Chaque élément de la réalité physique a une
   contrepartie dans la théorie, et chaque grandeur mesurée est dérivée.
4. **Correcte.** Les valeurs dérivées concordent avec la mesure, et les
   nouvelles prédictions sont figées avant d’être testées.

OPH ne prétend pas être la seule théorie dotée de ces propriétés.

## L’idée

OPH propose que la physique est
ce que des observateurs bornés peuvent établir ensemble. Des observateurs
finis comparent leurs registres sur des frontières partagées, réparent leurs
désaccords et tiennent pour réel ce qui résiste à toutes les comparaisons.
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
autre justification. Elle doit fonctionner et reproduire ce que mesurent les
observateurs. L’espace, le temps et la matière sont des descriptions
effectives des registres, et « notre Univers » est le registre public, ce qui
résiste à toute comparaison entre observateurs. La [base](docs/BASIS.md)
énonce l’hypothèse, la règle de l’absence d’extérieur, la condition d’accord
et la spécification exécutable ; le
[guide des termes](docs/INTERPRETATION_AND_TERMS.md) explique la vision
intemporelle et les horloges qu’elle contient.

## L’hypothèse exécutable

Le pari concret est un programme. Une seule loi locale simple de réparation et
de mémoire, exécutée sur un réseau de supports à douze ports, doit produire des
observateurs dont les registres partagés ressemblent à notre monde : des champs
quantiques en interaction en 3+1 dimensions, et la gravitation aux grandes
échelles. Les particules y seraient des nœuds du réseau, des défauts localisés
que les réparations locales ne peuvent pas défaire, et leur masse serait le
coût d’entretien supplémentaire qu’ils imposent.

C’est l’objectif. L’[hypothèse](docs/OPH_EXECUTABLE_HYPOTHESIS.md) dit ce qui
compterait comme un succès, et le [programme d’achèvement](docs/COMPLETION_PROGRAM.md)
découpe le travail en étapes avec des règles d’arrêt : soit les pièces se
rejoignent sur un même processus, soit une preuve montre que toute la classe
de lois échoue.

## Points forts

- **Un observateur qui fonctionne.** Un observateur à douze ports conserve ses
  registres et s’en sert pour choisir sa prochaine réparation. Des observateurs
  qui se chevauchent lisent une référence commune et son horloge, un processus
  quantique fini réunit flux électrique, matière chargée et champ scalaire, et
  les écrans raffinés donnent des limites de champ gaussiennes et des ondes de
  vitesse finie.
- **Un champ scalaire lorentzien à partir des registres.** Sur une extension
  radiale déclarée de l’écran raffiné, une seule action fixe corrélations,
  réponse et horloge modulaire, et une théorie lisse de champ scalaire en
  dimension (1+3) reproduit tout ensemble fini de lectures conjointes. La
  direction radiale et les unités physiques sont des entrées. Le
  [verdict intégré](model/CONSTRUCTION.md#integrated-effective-description)
  donne les détails.
- **Les particules sont la grande pièce manquante.** La théorie finie donne
  les charges d’une génération du Modèle standard et ses couplages de Higgs,
  mais aucune masse. Toutes les lois de masse testées sur le modèle exécutable
  manquent la hiérarchie mesurée des leptons, et le candidat de défaut évident
  sur le support natif se défait après quelques réparations. Une campagne
  finie exclut aussi la forme la plus simple de charge protégée : aucune
  grandeur cachée du support ne commute avec toutes les réparations. Trouver un défaut
  réellement stable est l’objet de
  l’[article sur les particules](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)
  et des [notes sur la réponse en masse](docs/research/OPH_PARTICLE_MASS_RESPONSE.md).
- **La physique standard à partir d’un accord fini.** Sous des prémisses
  énoncées, des théorèmes finis donnent à la spécification des structures
  familières : des registres publics indépendants de l’ordre des réparations ;
  la règle de Born et la borne de Tsirelson sur les registres achevés ; un
  support tridimensionnel dans sa réponse lente ; l’algèbre de Lie de jauge du
  Modèle standard et une représentation de matière sans anomalies
  conditionnelle ; un ordre causal de dimension (1+3) avec un temps propre lu
  dans les comptages ; et l’équation d’Einstein à partir d’un bilan nul.
- **Des mathématiques vérifiables.** La [bibliothèque Lean](Lean/) contient
  plus de 12700 théorèmes et lemmes publics sans preuve admise, et chaque reçu
  de calcul conservé est rejoué en intégration continue.
- **Des tests fixés à l’avance.** Les prédictions et leurs règles de rejet
  entrent dans l’[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md)
  avant les données de comparaison admissibles, avec des preuves OpenTimestamps,
  pour la plupart ancrées dans Bitcoin ; la prédiction de dispersion pour l’orbite des ports
  primitifs en fait partie. Le [registre des postdictions](docs/POSTDICTION_LEDGER.md)
  et le [programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md)
  séparent les comparaisons avec des mesures des tests des mathématiques.
  L’[audit inter-échelles](docs/EMPIRICAL_STATUS.md) teste des
  correspondances numériques partagées et consigne les identifications
  physiques nécessaires avant que des comparaisons séparées puissent compter
  comme éléments de preuve pour une même construction.

**L’accord avec l’ensemble des statistiques de mesure observées n’est pas
établi.** Les théorèmes valent sous des prémisses énoncées, et l’observateur
exécutable est un candidat jugé selon des critères d’acceptation déclarés.

## Lancer l’observateur

Depuis la racine du dépôt, avec Python 3.11 ou plus récent :

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

Cette commande exécute l’expérience minimale d’enregistrement et de
rétroaction et affiche les scores des prévisions et l’empreinte du résultat.
Le [guide du modèle](model/README.md) couvre les observateurs plus riches, les
horloges, les champs et les comparaisons aux données ; le
[guide de construction](model/CONSTRUCTION.md) énonce la portée vérifiée et
les entrées restantes.

## Fondements mathématiques

- **L’accord et les registres stables.** Sous des conditions explicites de
  terminaison et de cohérence, les réparations locales atteignent un même
  résultat public quel que soit leur ordre. D’autres résultats contrôlent
  quels registres survivent au raffinement.
  [Article sur le consensus](paper/reality_as_consensus_protocol.pdf).
- **Les observations quantiques.** Dans l’algèbre quantique finie énoncée, des
  résultats vérifiés couvrent les mesures répétables, les registres classiques
  que les opérations ultérieures ne peuvent perturber, l’incertitude et la
  borne de corrélation de Bell, avec un état qui l’atteint.
  [Algèbres finies d’événements quantiques](extra/machine_checked_finite_event_algebras.pdf).
- **Trois directions spatiales.** Le support fixe à douze ports et sa règle de
  réparation moyenne donnent une réponse lente exactement tridimensionnelle de
  géométrie euclidienne, que retrouvent les registres de courant d’un observateur.
  [Construction du support](code/a5_closure/README.md).
- **La symétrie de jauge.** Sous des prémisses énoncées de réponse réversible
  complète et de transport entre observateurs, l’algèbre de réponse à douze
  ports a le type de Lie de jauge du Modèle standard, et une représentation de
  matière spécifiée donne une génération sans anomalies.
  [Dérivation de jauge](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **La thermodynamique.** Une seule construction finie de réparation donne les
  quatre lois, les relations de fluctuation et la borne d’effacement de
  Landauer sous des prémisses énoncées de référence, de calibration de
  l’énergie et de basse température.
  [Exposé thermodynamique](paper/observers_are_all_you_need.pdf).
- **La capacité d’information.** Des résultats exacts comptent les registres
  qui restent distinguables après les opérations permises, bornent ce nombre
  par le support disponible et le suivent lors du raffinement.
  [Capacité des registres](code/capacity_readback/README.md).
- **Le candidat de couplage par paires intriquées.**
  [Alex Wissner-Gross](https://alexwg.org/publications/ChainTheoryAlpha.pdf)
  dérive une expression spectrale finie d’inverse `137.035999189566`. La
  [reconstruction OPH indépendante](model/observer_particles/README.md) lui
  attribue cette construction et teste ses hypothèses de résultats complets et
  de normalisation. Son identification avec le couplage
  électromagnétique physique et son rattachement à l’équation du grain d’écran
  sont conditionnels.
- **La matière et les champs.** Une action fournie de scalaire chargé et de
  Maxwell, avec sa géométrie et sa quantification, donne un hamiltonien en
  interaction après réduction de jauge, et un raffinement scalaire contrôlé
  relie les lectures classiques d’un détecteur à une équation de champ non
  linéaire. Un dictionnaire métrique commun est construit ; le produit radial
  statique testé échoue à son test gravitationnel de contrainte relative et
  nécessite une réponse géométrique.
  [Continuum de matière](paper/observers_are_all_you_need.pdf).

Ce sont des résultats mathématiques conditionnels. Faire tenir leurs prémisses
pour les observations d’un même modèle commun est le problème central de la
construction. La [carte des prémisses](model/PREMISES.md) indique quelles
prémisses la spécification exécutable fournit et lesquelles elle prend comme
entrées.

## Articles

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  est la synthèse principale : l’hypothèse de clôture, l’hypothèse du modèle
  exécutable et la tour des résultats, couche par couche, avec ce qu’il manque
  à chaque couche. Citez sa version publiée, la
  [fiche PhilPapers](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer for Observer Patch Holography](paper/executable_observer_construction.pdf)**
  énonce la spécification, les composants exécutés et les jonctions
  ordonnées qui restent à construire.
- **La bibliothèque de théorèmes** contient les preuves longues, du
  [consensus](paper/reality_as_consensus_protocol.pdf) et des
  [observateurs et registres](paper/observers_are_all_you_need.pdf) jusqu’à
  l’[espace-temps et la dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).
  L’[index des articles](paper/README.md) et l’[ensemble des articles](docs/PAPER_SET.md)
  donnent chaque article et son rôle.
- **Les articles compagnons** portent des résultats ciblés. Une relation entre
  masses de leptons, la relation de Koide, tient exactement sous une prémisse
  d’équilibre déclarée. Le programme de capacité demande si la capacité
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
un clone propre. Le vérificateur du modèle contrôle les empreintes des
artefacts et rejoue les reçus déclarés dans `model/verify.py` ; la suite de
régression vérifie aussi des résultats conservés séparément et des cas adverses :

```sh
python -m pip install -r model/requirements.txt
python -m model.verify
python -m pytest -q -W error model
```

Le guide donne aussi les contrôles autonomes des particules et des champs
sources utilisés par le [workflow de l’observateur](.github/workflows/observer-model.yml).
Ces contrôles ont des portées distinctes ; le vérificateur seul ne constitue
pas la validation complète du modèle.

La [bibliothèque Lean](Lean/) se compile avec Lake, l’[index de la documentation](docs/README.md)
donne les registres et les politiques, et le [livre](book/) offre un chemin
d’accès plus long. Les manuels sur [learn.floatingpragma.io](https://learn.floatingpragma.io/)
et la simulation interactive sur
[simulation.floatingpragma.io](https://simulation.floatingpragma.io/) sont
les introductions les plus accessibles.

## Contribuer

La tâche consiste à faire qu’un seul processus d’observation et de réparation
fournisse la physique effective : relier l’horloge de l’observateur, son
interprétation spatiale et ses mesures de champ sur les mêmes histoires, puis
comparer les probabilités de résultats, les corrélations et les réponses aux
interventions avec des observations indépendantes. Une équation de champ
fournie reste une entrée même lorsque sa simulation fonctionne. Les
contributions doivent résoudre une connexion nommée ou révéler un défaut du
modèle commun ; le [guide de construction](model/CONSTRUCTION.md), le
[guide des lois de réparation](model/CANDIDATES.md) et les
[tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
décrivent ce travail. Les publications d’articles suivent la
[politique de publication](model/CONSTRUCTION.md#research-programme-and-paper-release-hold) ;
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
[Pragma Research](https://floatingpragma.io/) réunit les projets.

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
