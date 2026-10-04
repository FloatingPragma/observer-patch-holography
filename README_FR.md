# Holographie des parcelles d’observateur

**La physique à partir des registres que les observateurs peuvent partager.**

[Article de physique](https://philpapers.org/rec/MUEFOC) · [Pragma Research](https://floatingpragma.io/research/) · [Cadence](https://floatingpragma.io/cadence/) · [English](README.md)

L’Holographie des parcelles d’observateur (OPH) étudie si la physique familière
peut être reconstruite à partir d’observateurs bornés qui parviennent à un
accord. Une parcelle d’observateur possède un état local, une frontière munie
de ports, une capacité de relecture, des registres et une rétroaction qui
répare les désaccords. Dans ce modèle, les faits deviennent publics lorsque
les registres résistent à la comparaison entre parcelles qui se recouvrent.

Ce dépôt contient le travail scientifique derrière cette idée : articles,
preuves vérifiées par machine et modèles exécutables. Il relie l’accord fini
aux probabilités quantiques, à la géométrie spatiale et à la structure de
symétrie des forces connues, avec une voie conditionnelle vers la dynamique
d’Einstein. Chaque construction précise ses hypothèses. Établir une
réalisation physique commune exige des mesures qui distinguent les modèles
proposés.

## Commencer ici

| Pour explorer… | Commencer par… |
| --- | --- |
| L’argument physique principal | [*Finite Observer Consensus as a Reconstruction Principle*](https://philpapers.org/rec/MUEFOC) |
| Les preuves et les calculs reproductibles | [Bibliothèque Lean](Lean/) et [guide de reproduction](REPRODUCE.md) |
| Les travaux associés sur les machines apprenantes | [Cadence](https://floatingpragma.io/cadence/) et [démonstrations interactives](https://floatingpragma.io/demos/) |
| Les explications pour le grand public | [Blog Pragma Research](https://blog.floatingpragma.io/) |

## Une seule architecture, toute la physique

Le programme de recherche relie ces domaines de la physique par une même
architecture d’observateurs, avec des hypothèses distinctes pour chaque résultat.

- **L’accord entre observateurs.** Sous les conditions de terminaison et de
  cohérence énoncées, les réparations locales aboutissent à un registre public
  protégé, indépendamment de leur ordre. Les théorèmes finis décrivent
  l’accord, la stabilité et le raffinement.
- **Les probabilités quantiques.** Les algèbres finies d’événements portent la
  règle de Born, le conditionnement de Lüders et des résultats exacts de
  Tsirelson, dont une construction de Bell qui atteint 2√2. Ces résultats
  mathématiques concernent la branche d’algèbre d’événements déclarée.
- **La géométrie et les forces connues.** Une construction de réponse à douze
  ports fournit une lecture spatiale tridimensionnelle. La réponse réversible
  complète et le transport interne entre observateurs imposent le type de Lie
  de jauge du Modèle standard sous les hypothèses énoncées. Une représentation
  de matière fournie donne une génération à quinze états sans anomalies ; sa
  réalisation physique et le groupe de jauge global exigent une structure
  supplémentaire.
- **L’espace-temps et la gravité.** Des lois de source et de lecture déclarées
  permettent des constructions contrôlées de géométrie causale. L’implication
  vers l’équation d’Einstein suppose une réalisation physique commune avec les
  données énoncées de contrainte mécanique, d’entropie, de continuum et
  d’échelle.
- **La matière classique et quantique.** Une action fournie de champ scalaire
  chargé couplé au champ de Maxwell permet un mouvement continu non linéaire
  contrôlé dans son secteur réel et un espace d’états quantiques en interaction
  sur un maillage fixé. Sélectionner cette action à partir des histoires des
  observateurs et l’identifier à la matière physique sont des exigences
  distinctes.
- **Les constantes comme problèmes de point fixe.** La relation de Koide tient
  exactement sous une prémisse d’équilibre déclarée. Le calcul de structure
  fine certifie une racine d’une application de fermeture déclarée et possède
  un statut diagnostique. Le programme de capacité demande si la capacité
  publique attribuée à l’univers s’accorde avec celle reconstruite depuis
  l’intérieur. Relier ces constructions aux constantes mesurées exige une
  identification physique.

L’[article principal](https://philpapers.org/rec/MUEFOC) développe ces liens
et leur portée précise. L’[index des articles](paper/) donne les exposés
spécialisés, notamment sur la thermodynamique, la géométrie des champs et les
constructions de points fixes pour les constantes.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Preuves et éléments de vérification

La [bibliothèque Lean](Lean/) contient plus de 12300 théorèmes et lemmes
publics sans preuve admise. Des certificats exacts et des simulations
reproductibles accompagnent les arguments mathématiques. Le
[simulateur de physique](https://github.com/muellerberndt/oph-physics-sim)
associé fournit une dynamique exécutable des observateurs et les éléments
de vérification conservés.

La [référence des axiomes](docs/AXIOM_REFERENCE.md) énonce les trois axiomes
fondamentaux ; le [registre des prémisses](docs/PREMISE_REGISTER_V3.md) et le
[registre des résultats](claims/claim_registry.yaml) consignent les prémisses
supplémentaires et les hypothèses propres à chaque résultat, y compris les
identifications physiques et les données empiriques.
Le [registre des postdictions](docs/POSTDICTION_LEDGER.md)
consigne les comparaisons avec des valeurs mesurées et l’origine de leurs
entrées ; l’[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md)
consigne les tests dont les conditions doivent être fixées avant l’examen
des données de comparaison. Le [programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md)
donne les observations qui réfuteraient des affirmations précises.

### Reproduire le noyau fini

Après avoir installé les dépendances décrites dans [REPRODUCE.md](REPRODUCE.md),
vérifiez le graphe des affirmations et une sélection de résultats sur
l’algèbre finie, la capacité des registres et le consensus :

```bash
python3 tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

Le guide de reproduction décrit les vérifications plus larges et la portée
de chaque famille de preuves.

## Carte de reconstruction

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="Chaîne de reconstruction OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La carte de reconstruction OPH relie registres d’observateurs, géométrie spatiale à trois dimensions, ordre causal, horloges, champs et états quantiques. Chaque flèche représente un lien mathématique ; les articles précisent les hypothèses qui les réunissent dans une description physique effective.</sub></p>

## Guide du dépôt

| Chemin | Contenu |
| --- | --- |
| [`flagship/`](flagship/) | Article principal autonome de physique, source TeX et PDF |
| [`paper/`](paper/) | Articles principaux et index des publications |
| [`Lean/`](Lean/) | Développement mathématique vérifié par machine |
| [`code/`](code/) et [`evidence/`](evidence/) | Modèles exécutables, certificats et éléments de reproduction |
| [`extra/`](extra/) et [`cosmology/`](cosmology/) | Recherches mathématiques et physiques spécialisées |
| [`book/`](book/) | *Reverse Engineering Reality*, source et livre téléchargeable |
| [`docs/`](docs/) | Politiques de lecture et registres scientifiques |

## Contribuer

OPH accueille les preuves, contre-exemples, simulations, revues indépendantes
et explications lisibles. Commencez par le [guide de reproduction](REPRODUCE.md)
et le [registre de sélection](docs/SELECTION_LEDGER.md), qui énonce les
prémisses et les frontières scientifiques utiles aux contributions.

## Des parcelles d’observateur aux machines apprenantes

[Cadence](https://github.com/muellerberndt/cadence) rend l’idée de stabilisation
exécutable sous forme d’architecture d’apprentissage. Ses parcelles logicielles
bornées portent un état local et une capacité de relecture, avec une
rétroaction qui répare les erreurs de prédiction. Son cerveau par
défaut, le Système 1, apprend par l’expérience comme les cerveaux animaux, et
non par rétropropagation. Des observateurs optionnels, le Système 2, ajoutent
une rétroaction récursive au sein de la même stabilisation.

L’[article Cadence](https://philpapers.org/rec/MUECAP-2) décrit l’architecture
et les expériences. Son évaluation repose sur le comportement
d’apprentissage et les ressources mesurées.
[Pragma Research](https://floatingpragma.io/) relie ce travail à l’IA incarnée.

## Licence

Le dépôt utilise des licences séparées par type d’artefact. Tout le logiciel, y compris la bibliothèque Lean, [`code/`](code), [`tools/`](tools) et les schémas de [`schemas/`](schemas), est publié sous [Apache-2.0](code/LICENSE). Les articles, le livre, la documentation, les figures, les données, les registres générés de [`tracking/`](tracking) et les données de particules empaquetées de [`pdg_data/`](pdg_data) sont publiés sous [CC BY-NC-SA 4.0](LICENSE). Les valeurs tabulées de `pdg_data/` conservent les conditions amont du Particle Data Group. Le fichier [LICENSE](LICENSE) donne la carte par répertoire.
