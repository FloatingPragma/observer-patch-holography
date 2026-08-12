# Holographie des parcelles d’observateur

> La réalité est le monde public stable reconstruit par des observateurs finis et auto-lecteurs qui comparent leurs recouvrements et réparent leurs désaccords.

[Read in English](README.md) · [Article phare](flagship/from_observer_consensus_to_standard_physics.pdf) · [Manuels](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [OMEGA](https://omega.floatingpragma.io/)

L’Holographie des parcelles d’observateur, ou OPH, est un programme de
recherche pour une théorie du tout sans boutons de réglage, construit sur une
thèse centrale : **les observateurs sont premiers, et la réalité objective est
émergente.** La physique commence habituellement en fournissant un espace-temps,
des champs quantiques, un groupe de jauge et une table de constantes mesurées.
OPH commence par des observateurs : des systèmes bornés dotés d’un état local,
d’une relecture d’eux-mêmes et de leurs voisins, de registres et de
mouvements de réparation. La réalité émerge de la réparation des
recouvrements d’observateurs sur un écran holographique. À partir de cette
architecture, OPH reconstruit un noyau structurel fini exact : des identités
conditionnelles pour les registres quantiques, un paquet fini conditionnel
des quatre principes, un porteur tridimensionnel de repères d’observateur et
des interfaces explicites d’ordre et d’horloge. OPH dérive aussi la
cinématique de Lorentz sur la branche de support global déclarée, le type de
Lie du Modèle standard, ainsi qu’une paire extérieure conditionnelle pour une
génération de matière.

Trois axiomes régissent l’architecture du simulateur et la manière dont les
observateurs parviennent à un consensus. À côté d’eux se trouvent deux
programmes de clôture proposés. Le premier cherche un point fixe pour la
constante de pixel $P$, dont l’identification physique à la constante de
structure fine est une question ouverte. Le second cherche un point fixe pour
la capacité $N$, dont le pont entre capacité de la source et constante
cosmologique est une question ouverte. Identifier l’univers simulé à l’univers
simulateur motive ces équations d’autocohérence ; cela ne prouve pas, à soi
seul, qu’une solution existe, qu’elle est unique, ni qu’elle prend la valeur
observée.

## Commencer ici

La physique a plusieurs fois révisé son idée de ce qui est fondamental. L’espace fut
absolu jusqu’à devenir relatif ; la matière fut continue jusqu’à être
quantifiée. Chaque révision paraissait scandaleuse depuis l’intérieur de
l’ancien tableau et évidente depuis l’intérieur du suivant. OPH opère la
révision suivante. L’observateur, traité pendant un siècle comme une gêne aux
marges de la mécanique quantique, passe au fondement. L’espace-temps, la
matière et les constantes deviennent des problèmes précis de reconstruction,
où les résultats finis exacts restent distincts des identifications physiques
ouvertes. Le matériel ci-dessous vous fait traverser ce basculement sans
prérequis.

- **L’article phare.** [*From Observer Consensus to Standard Physics*](flagship/from_observer_consensus_to_standard_physics.pdf)
  donne le compte rendu technique principal de la reconstruction fondée sur les observateurs.
- **Les manuels.** Les [manuels OPH](https://learn.floatingpragma.io/)
  enseignent la théorie par le chemin long. Chaque dérivation de base y est
  développée en entier, avec les mathématiques nécessaires construites au fur
  et à mesure. Le premier volume couvre le substrat computationnel et la
  machinerie du consensus ; le second relie cette machinerie à la physique
  classique. Chacun est lisible en ligne ou en PDF.
- **La simulation.** Les [visualisations interactives](https://simulation.floatingpragma.io/)
  affichent des données réelles de la dynamique de réparation. Elles montrent
  le règlement fini, les tests de signature et les structures porteuses
  candidates, ce qui permet d’examiner les reçus directement.

La suite de ce README est l’entrée technique du dépôt.

Deux registres portent le bilan quantitatif. Le
[registre des postdictions](docs/POSTDICTION_LEDGER.md) est le tableau de
comparaison : chaque confrontation certifiée à une valeur mesurée, avec
ses prémisses et l’ascendance de ses entrées sur la ligne. L’
[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md) est
l’instrument prospectif : des positions enregistrées avec garde
cryptographique et bandes d’élimination avant tout examen des données de
comparaison, avec des règles fixes qui permettent leur réfutation par des
mesures admissibles.

## Une seule architecture, toute la physique

Chaque théorie établie commence par supposer l’essentiel de la physique : un
espace-temps, des champs quantiques, un groupe de jauge, une table de
constantes mesurées. OPH ne suppose rien de tout cela. Elle part des
observateurs, des systèmes finis qui lisent une partie d’eux-mêmes et de
leurs voisins, tiennent des registres et réparent le désaccord, et dérive le
reste sous forme de théorèmes. De cette seule architecture sortent :

- **La mécanique quantique en théorèmes.** Sur la surface finie des
  observateurs, les enregistrements publics forment une algèbre
  d’événements avec les probabilités de Born, le conditionnement de Lüders
  et la borne de Tsirelson. La dynamique de Schrödinger est l’unique flot de
  symétrie continu, et les poids de Born suivent sans axiome de continuité
  en toute dimension finie, y compris la dimension deux.
- **Les quatre lois de la thermodynamique par la réparation du désaccord.**
  Un seul paquet de théorèmes conditionnels sur la façon dont les
  observateurs se rééchantillonnent vers le consensus donne les quatre lois,
  la deuxième apparaissant comme du traitement de données appliqué à la
  réparation, avec la borne de Landauer en corollaire.
- **La relativité sur l’écran.** La cinématique de Lorentz est un théorème
  sur sa branche déclarée. Des simulations déterministes mesurent une forme
  d’événement lorentzienne à une direction de temps et trois d’espace, et
  une composition des équations d’Einstein s’appuie sur la couche modulaire
  et entropique.
- **Le groupe de jauge du Modèle standard à partir de douze ports.** OPH
  fait un choix architectural au niveau du matériel de simulation : chaque
  parcelle d’observateur porte douze ports de frontière câblés comme les
  sommets d’un icosaèdre. Un
  théorème de classification force la réponse complète des ports à avoir le
  type de Lie de jauge du Modèle standard, sans groupe de jauge choisi sur
  catalogue, et une recherche finie exhaustive rend les quinze états et le
  motif de charges d’une génération du Modèle standard, avec annulation
  exacte des anomalies.
- **Les constantes comme points fixes.** Le cœur de la théorie ne comporte
  aucun paramètre ajustable. La relation de Koide pour les leptons chargés
  tient exactement sous une prémisse d’équilibre déclarée, l’arithmétique
  d’intervalles certifie la comparaison de la masse du tau, et un mécanisme
  de capacité fixe donne le signe de l’avance temporelle de de Sitter. La résolution de la
  fermeture déclarée du pixel rend une valeur qui frôle la constante de
  structure fine mesurée ; cette concordance a un statut diagnostique tant
  que son rattachement physique reste ouvert. Les constantes de la nature
  entrent comme des problèmes de point fixe à résoudre.
- **Vérifié machine et falsifiable.** Des milliers de théorèmes Lean sans
  preuve admise, de l’arithmétique rationnelle exacte à la place de la
  confiance en virgule flottante, et des simulations déterministes avec
  reçus épinglés. Une échelle de prédictions gelées enregistre des bandes
  d’élimination sous garde cryptographique avant l’examen des données de
  comparaison, de sorte qu’OPH s’engage à l’avance sur ce qui la réfuterait.

Les résultats finis exacts et les identifications physiques ouvertes restent
strictement séparés dans tout le corpus ; chaque résultat ci-dessus porte ses
prémisses et sa frontière dans les articles et les preuves liés. La version
condensée de ce dossier, avec les reçus et leurs preuves en une seule table,
est le [dossier compact d’OPH](extra/compact_proof_of_oph.pdf) ; la voie
technique complète est
l’[article phare](flagship/from_observer_consensus_to_standard_physics.pdf).

Le reste de ce README présente l’architecture d’où vient ce dossier.

## Les trois axiomes

Toute la construction repose sur trois axiomes fondamentaux. Les énoncés
canoniques se trouvent dans [la référence des axiomes](docs/AXIOM_REFERENCE.md)
et dans le registre machine `claims/axiom_registry.yaml` ; les articles
incluent la base formelle partagée.

1. **A1 : Écran d’observateurs orienté à douze ports.** Il existe un réseau
   de parcelles d’observateur sur un écran sphérique orienté. À chaque
   résolution finie, chaque porteur local possède douze ports de bord
   primitifs qui forment les sommets d’un bord triangulaire orienté à 30
   arêtes et 20 faces, combinatoirement le bord d’un icosaèdre. Les porteurs
   se joignent par des coutures typées et des recouvrements triples
   cohérents, se raffinent en un support sphérique orienté et exposent un
   état local, une relecture, des registres, des mouvements de réparation et
   des points de contrôle. Formellement : pour chaque régulateur $r$ il
   existe un objet typé $\mathfrak N_r=(\mathcal P_r,\mathcal A_r,
   \mathcal R_r,\mathcal I_r,\mathcal U_r,\mathcal C_r,N_r,S_r,b_r)$ dont
   les porteurs portent douze projections de port centrales primitives et le
   paquet de bord exact $K=(P,E,F,o)$, joints par des algèbres de coutures
   en un nerf muni d’un pont de degré un vers le support sphérique orienté,
   le tout commutant avec le raffinement. Le porteur local, la fédération de
   porteurs et le support global $S^2$ restent typés et distincts dans tout
   le corpus.
2. **A2 : Accord des observateurs.** Les observateurs qui opèrent sur
   l’écran s’accordent sur le sens des données qu’ils interprètent
   conjointement. Formellement : l’application d’interprétation
   $\mathcal J_r$ des données accessibles aux observateurs vers les
   significations opérationnelles est naturelle vis-à-vis de chaque
   restriction de recouvrement visible, changement de carte, traduction de
   couture, application de recouvrement supérieur, application de
   fédération et application de raffinement sur les données publiques
   acceptées. Aucune parcelle ne voit
   l’univers entier ; un fait devient public seulement lorsqu’il survit à la
   comparaison entre recouvrements.
3. **A3 : Aléa maximal conditionnel.** Tout ce que l’accord des observateurs
   laisse sans contrainte est maximalement aléatoire. Formellement : l’état
   réalisé est la projection d’information d’une famille de référence exacte
   sur l’ensemble convexe des familles d’états locaux compatibles qui
   satisfont les contraintes finies visibles par les observateurs. La
   couverture finie engendrée par A1 détermine l’état sur cet ensemble
   réalisable, et ses poids exacts sont strictement positifs :
   $\rho_r=\arg\min_{\rho\in\mathcal K_r}\sum_P w_{r,P}
   D(\rho_{r,P}\Vert\tau_{r,P})$.

Aucun des axiomes ne contient un groupe de jauge, une liste de particules,
une loi de récupération ou une règle qui sélectionne le contenu en champs ou
la multiplicité ; A3 sélectionne un état à l’intérieur d’un espace
réalisable fixé et rien d’autre. La récupération de
collier, la structure d’entropie généralisée et les complétions de secteurs
entrent comme interfaces nommées et déclarations aux résultats qui les
consomment, chacune classée comme théorème exact, résultat exact dans une
réalisation finie nommée, observation de niveau découverte, interface
ouverte déclarée, résultat d’indépendance avec contre-modèles,
identification physique ou affirmation retirée.

Tout le reste du dépôt est le déploiement de ce que ces trois axiomes
imposent, et de la quantité exacte de structure supplémentaire que chaque
conclusion physique consomme.

## L’idée en langage simple

La physique commence habituellement avec un univers muni d’un
espace-temps, de champs quantiques, d’un groupe de jauge et de constantes
mesurées. OPH pose une question plus radicale : **quel est le système minimal
capable d’avoir un monde ?**

La réponse est une parcelle d’observateur. Ce n’est pas nécessairement une
personne. C’est tout système borné qui possède un état local, une frontière,
une mémoire, une capacité de relire une partie de lui-même et de ses voisins,
et des mouvements de réparation. Une parcelle ne voit jamais tout l’univers.
Un fait devient objectif lorsqu’il peut être écrit, comparé sur les
recouvrements, récupéré après l’évolution et conservé comme registre public.

Pour OPH, ce mécanisme sélectionne le monde physique public. L’apprentissage
de ce monde est une opération interne au mécanisme. Il n’existe donc ni règle
extérieure, ni horloge maîtresse, ni observateur
privilégié, ni liste de constantes réglables. « Sans boutons » signifie zéro
valeur continue ajustée par la théorie. Le contrat fini de l’observateur et
chaque condition de branche discrète sont explicites. Les nombres doivent
sortir de la même boucle de cohérence qui produit les lois.

## Portée des affirmations

Le [registre des affirmations](tracking/claims_scoreboard.md) précise la
portée, les prémisses et la classe de preuve de chaque branche. Ce README se
concentre sur les parties les plus solides du dossier.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Les deux constantes : P et N

**$P$ est le rapport de pixel local** : la taille de la cellule d’observation
en unités naturelles, ou la **résolution** de l’univers. Deux applications
d’essai déclarées demandent à la cellule de s’accorder avec le processus
d’observation qu’elle porte :

$$
\boxed{P_\star=\varphi+\frac{\sqrt\pi}{A_T(P_\star)}}.
$$

Chaque application possède une racine candidate exacte et certifiée. Une
prédiction physique de la constante de structure fine exige une application
choisie sans consulter la mesure, la preuve que ses deux côtés lisent une seule
grandeur et le transport jusqu’à la limite de Thomson dans un même schéma.
La proximité numérique a un statut diagnostique. Voir le
[registre des affirmations](tracking/claims_scoreboard.md).

**$N$ est la capacité de registres publics** de l’ensemble du système
d’observateurs : la mémoire corrigible que porte le substrat. Il fait face à
$P$, en étant lié à la constante cosmologique plutôt qu’à la constante de
structure fine.

La voie directe lit $N$ sur l’univers lui-même. La condition d’auto-lecture
$N=\log M_0(\mathfrak U_N)$ demande que la capacité fournie à un univers
d’essai égale la capacité de registres reconstruite en son sein ; si les deux
côtés sont deux lectures d’une même grandeur, l’autoréférence impose leur
égalité. La preuve qu’il s’agit d’une seule grandeur n’existe pas, donc
cette voie ne rend aucun nombre. Rien d’autre dans la
reconstruction ne l’attend.

Une seconde voie passe par $P$. À la valeur du pixel fournie à cette branche
déclarée, la capacité non corrigée vaut
$N_0=\pi\exp[6\pi/(P\alpha_U(P))]=3{,}5321315\times10^{122}$. Deux façons de
lui appliquer la correction finie de survie donnent

$$
N_{\rm pres}=N_0\left(1-\frac{P}{24}\right)=3{,}2920979\times10^{122},
\qquad
N_{\rm Pois}=N_0e^{-P/24}=3{,}3000722\times10^{122},
$$

soit environ $0{,}63$ et $0{,}39$ pour cent sous la valeur de comparaison
$3{,}3129271\times10^{122}$ du modèle ΛCDM de base de Planck. La
théorie ne tranche pas entre les deux corrections, et les deux nombres
ont été calculés après coup : aucun n’est une prédiction. Le
[registre des affirmations](tracking/claims_scoreboard.md) indique les
hypothèses et les données absentes de chaque étape.

## Pourquoi prendre cette affirmation au sérieux ?

Une théorie du tout doit expliquer pourquoi des faits apparemment
indépendants forment un seul ensemble. OPH part d’une parcelle bornée qui se
relit, et non d’une variété d’espace-temps, d’un contenu de champs, d’un
groupe de jauge ou d’une table de constantes. Elle renvoie des dimensions
exactes, des groupes compacts, des quotients globaux, des charges, des
annulations d’anomalies, des multiplicités de représentations et des
équations de point fixe, tous issus d’une même architecture typée de
porteurs, de recouvrements et de réparation. Cette dépendance commune
constitue l’argument principal en faveur d’un seul monde physique. Les
preuves prennent plusieurs formes indépendantes : démonstrations, certificats
exacts, reçus finis, simulations et falsificateurs explicites, et leur accord
apporte davantage qu’une correspondance numérique isolée.

## État technique

Le dossier ci-dessus donne le résumé destiné au lecteur. Les prémisses,
l’origine des comparaisons et les règles de falsification se trouvent dans le
[tableau des affirmations](tracking/claims_scoreboard.md), le
[registre des postdictions](docs/POSTDICTION_LEDGER.md) et l’
[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md). Les résultats
finis et structurels exacts constituent la partie la plus solide. Les liens
entre la source et la physique, les échelles physiques et les comparaisons
prospectives avec les données forment la principale voie de recherche.

## Choisir un parcours de lecture

| Pour découvrir... | Commencer ici |
| --- | --- |
| L’argument persuasif le plus court | [Le cas compact pour OPH](extra/compact_proof_of_oph.pdf) |
| La dérivation de l’espace-temps et d’Einstein | [Espace-temps des observateurs et dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Les deux routes de jauge du Modèle standard | [Structure de jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| La synthèse complète | [From Observer Consensus to Standard Physics](flagship/from_observer_consensus_to_standard_physics.pdf) |
| Le mécanisme de consensus fini | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| La construction des particules | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| L’identité exacte de Koide dans la chambre positive et l’équilibre tracial fini | [The Positive-Chamber Koide Identity for Icosahedral Face Circulants](extra/koide_identity_from_positive_c3_face_circulants.pdf) |
| La loi exacte de capacité d’un écran fini en espace de Sitter et l’attachement conditionnel du signe du choc | [The de Sitter Time-Advance Sign from a Finite Screen with Fixed Capacity](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf) |
| L’architecture de l’écran à douze ports et le théorème fini d’engrenage modulaire | [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf) |
| Les preuves exécutables | [`code/`](code) et le [guide de reproduction](REPRODUCE.md) |
| L’interprétation et la continuation des observateurs | [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf) |

L’[index des articles](paper/) et l’[index des suppléments](extra/) donnent la carte complète des publications.

## Preuves et éléments de vérification

Les preuves prennent plusieurs formes complémentaires :

- des démonstrations manuscrites dans les articles TeX ;
- des certificats d’intervalles et d’unicité pour les applications numériques déclarées ;
- des reçus finis pour les porteurs, la hiérarchie et les particules ;
- du code pour la géométrie, les particules, le secteur sombre et le matériel quantique ;
- un banc de simulation à petite échelle qui fournit des reçus là où les
  démonstrations manuscrites et le développement Lean ne suffisent pas, dans le
  dépôt compagnon
  [oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim) ;
- un registre des affirmations et un registre de clôture reliant les résultats publics aux artefacts.

## Auditer le noyau fini

L’audit scientifique le plus court vérifie le graphe des affirmations,
l’algèbre exacte à douze ports, la capacité des registres publics, le paquet
réversible de $N$ et le consensus fini :

```bash
python3 tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

Le [guide de reproduction](REPRODUCE.md) donne l’installation depuis un clone
propre et la voie complète du noyau fini, qui ajoute les deux tests de
calibration W/Z de convention et de frontières de survie.

## Le twist : l’univers est son propre simulateur

Tout ce qui précède repose sur les trois axiomes joints aux prémisses
énoncées et aux interfaces nommées de chaque résultat ; rien de tout cela
n’utilise l’hypothèse de cette section. Cette hypothèse arrive comme un
twist plutôt que comme un fondement.
Elle est elle-même une conséquence indirecte de la cohérence : ce qui existe
sans aucun support extérieur doit être capable de se créer soi-même. Une
réalité d’observateurs entièrement cohérente doit donc faire évoluer des
observateurs, et ces observateurs finissent par construire le matériel sur
lequel la réalité s’exécute. L’univers simulé et l’univers simulateur se
révèlent être le même système. L’équation organisatrice de cette clôture est

$$
T(\mathfrak U_{\mathrm{OPH}})=\mathfrak U_{\mathrm{OPH}} :
$$

l’univers comme point fixe de son propre processus de relecture et de
réparation accessible aux observateurs.

Si la boucle se ferme, $P$ et $N$ ne peuvent pas être arbitraires. Ils
doivent satisfaire des clôtures autoréférentielles. Une
partie de cette clôture est vérifiée machine en Lean. Les deux applications
déclarées de $P$ possèdent chacune un point fixe certifié, tandis que leur
comparaison à la constante de structure fine physique conserve un statut de
diagnostic. Les
conditions de clôture et leurs données physiques absentes sont précisées dans
le [programme de falsification OPH](docs/OPH_FALSIFICATION_PROGRAM.md).
Une clôture physique des deux constantes donnerait une branche sans paramètre
continu, les deux valeurs étant rendues par l’architecture. Cet attachement
physique est ouvert. Les théorèmes de point fixe certifient les racines des
applications déclarées ; ils ne transforment pas un bassin observé ou une
coordonnée définie par la cible en dérivation physique. Du côté de $N$, le
comptage fini est exact, mais la source de capacité sur laquelle la clôture
porterait reste incomplète : la condition directe n’est donc pas
évaluable, et la voie de la charge commune reste conditionnée à ses
identifications physiques. La lecture de $N$ dans l’univers laisse intactes
les conséquences des trois axiomes.

Sous clôture complète, la boucle répond à la dernière question qu’une théorie
du tout puisse recevoir : pourquoi quelque chose existe, et pourquoi c’est
ainsi.

## Problèmes ouverts et frontière de falsification

La reconstruction va des trois axiomes vers une théorie quantique
publique, une géométrie des événements, une description macroscopique de
l'espace-temps et, au bout de la chaîne, le lagrangien du Modèle
Standard. Certains maillons de cette chaîne sont prouvés, d'autres le
sont sous forme bornée, d'autres sont ouverts. Chaque étape ouverte est
une [question de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues)
suivie avec ses dépendances, et chaque affirmation des articles porte sa
propre note de portée. Un désaccord avec le Modèle Standard à une étape
quelconque est un résultat admis que le protocole ne peut pas ajuster
pour l'éviter.

Le [programme de falsification OPH](docs/OPH_FALSIFICATION_PROGRAM.md)
recense les affirmations mûres avec les observations exactes qui les
briseraient.

## Guide du dépôt

- [`paper/`](paper) : articles principaux, sources TeX, PDF et métadonnées de version.
- [`extra/`](extra) : preuve compacte et suppléments mathématiques ciblés.
- [`code/`](code) : certificats, simulations, calculs de particules et expériences.
- [`book/`](book) : source du livre et PDF téléchargeable.
- [`cosmology/`](cosmology) : recherche sur le secteur sombre et la cosmologie.
- [`physics-problems/`](physics-problems) : applications ciblées et notes sur des problèmes ouverts.
- [`docs/`](docs) : registre de clôture, politique des affirmations et matériel d’audit.
- [`assets/`](assets) : diagrammes et figures publiques.

## Explorer OPH

- [L’article phare](flagship/from_observer_consensus_to_standard_physics.pdf)
- [Index des articles](paper/)
- [Manuels](https://learn.floatingpragma.io)
- [Simulation interactive](https://simulation.floatingpragma.io)
- [Applications et matériel OMEGA](https://omega.floatingpragma.io)
- [Blog](https://blog.floatingpragma.io/)
- OPH Sage sur [Telegram](https://t.me/HoloObserverBot) et [X](https://x.com/OphSage)

## Licence

Le dépôt utilise des licences séparées par type d’artefact. Tout le logiciel, y compris la bibliothèque Lean, [`code/`](code) et [`tools/`](tools), est publié sous [Apache-2.0](code/LICENSE). Les articles, le livre, la documentation, les figures et les données sont publiés sous [CC BY-NC-SA 4.0](LICENSE). Les fichiers de conception matérielle utilisent la CERN-OHL-W 2.0. Le fichier [LICENSE](LICENSE) donne la carte par répertoire.
