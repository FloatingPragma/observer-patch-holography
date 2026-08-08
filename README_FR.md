# Holographie des parcelles d’observateur

> La réalité est le monde public stable reconstruit par des observateurs finis et auto-lecteurs qui comparent leurs recouvrements et réparent leurs désaccords.

[Read in English](README.md) · [Livre](https://oph-book.floatingpragma.io/) · [Manuels](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [OMEGA](https://omega.floatingpragma.io/)

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
des interfaces explicites d’ordre et d’horloge qui ne constituent pas
un temps physique local d’observateur. OPH dérive aussi la cinématique de Lorentz sur la branche de support global
déclarée, le type de Lie du Modèle standard, ainsi qu’une paire extérieure
conditionnelle pour une génération de matière.

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

- **Le livre.** [*Reverse Engineering Reality*](https://oph-book.floatingpragma.io/),
  aussi disponible en [PDF qualité impression](https://cfxrbtseaimxxqsxlrku.supabase.co/storage/v1/object/public/books/reverse-engineering-reality.pdf),
  raconte toute l’histoire : ce que dit la théorie, comment elle a été
  découverte, et pourquoi le tournant observateur-d’abord est celui que la
  physique contourne depuis un siècle. Il est écrit pour divertir et la
  science y reste exacte.
- **Les manuels.** Les [manuels OPH](https://learn.floatingpragma.io/)
  enseignent la théorie par le chemin long. Chaque dérivation de base y est
  développée en entier, avec les mathématiques nécessaires construites au fur
  et à mesure. Des volumes couvrent la gravité, le Modèle Standard et
  l’unification, chacun lisible en ligne ou en PDF.
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

## Huit reçus de physique reproductibles

Ces résultats publics renvoient directement à leurs articles, preuves,
données et certificats :

1. **L’espace tridimensionnel émerge de l’algèbre des registres de
   réparation.** La réponse déclarée sur douze ports contient une complétion
   euclidienne abstraite exacte en trois dimensions. Additionner les registres
   de comparaison et compléter leur distance donne l’espace continu ordinaire,
   sans supposer de grille. L’identification physique de ces points et de leur
   échelle est en cours. Voir l’
   [article sur l’espace-temps et Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf),
   la [preuve Lean de complétion](Lean/Screen/PrimitivePortFrameQuotient.lean)
   et les [données finies de signature](evidence/einstein_convergence/).
2. **L’accord des observateurs reproduit les règles quantiques.** Dans le modèle
   fini déclaré, l’accord donne les probabilités et règles de mesure quantiques,
   ainsi que des limites exactes sur les corrélations et la copie. Le lien entre
   ces mathématiques et les appareils de mesure physiques est en cours. Voir l’
   [article principal](flagship/from_observer_consensus_to_standard_physics.pdf),
   la [preuve des registres publics](Lean/EventAlgebra/PublicRecordAlgebra.lean)
   et la [frontière de la règle de Born](Lean/EventAlgebra/FiniteEffectClosureBoundary.lean).
3. **La dynamique, les poids quantiques et l’action sont forcés, pas
   choisis.** L’évolution temporelle continue n’a aucune liberté au-delà
   d’un hamiltonien : la forme de Schrödinger est un théorème. Le poids
   quantique est la seule affectation additive sur les effets de mesure,
   jusqu’au cas du qubit. La dynamique réalisée fixe sa propre action à une
   jauge près : moindre action et histoire la plus probable sont deux
   lectures d’une même fonctionnelle, reliées au flot hamiltonien par un
   pont de Legendre. L’attachement des unités est en cours. Voir
   l’[article sur les observateurs](paper/observers_are_all_you_need.pdf),
   la [preuve de Born](Lean/EventAlgebra/FiniteBuschGleason.lean) et la
   [preuve de l’action dérivée](Lean/InformationProjection/LogTransitionAction.lean).
4. **Les quatre lois de la thermodynamique découlent du modèle fini des
   observateurs.** Avec un état de référence commun et l’information réparée
   visible par tous, la même règle finie donne l’équilibre, la croissance de
   l’entropie, le bilan entre chaleur et travail et la limite de basse
   température. Le travail restant relie les unités d’énergie et de temps à un
   système physique produit par la source. Voir l’
   [article sur les observateurs](paper/observers_are_all_you_need.pdf), la
   [preuve Lean de réparation conditionnelle](Lean/Thermodynamics/FiniteConditionalRepair.lean)
   et le [certificat exact](code/thermodynamics/conditional_repair_certificate.py).
5. **La structure de jauge du Modèle standard issue de douze ports.** La
   géométrie, la réponse réversible complète et l’accord des observateurs
   retrouvent les symétries des forces forte, faible et électromagnétique. Une
   structure de matière déclarée
   séparément donne leur forme globale familière
   $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. Les symétries supplémentaires
   responsables de la désintégration du proton dans l’unification minimale
   n’apparaissent pas. Dériver la matière et les champs physiques depuis la
   source est en cours. Voir l’
   [article sur la jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf),
   la [preuve Lean de jauge](Lean/Screen/A5OPH.lean) et la
   [preuve de la forme globale](Lean/Screen/Z6Descent.lean).
6. **Une génération de matière issue d’une recherche finie.** Dans l’ensemble
   des possibilités déclarées, le balayage exhaustif laisse les quinze états et
   les charges exactes d’une génération du Modèle standard, avec annulation de
   toutes les anomalies. Une sélection finie distincte donne un candidat de
   rang trois pour les familles. Le travail restant relie ces structures aux
   particules physiques et exclut les secteurs légers supplémentaires. Voir l’
   [article sur les particules](paper/deriving_the_particle_zoo_from_observer_consistency.pdf),
   la [preuve Lean de sélection](Lean/Screen/ExteriorSelection.lean) et la
   [preuve de la bande familiale](Lean/Screen/A5FamilyBand.lean).
7. **La relation de Koide devient un théorème.** Une réponse hermitienne $C_3$
   donne la relation exacte, dans la chambre positive, entre les masses de
   l’électron, du muon et du tau. Avec deux masses fournies, la formule fixe un
   intervalle large de 72 eV, centré à $1776{,}969027$ MeV et compatible avec
   la masse mesurée du tau. Il s’agit d’une postdiction conditionnelle informée
   par la cible, car sa prémisse d’équilibre vient de la structure connue des
   leptons. Sa dérivation depuis la source est en cours. Voir l’
   [article sur Koide](extra/koide_identity_from_positive_c3_face_circulants.pdf),
   la [preuve Lean](Lean/ObserverPatchHolography/KoideCirculant.lean) et le
   [registre des comparaisons](docs/POSTDICTION_LEDGER.md).
8. **Une empreinte gelée dans la façon dont les ondes se propagent.** Deux
   règles ondulatoires déclarées sur douze ports fixent des motifs directionnels
   distinctifs dont la première anisotropie apparaît au sixième ordre. Leurs
   rapports et règles de rejet sont sous garde cryptographique préalable à la
   comparaison. Le lien entre une règle sélectionnée par la source et un champ
   physique est en cours. Une mesure assez sensible peut alors éliminer cette
   branche. Voir l’
   [article sur la microphysique de l’écran](paper/screen_microphysics_and_observer_synchronization.pdf),
   les [reçus exacts](code/a5_fingerprint/runtime/) et l’
   [échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md).

Au-delà des huit reçus, la couche exacte porte une série de résultats
supplémentaires, chacun énoncé avec sa frontière au bout du lien :

- Un théorème de graphe signé prouve que l’écran n’a aucune excitation
  libre à coût nul : sur une capture de source sans donnée cible,
  l’opérateur signé déclaré vérifie $\lambda_{\min}\geq24^{-8661}>0$. Il
  ne fournit ni horloge physique, ni masse de particule. Voir l’
  [article sur la microphysique de l’écran](paper/screen_microphysics_and_observer_synchronization.pdf)
  et le
  [reçu figé de l’écart](https://github.com/muellerberndt/oph-physics-sim/blob/d99ca548a4853e83f819a3a2c9d813f7a3429bdb/data/local_domain/source_gap_receipt.json).
- Un théorème fini de capacité maximise l’entropie généralisée à $\log M$,
  donne le décalage exact du choc $\log(1-f)$ et fixe la relation de Sitter
  pure $\mu^2=d-2$ ; la lecture physique en avance temporelle attend ses
  dictionnaires déclarés. Voir l’
  [article ciblé sur de Sitter](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf)
  et sa [preuve Lean](Lean/ObserverPatchHolography/DeSitterCapacityShock.lean).
- Sur les branches déclarées séparément de Maxwell, de Yang--Mills et
  d’Einstein, les noyaux quadratiques ont un paramètre de masse dure nul et
  les modes classiques transverses attendus : des énoncés de porteurs
  classiques, pas des prédictions de pôles quantiques. Voir le
  [registre de la structure forcée](docs/POSTDICTION_LEDGER.md#forced-structure).
- La couche finie d’achèvement prouve un unique point terminal public sur
  les parcours achevés, un modèle de Lorentz interne avec contrat borné de
  soudure, et une interface de réseau régional munie de preuves, avec
  témoins exacts et théorèmes d’obstruction ; les attachements à la source
  restent ouverts et typés. Voir les
  [preuves du point terminal public](Lean/Tower/FixedPointEndpoint.lean), la
  [pile de théorèmes géométriques](Lean/Geometry.lean), l’
  [interface finie du réseau régional](Lean/QFT.lean) et les
  [notes de frontière Lean](Lean/docs/) pour les portées exactes.
- La bibliothèque définit ce qui compte comme un observateur à travers sept
  tests : fenêtre bornée, relecture, registres stables, action, prédiction
  vers l’avant et survie au raffinement. Un exemple concret passe les sept
  tests et trois systèmes volontairement défectueux échouent chacun à
  exactement un test. Voir la
  [preuve de l’observateur opérationnel](Lean/Tower/OperationalObserver.lean).

La bibliothèque Lean associée contient plus de 4200 théorèmes et lemmes et
aucune preuve admise. Des rapports d’axiomes explicites couvrent le
sous-ensemble audité. Vingt-trois preuves finies utilisent `native_decide` ;
leurs axiomes d’évaluation en code natif étendent la base de confiance au-delà
du noyau de Lean. Voir [Lean/](Lean/).

Le reste de ce README est l’architecture d’où viennent ces reçus.

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
concentre sur les reçus exacts et mesurés les plus forts.

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
$3{,}3129271\times10^{122}$ du modèle $\Lambda\mathrm{CDM}$ de base de Planck. La
théorie ne tranche pas entre les deux corrections, et les deux nombres
ont été calculés après coup : aucun n’est une prédiction. Le
[registre des affirmations](tracking/claims_scoreboard.md) indique les
hypothèses et les données absentes de chaque étape.

## Ce qui en sort

La relecture et la réparation finies transforment les états privés en
registres publics stables, et l’algèbre de ces registres donne les
probabilités quantiques et l’observation répétable. Sur la branche
géométrique certifiée, la géométrie conforme du support $S^2$ donne le
groupe de Lorentz connexe et exactement trois dimensions spatiales de
référentiels, et le flot modulaire avec la stationnarité de l’entropie donne
la relation de première variation d’Einstein.

La branche d’Einstein est instrumentée de bout en bout : chaque clause de
son antécédent est un théorème prouvé ou une quantité mesurée, munie d’un
instrument à fermeture sur échec et de contrôles adverses, jamais une
hypothèse. La mesure directe fournit le résultat empirique le plus fort du
corpus. Sur des configurations à 16 384, 65 536 et 262 144 porteurs, les
formes d’événements retenues portent la signature lorentzienne $(1,3)$ avec
des marges de cône de $-5{,}62$, $-3{,}22$ et $-1{,}41$ et une dispersion du
couplage décroissante, tandis qu’un contrôle de même taille à support étroit
retombe sur la signature $(2,2)$. Ces mesures établissent une sensibilité
reproductible à la structure du support ; elles n’établissent pas de loi de
convergence. Les données primaires sont archivées dans
[evidence/einstein_convergence](evidence/einstein_convergence/) et
reproductibles bit à bit depuis le
[dépôt de simulation](https://github.com/muellerberndt/oph-physics-sim).

La géométrie du porteur accomplit ensuite un travail exact surprenant.
L’incidence orientée seule dérive l’appariement antipodal, l’action propre
de $A_5$, le repère icosaédrique de rang trois et une loi de réponse inverse
sans donnée cible. La réponse réversible complète de l’axiome 1 et le
transport endogène de l’axiome 2 forcent alors le type de Lie local du
Modèle standard :

$$
P_{12}\cong_{A_5}\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5,
\qquad
(P_{12},[\ ,\ ]_\Theta)
\cong\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

Aucun groupe de jauge n’apparaît dans les prémisses ; l’alternative sans
centre est exclue par l’espace fixe de dimension un. Les matrices publiées
sont un témoin conditionnel exact, et les familles de crochets classifiées
portent désormais un secteur cinétique forcé : diagonal par blocs dans les
projecteurs certifiés, un rayon de couplage par facteur simple, et égal à
l’action dérivée de sa propre dynamique de source.

Dans l’algèbre extérieure de réponse déclarée, le balayage exhaustif des
1 024 sous-ensembles laisse exactement une sélection chirale sans anomalies :
quinze états portant les charges d’une génération du Modèle standard et un
noyau commun $\mathbb Z_6$, donc l’image fidèle maximale est
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. Sous deux prémisses
supplémentaires nommées, l’ordre exact des coûts $5-\sqrt5<6<5+\sqrt5$
sélectionne la bande familiale de rang trois.

Les frontières restent explicites d’un bout à l’autre : le courant matriciel,
l’action de matière et la sélection de la forme globale ne sont pas encore
reconstruits depuis les histoires de la source, et l’identification aux
courants de laboratoire, l’attachement des familles, la multiplicité scalaire
et les paquets de clôture restent ouverts. Le
[tableau des affirmations](tracking/claims_scoreboard.md) porte chaque
prémisse ; les reçus exacts de recherche de crochets et de Jacobi se trouvent
dans [code/a5_closure/](code/a5_closure/).

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

## Une reconstruction avec tronc commun et branches

```text
fédération de porteurs sélectionnée par la source
        ↓
parcelles avec registres, comparaison et réparation
        ↓
formes normales quotientées publiques
        ├─ reçus fédération-support → géométrie des calottes S2 et flot géométrique
        ├─ tour indépendante d’algèbres et d’états → flot modulaire
        │       composition sur la même tour → Lorentz et branche d’Einstein conditionnelle
        ├─ secteurs transportables → route indépendante de Tannaka
        └─ porteur local à douze ports → théorème exact de réponse inverse
                → type de Lie forcé par la réponse complète et le transport endogène
                → courant matriciel et matière conditionnels ; noyau tensoriel Z6 exact
                → catégorie complète de caractères et quotient global physique ouverts
                attachements du courant de laboratoire, du scalaire, du spectre et des familles ouverts
        ↓
tests quantitatifs de clôture et de lecture physique
```

Les hypothèses détaillées et les types de reçus sont énoncés dans les articles. La page d’accueil du dépôt est volontairement une carte du résultat positif, et non un substitut à ces énoncés de théorèmes.

## État technique

Les huit reçus ci-dessus donnent le résumé destiné au lecteur. Les prémisses,
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

Le bonus est quantitatif : si la boucle se ferme, $P$ et $N$ ne peuvent pas
être arbitraires. Ils doivent satisfaire des clôtures autoréférentielles. Une
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
ainsi. C’est le twist que le livre garde pour la fin de l’histoire, où
est sa place. Aucun des résultats ci-dessus n’en dépend.

## Obligations de preuve ouvertes et frontière de falsification

Le récit principal ci-dessus expose le résultat positif. Les obligations
ouvertes sont regroupées ici. L’objet fini central contient l’ensemble des
états accessibles, la politique qui décide quels registres sont publics, les
noyaux globaux de tous les points de contrôle, les projections du porteur et
les applications de raffinement.

Sur la branche réversible, ce paquet réduit la capacité à
$M_0=|X_{\rm reach}|$. Le problème mathématique ouvert consiste à
dériver la loi exacte du défaut fini $s(D)$ et à prouver qu’elle possède un
unique zéro physique. Viennent ensuite la saturation horizon-registre, le
porteur commun écran/faible-Higgs, les portes physiques qui forcent exactement
le Modèle standard, la tour commune de gravité, et les lectures quantitatives
des particules.

Le dépôt contient un
[paquet réversible de référence à douze ports](code/capacity_readback/reversible_public_checkpoint_packet.py)
et son [certificat lisible par machine](code/capacity_readback/runtime/reversible_public_checkpoint_packet_certificate.json).
Il ferme le schéma fini de la branche rapide sans se substituer au paquet
physique dérivé de la source ni à la loi d'unicité du défaut.

Le noyau structurel soutient aussi plusieurs prolongements actifs : géométrie
des neutrinos, cosmologie de capacité, spectre d’écran, gravité sombre,
transfert de Yang–Mills et systèmes matériels ou logiciels auto-lecteurs. Le
paquet radial démontre la non-identifiabilité à partir d’une seule coquille et
donne deux voies d’unicité distinctes : la dilatation physique de la source et
la tomographie par covariances radiales croisées.

Le dépôt contient aussi une calibration finie à 244 types du certificat de
colliers pour l’écart de Yang--Mills. Elle vérifie le modèle de données exact,
mais ne remplace pas un reçu physique compact-jauge dérivé de la source.

Tous partagent la même règle de conception : tout système physique proposé doit être représenté comme une parcelle bornée avec état local, frontières, relecture, registres, réparation et dossier public de preuves.

Le [programme de vérification OPH](docs/OPH_FALSIFICATION_PROGRAM.md) est volontairement limité aux affirmations mathématiques et aux branches réalisées suffisamment mûres. Il sert d’index de vérification, pas de récit principal du dépôt.

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

- [Le livre, édition web](https://oph-book.floatingpragma.io)
- [Le livre, PDF d’impression](https://cfxrbtseaimxxqsxlrku.supabase.co/storage/v1/object/public/books/reverse-engineering-reality.pdf)
- [Manuels](https://learn.floatingpragma.io)
- [Simulation interactive](https://simulation.floatingpragma.io)
- [Applications et matériel OMEGA](https://omega.floatingpragma.io)
- [Blog](https://blog.floatingpragma.io/)
- OPH Sage sur [Telegram](https://t.me/HoloObserverBot) et [X](https://x.com/OphSage)

## Licence

Le dépôt utilise des licences séparées par type d’artefact. Tout le logiciel, y compris la bibliothèque Lean, [`code/`](code) et [`tools/`](tools), est publié sous [Apache-2.0](code/LICENSE). Les articles, le livre, la documentation, les figures et les données sont publiés sous [CC BY-NC-SA 4.0](LICENSE). Les fichiers de conception matérielle utilisent la CERN-OHL-W 2.0. Le fichier [LICENSE](LICENSE) donne la carte par répertoire.
