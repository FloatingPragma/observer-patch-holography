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

## Sept reçus de physique reproductibles

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
3. **Les quatre lois de la thermodynamique découlent du modèle fini des
   observateurs.** Avec un état de référence commun et l’information réparée
   visible par tous, la même règle finie donne l’équilibre, la croissance de
   l’entropie, le bilan entre chaleur et travail et la limite de basse
   température. Le travail restant relie les unités d’énergie et de temps à un
   système physique produit par la source. Voir l’
   [article sur les observateurs](paper/observers_are_all_you_need.pdf), la
   [preuve Lean de réparation conditionnelle](Lean/Thermodynamics/FiniteConditionalRepair.lean)
   et le [certificat exact](code/thermodynamics/conditional_repair_certificate.py).
4. **La structure de jauge du Modèle standard issue de douze ports.** La
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
5. **Une génération de matière issue d’une recherche finie.** Dans l’ensemble
   des possibilités déclarées, le balayage exhaustif laisse les quinze états et
   les charges exactes d’une génération du Modèle standard, avec annulation de
   toutes les anomalies. Une sélection finie distincte donne un candidat de
   rang trois pour les familles. Le travail restant relie ces structures aux
   particules physiques et exclut les secteurs légers supplémentaires. Voir l’
   [article sur les particules](paper/deriving_the_particle_zoo_from_observer_consistency.pdf),
   la [preuve Lean de sélection](Lean/Screen/ExteriorSelection.lean) et la
   [preuve de la bande familiale](Lean/Screen/A5FamilyBand.lean).
6. **La relation de Koide devient un théorème.** Une réponse hermitienne $C_3$
   donne la relation exacte, dans la chambre positive, entre les masses de
   l’électron, du muon et du tau. Avec deux masses fournies, la formule fixe un
   intervalle large de 72 eV, centré à $1776{,}969027$ MeV et compatible avec
   la masse mesurée du tau. Il s’agit d’une postdiction conditionnelle informée
   par la cible, car sa prémisse d’équilibre vient de la structure connue des
   leptons. Sa dérivation depuis la source est en cours. Voir l’
   [article sur Koide](extra/koide_identity_from_positive_c3_face_circulants.pdf),
   la [preuve Lean](Lean/ObserverPatchHolography/KoideCirculant.lean) et le
   [registre des comparaisons](docs/POSTDICTION_LEDGER.md).
7. **Une empreinte gelée dans la façon dont les ondes se propagent.** Deux
   règles ondulatoires déclarées sur douze ports fixent des motifs directionnels
   distinctifs dont la première anisotropie apparaît au sixième ordre. Leurs
   rapports et règles de rejet sont sous garde cryptographique préalable à la
   comparaison. Le lien entre une règle sélectionnée par la source et un champ
   physique est en cours. Une mesure assez sensible peut alors éliminer cette
   branche. Voir l’
   [article sur la microphysique de l’écran](paper/screen_microphysics_and_observer_synchronization.pdf),
   les [reçus exacts](code/a5_fingerprint/runtime/) et l’
   [échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md).

Un théorème distinct de graphe signé prouve que l’écran n’a aucune
excitation libre à coût nul : une capture issue de la source et sans donnée
cible fixe l’ordre causal, la topologie de coutures, les sections typées et
38 triangles frustrés, et l’opérateur signé déclaré vérifie
$\lambda_{\min}\geq24^{-8661}>0$; un calcul numérique distinct donne
$0{,}1175367$ comme raffinement non certifié. Ce
résultat est distinct du spectre de réparation de jauge compacte et de
l’écart de masse continu de Yang–Mills, et il ne fournit ni horloge
physique, ni masse de particule. Voir l’
[article sur la microphysique de l’écran](paper/screen_microphysics_and_observer_synchronization.pdf)
et le
[reçu figé de l’écart](https://github.com/muellerberndt/oph-physics-sim/blob/d99ca548a4853e83f819a3a2c9d813f7a3429bdb/data/local_domain/source_gap_receipt.json).

Un théorème fini distinct maximise l’entropie généralisée à $\log M$,
donne le décalage exact du choc $\log(1-f)$ et fixe la relation de Sitter pure
$\mu^2=d-2$. Sa lecture physique comme avance temporelle exige les
dictionnaires déclarés de l’horizon, de la masse de l’observateur, de la
gravitation, des modes de jauge et de l’opérateur cinétique. Voir l’
[article ciblé sur de Sitter](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf)
et sa [preuve Lean](Lean/ObserverPatchHolography/DeSitterCapacityShock.lean).

Le registre structurel consigne trois autres conséquences au niveau des
actions. Sur les branches déclarées séparément de Maxwell non brisée, de
Yang--Mills perturbatif et d’Einstein pur, les noyaux quadratiques ont un
paramètre de masse dure nul et les modes classiques transverses, ou
transverses sans trace, attendus. Il s’agit d’énoncés sur des porteurs
classiques, pas de prédictions de pôles quantiques du photon, du gluon ou du
graviton. Voir le
[registre de la structure forcée](docs/POSTDICTION_LEDGER.md#forced-structure).

La bibliothèque Lean associée contient plus de 4000 théorèmes et lemmes et
aucune preuve admise. Des rapports d’axiomes explicites couvrent le
sous-ensemble audité. Vingt-trois preuves finies utilisent `native_decide` ;
leurs axiomes d’évaluation en code natif étendent la base de confiance au-delà
du noyau de Lean. Voir [Lean/](Lean/).

La couche finie d’achèvement V2 rend le quotient public littéral à un
régulateur : deux configurations brutes sont équivalentes lorsque leurs
signatures de relecture complètes sont égales. La terminaison finie,
l’achèvement terminal, la confluence, la complétude sémantique et la
congruence au quotient impliquent alors un unique point terminal public pour
tous les parcours achevés et tous les représentants. Ces prémisses sont
conditionnelles. L’achèvement terminal n’est pas une loi d’équité du
planificateur, et aucun monde physique choisi par la source ou cohérent entre
régulateurs n’en découle. Voir les
[preuves du point terminal public](Lean/Tower/FixedPointEndpoint.lean).

La même couche donne un modèle de Lorentz interne dans les coordonnées de
Pauli réelles des matrices hermitiennes $2\times2$. Le déterminant a l’inertie
$(1,3)$, les rayons nuls futurs positifs sont équivalents comme ensembles à
$S^2$, les repères unitaires de genre temps futurs ont des espaces de repos
euclidiens de dimension trois, et une carte linéaire exacte rejoint les
coordonnées du tenseur d’Einstein avec la convention de signe opposée. Cette
géométrie possède un contrat algébrique borné de soudure : descente
unique des lectures invariantes sur les classes de coïncidence. Séparément, à
partir d’un cocycle de recouvrement affine lorentzien orienté dans le temps et
d’une famille fournie de coordonnées de carte qui satisfait sa loi de
recouvrement, il prouve les lois d’identité, d’inverse et de cocycle, ainsi que la covariance des
déplacements, intervalles, directions célestes, repères et lectures dans
l’espace de repos local. La pile formelle n’identifie pas la lecture descendue
au quotient à cette famille de coordonnées. Le quotient de Gram de rang trois
est isométrique à
une fibre de repos interne standard comme lecture candidate, jamais comme
espace absolu. La famille de coordonnées et sa loi de recouvrement, l’atlas de
coïncidence issu de la source, la population des événements, la séparation
certifiée, les cartes ouvertes, le cône physique, le raffinement, la
causalité, la sélection du repère, l’horloge et l’attachement à l’espace-temps
physique restent des reçus explicites en aval. Voir la [pile de théorèmes
géométriques](Lean/Geometry.lean).

Une interface finie de réseau régional muni de ses preuves regroupe
l’isotonie, la commutation des régions déclarées disjointes, le raffinement,
la réparation locale idempotente, le recollement conditionnel des restrictions
sur des familles déclarées de sous-régions et la non-perturbation exacte des
observables distantes. Le type formel de famille n’impose aucune couverture
conjointe. Ses restrictions contravariantes et ses lois de recollement sont des
données fortes déclarées. Sa construction de cohérence paramétrée par une
partition et un état affecte la même algèbre commutative de registres à chaque
région ; elle ne constitue pas un témoin fermé issu de la source. Un réseau non
commutatif attaché à la source, une couverture véritable, une factorisation
régionale, des canaux
quantiques positifs, la localité d’un planificateur adaptatif, la causalité de
l’espace-temps et une théorie continue avec loi de tranche temporelle sont des
obligations ouvertes. La couche vérifie aussi la localité pour un mot de réparation fixé,
l’invariance marginale bipartite générique, la conservation et le transport
finis, ainsi que des auxiliaires conditionnels pour les histoires finies et
les variations réelles. Un espace de Gibbs fini ne peut pas fournir toutes
les variations réelles en un site. La loi d’histoire issue de la source, le
théorème de transfert, l’action physique et l’horloge n’ont aucune
construction. Il s’agit de résultats structurels, pas de prédictions physiques.
Voir l’[interface finie du réseau régional](Lean/QFT.lean), le
[registre des postdictions](docs/POSTDICTION_LEDGER.md#forced-structure) et
les [notes de frontière Lean](Lean/docs/) pour les portées exactes.

La bibliothèque définit ce qui compte comme un observateur : un système
doté d'une fenêtre bornée sur le monde, capable de relire ses propres
enregistrements, de les garder stables sous sa propre dynamique, d'agir
sur eux, de prédire vers l'avant plutôt que vers l'arrière, et de survivre
au raffinement. Un exemple concret passe les sept tests et trois systèmes
volontairement défectueux échouent chacun à exactement un test. Un
théorème compagnon relie l'enregistrement d'un événement observé, sa
localisation, l'accès de son propriétaire et ses coordonnées à un seul
paquet de données sous-jacent, et ce lien survit au raffinement.

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
$3{,}3129271\times10^{122}$ du modèle $\Lambda$CDM de base de Planck. La
théorie ne tranche pas entre les deux corrections, et les deux nombres
ont été calculés après coup : aucun n’est une prédiction. Le
[registre des affirmations](tracking/claims_scoreboard.md) indique les
hypothèses et les données absentes de chaque étape.

## Un univers complet imposé par la cohérence

OPH teste une proposition centrale : **une seule architecture de cohérence
entre observateurs peut reconstruire l’architecture de notre univers observé.**
Les résultats finis exacts et les producteurs physiques ouverts sont
distingués ci-dessous.

OPH est invariant sous les changements de présentation cachée qui préservent
le quotient visible aux observateurs. Il reste sensible à l’incidence des
ports, à la topologie, aux processus d’enregistrement et aux lois de réponse
lorsque ces données sont visibles. La géométrie du support physique peut donc
sélectionner la branche réalisée.

La relecture et la réparation finies transforment les états privés en
registres publics stables, et l’algèbre de ces registres donne les
probabilités quantiques et l’observation répétable. Sur la branche géométrique
certifiée, la géométrie conforme du support $S^2$ donne le groupe de Lorentz
connexe et exactement trois dimensions spatiales de référentiels, et le flot
modulaire avec la stationnarité de l’entropie donne la relation de première
variation d’Einstein.

La branche d’Einstein est instrumentée de bout en bout. Chaque
clause de son antécédent (normalisation modulaire géométrique, cyclicité GNS
et intersections modulaires, cône d’événements lorentzien, contrainte et
couplage de même source) possède un instrument certifié par machine, à
fermeture sur échec, avec contrôles négatifs adverses et contre-modèles
sémantiques : chaque clause est un théorème prouvé ou une quantité mesurée,
jamais une hypothèse. Deux clauses sont des théorèmes :
l’universalité du couplage vaut avec dispersion nulle pour toute loi source à
symétrie icosaédrique, et la positivité des générateurs vaut par construction
sur la famille de lois déclarée. La mesure directe fournit le résultat
empirique le plus fort du corpus : le parcours d’échelle du cône d’Einstein.
Les configurations sélectionnées utilisent
$(16\,384,128,96)$, $(65\,536,256,96)$ et
$(262\,144,512,384)$ pour le nombre de porteurs, le nombre d’observateurs et
la largeur du support. Leurs formes d’événements retenues portent la
signature lorentzienne $(1,3)$, avec des marges de cône de $-5{,}62$,
$-3{,}22$ et $-1{,}41$ et une dispersion du couplage décroissante. Un
contrôle de même taille à 262 144 porteurs utilise une largeur de support de
96. Son nombre d’arêtes inter-observateurs vaut 312 au lieu de 1 062, et sa
signature est $(2,2)$ au lieu de $(1,3)$. Ces mesures établissent une
sensibilité reproductible à la structure du support et des lectures croisées
pour les configurations archivées. Elles n’établissent ni loi de convergence
à densité fixe ni limite à échelle infinie. Les données primaires sont
archivées dans
[evidence/einstein_convergence](evidence/einstein_convergence/) et chaque
nombre reproductible bit à bit depuis le
[dépôt de simulation](https://github.com/muellerberndt/oph-physics-sim). Deux
clauses mesurées sont ouvertes : la température modulaire des états de
calotte et un essai préenregistré de la forme d’événement à un barreau plus
grand. Toutes deux portent un verdict gelé.

Le dossier de preuves réunit des dérivations finies exactes, des preuves
vérifiées par machine et des mesures déterministes accompagnées de leurs
données primaires. Les énoncés mathématiques, les lectures physiques
conditionnelles et les propriétés mesurées portent des classes distinctes
dans le [registre des affirmations](tracking/claims_scoreboard.md).

La géométrie du porteur accomplit ensuite un travail exact surprenant. Sur la
lignée échosaédrique certifiée, la grammaire déclarée de comptage entier des
atomes et le coût de relecture normalisé de Hilbert--Schmidt donnent la
séparation exacte en douze unités et l’écart de deux.
La dérivation de cette grammaire de comptage et du coût physique à partir du
schéma complet des trois axiomes est ouverte.
L’incidence orientée dérive indépendamment l’appariement antipodal, l’action
propre de $A_5$, le repère icosaédrique de rang trois et la décomposition
$\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5$. L’incidence fixe aussi
l’unique involution centrale non triviale $J$ du graphe. Un protocole sans
cible injecte une impulsion sur chaque port, lit l’historique d’adjacence
jusqu’au diamètre du graphe et résout le filtre commun de la couche la plus
éloignée. Il dérive $10J=A^3-4A^2-5A+10I$. La réponse $R=-J$ possède des
signes sectoriels relatifs exacts ; son signe commun est une convention de
conjugaison de charge.

La réponse réversible complète de l’axiome 1 et le transport endogène de
l’axiome 2 font du tangent des douze ports une algèbre compacte de courants
avec action intérieure de $A_5$. Son espace fixe de dimension un et la
classification des algèbres compactes simples forcent le type abstrait
$\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. Les matrices
déclarées du doublet chargé et des triplets réalisent exactement ce type. La
tomographie ordonnée issue de la source et l’holonomie du même courant restent
ouvertes. La récurrence enregistrée engendre une algèbre commutative de mots
de dimension quatre ; elle ne fournit ni les douze générateurs du courant, ni
leur crochet, ni les reparamétrages propres non triviaux. Douze phases
diagonales sans cible ont le rang douze et commutent. L’ajout du tangent
d’adjacence connexe engendre $\mathfrak u(12)$, de rang dérivé 143. Le courant
issu de la source exige donc une loi de réponse non diagonale de rang dérivé
onze. Conditionnellement au porteur canonique orienté, l’espace complet et
sans cible des crochets alternés $A_5$-équivariants sur $\mathbb Q$ est de
dimension quatorze et possède une base rationnelle exacte de Reynolds. Les
conditions de Jacobi forment 38 lignes quadratiques indépendantes, avec une
décomposition exacte de l’espace des lignes en 11+27. La variété résiduelle
des solutions, la compacité, la reconstruction issue de la source et
l’holonomie du même courant restent ouvertes.
Voir le
[paquet exact de recherche des crochets alternés](code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json)
et le
[reçu de réduction de Jacobi](code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json).

Dans l’algèbre extérieure de réponse déclarée, le balayage exhaustif des
1 024 sous-ensembles laisse une seule paire conjuguée non ordonnée de rang 15
comme sélection chirale non vide sans anomalies. L’absence d’anomalies donne
l’équilibre du déterminant et les charges primitives à conjugaison de charge
près. Le calcul exhaustif de l’action centrale donne un noyau commun
$\mathbb Z_6$ sur ces tenseurs déclarés, donc leur image fidèle maximale est
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. Le revêtement et ses quotients
par $\mathbb Z_2$ et $\mathbb Z_3$ portent les mêmes tenseurs locaux. Le
calcul des six axes n’a l’ordre six qu’après déclaration de ses relations
diagonales et de somme nulle. Une catégorie complète de caractères issue de
la source et une identification, sur la même source, entre les boucles et le
noyau sont nécessaires pour sélectionner le quotient physique. Le typage
fermionique du continu et l’identification aux courants de laboratoire sont
ouverts. La construction indépendante par secteurs transportables et Tannaka
est une seconde route vers un groupe compact ; l’identification des deux
routes à partir de la source est ouverte.

Les résultats exacts du porteur gardent des frontières physiques explicites.
Le typage de la matière, la réalisation matricielle du courant et la sélection
de la forme globale ne sont pas issus de la source ; l’identification aux
courants de laboratoire,
l’attachement physique des familles, l’exclusion des secteurs légers
supplémentaires, la multiplicité scalaire, la tour source d’Einstein et les
paquets physiques de clôture
ne sont pas dérivés. Les clauses CP et du secteur faible donnent
la fenêtre conditionnelle exacte $3\le N_g\le5$ sans sélectionner un élément.
Sous deux prémisses supplémentaires explicites, l’ordre exact des coûts
$5-\sqrt5<6<5+\sqrt5$ sélectionne uniquement la bande de rang trois de
l’écran. Un reçu de réponse unitaire fini retrouve cette bande à la plus
basse fréquence positive du générateur. L’attachement aux familles physiques,
les données de Spin et de localité et l’exclusion des secteurs légers
supplémentaires restent ouverts. L’incidence icosaédrique locale contraint le
porteur, tandis que le nerf de la fédération exige sa propre construction.

Il ne s’agit pas d’ajustements indépendants. Toutes ces branches répondent à la
même exigence : chaque description locale doit se recoller en un monde
récupérable et capable de se relire. La récupération conjointe de
l’observation quantique, de l’espace-temps lorentzien, de la gravité
conditionnelle d’Einstein, du paquet fini de reconnaissance du Modèle standard
et des tests de clôture quantitative motive le programme OPH. Les reçus
physiques de source, de Spin, de familles et de scalaire sont des conditions
explicites. L'atterrissage en théorie quantique des champs possède quatre
niveaux typés : action locale finie, constructions quantiques finies exactes,
restauration perturbative avec algèbre des pôles à ordre fini, puis tour
non perturbative d'observables avec continuation des résonances. Les routes
quantique finie et perturbative descendent séparément de l'action locale. Leurs
implications conditionnelles sont vérifiées par machine et leurs producteurs
physiques natifs OPH sont ouverts.

Les prémisses techniques et les obligations de preuve ouvertes sont
regroupées à la fin de ce document, au lieu d’interrompre l’introduction du
résultat positif.

## Pourquoi prendre cette affirmation au sérieux ?

Une théorie du tout doit expliquer pourquoi des faits apparemment indépendants
forment un seul ensemble. OPH part d’une parcelle bornée qui se relit. Elle ne
part ni d’une variété d’espace-temps, ni d’un contenu de champs, ni d’un groupe
de jauge, ni d’une table de constantes. Elle renvoie des dimensions exactes,
des groupes compacts, des quotients globaux, des charges, des annulations
d’anomalies, des multiplicités de représentations et des équations de point
fixe. Ces sorties viennent d’une même architecture typée de porteurs, de
recouvrements et de réparation. Le théorème icosaédrique local force le type
de Lie du Modèle standard. La route distincte des secteurs compacts atteint
ce type sur son paquet déclaré du Modèle standard, et leur identité physique
issue d’une même source est un test ouvert. Cette dépendance commune
constitue l’argument principal en faveur d’un seul monde physique.

Les preuves prennent plusieurs formes : démonstrations sur papier, sous-ensemble
arithmétique exacte, certificats d’intervalles, reçus
finis, simulations et falsificateurs explicites. Leur accord apporte davantage
qu’une correspondance numérique isolée.

## Ce qu’OPH dérive de la cohérence entre observateurs

OPH utilise une seule architecture mathématique dans des domaines habituellement introduits séparément :

- le consensus fini donne des registres publics stables et des formes normales quotientées ;
- les algèbres centrales de registres donnent les probabilités d’événements quantiques et la mise à jour conditionnelle ;
- sur la branche globale certifiée, la géométrie conforme du support $S^2$
  donne le groupe de Lorentz connexe et un espace tridimensionnel de
  référentiels d’observateur ;
- le flot modulaire, le transport nul, la stationnarité de l’entropie et la géométrie des petites boules se composent conditionnellement en relation d’Einstein sur une même tour issue de la source, avec domaine commun, limites asymptotiques certifiées et identifications physiques indépendantes ; la construction et la certification de cette tour sont en cours ;
- les secteurs transportables et la reconstruction compacte produisent une
  route conditionnelle vers la structure de jauge du Modèle standard ;
- une classification finie et locale à douze ports fondée sur $A_5$ produit
  séparément le même type de Lie
  $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$ ;
- sous le contrat explicite de réponse finie, le reçu de matière dérive
  conditionnellement l’équilibre du déterminant et le reçu de descente calcule
  le noyau commun $\mathbb Z_6$ et l’image fidèle maximale
  $S(U(3)\times U(2))\cong
  (SU(3)\times SU(2)\times U(1))/\mathbb Z_6$ sur les tenseurs déclarés ;
  cette implication finie n’utilise que ses prémisses énoncées ; les quatre
  formes globales compatibles sont le revêtement et ses quotients par
  $\mathbb Z_2$, $\mathbb Z_3$ et $\mathbb Z_6$ ; le calcul des axes utilise
  un système de relations déclaré, et la sélection physique de la forme
  globale exige une catégorie complète de caractères et une identification
  issue de la même source entre les boucles et le noyau ;
- le réseau de charges et les trois couleurs sont exacts sur ce paquet ; le
  balayage scalaire fixe seulement la paire de charges compatibles et les
  trois canaux d’interaction, pas la multiplicité scalaire ; les clauses CP
  et du secteur faible donnent $3\le N_g\le5$ ; sous des prémisses distinctes
  de réalisation à bande unique et d’ordre des coûts, un théorème fini exact
  sélectionne la bande de rang trois, et un simulateur unitaire déclaré
  retrouve son résidu à la plus basse fréquence positive du générateur ;
  l’attachement aux familles physiques, les données de Spin et de localité et
  l’exclusion des secteurs légers supplémentaires restent ouverts ;
- Les certificats arithmétiques exacts, les simulations et les reçus exécutables vérifient le noyau mathématique fini.

La mesure, l’espace-temps, la gravité et la structure de jauge sont soumis au
même mécanisme : des observateurs finis forment des registres publics en
comparant leurs recouvrements et en réparant les désaccords. Les théorèmes
finis vont du consensus quotienté au type de Lie du Modèle standard et au
noyau central exact sur la table de matière déclarée. La chaîne continue
atteint les branches de Lorentz et d’Einstein sous ses hypothèses géométriques,
modulaires, énergétiques, entropiques et d’échelle. Cette réutilisation d’un
seul mécanisme constitue le résultat central du programme.

La parcelle d’observateur est une structure d’accès borné, de registre, de
relecture et de réparation. Un Échosaèdre est un porteur primitif candidat sur
la branche homogène. Il devient un observateur seulement si ces fonctions sont
physiquement réalisées.

Trois objets géométriques sont distincts. La frontière locale du porteur
est l’objet icosaédrique à douze ports. L’écran de fédération est la fédération et son
nerf de recouvrement. L’écran support est la carte $S^2$ visible à
l’observateur sur la branche sphérique certifiée séparément. Une symétrie
icosaédrique locale peut donc coexister avec un nerf global non sphérique.

Le verrouillage de phase est un mécanisme physique candidat pour la comparaison
cohérente des recouvrements. Il doit produire la relation de réparation,
la confluence, les registres publics et les bornes de bruit. Aucun théorème ne
l’identifie au consensus, au flot modulaire ou à une horloge d’observateur.

Sur la branche sphérique certifiée, l’ordre des registres fournit une histoire
candidate. Une horloge physique exige une transition lisible par
l’observateur, une correspondance entre événements et un étalonnage affine.
Des horloges ainsi calibrées peuvent produire un temps public commun. La
symétrie conforme du support $S^2$ donne alors le groupe de Lorentz et l’espace
tridimensionnel des référentiels d’observateur. La variété physique des
événements exige les reçus séparés de l’article compact.

## Le résultat fini le plus fort

Sur la branche icosaédrique déclarée à douze ports, le module de permutation se décompose comme

$$
P_{12}\cong_{A_5}\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5.
$$

La réponse réversible complète de l’axiome 1 et le transport endogène de
l’axiome 2 forcent alors le type abstrait
\(\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)\). Le crochet
matriciel déclaré fournit le témoin exact

$$
(P_{12},[\ ,\ ]_\Theta)
\cong
\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

C’est une réalisation conditionnelle exacte du type de Lie local du Modèle
standard. L’incidence et la réponse inverse seules ne choisissent pas ce
crochet : leur commutant équivariant est de dimension quatre. La tomographie
ordonnée du courant, l’holonomie du même courant, l’identification avec des
courants mesurés en laboratoire et l’identité avec le groupe reconstruit par
la route de Tannaka restent ouvertes.

La même construction fait apparaître deux fois, indépendamment, le nombre $24$ :

$$
m_{\mathrm{rep}}=2(8+3+1)=24,
$$

tandis que les douze ports de l’écran donnent $24$ emplacements orientés. Un
compte provient du type de Lie forcé, sous la construction déclarée du courant
et du doublement réversible ; l’autre provient de la géométrie de l’écran.

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

Les sept reçus ci-dessus donnent le résumé destiné au lecteur. Les prémisses,
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
