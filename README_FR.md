# Holographie des parcelles d’observateur

> La réalité est le monde public stable reconstruit par des observateurs finis et auto-lecteurs qui comparent leurs recouvrements et réparent leurs désaccords.

[Read in English](README.md) · [Livre](https://oph-book.floatingpragma.io/) · [Manuels](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [OMEGA](https://omega.floatingpragma.io/)

L’Holographie des parcelles d’observateur, ou OPH, est une théorie du tout
sans boutons de réglage, construite sur une thèse centrale : **les
observateurs sont premiers, et la réalité objective est émergente.** La
physique commence habituellement en fournissant un espace-temps, des champs
quantiques, un groupe de jauge et une table de constantes mesurées. OPH
commence par des observateurs : des systèmes bornés dotés d’un état local,
d’une relecture d’eux-mêmes et de leurs voisins, de registres et de
mouvements de réparation. Elle en dérive le reste. La réalité émerge de la
réparation des recouvrements d’observateurs sur un écran holographique. De
trois axiomes et de deux constantes, $P$ et $N$, l’univers observé se déploie :
la mesure quantique, l’espace-temps lorentzien, la branche conditionnelle
d’Einstein, les symétries de jauge et la matière sont des relectures d’un
seul système fini de cohérence entre observateurs, sous leurs prémisses
déclarées.

## Commencer ici

La physique a plusieurs fois révisé son idée de ce qui est fondamental. L’espace fut
absolu jusqu’à devenir relatif ; la matière fut continue jusqu’à être
quantifiée. Chaque révision paraissait scandaleuse depuis l’intérieur de
l’ancien tableau et évidente depuis l’intérieur du suivant. OPH opère la
révision suivante. L’observateur, traité pendant un siècle comme une gêne aux
marges de la mécanique quantique, passe au fondement, et l’espace-temps, la
matière et les constantes en sortent comme des résultats. Le matériel
ci-dessous vous fait traverser ce basculement sans prérequis.

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
  affichent des données réelles de la dynamique de réparation. Vous regardez
  l’espace-temps et la matière émerger à l’écran au lieu de croire les
  articles sur parole.

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

## Six reçus reproductibles

Ces six artefacts publics renvoient directement à leurs preuves, données ou
certificats :

1. **Un espace-temps quadridimensionnel, mesuré en train d’émerger.** Un
   parcours ajusté par la largeur du support à 16k, 65k et 262k porteurs donne
   à la forme d’événement tenue à l’écart la signature lorentzienne $(1,3)$
   (un temps, trois espaces), avec des marges de cône de $-5{,}62$,
   $-3{,}22$ et $-1{,}41$. À 262k porteurs, réduire la largeur du support de
   384 à 96 fait passer le nombre d’arêtes inter-observateurs de 1 062 à 312
   et la signature de $(1,3)$ à $(2,2)$. Données brutes :
   [evidence/einstein_convergence](evidence/einstein_convergence/) ; chaque
   nombre se régénère bit à bit.
2. **Un noyau vérifié machine qui se surveille lui-même.** Une bibliothèque
   Lean 4 sans `sorry` de plus de 1000 théorèmes et lemmes couvre le noyau de
   consensus, le théorème d’identifiabilité de jauge, l’algèbre finie de
   l’écran et la composition de la branche d’Einstein. Chaque théorème public
   porte son rapport d’axiomes. [Lean/](Lean/)
3. **Une clôture sans dimension avec un statut arithmétique certifié.** La
   clôture de pixel $P=\varphi+\sqrt\pi/A_T(P)$ possède une racine unique
   certifiée machine pour chaque application déclarée, sans valeur continue
   ajustée. Son identification physique à la limite de Thomson exige un
   transport hadronique issu de la source. La comparaison enregistrée est un
   diagnostic, dont la portée et les preuves figurent dans le
   [registre des affirmations](tracking/claims_scoreboard.md).
4. **Un diagnostic des leptons chargés avec un test de clôture déclaré.** La
   surface de clôture empirique porte une cible de confirmation ou de rejet et
   une ascendance explicite des entrées. Elle n’établit aucune prédiction de
   masse issue de la source. Le [registre des postdictions](docs/POSTDICTION_LEDGER.md)
   donne la table complète.
5. **Un théorème exact de Koide dans la chambre positive.** Une réponse
   hermitienne $C_3$ sur la fibre d’une face icosaédrique obéit à
   $Q=1/3+(2/3)(|b|/a)^2$ ; ainsi $Q=2/3$ équivaut exactement à
   $|b|/a=1/\sqrt2$ dans la chambre aux valeurs propres non négatives. Des
   blocs de rang deux égaux et l’application finie de
   Gelfand–Naimark–Segal traciale donnent cet équilibre sous les prémisses du
   paquet d’événements déclaré. L’attachement physique aux familles chirales,
   la phase et les rapports numériques sont ouverts. L’[article autonome sur
   l’identité de Koide](extra/koide_identity_from_positive_c3_face_circulants.pdf)
   expose le théorème, sa preuve, la frontière de formalisation et la
   provenance du diagnostic numérique informé par la cible.
6. **Une loi exacte de capacité finie et une normalisation exacte du choc en
   espace de Sitter.** La maximisation de l’entropie généralisée finie sur les
   probabilités de secteur donne exactement $\log M$. Le transfert uniforme
   d’une fraction $f$ de la capacité de l’écran modifie l’entropie extrémale
   et la coordonnée logarithmique uniforme des secteurs de $\log(1-f)<0$.
   Dans l’espace-temps de Sitter pur, le coefficient du choc vérifie
   $\mu^2=d-2$, exactement la valeur propre $\ell=1$ du laplacien sphérique,
   indépendamment du rayon de l’horizon. Lean vérifie le noyau algébrique. La
   lecture du transfert comme
   choc physique produisant une avance temporelle exige les dictionnaires
   déclarés pour l’horizon, la masse de l’observateur, la gravitation, les
   modes de jauge et l’opérateur cinétique. L’[article ciblé sur de Sitter](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf)
   donne le théorème fini et sa frontière physique.

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

$P$ est le rapport de pixel local : la taille de la cellule élémentaire
d’observation en unités géométriques naturelles. OPH ne choisit pas ce grain
en ajustant la constante de structure fine. La cellule doit être cohérente avec
le processus d’observation qu’elle rend elle-même possible :

$$
\boxed{P_\star=\varphi+\frac{\sqrt\pi}{A_T(P_\star)}}.
$$

Intuitivement, $P$ est la **résolution** de l’univers. Le modifier à la main
changerait simultanément la géométrie cellulaire, le spectre de réparation, les
largeurs de jauge et la hiérarchie des particules. La clôture en fait une sortie
de l’architecture. Ici, $A_T(P)$ est l’inverse du couplage
électromagnétique dans la limite de Thomson produit par une cellule d’essai.

C’est le plus fort des deux résultats quantitatifs de clôture. Le théorème de
point fixe utilisé par le calcul établit qu’une application de
l’intervalle physique dans lui-même, dont la constante de contraction est
strictement inférieure à un, possède exactement un point fixe. Des certificats
d’intervalles à arrondi extérieur vérifient ces hypothèses pour chaque
application déclarée de $P$ et excluent une seconde racine sur son domaine
analytique complet. L’identification avec la limite physique de Thomson exige
le transport hadronique issu de la source dans le même schéma. Sa construction
est en cours.

Le [registre des affirmations](tracking/claims_scoreboard.md) donne la
racine, la comparaison externe, le résidu et la classe de l’énoncé. La
comparaison emploie $P_C$, défini à partir de la limite mesurée. Le transport
hadronique issu de la source dans le même schéma est absent. La comparaison
enregistrée a un statut de diagnostic, hors du
champ d’un énoncé sur la constante de structure fine physique.

**$N$ est la capacité de registres publics** de l’ensemble du système
d’observateurs, ou dans le langage de la simulation, la mémoire corrigible que
porte le substrat. $N$ est secondaire. L’univers observé peut simplement être
lu : $N$ se rétro-ingénie à partir de la mesure comme n’importe quel réglage
de machine se rétro-ingénie à partir du comportement de la machine, et aucun
résultat de la reconstruction centrale ne dépend de sa dérivation à partir de
premiers principes. Une condition conditionnelle d’auto-lecture,
$N=\log M_0(\mathfrak U_N)$, propose de la retrouver depuis la capacité des
registres publics corrigibles ; sa branche de comptage finie est exacte et son
attachement physique n’est pas établi.

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

À partir de ce contrat de réponse, un relèvement compact équivariant explicite
construit
$\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. Le paquet
extérieur à quinze états est sélectionné par le balayage d’anomalies
exhaustif des 1024 sous-ensembles : la paire conjuguée non ordonnée de rang
15 est l’unique sélection chirale non vide sans anomalies, avec la
graduation de parité fermionique en sortie. L’absence d’anomalies donne
l’équilibre du déterminant et les charges primitives à conjugaison de
charge près. Le calcul exhaustif de l’action centrale donne un noyau commun
$\mathbb Z_6$ sur ces tenseurs, donc leur image fidèle maximale est
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. Le revêtement et ses quotients
par $\mathbb Z_2$ et $\mathbb Z_3$ portent les mêmes tenseurs locaux, et les
données de secteurs de flux mesurées du certificat de descente sélectionnent
le quotient $\mathbb Z_6$ à portée de source finie. Cette
implication finie exacte n’utilise que ses prémisses énoncées.
Au-delà de cette portée de source finie, le typage fermionique du continu,
la forme globale du continu et
l’identification aux courants de laboratoire sont ouverts. La construction
indépendante par secteurs transportables et Tannaka est une seconde route
vers un groupe compact ; l’identification des deux routes à partir de la
source est ouverte.

Les résultats exacts du porteur gardent des frontières physiques explicites.
Le typage de la matière et la sélection de la forme globale ne sont mesurés
qu’à portée de source finie ; l’identification aux courants de laboratoire,
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
recouvrements et de réparation. La route icosaédrique locale et la route des
secteurs compacts se rencontrent au niveau du type de Lie du Modèle standard. Leur
identité physique issue de la source est un test ouvert. Cette dépendance
commune constitue l’argument principal en faveur d’un seul monde physique.

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
  $\mathbb Z_2$, $\mathbb Z_3$ et $\mathbb Z_6$, et la sélection physique de
  la forme globale est portée par les données de secteurs de flux mesurées
  du certificat de descente ;
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

La mesure, l’espace-temps, la gravité et la structure de jauge sont soumis au même mécanisme : des observateurs finis forment des registres publics en comparant leurs recouvrements et en réparant les désaccords. Les théorèmes finis vont du consensus quotienté à l’algèbre de jauge du Modèle standard et à sa forme globale. La chaîne continue atteint les branches de Lorentz et d’Einstein sous ses hypothèses géométriques, modulaires, énergétiques, entropiques et d’échelle. Cette réutilisation d’un seul mécanisme constitue le résultat central du programme.

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

Un rappel équivariant explicite du commutateur par blocs construit alors

$$
(P_{12},[\ ,\ ]_\Theta)
\cong
\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

C’est une construction exacte du type de Lie local du Modèle standard sous le
contrat explicite de réponse finie. L’incidence seule ne choisit pas une
réponse linéaire arbitraire : son commutant équivariant est de dimension
quatre. La sélection de cette réponse par une loi physique indépendante,
l’identification avec des courants mesurés en laboratoire et l’identité avec
le groupe reconstruit par la route de Tannaka restent ouvertes.

La même construction fait apparaître deux fois, indépendamment, le nombre $24$ :

$$
m_{\mathrm{rep}}=2(8+3+1)=24,
$$

tandis que les douze ports de l’écran donnent $24$ emplacements orientés. Un compte provient de l’algèbre de jauge reconstruite ; l’autre de la géométrie de l’écran.

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
                → constructions conditionnelles de courant et de matière sélectionnée par balayage
                → noyau tensoriel Z6 exact ; forme globale mesurée par flux à portée de source
                attachements du courant de laboratoire, du scalaire, du spectre et des familles ouverts
        ↓
tests quantitatifs de clôture et de lecture physique
```

Les hypothèses détaillées et les types de reçus sont énoncés dans les articles. La page d’accueil du dépôt est volontairement une carte du résultat positif, et non un substitut à ces énoncés de théorèmes.

## Résultats en un coup d’œil

| Résultat | Contribution d’OPH | Source principale |
| --- | --- | --- |
| Consensus fini | Réparation terminante, lecture protégée, formes normales quotientées indépendantes de l’ordonnancement et registres centraux | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| Surface d’événements quantiques | Probabilités de Born, conditionnement de Lüders et borne de Tsirelson sur la surface finie des registres | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| Relativité | Sur la branche globale certifiée avec une comparaison complète et indépendante des algèbres et états sur la même tour, $\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ et $H^3\cong\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$ | [Article espace-temps et Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Dynamique d’Einstein | Chaîne typée conditionnelle du transport modulaire et nul à $G_{ab}+\Lambda g_{ab}=8\pi G\langle T_{ab}\rangle$ sur une tour issue de la source et de domaine commun ; sa construction et sa certification sont en cours | [Article espace-temps et Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Théorème fini de courant $A_5$ | La réalisation déclarée avec comptage entier et coût normalisé de Hilbert--Schmidt donne la séparation exacte en douze unités. Indépendamment, l’incidence orientée donne l’appariement inverse, l’action propre de $A_5$, un repère de rang trois et l’unique involution centrale du graphe. La dérivation de la normalisation entière et du coût physique discret à partir du schéma complet des trois axiomes est ouverte. Sous le contrat explicite de réponse centrale signée et involutive, les réponses admissibles sont $\pm J$, avec des signes sectoriels relatifs exacts, et un relèvement compact explicite réalise $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. La sélection physique indépendante de la réponse et l’identification aux courants de laboratoire sont ouvertes ; aucune conclusion automatique sur un support global $S^2$ | [Article sur la jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Image fidèle maximale conditionnelle du Modèle standard | Sur la paire de modules extérieurs conjugués à quinze états sélectionnée par le balayage exhaustif, l’équilibre des anomalies fixe les charges primitives à conjugaison près. Le noyau commun exact est $\mathbb Z_6$, donc l’image fidèle maximale est $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. Le revêtement et ses quotients par $\mathbb Z_2$ et $\mathbb Z_3$ portent les mêmes tenseurs locaux ; la sélection physique de la forme globale est portée par les données de secteurs de flux mesurées du certificat de descente à portée de source finie. Cette implication finie n’utilise que ses prémisses énoncées | [Article sur la jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Structure de la matière | Modules extérieurs conditionnels exacts d’une génération, équilibre des hypercharges et des anomalies, porteur à trois couleurs, paire de charges scalaires compatibles et trois canaux d’interaction ; la multiplicité scalaire, la sélection physique de la matière, l’attachement des trois familles et l’exclusion de secteurs légers supplémentaires sont ouverts. Le nombre de générations est une complétion déclarée à l’intérieur de la fenêtre conditionnelle $3\le N_g\le5$ jusqu’à la dérivation de son attachement familial | [Article sur la jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Atterrissage en théorie quantique des champs | Invariance de l’action finie ; critères quantiques exacts de ligne déterminante et de hamiltonien ; restauration perturbative formelle et algèbre W/Z stricte à ordre fini ; reconstruction non perturbative et implications de résonance séparées. Les routes quantique finie et perturbative descendent en parallèle de l’action locale, avec leurs constructions issues de la source comme portes physiques explicites | [Article sur la jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Écran fini en espace de Sitter | Normalisation exacte du choc dans l’espace de Sitter pur, maximum d’entropie fini, loi uniforme de transfert de capacité pour la coordonnée logarithmique des secteurs et courbure analytique ; la lecture physique de l’avance temporelle est conditionnelle aux dictionnaires d’horizon et de choc énoncés dans l’article ciblé | [Article sur la capacité d’un écran fini en espace de Sitter](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf) |
| Pôles physiques W/Z | L’application stricte à une boucle qui transforme un paquet renormalisé complet en pôles complexes chargé et neutre est démontrée et vérifiée par machine ; les conventions de signe, de feuillet, d’ordre, de mélange neutre et la séparation entre coefficients stricts et racine carrée sont figées. Le jeu numérique est une régression de backend post-exposition ; l’appariement à la source, un moteur indépendant de symétrie de jauge, la covariance, les amplitudes de courant physique et l’horloge sont ouverts, donc aucun pôle natif d’OPH n’est promu | [Article sur les particules](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Clôture locale $P$ | $P=\varphi+\sqrt\pi/A_T(P)$ ; le schéma d’unicité du point fixe et les certificats d’intervalles donnent une racine par application déclarée ; le transport physique de Thomson est en cours | [Article sur la constante de structure fine](extra/fine_structure_constant_derivation.pdf) |
| Extension globale conditionnelle $N$ | $N=\log M_0(\mathfrak U_N)$, avec $M_0(q)=\alpha(G_q)$ et $M_0=\lvert X_{\rm reach}\rvert$ sur la branche réversible ; un paquet de simulateur issu de la source à coupure fixe $D=24$ est certifié, tandis que son attachement physique, la famille indexée et le zéro unique du défaut sont en cours | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| Pont $N$--Higgs | Relation conditionnelle $R_{\rm EW}=\alpha_U(P)\log(N/\pi)-6\pi/P$ issue du porteur de charge commun à l’écran et au secteur faible | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Vérification exacte | Certificats d’intervalles, reçus finis et simulations reproductibles | [`code/`](code) |

## Choisir un parcours de lecture

| Pour découvrir... | Commencer ici |
| --- | --- |
| L’argument persuasif le plus court | [Le cas compact pour OPH](extra/compact_proof_of_oph.pdf) |
| La dérivation de l’espace-temps et d’Einstein | [Espace-temps des observateurs et dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Les deux routes de jauge du Modèle standard | [Structure de jauge du Modèle standard](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| La synthèse complète | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
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
partie de cette clôture est vérifiée machine en Lean. L’application déclarée
de $P$ possède un point fixe certifié, tandis que sa comparaison à la
constante de structure fine physique conserve un statut de diagnostic. Les
conditions de clôture et leurs données physiques absentes sont précisées dans
le [programme de falsification OPH](docs/OPH_FALSIFICATION_PROGRAM.md).
Une clôture physique des deux constantes donnerait une branche sans paramètre
continu, les deux valeurs étant rendues par l’architecture. Cet attachement
physique est ouvert. Les théorèmes de point fixe certifient les racines des
applications déclarées ; ils ne transforment pas un bassin observé ou une
coordonnée définie par la cible en dérivation physique. La clôture de $N$ par
premiers principes est en cours. La lecture de $N$ dans l’univers laisse
intactes les conséquences des trois axiomes.

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
