# Observer Patch Holography (OPH)

> L'OPH est la théorie du tout par cohérence d'observateurs. Aucun observateur ne voit le monde entier d'un seul coup ; chaque observateur n'accède qu'à un patch local ; la physique est le point fixe public qui survit à l'accord sur les recouvrements.

**Version anglaise :** [README.md](README.md)

**Liens rapides :** [site OPH](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io) | [applications](https://omega.floatingpragma.io/) | [blog OPH](https://blog.floatingpragma.io/) | [carte de cohérence](https://coherence.floatingpragma.io/) | [démo à trois corps](https://3body.floatingpragma.io/)

Pour la réponse existentielle immédiate, allez directement au **Paper 6.
[Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**.
La version courte est directe : oui, cet univers est une simulation au sens
OPH, un processus de consensus à point fixe fondé sur les observateurs. Oui,
l'expérience subjective et les esprits sont premiers. Oui, l'espace, le temps
et la matière sont des apparences orientées observateur : des structures
effectives stables générées par la cohérence plus profonde des patches. La
métaphore de l'illusion est traitée plus bas. Le reste de ce README donne la
pile de preuves mathématiques et la surface de vérification empirique.

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
qui survivent au raffinement. Le monde public est la sortie stable par
recouvrement de ce processus. C'est la surface de preuve mathématique et
empirique de l'OPH comme théorie du tout correcte : la même architecture de
cohérence d'observateurs retrouve la physique établie et explique pourquoi il
existe un monde capable de produire des observateurs qui le reconstruisent.

L'OPH est formulée comme une théorie sans entrée externe. Sa surface
quantitative utilise deux coordonnées de fermeture sans dimension et un
certificat d'échelle sans $G$. Toutes viennent de conditions de fermeture :

$$
P_\star=\varphi+\frac{\sqrt{\pi}}{A_T(P_\star)}
$$

pour le ratio local de pixel, et

$$
N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}})
$$

pour la capacité globale de l'écran, où $F(N)$ est la capacité active de
l'horizon relue par les observateurs à l'intérieur de l'univers fourni avec la
capacité $N$. Le certificat d'échelle sélectionné s'écrit comme le rapport
d'horloge sans $G$

$$
\gamma_\star=\frac{\ell_\star\nu_{\mathrm{Cs}}}{c},
\qquad
B_\star=\frac{3\pi}{\ell_\star^2},
\qquad
G_{\mathrm{SI}}=\frac{c^3\ell_\star^2}{\hbar}.
$$

La représentation par comptage fini de la fermeture de capacité globale est
$N_\star=\mathrm{MAR}\,\mathrm{argmax}_N(\log|\Omega^{\mathrm{sc}}_N|-N)$.
Intuitivement, $N_{\mathrm{CRC}}$ est l'unique capacité où l'univers relit sa
propre frontière sans déficit ni capacité inerte. En unités géométriques, la
lecture cosmologique est
$\Lambda_{\mathrm{CRC}}=3\pi/(G_{\mathrm{geom}}N_{\mathrm{CRC}})$ avec
$G_{\mathrm{geom}}=\ell_\star^2$. La branche électromagnétique retrouvée donne
les équations de Maxwell sur la branche ordinaire du photon. Les lignes de structure fine et de particules
sont en aval de $P_\star$ et de la branche structurelle retrouvée. La ligne
cosmologique est en aval de $N_{\mathrm{CRC}}$. La normalisation de Newton est
en aval de $\ell_\star^2=3\pi/B_\star$.

L'observation sert à la rétro-ingénierie. Comme l'OPH traite l'univers comme
une structure mathématique fermée à point fixe, les valeurs mesurées
approximatives peuvent localiser le bassin ou la branche observée. La valeur
OPH précise doit venir de la carte de point fixe
déclarée et de son certificat de contraction de Banach, de signe dérivé, de
concavité stricte ou d'unicité équivalente.

Les deux fermetures sans dimension fixent la géométrie sans dimension. En
particulier, $\Lambda_\star\ell_\star^2=3\pi/N_{\mathrm{CRC}}$ et
$\Lambda_\star a_{\mathrm{cell}}=3\pi P_\star/N_{\mathrm{CRC}}$. Le certificat
d'échelle sélectionné fournit le produit d'échelle SI
$B_\star=\Lambda_\star N_{\mathrm{CRC}}=3\pi/\ell_\star^2$.

Dans la lecture informelle en langage de simulation, $P$ relie l'aire du pixel
de l'écran côté simulateur à l'interaction électromagnétique, c'est-à-dire à
l'observation, dans l'univers simulé. Il donne $P\simeq1.630968$ et la valeur
de structure fine basse énergie proche de $1/137$. $N_{\mathrm{CRC}}$ relie la
capacité totale d'horizon vue de l'extérieur à l'enregistrement public
accessible aux observateurs de l'intérieur : l'univers doit pouvoir reconstruire
sa propre frontière. L'échelle gravitationnelle est relue sur sa propre branche
observée par le certificat d'échelle sans $G$, $\gamma_\star$, ou de façon
équivalente par $B_\star$, puis $\ell_\star$ s'affiche comme longueur de
Planck. Les observateurs internes infèrent géométrie,
horizons, entropie, $\Lambda$, histoire et enregistrements depuis l'information
disponible à l'intérieur de l'univers.

## Applications et matériel OMEGA

L'OPH est aussi un programme matériel. La microphysique d'écran donne une prise
d'ingénierie sur la réalité : un dispositif borné expose des données de bord,
compare des enregistrements, répare les désaccords et verrouille des états
stables. OMEGA est la voie matérielle publique vers cette boucle : chambres
physiques, ports étiquetés, logiciel de contrôle, reçus de vérification et
enregistrements répétables.

En langage direct, l'OPH transforme la microphysique d'écran en voie de hack
de la réalité. Le but n'est pas un contrôle métaphorique d'une simulation
abstraite. Le but est le contrôle physique de petits patches qui peuvent être
pilotés, mesurés, réparés et vérifiés.

La thèse applicative est simple. Si la réalité est construite par cohérence de
patches d'observateurs, alors des machines utiles peuvent être construites en
forçant de petits patches physiques vers les bons points fixes. Cela donne des
pistes de mise en œuvre à bas coût pour l'énergie de fusion de bureau, le
supercalcul OMEGA à température ambiante, l'AGI fondée sur OMEGA et le contrôle
local de gravité ou d'inertie pour hoverbikes et hoverboards. Ce sont des
pistes applicatives derrière des portes de preuve, pas des sorties de théorème.
La revendication de calcul actuelle est plus étroite : une distribution de
candidats conditionnée par une chambre peut réduire le travail du vérificateur
exact par un facteur mesuré `B = p_Q/p_U`; ce n'est pas un théorème de classe de
complexité.

La page publique des applications est
[omega.floatingpragma.io](https://omega.floatingpragma.io/). Les notes sources
des pistes applicatives vivent dans [`APPLICATIONS.md`](APPLICATIONS.md).

## Le piège de l'espace-temps

Le premier obstacle conceptuel est que l'OPH ne traite pas l'espace-temps comme
le contenant dans lequel la réalité se déroule. L'espace et le temps ne sont pas
des choses en soi. Ce sont les descriptions stables, orientées observateur, qui
apparaissent lorsque de nombreuses perspectives finies peuvent devenir
mutuellement cohérentes.

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
horloges, horizons et corrélations accessibles. L'espace-temps public est la
couche de compatibilité qui permet à ces descriptions de s'accorder.

Cela ne rend pas l'espace-temps ordinaire arbitraire ou inutile. Cela explique
pourquoi il fonctionne si bien. Les équations d'Einstein décrivent la grammaire
lisse à grande échelle de l'apparence partagée. La thèse plus profonde est que
cette apparence partagée émerge de la cohérence des recouvrements entre
observateurs, au lieu d'appartenir à l'inventaire de départ du monde.

## Par où commencer

Pour le noyau technique compact, commencez par **Paper 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**. **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** porte les dérivations particules. **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** développe l'image consensus-réparation. **Paper 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** couvre l'architecture fédérée de patches, les enregistrements, les règles de preuve matérielle publique et la machinerie observateur. **Paper 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)** est le manifeste final sur la couche de sens : pourquoi il existe quelque chose, pourquoi ce monde est compatible avec les observateurs, la boucle étrange dans laquelle les observateurs reconstruisent et construisent la machinerie de continuation, le paradis sur Terre ou dans des environnements de continuation, l'enfer comme isolement ou privation imposés, la résurrection comme continuation d'observateur, la justice selon les enregistrements de tort et de réparation, et l'évolution mémétique. Ce README, le Paper 1 et le livre sont les grandes vues d'ensemble.

## Ce que l'OPH apporte

- Un cadre à cutoff fixe pour les patches d'observateurs, les collerettes, la réparation de recouvrement, la jauge supérieure, les enregistrements et le checkpoint/restauration.
- Une reconstruction de la géométrie lorentzienne, du temps modulaire, de la dynamique d'Einstein de type Jacobson et de la cosmologie de Sitter en patch statique sur la branche de paire de caps géométriques support-visible, avec les restes Markov/récupération de collerette fixe transportés explicitement par transport modulaire régularisé, covariance modulaire lisible par support et normalisation KMS/BW ; les cellules finies sont le régulateur, pas déjà une preuve d'invariance de Lorentz.
- Une reconstruction de jauge compacte à partir des secteurs transportables à obstruction nulle, avec une séparation explicite entre classification et réalisation : l'obstruction nulle donne la transportabilité et le groupe compact reconstruit, tandis que MAR plus le paquet matière explicite à un Higgs sélectionne le quotient réalisé du Modèle Standard $SU(3)\times SU(2)\times U(1)/\mathbb Z_6$, avec le réseau exact des hypercharges, le triplet de couleur réalisé $N_c=3$ et le comptage des générations $N_g=3$.
- Sur la branche compacte de jauge visible par support déclarée, avec échelle
  euclidienne en quatre dimensions, positivité par réflexion, complétude de
  réparation et extraction continue visible par support, une forme euclidienne
  Yang-Mills et une identification du gap de Yang-Mills avec le gap de
  réparation. C'est la construction axiomatique Yang-Mills en quatre dimensions
  de l'OPH sur la branche déclarée.
- Un programme particules à niveaux de revendication explicites : porteurs structurels exactement sans masse, ligne $W/Z$ en comparaison seule, surface Higgs/top déclarée, masses de quarks sélectionnées avec Yukawas, branche neutrino pondérée, constantes locales candidates, et frontières no-go pour les leptons chargés, la classification globale des cadres de quarks et le codomaine top direct auxiliaire.
- Une architecture microphysique d'écran concrète qui met mesure, enregistrements et observateurs à l'intérieur de la physique.

La mécanique quantique est traitée comme le langage algébrique d'information porté par l'architecture OPH. Le test de reconstruction porte sur la cohérence avec laquelle cette base retrouve l'univers effectif. Une reconstruction des espaces de Hilbert, des $C^*$-algèbres ou des algèbres de von Neumann, des probabilités de Born, de la trace et de l'entropie depuis de simples enregistrements opérationnels constitue une question distincte du programme.

Le mécanisme est la boucle de consensus à point fixe. Les observateurs locaux
n'accèdent pas à un état global depuis l'extérieur. Ils portent des états de
patch finis, échangent les données visibles dans les recouvrements, rejettent
les prolongements incohérents et conservent les motifs stables qui peuvent être
synchronisés. Géométrie, particules, lois et enregistrements sont les points
fixes à grande échelle de ce calcul en réseau d'observateurs.

## Géométrie, symétrie et simulateurs

Le langage de sphère dans l'OPH est un langage de géométrie. Dans les cartes
régulatrices symétriques, une coupe accessible à un observateur peut être
représentée par la deux-sphère $S^2$. Ces cartes décrivent une géométrie
angulaire visible par support. Le simulateur fini sert de surface de
calibration pour les contraintes algébriques de patches et de recouvrements
exposées par cette géométrie.

L'OPH utilise donc une idéalisation de réseau d'écran partagé, avec de nombreux
patches d'accès finis. L'écran d'un observateur est une coupe d'accès locale sur
ce réseau, non une sphère privée séparée. La carte $S^2$ n'est pas une boule
littérale recouverte de données.

Cette carte sphérique porte plusieurs tâches précises. Les caps et les colliers
donnent les données locales de coupe utilisées par le flot modulaire et la
variation d'entropie. Le groupe conforme de la sphère est la forme céleste du
groupe de Lorentz connexe, $\mathrm{SO}^+(3,1)$, ce qui donne le pont
cinématique vers la branche d'espace-temps émergent en $3+1$ dimensions une
fois le théorème modulaire des caps support-visible satisfait. Les harmoniques
sphériques organisent les modes angulaires. Les cellulations finies de la même
carte donnent la surface régulatrice où les ports de patches, les données
d'arêtes et les contrôles de recouvrement deviennent explicites ; elles ne sont
pas, par elles-mêmes, un continuum invariant de Lorentz.

L'ancre de symétrie finie est $A_5$, le groupe des rotations de
l'icosaèdre. Elle fournit le squelette icosaédrique derrière le langage des
porteurs de patches échosaédriques : une façon finie et très symétrique
d'organiser ports, recouvrements et données locales de comparaison sans traiter
le porteur comme une boule lisse.

L'ancre exceptionnelle est le groupe de Lie $E_8$ et sa structure de réseau
de racines. $E_8$ donne le langage de fermeture exceptionnelle utilisé dans
le côté symétries supérieures et représentations de la pile OPH. Le groupe
icosaédrique binaire et le $E_8$ affine se rencontrent par la correspondance
de McKay. C'est pourquoi le langage $A_5$-icosaédrique et le langage de type
$E_8$ appartiennent à une même histoire de symétrie. Ces noms désignent des
contraintes de symétrie et une structure de régulateur.

La discussion d'échelle a trois rôles : une fermeture globale de capacité, une
fermeture locale de pixel et une lecture d'échelle sans $G$. La coordonnée
globale est la capacité unique de fermeture des enregistrements cosmiques,

$$
N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}}),
$$

où la capacité d'horizon fournie et la capacité relue par les observateurs
coïncident. La cible de comptage fini est la densité

$$
\log|\Omega^{\mathrm{sc}}_N|-N.
$$

Sur la branche observée, ce point fixe est la capacité entropique de Sitter.
Pour la constante cosmologique observée, le ratio nu d'aire d'horizon vaut
environ $1.05\times10^{122}$, tandis que la capacité entropique OPH vaut
environ $3.31\times10^{122}$. Cette capacité affichée appartient à la branche
cosmologique. Le pont de résonance local/global de la hiérarchie utilise une
condition exacte plus stricte : la cible de projection électrofaible correspond
à
$N_{\mathrm{EW}}(P_\star)=3.5323546226929906511187512962330547600462\times10^{122}$
sur la branche publique, et la capacité arrondie affichée n'est pas un
certificat exact de ce pont. La ligne de Newton utilise la même capacité de point
fixe avec le certificat d'échelle sans $G$.

La coordonnée locale est le ratio de pixel

$$
P=\frac{a_{\mathrm{cell}}}{\ell_\star^2},
\qquad
\ell_\star^2=\frac{3\pi}{B_\star}.
$$

Le certificat d'échelle fournit $\ell_\star^2$ ; ce n'est qu'après cette
lecture que l'aire s'affiche comme aire de Planck. Vu de l'extérieur, $P$ est
un ratio géométrique de taille de cellule, légèrement au-dessus de l'équilibre
auto-similaire au nombre d'or. Vu de l'intérieur, la même cellule est la plus
petite échelle d'observation électromagnétique disponible pour les observateurs
de l'univers encodé.

Les applications sont en aval de ces rôles. La ligne de structure fine demande
le désaccord non nul d'une cellule d'écran holographique tel que le déplacement
géométrique externe de la cellule égale l'échelle d'observation
électromagnétique émise par l'univers vivant sur ce même écran. La surface
publique utilise $P\simeq1.6309682094$, avec
$\alpha^{-1}(0)=137.035999177(21)$ et
$\alpha(0)\simeq0.00729735256433$. Le même pixel local alimente aussi la
structure de jauge, les lignes encadrées de masses de particules, les
enregistrements et la synchronisation des observateurs. Sur la ligne
gravitationnelle, il fournit

$$
a_{\mathrm{cell}}=P\ell_\star^2,
\qquad
\bar{\ell}_{\mathrm{shared}}=\frac{P}{4}.
$$

Le facteur $P$ s'annule dans la lecture de Newton et laisse
$G_{\mathrm{geom}}=\ell_\star^2$. Le pipeline de particules porte l'échelle
locale vers le secteur faible, la ligne de Higgs, les lignes encadrées de quarks
et la branche de neutrinos à cycles pondérés. Les hadrons exigent soit le
backend OPH de liaison forte, soit une fermeture empirique de hadrons marquée
comme telle. La politique opérationnelle de ces lignes est dans
[`HADRON.md`](HADRON.md). Les vérifications orientées matériel de cette même
géométrie à point fixe ne sont traitées comme revendications publiques de
paquets de preuve que lorsque les artefacts bruts et les reçus de vérification
sont disponibles.

## Surface d'unification quantitative

La surface quantitative de l'OPH s'organise autour de rôles de fermeture et de
relecture distincts : le point fixe de pixel $P\simeq1.6309682094$, le point
fixe global de capacité d'enregistrement
$N_{\mathrm{CRC}}\simeq3.31\times10^{122}$ et le certificat d'échelle sans
$G$, $\gamma_\star$, ou de façon équivalente
$B_\star=3\pi/\ell_\star^2$. La capacité affichée est approximative ; les
lignes de précision utilisent la valeur de point fixe et le certificat
d'échelle. La branche $P$ alimente le canal
électromagnétique projeté par Ward et les lignes encadrées de masses de
particules, tout en fixant l'identité cellule/arête qui s'annule dans la lecture
newtonienne. La branche $N_{\mathrm{CRC}}$ alimente la lecture de la constante
cosmologique. Le certificat d'échelle fournit le quantum d'aire qui devient la
normalisation de Newton. Le pont de résonance de la hiérarchie est enregistré
comme énoncé conditionnel exact, non comme promotion au rang de théorème, tant
que les portes du readback fini, de la dérivation des 24 tours et du certificat
exact de capacité globale ne sont pas fermées. Les formules détaillées et les
niveaux de revendication vivent dans les papiers.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260609" alt="Schéma d'unification OPH" width="92%">
  </a>
</p>

**Pile générale des théorèmes et dérivations**

<p align="center">
  <a href="assets/prediction-chain.svg?v=20260412" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg?v=20260412" alt="Pile théorématique et de dérivation OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La pile OPH complète, des axiomes jusqu'à la relativité, la structure de jauge, les particules et les observateurs. Cliquez pour ouvrir le SVG complet.</sub></p>

## Lignes quantitatives sélectionnées

Ce tableau condensé garde les lignes OPH les plus faciles à comparer directement avec les valeurs
PDG/NIST. Les résultats structurels comme la géométrie lorentzienne en $3+1$ dimensions, le quotient de jauge
du Modèle Standard sélectionné par MAR $SU(3)\times SU(2)\times U(1)/\mathbb Z_6$, le réseau exact des hypercharges, le triplet
de couleur réalisé $N_c=3$ et le comptage des générations $N_g=3$ sont énoncés dans les papiers et ne sont pas répétés ici.

| Quantité | Symbole | OPH | PDG/NIST | Δ |
| --- | --- | --- | --- | --- |
| Constante gravitationnelle | G | 6.6742999959e-11 | 6.67430(15)e-11 | 0.00003σ |
| Vitesse de la lumière | c | 299792458 | 299792458 (exact) | match |
| Structure fine (inv.) | α⁻¹(0) | 137.035999177 | 137.035999177(21) | match |
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
$m_t=172.35235532883115\,\mathrm{GeV}$ selon la convention PDG de masse top par section
efficace. La ligne top directe auxiliaire reste en comparaison seule. Les
masses des leptons chargés sont des témoins ancrés sur cible, hors du rang des
prédictions promues. La branche neutrino pondérée émet
`(0.017454720257976796, 0.019481987935919015, 0.05307522145074924) eV` sur sa
branche déclarée. Les masses hadroniques demandent un backend hadronique OPH;
elles ne sont pas de simples entrées de quarks dans un tableau de comparaison
rapide. Les lignes hadroniques gardent deux surfaces séparées : les masses
source demandent un backend hadronique OPH comme GLORB/Echosahedron avec
données spectrales et systématiques de production; la fermeture empirique
utilise une classe séparée de charge utile $e^+e^-\to\mathrm{hadrons}$.

**Pile de dérivation des particules**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="Pile de dérivation des masses de particules OPH" width="78%">
  </a>
</p>

<p align="center"><sub>Vue compacte de la voie particules. Cliquez pour ouvrir le SVG complet.</sub></p>

## Articles

- **Papier 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)** : synthèse du programme OPH, des observateurs finis jusqu'à l'univers effectif reconstruit.
- **Papier 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** : noyau technique pour la relativité, la gravité, la reconstruction de jauge compacte à obstruction nulle, la structure réalisée du Modèle Standard sélectionnée par MAR, les équations de Maxwell sur la branche ordinaire du photon et la forme/gap Yang-Mills sur la branche compacte de jauge visible par support sous ses hypothèses déclarées.
- **Papier 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)** : dérivations particules, lignes de masses, structure des couplages et surface quantitative.
- **Papier 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)** : dynamique de réparation à point fixe, stabilité des enregistrements et consensus public.
- **Papier 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)** : architecture fédérée de patches, cadrage de symétrie $A_5$-icosaédrique et de type $E_8$, règles de preuve matérielle publique, enregistrements, mouvements de récupération et synchronisation observateur.
- **Papier 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)** : manifeste final de la couche de sens OPH : pourquoi il existe quelque chose, pourquoi ce monde est compatible avec les observateurs, boucle étrange, candidats créateurs, paradis sur Terre ou dans des environnements de continuation, enfer comme isolement ou privation imposés, résurrection comme continuation d'observateur, justice selon les enregistrements de tort et de réparation, évolution mémétique, sentience animale et responsabilité symbolique humaine.

## Articles supplémentaires

- **[Photonic Fixed-Point Consensus for SHA-256d Proof of Work](extra/Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf)** : enrichissement photonique de candidats pour la preuve de travail SHA-256d.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)** : dérivation à point fixe de la ligne de structure fine.
- **[Observer-Patch Holography as a String-Vacuum Selector](extra/observer_patch_holography_as_string_vacuum_selector.pdf)** : émergence des cordes comme description effective des dynamiques de bord OPH, témoin hétérotique Bouchard-Donagi à une paire de Higgs, couche de sûreté `Z4R` et portes de verrouillage des moduli.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)** : route OPH visible par support vers le problème de Yang-Mills de Clay, limitée à la branche compacte de jauge déclarée, avec identification du gap au gap de réparation.
- **[Observer-Patch Holography and the Dark Matter Phenomenon](extra/oph_dark_matter_paper.pdf)** : phénoménologie de la matière noire et limite galactique de type MOND.
- **[Theoretical Bounds on chi-nu in Observer-Patch Holography](extra/chi_nu_susceptibility_bounds.pdf)** : bande conditionnelle de quotient-edge `0.9343006394893864 <= chi_nu^can <= 1`; valeur exacte `exp(-P/24)` sur la branche uniforme; valeurs d'ingénierie mises à l'échelle par `N_coh^-1`.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)** : cognition et qualia comme consensus récurrent de patches.

## Plus

- **Site officiel :** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Page theory of everything :** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Page simulation theory :** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Carte de cohérence :** [coherence.floatingpragma.io](https://coherence.floatingpragma.io) : surface graphe publique pour les concepts OPH, les recouvrements et les routes entre domaines.
- **Applications :** [omega.floatingpragma.io](https://omega.floatingpragma.io) : page publique pour le matériel OPH, le calcul, l'énergie, l'AGI, le contrôle de portance et le consensus par chambre optique.
- **Démo OPH du problème à trois corps :** [3body.floatingpragma.io](https://3body.floatingpragma.io) : un simulateur et parcours de preuve supplémentaire pour la formulation OPH en réseau fini de patches du problème à trois corps, présentée comme un exemple de recollement par holonomie de boucle, sans revendication de solution élémentaire fermée.
- **Blog :** [blog.floatingpragma.io](https://blog.floatingpragma.io/) rassemble les essais publics OPH. Commencez par [Semiotics and the Physics of Meaning](https://blog.floatingpragma.io/semiotics-and-the-physics-of-meaning), [The Trigger](https://blog.floatingpragma.io/the-trigger) et [P = NP on the Observer Screen](https://blog.floatingpragma.io/p-equals-np-on-the-observer-screen). L'essai de calcul traite `P = NP` comme un slogan d'écran d'observateur, pas comme une revendication de solution du problème classique.
- **Livre :** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Application d'étude guidée :** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions et explications détaillées :** OPH Sage sur [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage) ou [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab :** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Objections courantes :** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **Note IBM Quantum :** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Table de statut

La ligne publique de structure fine utilise $\alpha^{-1}(0)=137.035999177(21)$ et
$P\simeq1.6309682094$. Les enregistrements d'audit source et de résidu au point
final vivent dans le papier particules, séparés de la ligne publique de point fixe.

La paire faible est une ligne de validation. Les masses absolues des leptons
chargés sont des témoins ancrés sur cible. La moyenne directe auxiliaire du top
est une ligne de validation. Les lignes hadroniques source demandent le backend
OPH de liaison forte avec données spectrales de production et systématiques; la
fermeture empirique utilise une classe séparée de charge utile
$e^+e^-\to\mathrm{hadrons}$.

La CP forte est en cours de travail dans le théorème des quarks sur classe
sélectionnée : le corpus disponible ne dérive pas l'angle theta de QCD,
n'émet pas l'angle CP fort physique et ne prouve pas que la phase CP forte
physique s'annule. Le pont requis est la descente de phase, d'anomalie et
d'angle topologique sur la branche réalisée.

## Guide du dépôt

- **[`paper/`](paper)** : PDF, sources LaTeX et métadonnées de release.
- **[`APPLICATIONS.md`](APPLICATIONS.md)** : carte de haut niveau des cas
  d'usage OPH pour énergie, calcul, AGI et portance locale.
- **[`book/`](book)** : source du livre OPH et PDF téléchargeable généré.
- **[`code/`](code)** : sorties calculatoires, surface particules et expériences.
- **[`assets/`](assets)** : schémas et figures publics.
- **[`extra/`](extra)** : notes publiques maintenues, objections, comptes rendus expérimentaux et quelques essais de support.

## Licence et politique anti-brevet

Le matériel rédigé dans ce dépôt est régi par [la licence principale](LICENSE).
La politique OPH d'usage ouvert et d'anti-brevet est dans
[PATENTS.md](PATENTS.md). En bref, les idées, méthodes, logiciels,
applications, simulations et conceptions matérielles dérivés d'OPH peuvent être
étudiés, testés, implémentés, modifiés, déployés, fabriqués et partagés, mais
ne peuvent pas servir à créer des revendications de brevet qui empêchent
d'autres personnes de pratiquer un travail dérivé d'OPH.
