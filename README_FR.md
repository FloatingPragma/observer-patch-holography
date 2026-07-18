# Holographie des parcelles d’observateur

> La réalité est le monde public stable reconstruit par des observateurs finis et auto-lecteurs qui comparent leurs recouvrements et réparent leurs désaccords.

[Read in English](README.md) · [Site OPH](https://floatingpragma.io/oph/) · [Livre](https://oph-book.floatingpragma.io/) · [Étude guidée](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [OMEGA](https://omega.floatingpragma.io/)

L’Holographie des parcelles d’observateur, ou OPH, est un programme de reconstruction de la physique fondamentale qui part de l’observateur. Il commence par des parcelles finies dotées d’un état local, de frontières, de registres, d’une relecture et de mouvements de réparation. Un fait public est un fait qui survit lorsque les parcelles qui se recouvrent comparent ce qu’elles peuvent voir et convergent vers une forme normale commune.

L’équation organisatrice est simple :

$$
T(\mathfrak U_{\mathrm{OPH}})=\mathfrak U_{\mathrm{OPH}}.
$$

L’univers est modélisé comme un point fixe de son propre processus de relecture et de réparation accessible aux observateurs.

## Pourquoi OPH est intéressant

OPH utilise une seule architecture mathématique dans des domaines habituellement introduits séparément :

- le consensus fini donne des registres publics stables et des formes normales quotientées ;
- les algèbres centrales de registres donnent les probabilités d’événements quantiques et la mise à jour conditionnelle ;
- la géométrie conforme d’un écran d’observateur donne le groupe de Lorentz connexe et un espace tridimensionnel de référentiels d’observateur ;
- le flot modulaire, le transport nul, la stationnarité de l’entropie et la géométrie des petites boules se composent conditionnellement en relation d’Einstein sur une même tour issue de la source, avec domaine commun, limites asymptotiques certifiées et identifications physiques indépendantes ; la construction et la certification de cette tour sont en cours ;
- les charges transportables et la reconstruction compacte donnent la structure de jauge du Modèle standard ;
- une construction finie à douze ports fondée sur $A_5$ produit l’algèbre de Lie
  $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$ ;
- l’équilibre des traces intègre cette algèbre en
  $S(U(3)\times U(2))\cong(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$ ;
- la branche de réalisation admissible minimale sélectionne le réseau de charges du Modèle standard, trois couleurs, trois générations et un doublet de Higgs ;
- Lean, les certificats arithmétiques exacts, les simulations et les reçus exécutables vérifient le noyau mathématique fini.

L’importance vient de la convergence. OPH n’introduit pas un mécanisme distinct pour la mesure, l’espace-temps, la gravité et la structure de jauge. Il demande quelle part de ces quatre domaines découle d’une seule exigence : des observateurs finis doivent pouvoir former un monde public stable.

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

C’est l’algèbre de Lie locale des forces de jauge du Modèle standard. Elle est obtenue à partir de la géométrie finie des coefficients plutôt que posée comme symétrie initiale.

La même construction fait apparaître deux fois, indépendamment, le nombre $24$ :

$$
m_{\mathrm{rep}}=2(8+3+1)=24,
$$

tandis que les douze ports de l’écran donnent $24$ emplacements orientés. Un compte provient de l’algèbre de jauge reconstruite ; l’autre de la géométrie de l’écran.

## Une chaîne de reconstruction

```text
parcelles auto-lectrices
        ↓
registres, comparaison des recouvrements, réparation
        ↓
formes normales quotientées publiques
        ↓
géométrie conforme de l’écran et flot modulaire
        ↓
cinématique de Lorentz, temps d’observateur, dynamique d’Einstein
        ↓
charges transportables et reconstruction compacte
        ↓
SU(3) × SU(2) × U(1) / Z6 et matière du Modèle standard
        ↓
programmes de points fixes quantitatifs et de lecture physique
```

Les hypothèses détaillées et les types de reçus sont énoncés dans les articles. La page d’accueil du dépôt est volontairement une carte du résultat positif, et non un substitut à ces énoncés de théorèmes.

## Résultats en un coup d’œil

| Résultat | Contribution d’OPH | Source principale |
| --- | --- | --- |
| Consensus fini | Réparation terminante, lecture protégée, formes normales quotientées indépendantes de l’ordonnancement et registres centraux | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| Surface d’événements quantiques | Probabilités de Born, conditionnement de Lüders et borne de Tsirelson sur la surface finie des registres | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| Relativité | $\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ et $H^3\cong\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$ | [Article compact de reconstruction](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Dynamique d’Einstein | Chaîne typée conditionnelle du transport modulaire et nul à $G_{ab}+\Lambda g_{ab}=8\pi G\langle T_{ab}\rangle$ sur une tour issue de la source et de domaine commun ; sa construction et sa certification sont en cours | [Article compact de reconstruction](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Algèbre de jauge finie $A_5$ | Construction exacte de $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$ sur douze ports | [Article compact de reconstruction](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Forme globale du Modèle standard | $S(U(3)\times U(2))$ et quotient par le centre commun $\mathbb Z_6$ | [Article compact de reconstruction](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Structure de la matière | Réseau d’hypercharges, trois couleurs, trois générations et un doublet de Higgs sur la branche de réalisation admissible minimale (MAR) | [Article compact de reconstruction](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Vérification exacte | Sous-ensemble Lean, certificats d’intervalles, reçus finis et simulations reproductibles | [`Lean/`](Lean) et [`code/`](code) |

## Choisir un parcours de lecture

| Pour découvrir... | Commencer ici |
| --- | --- |
| L’argument persuasif le plus court | [Le cas compact pour OPH](extra/compact_proof_of_oph.pdf) |
| Le centre technique | [Reconstruction de la relativité et du Modèle standard](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| La synthèse complète | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| Le mécanisme de consensus fini | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| La construction des particules | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| L’architecture de l’écran à douze ports | [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf) |
| Les preuves exécutables | [`Lean/`](Lean), [`code/`](code) et le [registre de clôture](docs/CLOSURE_LEDGER.md) |
| L’interprétation et la continuation des observateurs | [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf) |

L’[index des articles](paper/) et l’[index des suppléments](extra/) donnent la carte complète des publications.

## Preuves et éléments de vérification

Le dépôt réunit plusieurs formes complémentaires de preuve :

- des démonstrations manuscrites dans les articles TeX ;
- un sous-ensemble Lean sans `sorry` de 111 théorèmes ;
- des certificats d’intervalles et d’unicité pour les applications numériques déclarées ;
- des reçus finis pour les porteurs, la hiérarchie et les particules ;
- du code pour la géométrie, les particules, le secteur sombre et le matériel quantique ;
- des simulations dépassant un million de parcelles ;
- un registre des affirmations et un registre de clôture reliant les résultats publics aux artefacts.

## Frontière de recherche

Le noyau structurel d’OPH soutient plusieurs prolongements actifs : lecture quantitative des particules, géométrie des neutrinos, cosmologie de capacité, spectre d’écran conditionnel avec amplitude issue de la source et pente de bord, gravité sombre comme condensat de charge de réparation, transfert de Yang–Mills et systèmes matériels ou logiciels auto-lecteurs. Le paquet radial démontre la non-identifiabilité à partir d’une seule coquille et donne deux voies d’unicité distinctes : la dilatation physique de la source et la tomographie par covariances radiales croisées. Un dossier fini provenant d’une seule source et satisfaisant tous les reçus du spectre est un travail en cours.

Tous partagent la même règle de conception : tout système physique proposé doit être représenté comme une parcelle bornée avec état local, frontières, relecture, registres, réparation et dossier public de preuves.

Le [programme de vérification OPH](docs/OPH_FALSIFICATION_PROGRAM.md) est volontairement limité aux affirmations mathématiques et aux branches réalisées suffisamment mûres. Il sert d’index de vérification, pas de récit principal du dépôt.

## Guide du dépôt

- [`paper/`](paper) : articles principaux, sources TeX, PDF et métadonnées de version.
- [`extra/`](extra) : preuve compacte et suppléments mathématiques ciblés.
- [`Lean/`](Lean) : développement de théorèmes vérifiés par machine.
- [`code/`](code) : certificats, simulations, calculs de particules et expériences.
- [`book/`](book) : source du livre et PDF téléchargeable.
- [`cosmology/`](cosmology) : recherche sur le secteur sombre et la cosmologie.
- [`physics-problems/`](physics-problems) : applications ciblées et notes sur des problèmes ouverts.
- [`docs/`](docs) : registre de clôture, politique des affirmations et matériel d’audit.
- [`assets/`](assets) : diagrammes et figures publiques.

## Explorer OPH

- [Présentation de la théorie](https://floatingpragma.io/oph/theory-of-everything)
- [Simulation interactive](https://simulation.floatingpragma.io)
- [Applications et matériel OMEGA](https://omega.floatingpragma.io)
- [Livre](https://oph-book.floatingpragma.io)
- [Étude guidée](https://learn.floatingpragma.io)
- [Blog](https://blog.floatingpragma.io/)
- OPH Sage sur [Telegram](https://t.me/HoloObserverBot) et [X](https://x.com/OphSage)

## Licence

Les contenus sont publiés sous [CC BY-NC-SA 4.0](LICENSE). Le [pacte OPH d’usage ouvert et anti-brevet](PATENTS.md) maintient les idées, logiciels, méthodes, dispositifs, simulations et matériels dérivés d’OPH ouverts à l’étude, à l’implémentation, à la modification et au partage, sans monopoles privés par brevet.
