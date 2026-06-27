# NOTICE BERGE v0.7

## Sommaire

1. [Présentation](#1-présentation)
2. [Ce qui change par rapport à BERGE v0.6](#2-ce-qui-change-par-rapport-à-berge-v0.6)
3. [Installation dans QGIS](#3-installation-dans-qgis)
4. [Paramètres du formulaire](#4-paramètres-du-formulaire)
5. [Réglages de départ recommandés](#5-réglages-de-départ-recommandés)
6. [Comprendre les garde-fous](#6-comprendre-les-garde-fous)
7. [Guide de diagnostic et de réglage](#7-guide-de-diagnostic-et-de-réglage)
8. [Sorties produites](#8-sorties-produites)
9. [Reproductibilité entre campagnes](#9-reproductibilité-entre-campagnes)
10. [Limites connues](#10-limites-connues)
11. [Référence scientifique des indices](#11-référence-scientifique-des-indices)

---

## 1. Présentation

BERGE v0.7 mesure automatiquement la couverture végétale de berges, fossés et
micro-habitats de zone humide (emprise de 3×2 m à 25×15 m) à partir d'une
orthophoto drone **RGB uniquement**. Aucune bande NIR n'est requise.

### Classes produites

| Code | Classe | Couleur symbologie |
|------|--------|--------------------|
| 0 | NoData / hors emprise / eau libre | transparent |
| 1 | Sol nu / substrat homogène | brun |
| 2 | Végétation sèche / mixte / transition texturée | ocre |
| 3 | Végétation verte | vert clair |
| 4 | Végétation dense | vert foncé |

### Indicateurs exportés

| Indicateur | Classes | Usage recommandé |
|------------|---------|------------------|
| `couverture_vegetale_ecologique_pct` | 2+3+4 | Suivi écologique général |
| `couverture_vegetale_stricte_pct` | 3+4 | Végétation photosynthétiquement active |
| `couverture_vegetale_dense_pct` | 4 | Zones de forte biomasse |
| `substrat_ouvert_pct` | 1 | Proportion de sol nu |

<img width="1535" height="869" alt="image" src="https://github.com/user-attachments/assets/9144dd78-0530-4726-98ec-ea9420f55095" />

<img width="1098" height="330" alt="image" src="https://github.com/user-attachments/assets/0946bcfa-96c9-4c38-9061-ab26e6a619e0" />


### Compatibilité

- QGIS 3.36 LTR ou plus récent
- Windows 10/11 et Ubuntu 22.04 / Debian 12
- Python 3.10 + NumPy + GDAL 3.8 (fournis avec QGIS)
- **Aucun plugin externe requis**

---

## 2. Ce qui change par rapport à BERGE v0.6

### 2.1 Score spectral : 3 indices peu corrélés remplacent 4 indices redondants

BERGE v0.6 utilisait VARI, ExG, GLI et NGRDI. Ces quatre indices dérivent tous
de la même quantité de base `G − R`, avec une corrélation mutuelle supérieure
à 0,85 sur des données drone. Ils ne fournissaient qu'une seule dimension
d'information utile.

BERGE v0.7 les remplace par trois indices choisis pour leur complémentarité :

| Indice | Formule (valeurs normalisées 0–1) | Apport principal |
|--------|-----------------------------------|------------------|
| **CIVE** | `0.441·R − 0.811·G + 0.385·B` | Robustesse aux variations d'illumination (ombre, nuage) |
| **ExG** | `2G − R − B` | Sensible à la végétation verte et mixte |
| **VEG** | `G / (R^0.667 · B^0.333)` | Non-linéaire, moins saturé sur végétation dense ; complémentaire d'ExG |

**Attention** : CIVE est un indice « sol », c'est-à-dire qu'une valeur élevée
signale du sol nu. Il est automatiquement **inversé** avant d'être combiné
(`1 − CIVE_normalisé`) pour que le score final reste croissant avec la
végétation.

### 2.2 Texture sèche simplifiée

BERGE v5/v0.6 calculait `dry_texture = 0.55·var(Lum) + 0.45·var(VARI)`.
La variance de VARI est très corrélée à celle de la luminance, et diverge
là où `G + R − B ≈ 0`. BERGE v0.7 utilise uniquement la **variance locale de
luminance** sur 3 fenêtres (7, 15, 31 px). Le résultat est plus stable et
l'interprétation plus directe.

### 2.3 Richesse texturale : deux nouveaux discriminants indépendants

| Métrique | Calcul | Discriminant principal |
|----------|--------|----------------------|
| **Saturation chromatique** | `max(R,G,B) − min(R,G,B)` | Sol nu → saturation faible et homogène |
| **Entropie locale de luminance** | Shannon sur 8 niveaux, fenêtre 15 px | Sol nu → entropie basse même si rugueux |

Ces deux métriques sont indépendantes de la teinte (elles fonctionnent même
si la végétation sèche est aussi beige que le sol nu).

### 2.4 Logique micro-pixel assouplie

BERGE v0.6 imposait `g_ch > MICRO_G ET (ExG_ch > s OU NGRDI > s)`.
Le `ET` sur `g_ch` excluait la végétation sèche claire dont la chrominance
verte est modeste.

BERGE v0.7 utilise un `OU` à double seuil :

```
pixel compté comme micro-vert si :
  g_ch > MICRO_G_STRICT  (suffisant seul : pixel nettement vert)
  OU
  (g_ch > MICRO_G ET (ExG_ch > MICRO_EXG OU NGRDI > MICRO_NGRDI))
```

`MICRO_G_STRICT` (défaut 0.40) est le seuil "sans doute".
`MICRO_G` (défaut 0.34) est le seuil assoupli, valable seulement combiné.

### 2.5 Quatrième garde-fou : eau libre et pixels très sombres

Les bords inondés d'un fossé ou d'un marais produisent des pixels très
sombres et peu saturés qui peuvent perturber les statistiques. BERGE v0.7 les
détecte et les reclasse en **NoData** (classe 0) avant le calcul des %.

### 2.6 Poids spectraux internes exposés dans l'interface

BERGE v0.6 codait en dur `0.40·VARI + 0.30·ExG + 0.20·GLI + 0.10·NGRDI`.
BERGE v0.7 expose les poids `WS_CIVE`, `WS_EXG`, `WS_VEG` dans le formulaire
Processing pour permettre leur ajustement sans modifier le code.

---

## 3. Installation dans QGIS

### 3.1 Méthode recommandée : Éditeur de scripts Processing

1. [⬇️ Télécharger Berge v0-7](https://github.com/jancelin/berge/releases/download/0.7/berge_v7_vegetation_rgb_texture.py)
2. Dans QGIS, ouvrir le **Panneau Processing** (menu *Traitement → Boîte à
   outils*).
3. En haut du panneau, cliquer sur l'icône **Python** puis
   *Ouvrir l'éditeur de scripts Python…*
4. Dans l'éditeur, cliquer sur **Ouvrir un script…** et sélectionner le
   fichier `berge_v7_vegetation_rgb_texture.py`.
5. Cliquer sur **Enregistrer sous…** et placer le script dans le dossier
   des scripts utilisateur QGIS :
   - Windows : `C:\Users\<nom>\AppData\Roaming\QGIS\QGIS3\profiles\default\processing\scripts\`
   - Linux   : `~/.local/share/QGIS/QGIS3/profiles/default/processing/scripts/`
6. Cliquer sur **Exécuter dans l'éditeur** une première fois pour valider
   l'absence d'erreur de syntaxe.
7. Fermer l'éditeur. Dans la boîte à outils Processing, le groupe **BERGE**
   doit apparaître avec l'algorithme
   *BERGE v7 - RGB · CIVE+ExG+VEG · entropie · saturation · 4 garde-fous*.

> **Si le groupe BERGE n'apparaît pas** : cliquer sur l'icône de rafraîchissement
> (flèche circulaire) en haut de la boîte à outils Processing, ou fermer et
> rouvrir QGIS.

### 3.2 Méthode alternative : dépôt dans le dossier scripts

Copier directement `berge_v7_vegetation_rgb_texture.py` dans le dossier
scripts utilisateur indiqué ci-dessus, puis relancer QGIS.

### 3.3 Vérifier les dépendances

BERGE v0.7 n'utilise que des modules inclus dans QGIS :

```python
import numpy       # fourni avec QGIS
from osgeo import gdal  # fourni avec QGIS
import sqlite3     # module standard Python
import json, csv, os, math, tempfile, shutil  # modules standard Python
```

Aucun `pip install` n'est nécessaire.

---

## 4. Paramètres du formulaire

Les paramètres sont organisés en groupes dans le formulaire Processing.

### 4.1 Entrées obligatoires

| Paramètre | Description |
|-----------|-------------|
| `Orthophoto RGB` | Raster d'entrée (GeoTIFF 3 bandes RGB, 8 ou 16 bits) |
| `Polygone d'emprise` | Vecteur délimitant la zone à analyser (GeoPackage, Shapefile…) |
| `SCR` | Système de coordonnées si le raster n'en contient pas (ex. `EPSG:2154`) |
| `Dossier de sortie` | Dossier où seront créés le GeoPackage, le CSV et le JSON |
| `Nom du site` | Préfixe pour nommer les fichiers de sortie |
| `Date / code campagne` | Deuxième partie du préfixe (ex. `2025-11`) |

### 4.2 Reproductibilité

| Paramètre | Description |
|-----------|-------------|
| `metadata.json à importer` | JSON d'une campagne précédente (optionnel) |
| `Utiliser les paramètres du JSON` | Si coché, les valeurs du formulaire sont ignorées |
| `Réutiliser les bornes de normalisation` | Garantit des scores comparables entre campagnes |

### 4.3 Normalisation robuste

| Paramètre | Défaut | Description |
|-----------|--------|-------------|
| `P_LOW` | 2 | Percentile bas (écrête les valeurs aberrantes basses) |
| `P_HIGH` | 98 | Percentile haut (écrête les valeurs aberrantes hautes) |

### 4.4 Poids globaux du score BERGE (4 composantes)

La somme est automatiquement normalisée à 1.

| Paramètre | Défaut | Composante |
|-----------|--------|-----------|
| `W_SPECTRAL` | 0.35 | Score spectral composite (CIVE+ExG+VEG) |
| `W_GREEN_DENSITY` | 0.25 | Densité locale de micro-pixels verts |
| `W_DRY_TEXTURE` | 0.25 | Variance Lum multi-échelle (végétation sèche) |
| `W_RICHNESS` | 0.15 | Richesse texturale (saturation + entropie) |

### 4.5 Poids internes du score spectral (3 indices)

La somme est automatiquement normalisée à 1.

| Paramètre | Défaut | Indice |
|-----------|--------|--------|
| `WS_CIVE` | 0.40 | CIVE (robustesse illumination) — **inversé automatiquement** |
| `WS_EXG` | 0.35 | ExG (végétation verte) |
| `WS_VEG` | 0.25 | VEG (non-linéaire, complémentaire) |

### 4.6 Densité locale de micro-pixels verts

| Paramètre | Défaut | Description |
|-----------|--------|-------------|
| `DENSITY_WINDOW` | 25 | Fenêtre de comptage en pixels (impaire) |
| `MICRO_EXG` | 0.03 | Seuil ExG chromatique (mode combiné) |
| `MICRO_NGRDI` | -0.02 | Seuil NGRDI (mode combiné) |
| `MICRO_G` | 0.34 | Seuil g_ch mode souple (valable combiné avec ExG ou NGRDI) |
| `MICRO_G_STRICT` | 0.40 | Seuil g_ch mode strict (suffisant seul) |

### 4.7 Texture sèche

| Paramètre | Défaut | Description |
|-----------|--------|-------------|
| `MULTISCALE_TEXTURE` | Oui | Calcul sur plusieurs fenêtres (recommandé) |
| `TEXTURE_WINDOWS` | `7,15,31` | Fenêtres en pixels, séparées par des virgules |

### 4.8 Richesse texturale

| Paramètre | Défaut | Description |
|-----------|--------|-------------|
| `WR_SAT` | 0.45 | Poids saturation dans richesse |
| `WR_ENT` | 0.55 | Poids entropie dans richesse |
| `ENTROPY_WINDOW` | 15 | Fenêtre pour l'entropie en pixels |
| `ENTROPY_BINS` | 8 | Niveaux de quantification (4–16) |

### 4.9 Seuils de classification

| Paramètre | Défaut | Frontière |
|-----------|--------|-----------|
| `T_SOL` | 0.28 | Sol nu / végétation sèche-transition |
| `T_VEG` | 0.52 | Végétation sèche-transition / végétation verte |
| `T_DENSE` | 0.74 | Végétation verte / végétation dense |

Contrainte vérifiée : `0 ≤ T_SOL < T_VEG < T_DENSE ≤ 1`.

### 4.10 Garde-fous

| Paramètre | Défaut | Description |
|-----------|--------|-------------|
| `T_DRY_RECOVERY` | 0.48 | Texture min pour passer cl.1 → cl.2 (récupération vég. sèche) |
| `T_SOIL_HOMOGENEITY` | 0.22 | Texture max pour forcer cl.2 → cl.1 (sol homogène) |
| `T_CRACK_SOIL_SPECTRAL_MAX` | 0.32 | Score spectral max pour le garde-fou sol craquelé |
| `T_CRACK_SOIL_CHROMA_MAX` | 0.14 | Texture chromatique max (paramètre clé du garde-fou craquelé) |
| `T_CRACK_SOIL_TEXTURE_MIN` | 0.18 | Texture sèche min pour cibler le sol craquelé |
| `T_DARK_WATER_LUM` | 0.08 | Luminance max pour le garde-fou eau libre |
| `T_DARK_WATER_SAT` | 0.06 | Saturation max pour le garde-fou eau libre |

### 4.11 Options générales

| Paramètre | Défaut | Description |
|-----------|--------|-------------|
| `MIN_PIXELS` | 16 | Taille minimale des objets après SieveFilter. 0 = désactivé |
| `Pixels noirs = NoData` | Oui | Traiter RGB 0,0,0 comme hors emprise |
| `Conserver les intermédiaires` | Oui | Exporter les rasters de diagnostic dans le GeoPackage |

---

## 5. Réglages de départ recommandés

### 5.1 Point de départ général (valeurs par défaut)

Les valeurs par défaut sont calibrées pour les micro-habitats de zone humide
du littoral atlantique français (berges et fossés de marais, GSD 5–30 mm).

```text
W_SPECTRAL      = 0.35
W_GREEN_DENSITY = 0.25
W_DRY_TEXTURE   = 0.25
W_RICHNESS      = 0.15

WS_CIVE = 0.40  WS_EXG = 0.35  WS_VEG = 0.25

T_SOL   = 0.28
T_VEG   = 0.52
T_DENSE = 0.74

T_DRY_RECOVERY            = 0.48
T_SOIL_HOMOGENEITY        = 0.22
T_CRACK_SOIL_SPECTRAL_MAX = 0.32
T_CRACK_SOIL_CHROMA_MAX   = 0.14
T_CRACK_SOIL_TEXTURE_MIN  = 0.18
T_DARK_WATER_LUM          = 0.08
T_DARK_WATER_SAT          = 0.06

MICRO_G        = 0.34
MICRO_G_STRICT = 0.40
MICRO_EXG      = 0.03
MICRO_NGRDI    = -0.02
DENSITY_WINDOW = 25

TEXTURE_WINDOWS = 7,15,31
ENTROPY_WINDOW  = 15
ENTROPY_BINS    = 8
MIN_PIXELS      = 16
```

### 5.2 Si la résolution est très fine (GSD < 5 mm)

Augmenter les fenêtres de texture pour couvrir des structures physiques
réalistes (brin d'herbe = 1–5 mm, touffe = 5–50 mm).

```text
DENSITY_WINDOW  = 31
TEXTURE_WINDOWS = 9,21,41
ENTROPY_WINDOW  = 21
```

### 5.3 Si la résolution est plus grossière (GSD > 3 cm)

Réduire les fenêtres pour éviter les effets de bord.

```text
DENSITY_WINDOW  = 15
TEXTURE_WINDOWS = 5,11,21
ENTROPY_WINDOW  = 11
```

### 5.4 Si la végétation sèche est abondante et mal détectée

Donner plus de poids à la texture et réduire le seuil `T_SOL` :

```text
W_DRY_TEXTURE   = 0.35
W_GREEN_DENSITY = 0.20
T_SOL           = 0.22
T_DRY_RECOVERY  = 0.42
```

### 5.5 Si le sol nu est très structuré (argile craquelée, sédiment ridé)

Renforcer le garde-fou sol craquelé :

```text
T_CRACK_SOIL_CHROMA_MAX   = 0.16   (remonter si sol craquelé encore en cl.2)
T_CRACK_SOIL_SPECTRAL_MAX = 0.35
T_CRACK_SOIL_TEXTURE_MIN  = 0.15
```

---

## 6. Comprendre les garde-fous

BERGE v0.7 applique quatre garde-fous dans l'ordre suivant, après la
classification initiale par le score.

### GF1 — Récupération végétation sèche (cl.1 → cl.2)

**Quand :** un pixel est en classe 1 (sol nu) mais sa texture de luminance
est forte et il est entouré de micro-pixels verts.

**Condition :** `dry_texture ≥ T_DRY_RECOVERY ET green_density ≥ 0.05`

**Interpétation :** un pixel isolé classé sol nu dans une zone texturée avec
au moins 5 % de micro-vert voisin est probablement de la végétation sèche mal
détectée par le score spectral.

**Régler si** la végétation sèche est toujours en cl.1 :
→ baisser `T_DRY_RECOVERY` de 0.48 à 0.40–0.44.

### GF2 — Sol homogène (cl.2 → cl.1)

**Quand :** un pixel est en classe 2 mais sa texture est très faible et son
score spectral faible.

**Condition :** `dry_texture ≤ T_SOIL_HOMOGENEITY ET spectral < 0.35`

**Interprétation :** une surface peu texturée et spectralement neutre ne peut
pas être de la végétation sèche — c'est du sol nu (limon, sable, béton).

**Régler si** du sol nu reste en cl.2 sans être craquelé :
→ remonter `T_SOIL_HOMOGENEITY` de 0.22 à 0.26.

### GF3 — Sol craquelé (cl.≥2 → cl.1)

**Quand :** un sol craquelé a une texture de luminance forte mais
une texture chromatique faible et un score spectral faible.

**Conditions (les 3 doivent être vraies) :**
```
spectral     ≤ T_CRACK_SOIL_SPECTRAL_MAX   (peu végétalisé spectralement)
chroma_tex   ≤ T_CRACK_SOIL_CHROMA_MAX     (faible variation de chrominance)
dry_texture  ≥ T_CRACK_SOIL_TEXTURE_MIN    (forte texture de luminance)
```

**Paramètre clé :** `T_CRACK_SOIL_CHROMA_MAX` (défaut 0.14).

| Problème | Action |
|----------|--------|
| Sol craquelé encore en cl.2 | Augmenter `T_CRACK_SOIL_CHROMA_MAX` de 0.02 |
| Végétation sèche reclassée en cl.1 | Diminuer `T_CRACK_SOIL_CHROMA_MAX` de 0.02 |

### GF4 — Eau libre / pixels très sombres (cl.≥1 → NoData)

**Quand :** un pixel valide est très sombre et très peu saturé.

**Conditions :**
```
luminance  ≤ T_DARK_WATER_LUM  (pixel très sombre)
saturation ≤ T_DARK_WATER_SAT  (peu de chrominance)
```

**Interprétation :** eau libre, ombre profonde ou bâche noire.

**Régler si** des zones sombres légitimes (sols très sombres, mousse noire)
disparaissent :
→ baisser `T_DARK_WATER_LUM` à 0.05 ou désactiver avec `T_DARK_WATER_LUM = 0`.

---

## 7. Guide de diagnostic et de réglage

Activer `Conserver les rasters intermédiaires` pour accéder aux couches de
diagnostic dans le GeoPackage. Toutes les couches sont accessibles via
*Couche → Ajouter une couche → Raster…* en sélectionnant le `.gpkg`.

### 7.1 Couches de diagnostic disponibles

| Couche | Valeurs | Usage |
|--------|---------|-------|
| `spectral_index` | 0–1 | Score spectral seul. Doit être élevé sur végétation verte, bas sur sol nu. |
| `green_density` | 0–1 | Densité micro-pixels verts. Élevé sur toute végétation, bas sur sol nu pur. |
| `dry_texture` | 0–1 | Variance Lum. Élevé sur végétation sèche ET sol craquelé. |
| `saturation` | 0–1 | Saturation chromatique normalisée. |
| `entropy_local` | 0–1 | Entropie locale de luminance normalisée. |
| `richness` | 0–1 | Combinaison saturation + entropie. |
| `chroma_texture` | 0–1 | Variance ExG_ch. Clé du garde-fou sol craquelé. |
| `berge_score` | 0–1 | Score final avant classification. |
| `guard_dry_recovery` | 0/1 | Pixels passés cl.1 → cl.2 par GF1. |
| `guard_soil_homogeneity` | 0/1 | Pixels passés cl.2 → cl.1 par GF2. |
| `guard_cracked_soil` | 0/1 | Pixels forcés cl.1 par GF3 (sol craquelé). |
| `guard_dark_water` | 0/1 | Pixels exclus par GF4 (eau libre). |

### 7.2 Arbre de décision de réglage

```
La végétation sèche est en classe 1 ?
├── dry_texture élevé sur ces zones ?
│   ├── OUI → baisser T_SOL de 0.04
│   │         et/ou baisser T_DRY_RECOVERY de 0.04
│   └── NON → dry_texture trop faible
│             → baisser MICRO_G de 0.02
│             → ou augmenter W_DRY_TEXTURE de 0.05
│
Le sol nu craquelé est en classe 2 ?
├── guard_cracked_soil montre peu de pixels ?
│   ├── chroma_texture élevé sur ces zones ?
│   │   └── OUI → ce n'est pas du sol craquelé — c'est de la vég. sèche
│   └── NON → augmenter T_CRACK_SOIL_CHROMA_MAX de 0.02
│             et/ou augmenter T_CRACK_SOIL_SPECTRAL_MAX de 0.03
│
La végétation verte est en classe 2 ?
├── spectral_index faible ?
│   └── OUI → augmenter WS_EXG ou WS_VEG
│             et/ou baisser T_VEG de 0.04
└── green_density faible ?
    └── OUI → baisser MICRO_G de 0.02
              ou baisser DENSITY_WINDOW de 6

Des pixels sombres légitimes sont en NoData ?
└── Baisser T_DARK_WATER_LUM (ex. 0.05 au lieu de 0.08)
    ou T_DARK_WATER_SAT (ex. 0.04)
```

### 7.3 Interpréter le JSON de diagnostic (`diagnostic_stats`)

Le fichier `metadata.json` contient les percentiles p10/p25/p50/p75/p90 de
chaque composante sur l'ensemble des pixels valides.

Valeurs de référence attendues sur une image bien équilibrée :

| Composante | p50 attendu | Interprétation si trop bas/haut |
|------------|-------------|--------------------------------|
| `score` | 0.30–0.55 | < 0.25 → seuils trop hauts ou image très peu végétalisée |
| `spectral` | 0.30–0.55 | < 0.20 → WS_EXG ou WS_VEG insuffisants |
| `green_density` | 0.10–0.40 | < 0.05 → MICRO_G trop strict |
| `dry_texture` | 0.20–0.50 | < 0.10 → fenêtres trop petites pour la GSD |
| `richness` | 0.15–0.45 | < 0.10 → ENTROPY_BINS trop élevé ou image peu contrastée |

---

## 8. Sorties produites

Pour un site `exclos3` et la date `2025-11`, les fichiers produits sont :

```
dossier_de_sortie/
├── exclos3_2025-11_BERGE_v0.7.gpkg      ← GeoPackage principal
├── exclos3_2025-11_stats_pct.csv      ← Statistiques tabulaires
└── exclos3_2025-11_metadata.json      ← Paramètres + diagnostics
```

### Contenu du GeoPackage

**Rasters toujours présents :**

| Couche | Description |
|--------|-------------|
| `rgb_clip` | Orthophoto découpée par l'emprise |
| `classes` | Classification finale avec table de couleurs |
| `vegetation_stricte` | Masque binaire classes 3+4 |
| `vegetation_ecologique` | Masque binaire classes 2+3+4 |
| `berge_score` | Score continu BERGE v0.7 [0–1] |

**Rasters de diagnostic (si `Conserver intermédiaires = Oui`) :**
`cive`, `exg`, `veg_index`, `spectral_index`, `green_density`,
`dry_texture`, `saturation`, `entropy_local`, `richness`,
`chroma_texture`, `guard_dry_recovery`, `guard_soil_homogeneity`,
`guard_cracked_soil`, `guard_dark_water`

**Tables SQLite :**

| Table | Contenu |
|-------|---------|
| `stats_pct` | Statistiques par classe et par indicateur |
| `berge_parameters` | Tous les paramètres effectifs utilisés |
| `berge_metadata_json` | JSON complet de métadonnées (1 ligne) |
| `emprise` | Polygone d'emprise copié en vecteur |

### Charger les couches dans QGIS après traitement

```
Couche → Ajouter une couche → Raster → Parcourir → sélectionner le .gpkg
→ choisir la table dans la liste déroulante
```

Ou glisser-déposer le `.gpkg` dans le panneau de couches pour charger toutes
les tables d'un coup.

---

## 9. Reproductibilité entre campagnes

### 9.1 Principe

La normalisation percentile recalcule les bornes sur chaque image
indépendamment. Si la composition végétale change entre campagnes, les bornes
changent aussi et les scores ne sont plus comparables.

La solution est de **fixer les bornes** sur une campagne de référence.

### 9.2 Procédure

**Campagne 1 (référence) :**

1. Lancer BERGE v0.7 normalement (sans JSON de référence).
2. Conserver le fichier `*_metadata.json` produit dans le dossier de sortie.

**Campagnes suivantes :**

1. Dans le formulaire, renseigner le champ `metadata.json à importer` avec
   le JSON de la campagne de référence.
2. Cocher `Réutiliser les bornes de normalisation`.
3. Laisser décoché `Utiliser les paramètres du JSON` si on veut conserver les
   paramètres actuels du formulaire.

Les bornes de normalisation des 7 composantes
(`CIVE`, `ExG`, `VEG`, `LUM_VAR`, `SAT`, `ENTROPY`, `CHROMA_VAR`)
seront lues dans le JSON de référence.

### 9.3 Vérification dans les logs

Le log Processing indique pour chaque indice le mode de normalisation utilisé :

```
CIVE : normalisation référence [-0.123456, 0.234567]
ExG  : normalisation référence [-0.456789, 0.678901]
```

Si la mention est `locale` au lieu de `référence`, le JSON n'a pas été lu.

### 9.4 Recommandation

Conserver la campagne réalisée en conditions de végétation maximale (été /
début automne) comme référence. Elle contient la plus grande plage spectrale
et garantit que les bornes ne seront pas saturées sur les campagnes suivantes.

---

## 10. Limites connues

### 10.1 Classification pixel à pixel

BERGE v0.7 reste essentiellement une classification pixel par pixel avec un
lissage par fenêtres glissantes. Elle produit des contours irréguliers et de
petits îlots de bruit. Le paramètre `MIN_PIXELS` (SieveFilter GDAL) atténue
ce problème.

Une approche OBIA (segmentation + classification Random Forest des objets) est
envisagée pour BERGE v8 et devrait améliorer significativement la qualité des
contours et la robustesse sur les textures complexes.

### 10.2 CIVE calibré pour des DN 0–255

L'indice CIVE a été défini dans la littérature pour des images en niveaux de
gris 0–255. BERGE v0.7 le calcule sur des valeurs normalisées [0–1], ce qui
déplace l'espace des valeurs. La normalisation percentile corrige cet effet,
mais les valeurs brutes de CIVE exportées dans le GeoPackage ne correspondent
pas aux valeurs tabulées dans la littérature originale.

### 10.3 VEG et les pixels très sombres

L'indice VEG (`G / R^0.667 · B^0.333`) diverge sur les pixels sombres
(R ≈ 0 ou B ≈ 0). Un garde minimal `eps = 1e-6` est appliqué, mais les
valeurs de VEG restent instables dans les zones très sombres. Le garde-fou
eau libre (GF4) écarte la plupart de ces pixels en amont.

### 10.4 Variations d'illumination

Même si CIVE est plus robuste que VARI aux variations d'illumination, des
orthophotos acquises à des heures très différentes ou avec des conditions
nuageuses inégales peuvent produire des scores non comparables. Préférer des
acquisitions en lumière diffuse et à heure fixe pour le suivi temporel.

### 10.5 Entropie locale et temps de calcul

Le calcul de l'entropie locale sur 8 niveaux et une fenêtre de 15 px est le
plus lent de l'algorithme (8 passes box_sum). Sur une image de 5000×5000 px,
compter environ 30–60 s de plus par rapport à BERGE v0.6. Réduire
`ENTROPY_BINS` à 6 ou `ENTROPY_WINDOW` à 11 pour accélérer.

---

## 11. Référence scientifique des indices

| Indice | Référence |
|--------|-----------|
| CIVE | Kataoka T. et al. (2003). *Crop growth estimation system using convolutional neural networks*. Proceedings of the ICCV Workshop on Computer Vision for Biomedical Image Applications. |
| ExG | Woebbecke D.M. et al. (1995). *Color indices for weed identification under various soil, residue, and lighting conditions*. Transactions of the ASAE, 38(1), 259–269. |
| VEG | Hague T. et al. (2006). *Automated crop and weed monitoring in widely spaced cereals*. Precision Agriculture, 7(1), 21–32. |
| NGRDI | Tucker C.J. (1979). *Red and photographic infrared linear combinations for monitoring vegetation*. Remote Sensing of Environment, 8(2), 127–150. |
| Entropie locale de Shannon | Shannon C.E. (1948). *A mathematical theory of communication*. Bell System Technical Journal, 27(3), 379–423. |
| SieveFilter | GDAL Documentation. *gdal_sieve : removes small raster polygons*. https://gdal.org/en/stable/programs/gdal_sieve.html |

---

*BERGE v0.7 — Julien Ancelin, INRAE — 2026*
*Algorithme reproductible sous licence libre. Citer comme : Ancelin J. (2026).
BERGE v0.7 — Algorithme de suivi de la couverture végétale RGB pour les berges
et fossés de zone humide. INRAE DSLP, Saint Laurent de la Prée.*
