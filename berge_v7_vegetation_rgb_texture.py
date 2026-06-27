# -*- coding: utf-8 -*-
"""
BERGE v7 - Couverture végétale RGB · indices non redondants · texture enrichie · GeoPackage
Compatible QGIS 3.36.x / Python 3.10 / GDAL 3.8

Objectif scientifique
---------------------
Mesurer en pourcentage la couverture végétale de berges, fossés et micro-habitats de zone humide
(3×2 m à 25×15 m) à partir d'une orthophoto drone RGB uniquement. Aucune bande NIR requise.

Évolutions clés par rapport à BERGE v6
---------------------------------------
1. Remplacement du score spectral par trois indices faiblement corrélés :
   • CIVE  = 0.441·R − 0.811·G + 0.385·B + 18.787  (robuste aux variations d'illumination)
   • ExG   = 2G − R − B                              (sensibilité végétation verte et mixte)
   • VEG   = G / (R^0.667 · B^0.333)                (non-linéaire, peu saturé en dense)
   Les 4 anciens indices (VARI, ExG, GLI, NGRDI) étaient corrélés > 0.85 entre eux sur données
   drone — ils ne représentaient qu'une seule dimension d'information.

2. Simplification de la texture sèche :
   dry_texture = variance locale Lum seule, multi-échelle (7/15/31).
   La variance de VARI ajoutée en v5/v6 est très corrélée à lum_var et diverge là où G+R-B ≈ 0.
   La lum_var est plus stable et plus interprétable.

3. Ajout de deux nouvelles métriques de texture indépendantes :
   • Saturation chromatique : sat = max(R,G,B) − min(R,G,B)
     Sol nu → sat faible et homogène. Vég. sèche → sat modérée, variable. Vég. verte → sat élevée.
   • Entropie locale de luminance (8 niveaux, fenêtre 15 px) :
     H = −Σ p·log2(p). Sol nu → entropie basse. Vég. sèche → entropie élevée.
     Calculée en NumPy pur sans SciPy ni scikit-image.

4. Correction logique micro_green :
   v6 : g_ch > MICRO_G AND (exg_ch > s OR ngrdi > s)  → exclut trop de végétation sèche claire
   v7 : OR logique assoupli avec seuil séparé pour g_ch (voir MICRO_G_STRICT)

5. Troisième garde-fou : eau libre / pixels très sombres (bords de fossé inondés)
   Pixels très sombres et sans chrominance → forcés NoData pour éviter de polluer les %.

6. Diagnostic enrichi : pourcentages de pixels affectés par chaque garde-fou dans le JSON.

Classes
-------
0 = NoData / hors emprise
1 = sol nu / substrat homogène
2 = végétation sèche / mixte / transition texturée
3 = végétation verte
4 = végétation dense

Indicateur principal : couverture_vegetale_ecologique_pct (classes 2+3+4)
Indicateur strict    : couverture_vegetale_stricte_pct   (classes 3+4)

Dépendances : QGIS Processing, GDAL, NumPy, sqlite3 — aucun plugin externe.
"""

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterFile,
    QgsProcessingParameterString,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsRasterShader,
    QgsColorRampShader,
    QgsSingleBandPseudoColorRenderer,
    QgsPalettedRasterRenderer,
    QgsFillSymbol,
)
from qgis.PyQt.QtGui import QColor

import os
import csv
import json
import math
import sqlite3
import tempfile
import shutil
import numpy as np
from osgeo import gdal
import processing


class BergeV7Algorithm(QgsProcessingAlgorithm):

    # ── identifiants des paramètres ────────────────────────────────────────────
    INPUT_RGB    = "INPUT_RGB"
    MASK         = "MASK"
    SOURCE_CRS   = "SOURCE_CRS"
    OUTPUT_FOLDER= "OUTPUT_FOLDER"

    PARAM_METADATA_JSON        = "PARAM_METADATA_JSON"
    APPLY_METADATA_PARAMS      = "APPLY_METADATA_PARAMS"
    APPLY_METADATA_NORMALIZATION = "APPLY_METADATA_NORMALIZATION"

    SITE = "SITE"
    DATE = "DATE"

    P_LOW  = "P_LOW"
    P_HIGH = "P_HIGH"

    # poids des 4 composantes du score final
    W_SPECTRAL      = "W_SPECTRAL"
    W_GREEN_DENSITY = "W_GREEN_DENSITY"
    W_DRY_TEXTURE   = "W_DRY_TEXTURE"
    W_RICHNESS      = "W_RICHNESS"       # saturation + entropie

    # score spectral : poids des 3 indices (CIVE, ExG, VEG)
    WS_CIVE = "WS_CIVE"
    WS_EXG  = "WS_EXG"
    WS_VEG  = "WS_VEG"

    # densité micro-pixels verts
    DENSITY_WINDOW  = "DENSITY_WINDOW"
    MICRO_EXG       = "MICRO_EXG"
    MICRO_NGRDI     = "MICRO_NGRDI"
    MICRO_G         = "MICRO_G"
    MICRO_G_STRICT  = "MICRO_G_STRICT"

    # texture luminance multi-échelle
    MULTISCALE_TEXTURE = "MULTISCALE_TEXTURE"
    TEXTURE_WINDOWS    = "TEXTURE_WINDOWS"

    # richesse : saturation + entropie
    WR_SAT  = "WR_SAT"
    WR_ENT  = "WR_ENT"
    ENTROPY_WINDOW = "ENTROPY_WINDOW"
    ENTROPY_BINS   = "ENTROPY_BINS"

    # seuils de classification
    T_SOL            = "T_SOL"
    T_VEG            = "T_VEG"
    T_DENSE          = "T_DENSE"

    # garde-fous
    T_DRY_RECOVERY             = "T_DRY_RECOVERY"
    T_SOIL_HOMOGENEITY         = "T_SOIL_HOMOGENEITY"
    T_CRACK_SOIL_SPECTRAL_MAX  = "T_CRACK_SOIL_SPECTRAL_MAX"
    T_CRACK_SOIL_CHROMA_MAX    = "T_CRACK_SOIL_CHROMA_MAX"
    T_CRACK_SOIL_TEXTURE_MIN   = "T_CRACK_SOIL_TEXTURE_MIN"
    T_DARK_WATER_LUM           = "T_DARK_WATER_LUM"
    T_DARK_WATER_SAT           = "T_DARK_WATER_SAT"

    # options générales
    MIN_PIXELS          = "MIN_PIXELS"
    TREAT_BLACK_AS_NODATA = "TREAT_BLACK_AS_NODATA"
    KEEP_INTERMEDIATE   = "KEEP_INTERMEDIATE"

    # ──────────────────────────────────────────────────────────────────────────

    def initAlgorithm(self, config=None):

        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT_RGB, "Orthophoto RGB géoréférencée (GeoTIFF recommandé)"
        ))
        self.addParameter(QgsProcessingParameterVectorLayer(
            self.MASK, "Polygone d'emprise de l'analyse",
            [QgsProcessing.TypeVectorPolygon]
        ))
        self.addParameter(QgsProcessingParameterString(
            self.SOURCE_CRS,
            "SCR à assigner si l'image n'en contient pas (ex. EPSG:2154)",
            defaultValue="EPSG:2154"
        ))
        self.addParameter(QgsProcessingParameterFolderDestination(
            self.OUTPUT_FOLDER, "Dossier de sortie"
        ))
        self.addParameter(QgsProcessingParameterFile(
            self.PARAM_METADATA_JSON,
            "metadata.json BERGE à importer (optionnel — pour rejouer paramètres/normalisation)",
            behavior=QgsProcessingParameterFile.File,
            extension="json", optional=True, defaultValue=None
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.APPLY_METADATA_PARAMS,
            "Utiliser les paramètres du metadata.json au lieu des valeurs du formulaire",
            defaultValue=False
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.APPLY_METADATA_NORMALIZATION,
            "Réutiliser les bornes de normalisation du metadata.json",
            defaultValue=True
        ))

        self.addParameter(QgsProcessingParameterString(
            self.SITE, "Nom du site", defaultValue="berge"
        ))
        self.addParameter(QgsProcessingParameterString(
            self.DATE, "Date ou code campagne", defaultValue="2025-11"
        ))

        # Normalisation
        self.addParameter(QgsProcessingParameterNumber(
            self.P_LOW, "Percentile bas de normalisation robuste",
            QgsProcessingParameterNumber.Double,
            defaultValue=2.0, minValue=0.0, maxValue=20.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.P_HIGH, "Percentile haut de normalisation robuste",
            QgsProcessingParameterNumber.Double,
            defaultValue=98.0, minValue=80.0, maxValue=100.0
        ))

        # Poids globaux (4 composantes)
        self.addParameter(QgsProcessingParameterNumber(
            self.W_SPECTRAL, "Poids global du score spectral (CIVE+ExG+VEG)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.35, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.W_GREEN_DENSITY, "Poids densité locale de micro-pixels verts",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.25, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.W_DRY_TEXTURE, "Poids texture sèche (variance Lum multi-échelle)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.25, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.W_RICHNESS, "Poids richesse texturale (saturation + entropie)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.15, minValue=0.0, maxValue=1.0
        ))

        # Poids internes score spectral
        self.addParameter(QgsProcessingParameterNumber(
            self.WS_CIVE, "Poids spectral CIVE (robuste illumination)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.40, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.WS_EXG, "Poids spectral ExG (végétation verte/mixte)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.35, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.WS_VEG, "Poids spectral VEG (non-linéaire, peu saturé en dense)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.25, minValue=0.0, maxValue=1.0
        ))

        # Densité micro-pixels
        self.addParameter(QgsProcessingParameterNumber(
            self.DENSITY_WINDOW,
            "Fenêtre de densité locale des micro-pixels verts (pixels, impaire)",
            QgsProcessingParameterNumber.Integer,
            defaultValue=25, minValue=3, maxValue=101
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.MICRO_EXG, "Seuil micro-pixel : ExG chromatique",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.03, minValue=-1.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.MICRO_NGRDI, "Seuil micro-pixel : NGRDI",
            QgsProcessingParameterNumber.Double,
            defaultValue=-0.02, minValue=-1.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.MICRO_G,
            "Seuil micro-pixel g_ch (mode souple) — complément OU si ExG ou NGRDI dépassés",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.34, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.MICRO_G_STRICT,
            "Seuil micro-pixel g_ch (mode strict) — suffisant seul si très élevé",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.40, minValue=0.0, maxValue=1.0
        ))

        # Texture sèche
        self.addParameter(QgsProcessingParameterBoolean(
            self.MULTISCALE_TEXTURE,
            "Texture multi-échelle (7, 15, 31 px). Recommandé pour végétation sèche",
            defaultValue=True
        ))
        self.addParameter(QgsProcessingParameterString(
            self.TEXTURE_WINDOWS,
            "Fenêtres de texture si multi-échelle, séparées par virgules",
            defaultValue="7,15,31"
        ))

        # Richesse texturale
        self.addParameter(QgsProcessingParameterNumber(
            self.WR_SAT, "Poids saturation chromatique dans richesse texturale",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.45, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.WR_ENT, "Poids entropie locale dans richesse texturale",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.55, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.ENTROPY_WINDOW,
            "Fenêtre pour l'entropie locale (pixels, impaire, 9–31 recommandé)",
            QgsProcessingParameterNumber.Integer,
            defaultValue=15, minValue=5, maxValue=61
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.ENTROPY_BINS,
            "Nombre de niveaux de quantification pour l'entropie (4–16)",
            QgsProcessingParameterNumber.Integer,
            defaultValue=8, minValue=4, maxValue=16
        ))

        # Seuils de classification
        self.addParameter(QgsProcessingParameterNumber(
            self.T_SOL, "Seuil score : sol nu / végétation sèche-transition",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.28, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_VEG, "Seuil score : végétation sèche-transition / végétation verte",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.52, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_DENSE, "Seuil score : végétation verte / végétation dense",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.74, minValue=0.0, maxValue=1.0
        ))

        # Garde-fous
        self.addParameter(QgsProcessingParameterNumber(
            self.T_DRY_RECOVERY,
            "Récupération végétation sèche : seuil dry_texture pour reclasser classe 1 → 2",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.48, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_SOIL_HOMOGENEITY,
            "Garde-fou sol homogène : dry_texture max pour forcer classe 2 → 1",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.22, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_CRACK_SOIL_SPECTRAL_MAX,
            "Garde-fou sol craquelé : score spectral max pour forcer classe ≥2 → 1",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.32, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_CRACK_SOIL_CHROMA_MAX,
            "Garde-fou sol craquelé : chroma_texture max (clé principale — ajuster ±0.02)",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.14, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_CRACK_SOIL_TEXTURE_MIN,
            "Garde-fou sol craquelé : dry_texture min pour cibler les zones craquelées",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.18, minValue=0.0, maxValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_DARK_WATER_LUM,
            "Garde-fou eau libre / sombre : luminance max pour reclasser en NoData",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.08, minValue=0.0, maxValue=0.3
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.T_DARK_WATER_SAT,
            "Garde-fou eau libre / sombre : saturation max pour reclasser en NoData",
            QgsProcessingParameterNumber.Double,
            defaultValue=0.06, minValue=0.0, maxValue=0.2
        ))

        # Options générales
        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_PIXELS,
            "Nettoyage spatial : taille minimale des objets en pixels. 0 = désactivé",
            QgsProcessingParameterNumber.Integer,
            defaultValue=16, minValue=0, maxValue=10000
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.TREAT_BLACK_AS_NODATA,
            "Considérer les pixels RGB 0,0,0 comme NoData",
            defaultValue=True
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.KEEP_INTERMEDIATE,
            "Conserver les rasters intermédiaires de diagnostic dans le GeoPackage",
            defaultValue=False
        ))

    # ── utilitaires numériques ────────────────────────────────────────────────

    def _safe_prefix(self, site, date):
        p = f"{site}_{date}"
        for c in ' /\\:;,.()[]':
            p = p.replace(c, '_')
        while '__' in p:
            p = p.replace('__', '_')
        return p.strip('_') or 'berge_campagne'

    def _odd(self, k, minimum=3):
        k = int(k)
        k = max(k, minimum)
        return k if k % 2 == 1 else k + 1

    def _parse_windows(self, s):
        vals = []
        for part in str(s).replace(';', ',').split(','):
            part = part.strip()
            if part:
                try:
                    vals.append(self._odd(int(float(part))))
                except ValueError:
                    pass
        vals = sorted(set(vals))
        return vals if vals else [15]

    def _box_sum(self, arr, k):
        """Somme glissante k×k via image intégrale. Pure NumPy."""
        pad = k // 2
        arr_pad = np.pad(arr, pad, mode='reflect')
        integral = np.pad(arr_pad, ((1, 0), (1, 0)), mode='constant').cumsum(0).cumsum(1)
        return integral[k:, k:] - integral[:-k, k:] - integral[k:, :-k] + integral[:-k, :-k]

    def _local_mean(self, arr, valid, k):
        num = self._box_sum(np.where(valid, arr, 0.0).astype(np.float64), k)
        den = self._box_sum(valid.astype(np.float64), k)
        return np.divide(num, den, out=np.zeros_like(num), where=den > 0).astype(np.float32)

    def _local_variance(self, arr, valid, k):
        """E[X²] − E[X]² sur les pixels valides uniquement."""
        a = arr.astype(np.float32)
        m1 = self._local_mean(a, valid, k).astype(np.float64)
        m2 = self._local_mean((a * a).astype(np.float32), valid, k).astype(np.float64)
        return np.maximum(m2 - m1 * m1, 0.0).astype(np.float32)

    def _local_fraction(self, binary, valid, k):
        """Densité locale d'un masque binaire (fraction de pixels valides)."""
        num = self._box_sum(np.where(valid, binary, 0).astype(np.float64), k)
        den = self._box_sum(valid.astype(np.float64), k)
        return np.divide(num, den, out=np.zeros_like(num), where=den > 0).astype(np.float32)

    def _local_entropy(self, arr, valid, k, n_bins):
        """
        Entropie locale de Shannon en bits sur une fenêtre k×k.
        arr     : tableau 2D float32 normalisé [0, 1]
        valid   : masque booléen des pixels exploitables
        k       : taille de la fenêtre (pixels, impaire)
        n_bins  : nombre de niveaux de quantification (typiquement 8)

        Méthode :
        1. Quantifier arr en n_bins niveaux entiers.
        2. Pour chaque niveau b, calculer la densité locale p_b via _local_fraction.
        3. H = −Σ p_b · log2(p_b + ε), sommé sur les niveaux non vides.
        H_max = log2(n_bins) ≈ 3 bits pour 8 niveaux.
        """
        lum_q = np.clip(np.floor(arr * n_bins).astype(np.int32), 0, n_bins - 1)
        entropy = np.zeros(arr.shape, dtype=np.float32)
        eps = np.float32(1e-9)
        for b in range(n_bins):
            is_b = valid & (lum_q == b)
            p = self._local_fraction(is_b, valid, k).astype(np.float32)
            # p · log2(p) seulement là où p > eps
            nonzero = p > eps
            entropy[nonzero] -= p[nonzero] * np.log2(p[nonzero] + eps)
        # Normaliser par H_max = log2(n_bins) pour obtenir [0, 1]
        h_max = math.log2(n_bins) if n_bins > 1 else 1.0
        entropy = np.clip(entropy / h_max, 0.0, 1.0).astype(np.float32)
        entropy[~valid] = 0.0
        return entropy

    def _robust_normalize(self, arr, mask, low_pct, high_pct, ref=None):
        """Normalisation [0,1] par percentiles. ref=(lo,hi) pour mode référence."""
        out = np.full(arr.shape, -9999.0, dtype=np.float32)
        if ref is not None:
            lo, hi = float(ref[0]), float(ref[1])
        else:
            vals = arr[mask]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                out[mask] = 0.0
                return out, 0.0, 1.0
            lo = float(np.nanpercentile(vals, low_pct))
            hi = float(np.nanpercentile(vals, high_pct))
        if not np.isfinite(lo) or not np.isfinite(hi) or abs(hi - lo) < 1e-12:
            out[mask] = 0.0
        else:
            out[mask] = np.clip((arr[mask] - lo) / (hi - lo), 0.0, 1.0)
        return out, lo, hi

    def _class_color_table(self):
        ct = gdal.ColorTable()
        ct.SetColorEntry(0, (0,   0,   0,   0))   # NoData – transparent
        ct.SetColorEntry(1, (180, 110,  60, 255))  # sol nu – brun
        ct.SetColorEntry(2, (210, 190,  85, 255))  # vég. sèche – ocre
        ct.SetColorEntry(3, ( 70, 185,  95, 255))  # vég. verte – vert clair
        ct.SetColorEntry(4, (  0, 120,  40, 255))  # vég. dense – vert foncé
        return ct

    def _write_tif(self, path, array, dtype, nodata, gt, proj, color_table=None):
        driver = gdal.GetDriverByName('GTiff')
        ysize, xsize = array.shape
        ds = driver.Create(path, xsize, ysize, 1, dtype,
                           options=['COMPRESS=LZW', 'TILED=YES'])
        if ds is None:
            raise QgsProcessingException('Impossible de créer : ' + path)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        b = ds.GetRasterBand(1)
        b.WriteArray(array)
        b.SetNoDataValue(nodata)
        if color_table:
            b.SetRasterColorTable(color_table)
            b.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
        try:
            b.ComputeStatistics(False)
        except Exception:
            pass
        b.FlushCache()
        ds.FlushCache()
        ds = None

    # ── gestion metadata JSON ─────────────────────────────────────────────────

    def _load_metadata(self, path, feedback):
        if not path or not os.path.isfile(path):
            return {}
        try:
            with open(path, 'r', encoding='utf-8') as f:
                meta = json.load(f)
            feedback.pushInfo('metadata.json chargé : ' + os.path.basename(path))
            return meta
        except Exception as exc:
            feedback.reportError(
                f'metadata.json illisible : {exc}. Paramètres du formulaire utilisés.',
                fatalError=False
            )
            return {}

    def _get_meta_param(self, meta, key, default):
        """Lecture robuste compatible BERGE v5/v6/v7."""
        params = meta.get('parameters', {}) if isinstance(meta, dict) else {}
        if key in params:
            return params[key]
        # Compatibilité v6 : normalization_bounds → bornes, pas paramètres
        return default

    def _normalization_bounds_from_metadata(self, meta):
        if not isinstance(meta, dict):
            return {}
        if 'normalization_bounds' in meta:
            return meta['normalization_bounds']
        # Compat v4 : index_percentiles
        if 'index_percentiles' in meta:
            b = {}
            for k, v in meta['index_percentiles'].items():
                if 'p_low' in v and 'p_high' in v:
                    b[k] = {'low': v['p_low'], 'high': v['p_high']}
            return b
        return {}

    def _ref_pair(self, bounds, name):
        if name in bounds:
            bv = bounds[name]
            lo = bv.get('low', bv.get('p_low'))
            hi = bv.get('high', bv.get('p_high'))
            if lo is not None and hi is not None:
                return (lo, hi)
        return None

    # ── GeoPackage ────────────────────────────────────────────────────────────

    def _add_raster_to_gpkg(self, tif_path, gpkg_path, table_name, first, feedback):
        co = [f'RASTER_TABLE={table_name}', 'TILE_FORMAT=PNG_JPEG', 'BLOCKSIZE=256']
        if not first:
            co.append('APPEND_SUBDATASET=YES')
        ds = gdal.Translate(gpkg_path, tif_path, format='GPKG', creationOptions=co)
        if ds is None:
            feedback.reportError(
                f'Export raster GeoPackage non réalisé pour : {table_name}', fatalError=False
            )
        else:
            ds = None
            feedback.pushInfo(f'  raster GeoPackage : {table_name}')

    def _copy_vector_to_gpkg(self, vector_source, gpkg_path, layer_name, feedback):
        try:
            src_path = vector_source.split('|')[0]
            opts = gdal.VectorTranslateOptions(
                format='GPKG', accessMode='update', layerName=layer_name,
                layerCreationOptions=['OVERWRITE=YES']
            )
            out = gdal.VectorTranslate(gpkg_path, src_path, options=opts)
            if out is None:
                feedback.reportError('Copie du polygone d\'emprise non réalisée.', fatalError=False)
            else:
                out = None
                feedback.pushInfo(f'  vecteur GeoPackage : {layer_name}')
        except Exception as exc:
            feedback.reportError(f'Copie du polygone non réalisée : {exc}', fatalError=False)

    def _write_tables_to_gpkg(self, gpkg_path, rows, metadata, parameters, feedback):
        conn = sqlite3.connect(gpkg_path)
        cur = conn.cursor()

        cur.execute('DROP TABLE IF EXISTS stats_pct')
        cur.execute('''CREATE TABLE stats_pct (
            site TEXT, date TEXT, algorithme TEXT, type TEXT,
            code INTEGER, classe_ou_indicateur TEXT, pixels INTEGER, pourcentage REAL
        )''')
        cur.executemany('INSERT INTO stats_pct VALUES (?,?,?,?,?,?,?,?)', rows)

        cur.execute('DROP TABLE IF EXISTS berge_parameters')
        cur.execute('CREATE TABLE berge_parameters (key TEXT PRIMARY KEY, value TEXT)')
        cur.executemany(
            'INSERT INTO berge_parameters(key, value) VALUES (?, ?)',
            [(k, json.dumps(v, ensure_ascii=False)) for k, v in sorted(parameters.items())]
        )

        cur.execute('DROP TABLE IF EXISTS berge_metadata_json')
        cur.execute('CREATE TABLE berge_metadata_json (id INTEGER PRIMARY KEY, metadata TEXT)')
        cur.execute('INSERT INTO berge_metadata_json(id, metadata) VALUES (1, ?)',
                    (json.dumps(metadata, ensure_ascii=False, indent=2),))

        try:
            for table in ('stats_pct', 'berge_parameters', 'berge_metadata_json'):
                cur.execute('DELETE FROM gpkg_contents WHERE table_name=?', (table,))
                cur.execute(
                    "INSERT INTO gpkg_contents(table_name, data_type, identifier, description, last_change) "
                    "VALUES (?, 'attributes', ?, ?, strftime('%Y-%m-%dT%H:%M:%fZ','now'))",
                    (table, table, f'Table BERGE v7 : {table}')
                )
        except Exception:
            pass

        conn.commit()
        conn.close()
        feedback.pushInfo('  tables GeoPackage : stats_pct, berge_parameters, berge_metadata_json')

    # ── Styles QGIS intégrés au GeoPackage ────────────────────────────────────

    def _berge_continuous_colors(self):
        """Rampe issue du style QGIS manuel BERGE : violet → bleu → cyan → vert → jaune → rouge."""
        return [
            '#30123b', '#4147ad', '#4777ef', '#38a5fb', '#1bd0d5',
            '#26eda6', '#64fd6a', '#a4fc3c', '#d3e835', '#f5c63a',
            '#fe992c', '#f36315', '#d93807', '#b01901', '#7a0403',
        ]

    def _format_style_label(self, value):
        try:
            return f'{float(value):.4f}'.replace('.', ',')
        except Exception:
            return str(value)

    def _quantile_breaks(self, array, nodata=None, n_classes=15):
        """
        Calcule les ruptures de style comme QGIS en mode Quantile.
        Les NoData, NaN et infinis sont exclus.
        """
        vals = np.asarray(array, dtype=np.float64).ravel()
        vals = vals[np.isfinite(vals)]
        if nodata is not None:
            vals = vals[vals != float(nodata)]
        # Exclure explicitement la sentinelle BERGE des rasters Float32.
        vals = vals[vals != -9999.0]
        if vals.size == 0:
            return None

        qs = np.linspace(0.0, 100.0, int(n_classes))
        breaks = np.nanpercentile(vals, qs).astype(float)

        # QGIS tolère mal les listes où beaucoup de quantiles sont identiques
        # (cas des rasters normalisés saturés à 0/1). On force une croissance
        # stricte très faible sans modifier visuellement la classification.
        vmin = float(np.nanmin(vals))
        vmax = float(np.nanmax(vals))
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            return None
        if abs(vmax - vmin) < 1e-12:
            vmax = vmin + 1.0
            breaks = np.linspace(vmin, vmax, int(n_classes))
        else:
            eps = max(abs(vmax - vmin) * 1e-9, 1e-9)
            for i in range(1, len(breaks)):
                if breaks[i] <= breaks[i - 1]:
                    breaks[i] = breaks[i - 1] + eps
        return breaks

    def _set_continuous_raster_style(self, layer, array, nodata, opacity=0.60, n_classes=15):
        breaks = self._quantile_breaks(array, nodata=nodata, n_classes=n_classes)
        if breaks is None:
            return False

        colors = self._berge_continuous_colors()
        if len(colors) != len(breaks):
            raise QgsProcessingException('Nombre de couleurs BERGE incohérent avec les classes de style.')

        ramp_items = []
        for value, color in zip(breaks, colors):
            ramp_items.append(
                QgsColorRampShader.ColorRampItem(
                    float(value), QColor(color), self._format_style_label(value)
                )
            )

        color_shader = QgsColorRampShader()
        color_shader.setColorRampType(QgsColorRampShader.Interpolated)
        try:
            # 3 = Quantile dans le XML QGIS ; l'enum existe dans QGIS 3.36.
            color_shader.setClassificationMode(QgsColorRampShader.Quantile)
        except Exception:
            pass
        color_shader.setColorRampItemList(ramp_items)

        raster_shader = QgsRasterShader()
        raster_shader.setRasterShaderFunction(color_shader)

        renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, raster_shader)
        layer.setRenderer(renderer)
        layer.setOpacity(float(opacity))
        return True

    def _set_classes_style(self, layer):
        entries = [
            QgsPalettedRasterRenderer.Class(0, QColor(0, 0, 0, 0), 'NoData / hors emprise / eau'),
            QgsPalettedRasterRenderer.Class(1, QColor('#b46e3c'), '1 · sol nu / substrat homogène'),
            QgsPalettedRasterRenderer.Class(2, QColor('#d2be55'), '2 · végétation sèche / mixte'),
            QgsPalettedRasterRenderer.Class(3, QColor('#46b95f'), '3 · végétation verte'),
            QgsPalettedRasterRenderer.Class(4, QColor('#007828'), '4 · végétation dense'),
        ]
        layer.setRenderer(QgsPalettedRasterRenderer(layer.dataProvider(), 1, entries))
        layer.setOpacity(0.60)
        return True

    def _set_binary_vegetation_style(self, layer, label):
        entries = [
            QgsPalettedRasterRenderer.Class(0, QColor(0, 0, 0, 0), '0 · hors indicateur'),
            QgsPalettedRasterRenderer.Class(1, QColor('#f7fcf5'), label),
        ]
        layer.setRenderer(QgsPalettedRasterRenderer(layer.dataProvider(), 1, entries))
        layer.setOpacity(0.60)
        return True

    def _set_binary_guard_style(self, layer, label):
        entries = [
            QgsPalettedRasterRenderer.Class(0, QColor(0, 0, 0, 0), '0 · non activé'),
            QgsPalettedRasterRenderer.Class(1, QColor('#d7301f'), label),
        ]
        layer.setRenderer(QgsPalettedRasterRenderer(layer.dataProvider(), 1, entries))
        layer.setOpacity(0.70)
        return True

    def _set_emprise_style(self, layer):
        symbol = QgsFillSymbol.createSimple({
            'style': 'no',
            'outline_style': 'solid',
            'outline_color': '35,35,35,255',
            'outline_width': '0.35',
            'outline_width_unit': 'MM',
        })
        layer.renderer().setSymbol(symbol)
        return True

    def _gpkg_style_exists(self, gpkg_path, table_name, style_name):
        try:
            conn = sqlite3.connect(gpkg_path)
            cur = conn.cursor()
            cur.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='layer_styles'"
            )
            if cur.fetchone() is None:
                conn.close()
                return False
            cur.execute(
                "SELECT COUNT(*) FROM layer_styles WHERE f_table_name=? AND styleName=?",
                (table_name, style_name)
            )
            count = int(cur.fetchone()[0])
            conn.close()
            return count > 0
        except Exception:
            return False

    def _save_layer_style_to_gpkg(self, layer, gpkg_path, table_name, feedback):
        """
        Enregistre le style dans la table QGIS `layer_styles` du GeoPackage.
        Méthode principale : API QGIS. Fallback : insertion SQLite du QML généré.
        """
        style_name = 'BERGE défaut'
        description = f'Style par défaut généré automatiquement par BERGE v7 pour {table_name}'

        try:
            res = layer.saveStyleToDatabase(style_name, description, True, '')
            # Selon les versions QGIS, le retour peut être None, bool, ou tuple.
            api_ok = False
            if isinstance(res, tuple):
                api_ok = any(isinstance(x, bool) and x for x in res)
                if len(res) >= 2 and res[1] is True:
                    api_ok = True
            elif res is True or res is None:
                api_ok = True
            if api_ok and self._gpkg_style_exists(gpkg_path, table_name, style_name):
                return True
        except Exception as exc:
            feedback.reportError(
                f'API QGIS saveStyleToDatabase indisponible pour {table_name} : {exc}. '
                'Tentative fallback SQLite.', fatalError=False
            )

        tmp_qml = os.path.join(tempfile.gettempdir(), f'berge_style_{table_name}.qml')
        try:
            layer.saveNamedStyle(tmp_qml)
            if not os.path.exists(tmp_qml):
                return False
            with open(tmp_qml, 'r', encoding='utf-8') as f:
                qml = f.read()

            conn = sqlite3.connect(gpkg_path)
            cur = conn.cursor()
            cur.execute("""CREATE TABLE IF NOT EXISTS layer_styles (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                f_table_catalog TEXT,
                f_table_schema TEXT,
                f_table_name TEXT,
                f_geometry_column TEXT,
                styleName TEXT,
                styleQML TEXT,
                styleSLD TEXT,
                useAsDefault INTEGER,
                description TEXT,
                owner TEXT,
                ui TEXT,
                update_time TEXT
            )""")
            cur.execute(
                "DELETE FROM layer_styles WHERE f_table_name=? AND styleName=?",
                (table_name, style_name)
            )
            cur.execute("""INSERT INTO layer_styles(
                f_table_catalog, f_table_schema, f_table_name, f_geometry_column,
                styleName, styleQML, styleSLD, useAsDefault, description, owner, ui, update_time
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, strftime('%Y-%m-%dT%H:%M:%fZ','now'))""",
                ('', '', table_name, '', style_name, qml, '', 1, description, 'BERGE', '')
            )
            conn.commit()
            conn.close()
            return True
        except Exception as exc:
            feedback.reportError(f'Fallback SQLite style échoué pour {table_name} : {exc}', fatalError=False)
            return False

    def _apply_gpkg_styles(self, gpkg_path, raster_style_specs, feedback):
        """Applique et enregistre les styles QGIS par défaut pour les couches du GeoPackage."""
        feedback.pushInfo('Écriture des styles QGIS intégrés au GeoPackage...')

        ok_count = 0
        for table_name, spec in raster_style_specs.items():
            uri = f'GPKG:{gpkg_path}:{table_name}'
            layer = QgsRasterLayer(uri, table_name, 'gdal')
            if not layer.isValid():
                feedback.reportError(f'Style ignoré : couche raster invalide {table_name}', fatalError=False)
                continue

            kind = spec.get('kind')
            styled = False
            if kind == 'continuous':
                styled = self._set_continuous_raster_style(
                    layer,
                    spec.get('array'),
                    spec.get('nodata', -9999),
                    opacity=spec.get('opacity', 0.60),
                    n_classes=spec.get('classes', 15),
                )
            elif kind == 'classes':
                styled = self._set_classes_style(layer)
            elif kind == 'binary_vegetation':
                styled = self._set_binary_vegetation_style(layer, spec.get('label', '1 · végétation'))
            elif kind == 'binary_guard':
                styled = self._set_binary_guard_style(layer, spec.get('label', '1 · garde-fou activé'))
            elif kind == 'rgb':
                # Le fournisseur GDAL ouvre déjà le RGB en multibande couleur ; on enregistre le style natif.
                layer.setOpacity(1.0)
                styled = True

            if styled and self._save_layer_style_to_gpkg(layer, gpkg_path, table_name, feedback):
                ok_count += 1
                feedback.pushInfo(f'  style GeoPackage : {table_name}')
            else:
                feedback.reportError(f'Style non enregistré pour : {table_name}', fatalError=False)

        # Style du vecteur d'emprise.
        v_uri = f'{gpkg_path}|layername=emprise'
        v_layer = QgsVectorLayer(v_uri, 'emprise', 'ogr')
        if v_layer.isValid():
            try:
                self._set_emprise_style(v_layer)
                if self._save_layer_style_to_gpkg(v_layer, gpkg_path, 'emprise', feedback):
                    ok_count += 1
                    feedback.pushInfo('  style GeoPackage : emprise')
            except Exception as exc:
                feedback.reportError(f'Style emprise non enregistré : {exc}', fatalError=False)

        feedback.pushInfo(f'  styles enregistrés : {ok_count}')

    # ── processAlgorithm ──────────────────────────────────────────────────────

    def processAlgorithm(self, parameters, context, feedback):

        rgb_layer  = self.parameterAsRasterLayer(parameters, self.INPUT_RGB, context)
        mask_layer = self.parameterAsVectorLayer(parameters, self.MASK, context)
        source_crs = self.parameterAsString(parameters, self.SOURCE_CRS, context).strip()
        out_dir    = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        os.makedirs(out_dir, exist_ok=True)

        # Metadata JSON import
        meta_path_in  = self.parameterAsFile(parameters, self.PARAM_METADATA_JSON, context)
        apply_params  = self.parameterAsBool(parameters, self.APPLY_METADATA_PARAMS, context)
        apply_norm    = self.parameterAsBool(parameters, self.APPLY_METADATA_NORMALIZATION, context)
        imported_meta = self._load_metadata(meta_path_in, feedback)
        norm_bounds   = self._normalization_bounds_from_metadata(imported_meta) if apply_norm else {}

        def pf(key):
            val = self.parameterAsDouble(parameters, key, context)
            return float(self._get_meta_param(imported_meta, key, val)) if apply_params else val
        def pi(key):
            val = self.parameterAsInt(parameters, key, context)
            return int(float(self._get_meta_param(imported_meta, key, val))) if apply_params else val
        def pb(key):
            val = self.parameterAsBool(parameters, key, context)
            if apply_params:
                mv = self._get_meta_param(imported_meta, key, val)
                if isinstance(mv, str):
                    return mv.lower() in ('1', 'true', 'yes', 'oui')
                return bool(mv)
            return val
        def ps(key, default=''):
            val = self.parameterAsString(parameters, key, context)
            return str(self._get_meta_param(imported_meta, key, val)) if apply_params else val

        site   = ps(self.SITE, 'berge').strip() or 'berge'
        date   = ps(self.DATE, 'campagne').strip() or 'campagne'
        prefix = self._safe_prefix(site, date)

        p_low  = pf(self.P_LOW)
        p_high = pf(self.P_HIGH)
        if not (0.0 <= p_low < p_high <= 100.0):
            raise QgsProcessingException('Percentiles invalides : vérifier P_LOW < P_HIGH.')

        # Poids globaux
        w_spectral     = pf(self.W_SPECTRAL)
        w_green_dens   = pf(self.W_GREEN_DENSITY)
        w_dry_tex      = pf(self.W_DRY_TEXTURE)
        w_richness     = pf(self.W_RICHNESS)
        w_sum = w_spectral + w_green_dens + w_dry_tex + w_richness
        if w_sum <= 0:
            raise QgsProcessingException('La somme des poids globaux doit être > 0.')
        w_spectral   /= w_sum
        w_green_dens /= w_sum
        w_dry_tex    /= w_sum
        w_richness   /= w_sum

        # Poids internes score spectral
        ws_cive = pf(self.WS_CIVE)
        ws_exg  = pf(self.WS_EXG)
        ws_veg  = pf(self.WS_VEG)
        ws_sum  = ws_cive + ws_exg + ws_veg
        if ws_sum <= 0:
            raise QgsProcessingException('La somme des poids spectraux (CIVE+ExG+VEG) doit être > 0.')
        ws_cive /= ws_sum
        ws_exg  /= ws_sum
        ws_veg  /= ws_sum

        # Richesse texturale
        wr_sat = pf(self.WR_SAT)
        wr_ent = pf(self.WR_ENT)
        wr_sum = wr_sat + wr_ent
        if wr_sum <= 0:
            wr_sat, wr_ent = 0.45, 0.55
        else:
            wr_sat /= wr_sum
            wr_ent /= wr_sum

        entropy_win  = self._odd(pi(self.ENTROPY_WINDOW))
        entropy_bins = max(4, min(16, pi(self.ENTROPY_BINS)))

        density_win  = self._odd(pi(self.DENSITY_WINDOW))
        micro_exg    = pf(self.MICRO_EXG)
        micro_ngrdi  = pf(self.MICRO_NGRDI)
        micro_g      = pf(self.MICRO_G)
        micro_g_strict = pf(self.MICRO_G_STRICT)

        multi   = pb(self.MULTISCALE_TEXTURE)
        windows = self._parse_windows(ps(self.TEXTURE_WINDOWS, '7,15,31'))
        if not multi:
            windows = [15]

        t_sol   = pf(self.T_SOL)
        t_veg   = pf(self.T_VEG)
        t_dense = pf(self.T_DENSE)
        if not (0 <= t_sol < t_veg < t_dense <= 1):
            raise QgsProcessingException('Seuils invalides : 0 ≤ T_SOL < T_VEG < T_DENSE ≤ 1.')

        t_dry_rec      = pf(self.T_DRY_RECOVERY)
        t_soil_hom     = pf(self.T_SOIL_HOMOGENEITY)
        t_crack_spec   = pf(self.T_CRACK_SOIL_SPECTRAL_MAX)
        t_crack_chroma = pf(self.T_CRACK_SOIL_CHROMA_MAX)
        t_crack_tex    = pf(self.T_CRACK_SOIL_TEXTURE_MIN)
        t_dark_lum     = pf(self.T_DARK_WATER_LUM)
        t_dark_sat     = pf(self.T_DARK_WATER_SAT)

        min_pixels      = pi(self.MIN_PIXELS)
        black_as_nodata = pb(self.TREAT_BLACK_AS_NODATA)
        keep_intermed   = pb(self.KEEP_INTERMEDIATE)

        gpkg_out     = os.path.join(out_dir, f'{prefix}_BERGE_v7.gpkg')
        csv_out      = os.path.join(out_dir, f'{prefix}_stats_pct.csv')
        metadata_out = os.path.join(out_dir, f'{prefix}_metadata.json')
        if os.path.exists(gpkg_out):
            os.remove(gpkg_out)

        tmp_dir = tempfile.mkdtemp(prefix='berge_v7_')
        feedback.pushInfo('Dossier temporaire : ' + tmp_dir)

        try:
            rgb_clip      = os.path.join(tmp_dir, f'{prefix}_rgb_clip.tif')
            classes_raw   = os.path.join(tmp_dir, f'{prefix}_classes_raw.tif')
            classes_clean = os.path.join(tmp_dir, f'{prefix}_classes.tif')

            # ── 1. Découpage ────────────────────────────────────────────────
            feedback.pushInfo('BERGE v7 : découpage de l\'orthophoto par l\'emprise...')
            try:
                processing.run('gdal:cliprasterbymasklayer', {
                    'INPUT':            rgb_layer.source(),
                    'MASK':             mask_layer.source(),
                    'SOURCE_CRS':       source_crs or None,
                    'TARGET_CRS':       source_crs or None,
                    'NODATA':           0,
                    'ALPHA_BAND':       False,
                    'CROP_TO_CUTLINE':  True,
                    'KEEP_RESOLUTION':  True,
                    'SET_RESOLUTION':   False,
                    'X_RESOLUTION':     None,
                    'Y_RESOLUTION':     None,
                    'MULTITHREADING':   True,
                    'OPTIONS':          '',
                    'DATA_TYPE':        0,
                    'EXTRA':            '',
                    'OUTPUT':           rgb_clip,
                }, context=context, feedback=feedback)
            except Exception as exc:
                raise QgsProcessingException(
                    'Échec du découpage. Vérifier que le raster est géoréférencé, que SOURCE_CRS '
                    'est correct et que le polygone d\'emprise est dans le même SCR. Détail : ' + str(exc)
                )

            # ── 2. Lecture RGB ───────────────────────────────────────────────
            ds = gdal.Open(rgb_clip)
            if ds is None:
                raise QgsProcessingException('Impossible d\'ouvrir le raster RGB découpé.')
            if ds.RasterCount < 3:
                raise QgsProcessingException('Le raster doit avoir au moins 3 bandes RGB.')

            gt    = ds.GetGeoTransform()
            proj  = ds.GetProjection()
            xsize = ds.RasterXSize
            ysize = ds.RasterYSize
            bands   = [ds.GetRasterBand(i) for i in (1, 2, 3)]
            R0, G0, B0 = [b.ReadAsArray().astype(np.float32) for b in bands]
            nodatas = [b.GetNoDataValue() for b in bands]
            ds = None

            nodata_mask = np.zeros((ysize, xsize), dtype=bool)
            for arr, nd in zip((R0, G0, B0), nodatas):
                if nd is not None:
                    nodata_mask |= (arr == float(nd))
            if black_as_nodata:
                nodata_mask |= ((R0 == 0) & (G0 == 0) & (B0 == 0))

            max_val = float(np.nanmax([np.nanmax(R0), np.nanmax(G0), np.nanmax(B0)]))
            if not np.isfinite(max_val) or max_val <= 0:
                raise QgsProcessingException('Valeurs RGB invalides : maximum ≤ 0.')
            scale = 255.0 if max_val <= 255.0 else (65535.0 if max_val <= 65535.0 else max_val)

            R = (R0 / scale).astype(np.float32)
            G = (G0 / scale).astype(np.float32)
            B = (B0 / scale).astype(np.float32)
            del R0, G0, B0

            eps   = np.float32(1e-6)
            valid = (~nodata_mask) & np.isfinite(R) & np.isfinite(G) & np.isfinite(B)
            if int(np.sum(valid)) < 100:
                raise QgsProcessingException(
                    'Moins de 100 pixels valides. Vérifier l\'emprise, le NoData et le SCR.'
                )

            sum_rgb = R + G + B + eps
            r_ch = (R / sum_rgb).astype(np.float32)
            g_ch = (G / sum_rgb).astype(np.float32)
            b_ch = (B / sum_rgb).astype(np.float32)

            # ── 3. Indices spectraux (CIVE, ExG, VEG) ───────────────────────
            feedback.pushInfo('Calcul des 3 indices spectraux (CIVE, ExG, VEG)...')
            NODATA_F = np.float32(-9999.0)

            # CIVE : Colour Index of Vegetation Extraction
            # CIVE = 0.441·R − 0.811·G + 0.385·B + 18.787
            # Défini sur DN 0-255 → adapter à [0,1] en retirant le décalage constant
            # puis renormaliser. On travaille en valeurs normalisées.
            cive = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            cive[valid] = (0.441 * R[valid] - 0.811 * G[valid] + 0.385 * B[valid])
            # Note : le terme +18.787 est un biais constant sur DN 0-255.
            # Sur valeurs [0,1] il disparaît dans la normalisation percentile.

            # ExG : Excess Green
            exg = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            exg[valid] = (2.0 * G[valid] - R[valid] - B[valid])

            # VEG : Vegetative Index (Hague 2006)
            # VEG = G / (R^a · B^(1-a))  avec a = 0.667
            # Éviter division par zéro
            veg = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            r_pow = np.where(valid & (R > eps), np.power(R, 0.667, where=R > 0), eps)
            b_pow = np.where(valid & (B > eps), np.power(B, 0.333, where=B > 0), eps)
            denom_veg = r_pow * b_pow
            m_veg = valid & (denom_veg > eps) & (G > eps)
            veg[m_veg] = G[m_veg] / denom_veg[m_veg]

            # NGRDI pour le critère micro-pixel (conservé en interne)
            ngrdi = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            denom_ngrdi = G + R
            m_ngrdi = valid & (np.abs(denom_ngrdi) > eps)
            ngrdi[m_ngrdi] = (G[m_ngrdi] - R[m_ngrdi]) / denom_ngrdi[m_ngrdi]

            index_valid = valid & (cive != NODATA_F) & (exg != NODATA_F) & \
                          (veg != NODATA_F) & (ngrdi != NODATA_F)
            n_valid = int(np.sum(index_valid))
            if n_valid < 100:
                raise QgsProcessingException('Indices invalides : moins de 100 pixels exploitables.')
            feedback.pushInfo(f'Pixels valides : {n_valid:,}')

            # ── 4. Normalisation ─────────────────────────────────────────────
            feedback.pushInfo('Normalisation des indices...')

            def norm(name, arr):
                ref = self._ref_pair(norm_bounds, name)
                n, lo, hi = self._robust_normalize(arr, index_valid, p_low, p_high, ref=ref)
                mode = 'référence' if ref else 'locale'
                feedback.pushInfo(f'  {name} : normalisation {mode} [{lo:.6f}, {hi:.6f}]')
                return n, lo, hi

            cive_n, cive_lo, cive_hi = norm('CIVE', cive)
            exg_n,  exg_lo,  exg_hi  = norm('ExG',  exg)
            veg_n,  veg_lo,  veg_hi  = norm('VEG',  veg)

            # Score spectral composite : CIVE est inversé (plus négatif = plus vert)
            # Après normalisation, CIVE_n élevé = peu végétalisé → inverser
            spectral = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            spectral[index_valid] = np.clip(
                ws_cive * (1.0 - cive_n[index_valid]) +   # CIVE inversé
                ws_exg  *        exg_n[index_valid]   +
                ws_veg  *        veg_n[index_valid],
                0.0, 1.0
            )

            # ── 5. Densité locale de micro-pixels verts ──────────────────────
            feedback.pushInfo('Calcul de la densité locale de micro-pixels verts...')
            exg_ch = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            exg_ch[valid] = (2.0 * g_ch[valid] - r_ch[valid] - b_ch[valid])

            # Logique v7 (assouplie) :
            # Un pixel est compté micro-vert si :
            #   g_ch > micro_g_strict (très vert, condition suffisante)
            #   OU (g_ch > micro_g ET (exg_ch > seuil OU ngrdi > seuil))
            micro_green = index_valid & (
                (g_ch > micro_g_strict) |
                (
                    (g_ch > micro_g) &
                    ((exg_ch > micro_exg) | (ngrdi > micro_ngrdi))
                )
            )
            green_density = self._local_fraction(micro_green, index_valid, density_win)

            # ── 6. Texture sèche : variance Lum multi-échelle ────────────────
            feedback.pushInfo('Calcul de la texture sèche (variance Lum multi-échelle)...')
            lum = (0.299 * R + 0.587 * G + 0.114 * B).astype(np.float32)
            lum[~index_valid] = 0.0

            lum_vars = []
            for k in windows:
                feedback.pushInfo(f'  variance Lum fenêtre {k}×{k}')
                lum_vars.append(self._local_variance(lum, index_valid, k))
            lum_var = np.maximum.reduce(lum_vars)
            del lum_vars

            lum_var_n, lumv_lo, lumv_hi = norm('LUM_VAR', lum_var)
            dry_texture = lum_var_n  # alias sémantique

            # ── 7. Richesse texturale : saturation + entropie ────────────────
            feedback.pushInfo('Calcul de la richesse texturale (saturation + entropie)...')

            # 7a. Saturation chromatique : sat = max(R,G,B) − min(R,G,B)
            sat = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            sat[valid] = (
                np.maximum(np.maximum(R[valid], G[valid]), B[valid]) -
                np.minimum(np.minimum(R[valid], G[valid]), B[valid])
            )
            sat_n, sat_lo, sat_hi = norm('SAT', sat)

            # 7b. Entropie locale de Shannon de la luminance
            lum_norm = np.clip(lum, 0.0, 1.0).astype(np.float32)   # déjà [0,1]
            feedback.pushInfo(f'  entropie locale Lum (fenêtre {entropy_win}×{entropy_win}, {entropy_bins} niveaux)...')
            entropy = self._local_entropy(lum_norm, index_valid, entropy_win, entropy_bins)
            ent_n, ent_lo, ent_hi = norm('ENTROPY', entropy)

            # Richesse = combinaison pondérée saturation + entropie
            richness = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            richness[index_valid] = np.clip(
                wr_sat * sat_n[index_valid] +
                wr_ent * ent_n[index_valid],
                0.0, 1.0
            )

            # Texture chromatique (variance ExG_ch) pour le garde-fou sol craquelé
            exgch_safe = np.where(index_valid, exg_ch, 0.0).astype(np.float32)
            chroma_var_list = []
            for k in windows:
                chroma_var_list.append(self._local_variance(exgch_safe, index_valid, k))
            chroma_var = np.maximum.reduce(chroma_var_list)
            del chroma_var_list
            chroma_tex_n, chromav_lo, chromav_hi = norm('CHROMA_VAR', chroma_var)

            # ── 8. Score BERGE v7 ────────────────────────────────────────────
            feedback.pushInfo('Calcul du score BERGE v7...')
            berge_score = np.full((ysize, xsize), NODATA_F, dtype=np.float32)
            berge_score[index_valid] = np.clip(
                w_spectral   * spectral[index_valid]      +
                w_green_dens * green_density[index_valid] +
                w_dry_tex    * dry_texture[index_valid]   +
                w_richness   * richness[index_valid],
                0.0, 1.0
            )

            feedback.pushInfo(
                f'  Poids effectifs — spectral:{w_spectral:.2f}  '
                f'densité:{w_green_dens:.2f}  texture:{w_dry_tex:.2f}  richesse:{w_richness:.2f}'
            )
            feedback.pushInfo(
                f'  Poids spectraux internes — CIVE:{ws_cive:.2f}  '
                f'ExG:{ws_exg:.2f}  VEG:{ws_veg:.2f}'
            )

            # ── 9. Classification ────────────────────────────────────────────
            feedback.pushInfo('Classification du score BERGE v7...')
            classes = np.zeros((ysize, xsize), dtype=np.uint8)
            classes[index_valid & (berge_score <  t_sol)]                              = 1
            classes[index_valid & (berge_score >= t_sol)  & (berge_score < t_veg)]    = 2
            classes[index_valid & (berge_score >= t_veg)  & (berge_score < t_dense)]  = 3
            classes[index_valid & (berge_score >= t_dense)]                            = 4

            # ── Garde-fous ───────────────────────────────────────────────────

            # GF1 : récupération végétation sèche (classe 1 → 2 si texture forte)
            dry_recovery = (
                index_valid &
                (classes == 1) &
                (dry_texture >= t_dry_rec) &
                (green_density >= 0.05)
            )
            classes[dry_recovery] = 2

            # GF2 : sol homogène (classe 2 → 1 si texture faible et score spectral faible)
            soil_guard = (
                index_valid &
                (classes == 2) &
                (dry_texture <= t_soil_hom) &
                (spectral < 0.35)
            )
            classes[soil_guard] = 1

            # GF3 : sol craquelé (classe ≥2 → 1 si spectral faible + chroma faible + texture forte)
            cracked_soil_guard = (
                index_valid &
                (classes >= 2) &
                (spectral     <= t_crack_spec) &
                (chroma_tex_n <= t_crack_chroma) &
                (dry_texture  >= t_crack_tex)
            )
            classes[cracked_soil_guard] = 1

            # GF4 : eau libre / pixels très sombres (forcer NoData)
            dark_water_guard = (
                valid &
                (lum      <= t_dark_lum) &
                (sat      <= t_dark_sat) &
                (classes  >= 1)
            )
            classes[dark_water_guard] = 0   # NoData

            # ── 10. Nettoyage spatial ────────────────────────────────────────
            self._write_tif(classes_raw, classes, gdal.GDT_Byte, 0, gt, proj,
                            self._class_color_table())
            final_classes = classes.copy()

            if min_pixels > 0:
                feedback.pushInfo(
                    f'Nettoyage spatial GDAL SieveFilter : objets < {min_pixels} pixels'
                )
                raw_ds  = gdal.Open(classes_raw)
                driver  = gdal.GetDriverByName('GTiff')
                out_ds  = driver.Create(classes_clean, xsize, ysize, 1, gdal.GDT_Byte,
                                        options=['COMPRESS=LZW', 'TILED=YES'])
                out_ds.SetGeoTransform(gt)
                out_ds.SetProjection(proj)
                src_b = raw_ds.GetRasterBand(1)
                dst_b = out_ds.GetRasterBand(1)
                dst_b.SetNoDataValue(0)
                dst_b.SetRasterColorTable(self._class_color_table())
                dst_b.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
                gdal.SieveFilter(src_b, src_b.GetMaskBand(), dst_b, int(min_pixels), 8)
                dst_b.FlushCache()
                out_ds.FlushCache()
                out_ds = None
                raw_ds = None
                clean_ds = gdal.Open(classes_clean)
                final_classes = clean_ds.GetRasterBand(1).ReadAsArray().astype(np.uint8)
                clean_ds = None
            else:
                self._write_tif(classes_clean, final_classes, gdal.GDT_Byte, 0, gt, proj,
                                self._class_color_table())

            # Rasters binaires dérivés
            veg_strict = np.zeros((ysize, xsize), dtype=np.uint8)
            veg_strict[(final_classes == 3) | (final_classes == 4)] = 1
            veg_eco = np.zeros((ysize, xsize), dtype=np.uint8)
            veg_eco[(final_classes == 2) | (final_classes == 3) | (final_classes == 4)] = 1

            # ── 11. GeoPackage ───────────────────────────────────────────────
            raster_items = []
            raster_style_specs = {}

            def add_raster(table, arr, dtype, nodata, color=None, style_kind=None, style_label=None):
                path = os.path.join(tmp_dir, f'{prefix}_{table}.tif')
                self._write_tif(path, arr, dtype, nodata, gt, proj, color)
                raster_items.append((table, path))
                if style_kind:
                    raster_style_specs[table] = {
                        'kind': style_kind,
                        'array': arr,
                        'nodata': nodata,
                        'label': style_label,
                    }
                return path

            raster_items.append(('rgb_clip', rgb_clip))
            raster_style_specs['rgb_clip'] = {'kind': 'rgb'}

            add_raster('classes', final_classes, gdal.GDT_Byte, 0,
                       self._class_color_table(), style_kind='classes')
            add_raster('vegetation_stricte', veg_strict, gdal.GDT_Byte, 0,
                       style_kind='binary_vegetation',
                       style_label='1 · végétation stricte classes 3+4')
            add_raster('vegetation_ecologique', veg_eco, gdal.GDT_Byte, 0,
                       style_kind='binary_vegetation',
                       style_label='1 · végétation écologique classes 2+3+4')
            add_raster('berge_score', berge_score, gdal.GDT_Float32, -9999,
                       style_kind='continuous')

            if keep_intermed:
                add_raster('cive', cive, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('exg', exg, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('veg_index', veg, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('spectral_index', spectral, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('green_density', green_density, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('dry_texture', dry_texture, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('saturation', sat, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('entropy_local', entropy, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('richness', richness, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('chroma_texture', chroma_tex_n, gdal.GDT_Float32, -9999,
                           style_kind='continuous')
                add_raster('guard_dry_recovery', dry_recovery.astype(np.uint8), gdal.GDT_Byte, 0,
                           style_kind='binary_guard', style_label='1 · récupération végétation sèche')
                add_raster('guard_soil_homogeneity', soil_guard.astype(np.uint8), gdal.GDT_Byte, 0,
                           style_kind='binary_guard', style_label='1 · correction sol homogène')
                add_raster('guard_cracked_soil', cracked_soil_guard.astype(np.uint8), gdal.GDT_Byte, 0,
                           style_kind='binary_guard', style_label='1 · correction sol craquelé')
                add_raster('guard_dark_water', dark_water_guard.astype(np.uint8), gdal.GDT_Byte, 0,
                           style_kind='binary_guard', style_label='1 · correction eau sombre')

            # ── 12. Statistiques ─────────────────────────────────────────────
            labels = {
                1: 'sol_nu_substrat_homogene',
                2: 'vegetation_seche_mixte_transition_texturee',
                3: 'vegetation_verte',
                4: 'vegetation_dense',
            }
            total = int(np.sum(final_classes > 0))
            if total == 0:
                raise QgsProcessingException('Aucun pixel classé.')
            counts = {code: int(np.sum(final_classes == code)) for code in labels}
            pcts   = {code: counts[code] / total * 100.0 for code in labels}
            pct_eco    = pcts[2] + pcts[3] + pcts[4]
            pct_strict = pcts[3] + pcts[4]
            pct_dense  = pcts[4]

            rows = []
            for code, label in labels.items():
                rows.append((site, date, 'BERGE_v7', 'classe', code, label,
                             counts[code], round(pcts[code], 3)))
            rows.extend([
                (site, date, 'BERGE_v7', 'indicateur', None,
                 'couverture_vegetale_ecologique_pct_classes_2_3_4',
                 counts[2]+counts[3]+counts[4], round(pct_eco, 3)),
                (site, date, 'BERGE_v7', 'indicateur', None,
                 'couverture_vegetale_stricte_pct_classes_3_4',
                 counts[3]+counts[4], round(pct_strict, 3)),
                (site, date, 'BERGE_v7', 'indicateur', None,
                 'couverture_vegetale_dense_pct_classe_4',
                 counts[4], round(pct_dense, 3)),
                (site, date, 'BERGE_v7', 'indicateur', None,
                 'substrat_ouvert_pct_classe_1',
                 counts[1], round(pcts[1], 3)),
            ])

            with open(csv_out, 'w', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                w.writerow(['site','date','algorithme','type','code',
                            'classe_ou_indicateur','pixels','pourcentage'])
                w.writerows(rows)

            # ── 13. Metadata JSON ────────────────────────────────────────────
            def pstats(arr):
                v = arr[index_valid]
                return {
                    'p10': float(np.nanpercentile(v, 10)),
                    'p25': float(np.nanpercentile(v, 25)),
                    'p50': float(np.nanpercentile(v, 50)),
                    'p75': float(np.nanpercentile(v, 75)),
                    'p90': float(np.nanpercentile(v, 90)),
                }

            effective_parameters = {
                'SITE': site, 'DATE': date, 'SOURCE_CRS': source_crs,
                'P_LOW': p_low, 'P_HIGH': p_high,
                'W_SPECTRAL': w_spectral, 'W_GREEN_DENSITY': w_green_dens,
                'W_DRY_TEXTURE': w_dry_tex, 'W_RICHNESS': w_richness,
                'WS_CIVE': ws_cive, 'WS_EXG': ws_exg, 'WS_VEG': ws_veg,
                'DENSITY_WINDOW': density_win,
                'MICRO_EXG': micro_exg, 'MICRO_NGRDI': micro_ngrdi,
                'MICRO_G': micro_g, 'MICRO_G_STRICT': micro_g_strict,
                'MULTISCALE_TEXTURE': multi,
                'TEXTURE_WINDOWS': ','.join(str(k) for k in windows),
                'WR_SAT': wr_sat, 'WR_ENT': wr_ent,
                'ENTROPY_WINDOW': entropy_win, 'ENTROPY_BINS': entropy_bins,
                'T_SOL': t_sol, 'T_VEG': t_veg, 'T_DENSE': t_dense,
                'T_DRY_RECOVERY': t_dry_rec,
                'T_SOIL_HOMOGENEITY': t_soil_hom,
                'T_CRACK_SOIL_SPECTRAL_MAX': t_crack_spec,
                'T_CRACK_SOIL_CHROMA_MAX': t_crack_chroma,
                'T_CRACK_SOIL_TEXTURE_MIN': t_crack_tex,
                'T_DARK_WATER_LUM': t_dark_lum,
                'T_DARK_WATER_SAT': t_dark_sat,
                'MIN_PIXELS': min_pixels,
                'TREAT_BLACK_AS_NODATA': black_as_nodata,
                'KEEP_INTERMEDIATE': keep_intermed,
                'APPLY_METADATA_PARAMS': apply_params,
                'APPLY_METADATA_NORMALIZATION': apply_norm,
                'PARAM_METADATA_JSON': meta_path_in or '',
            }

            metadata = {
                'algorithm':     'BERGE v7',
                'qgis_tested':   '3.36.3',
                'site':          site,
                'date':          date,
                'input_rgb':     rgb_layer.source(),
                'mask':          mask_layer.source(),
                'output_gpkg':   gpkg_out,
                'output_csv':    csv_out,
                'parameters':    effective_parameters,
                'metadata_import': {
                    'source_json':        meta_path_in or None,
                    'apply_parameters':   apply_params,
                    'apply_normalization': apply_norm,
                },
                'principle': (
                    'Score spectral RGB (CIVE+ExG+VEG, indices peu corrélés) '
                    '+ densité locale micro-pixels verts (logique OR assouplie) '
                    '+ variance Lum multi-échelle (végétation sèche) '
                    '+ richesse texturale (saturation chromatique + entropie locale)'
                ),
                'normalization_mode': 'metadata_bounds' if norm_bounds else 'local_percentiles',
                'normalization_bounds': {
                    'CIVE':       {'low': cive_lo,    'high': cive_hi},
                    'ExG':        {'low': exg_lo,     'high': exg_hi},
                    'VEG':        {'low': veg_lo,     'high': veg_hi},
                    'LUM_VAR':    {'low': lumv_lo,    'high': lumv_hi},
                    'SAT':        {'low': sat_lo,     'high': sat_hi},
                    'ENTROPY':    {'low': ent_lo,     'high': ent_hi},
                    'CHROMA_VAR': {'low': chromav_lo, 'high': chromav_hi},
                },
                'diagnostic_stats': {
                    'score':         pstats(berge_score),
                    'spectral':      pstats(spectral),
                    'green_density': pstats(green_density),
                    'dry_texture':   pstats(dry_texture),
                    'saturation':    pstats(sat),
                    'entropy':       pstats(entropy),
                    'richness':      pstats(richness),
                    'chroma_texture':pstats(chroma_tex_n),
                },
                'guard_pixels': {
                    'dry_recovery':     int(np.sum(dry_recovery)),
                    'soil_homogeneity': int(np.sum(soil_guard)),
                    'cracked_soil':     int(np.sum(cracked_soil_guard)),
                    'dark_water':       int(np.sum(dark_water_guard)),
                },
                'classes':     labels,
                'counts':      counts,
                'percentages': {str(k): round(v, 3) for k, v in pcts.items()},
                'indicators': {
                    'couverture_vegetale_ecologique_pct_classes_2_3_4': round(pct_eco, 3),
                    'couverture_vegetale_stricte_pct_classes_3_4':      round(pct_strict, 3),
                    'couverture_vegetale_dense_pct_classe_4':           round(pct_dense, 3),
                    'substrat_ouvert_pct_classe_1':                     round(pcts[1], 3),
                },
                'valid_pixels': total,
                'rasters_in_geopackage': [name for name, _ in raster_items],
                'tables_in_geopackage':  ['stats_pct', 'berge_parameters', 'berge_metadata_json'],
            }

            with open(metadata_out, 'w', encoding='utf-8') as f:
                json.dump(metadata, f, ensure_ascii=False, indent=2)

            # ── 14. Assemblage GeoPackage ────────────────────────────────────
            feedback.pushInfo('Création du GeoPackage unique de résultats...')
            first = True
            for table, tif in raster_items:
                self._add_raster_to_gpkg(tif, gpkg_out, table, first, feedback)
                first = False

            self._copy_vector_to_gpkg(mask_layer.source(), gpkg_out, 'emprise', feedback)
            self._write_tables_to_gpkg(gpkg_out, rows, metadata, effective_parameters, feedback)
            self._apply_gpkg_styles(gpkg_out, raster_style_specs, feedback)

            # ── 15. Résumé ───────────────────────────────────────────────────
            feedback.pushInfo('=' * 64)
            feedback.pushInfo('BERGE v7 terminé')
            feedback.pushInfo(f'GeoPackage : {gpkg_out}')
            feedback.pushInfo(f'CSV stats  : {csv_out}')
            feedback.pushInfo(f'metadata   : {metadata_out}')
            feedback.pushInfo(f'Couverture écologique (cl.2+3+4) : {pct_eco:.2f} %')
            feedback.pushInfo(f'Couverture stricte    (cl.3+4)   : {pct_strict:.2f} %')
            feedback.pushInfo(f'Végétation dense      (cl.4)     : {pct_dense:.2f} %')
            feedback.pushInfo(f'Substrat ouvert       (cl.1)     : {pcts[1]:.2f} %')
            feedback.pushInfo(
                f'Garde-fous activés — '
                f'dry_recovery:{int(np.sum(dry_recovery)):,}  '
                f'soil_hom:{int(np.sum(soil_guard)):,}  '
                f'cracked:{int(np.sum(cracked_soil_guard)):,}  '
                f'water:{int(np.sum(dark_water_guard)):,}'
            )
            feedback.pushInfo('=' * 64)
            if not norm_bounds:
                feedback.pushInfo(
                    f'Pour la prochaine campagne, fournir : {metadata_out}\n'
                    f'en paramètre PARAM_METADATA_JSON pour garantir la comparabilité temporelle.'
                )

            return {
                'GEOPACKAGE': gpkg_out,
                'CSV':        csv_out,
                'METADATA':   metadata_out,
            }

        finally:
            try:
                shutil.rmtree(tmp_dir)
            except Exception:
                pass

    # ── métadonnées QGIS ─────────────────────────────────────────────────────

    def name(self):
        return 'berge_v7_vegetation_rgb_texture'

    def displayName(self):
        return 'BERGE v7 - RGB · CIVE+ExG+VEG · entropie · saturation · 4 garde-fous'

    def group(self):
        return 'BERGE'

    def groupId(self):
        return 'berge'

    def shortHelpString(self):
        return (
            'BERGE v7 mesure la couverture végétale de berges et fossés à partir d\'une orthophoto '
            'drone RGB.\n\n'
            'Améliorations vs BERGE v6 :\n'
            '• Score spectral basé sur 3 indices peu corrélés : CIVE (robuste illumination), '
            'ExG (végétation verte/mixte), VEG (non-linéaire).\n'
            '• Texture sèche simplifiée : variance Lum seule, multi-échelle 7/15/31 px.\n'
            '• Richesse texturale : saturation chromatique + entropie locale (NumPy pur).\n'
            '• Logique micro-pixel assouplie : OR strict (g_ch seul si très élevé) '
            'OU OR combiné (g_ch modéré + ExG ou NGRDI).\n'
            '• 4ᵉ garde-fou : eau libre / pixels très sombres → exclus des statistiques.\n'
            '• Poids spectraux internes (CIVE/ExG/VEG) exposés dans l\'interface.\n\n'
            'Indicateur principal : couverture_vegetale_ecologique_pct (classes 2+3+4).\n'
            'Indicateur strict    : couverture_vegetale_stricte_pct   (classes 3+4).'
        )

    def createInstance(self):
        return BergeV7Algorithm()
