#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistik_deppSeq-ver21.py
MLVA/VNTR contamination analysis — output module v21.

========================================================================
MATHEMATICAL OVERHAUL per specification v21
========================================================================

IDEALIZED MODEL (single scale, no per-locus renormalization of Q_i):
    N_filt   = Σ Raw_Count          over all detected alleles of a sample
    G_mini   = Σ Fragment_Length    over all detected alleles of a sample
    C_mini   = N_filt · R / G_mini  where R = Mean_Read_Length_bp
    p_i      = p_star  (fallback: max(0, (R - L_i + 1) / R); 0 if L_i > R)
    E_ideal  = C_mini · p_i                         ← NO length multiplier
    Q_i      = Raw_Count / E_ideal                  ← direct ratio
    Idealized_Share_% (per locus) = Q_i / Σ(Q_j in locus) · 100

EXPERIMENTAL MODEL (two-step; pure loci / pure alleles):
    STEP 1 — sample baseline over clean monoallelic loci of the sample:
        B_sample = median( Q_i over clean monoallelic loci )
        Q_i_sample_adjusted = Q_i / B_sample

    STEP 2 — cohort locus baseline over clean monoallelic observations
             of the same locus across samples:
        B_locus_cohort = median( Q_i_sample_adjusted )
        Q_i_experimental = Q_i_sample_adjusted / B_locus_cohort
        (fallback: Q_i_experimental = Q_i_sample_adjusted if no cohort baseline)

    Experimental_Share_% (per locus) = Q_i_experimental / Σ(..) · 100

CLEAN MONOALLELIC LOCUS:
    exactly one detected allele AND Raw_Count > 0 AND p_i > 0
    AND Detectability_Class != 'Unobservable'
    AND Dropout_Adjusted_Status ∉ EXCLUDE_DROPOUT_STATUSES

REMOVED (vs. v20):
    • Q_e, old per-locus B_locus (median within locus used as bias)
    • Integrated_Mass, Integrated_Percent, Sample_Percent
    • Major_Allele_Reliability_Pct, Minor_to_Major_Ratio
    • Legacy Idealized_Support = Raw / p_star
    • All plot code (build_plots_html) and Plot_* columns
    • classify_locus_status_* (all model-based status calls)

OUTPUT FILES:
    • allele_level.xlsx, locus_level, sample_level  — master workbook
    • publication_tables.xlsx → primary_allele_table, primary_locus_table
    • output_dictionary.tsv
    • run_metadata.txt
    • interpretation_notes.txt
    • plots_idealized_share.html          (Plotly; per-sample Idealized_Share_%)
    • all_samples_idealized_alleles.html  (Plotly; all samples, scrollable)
    • plots_experimental_share.html       (Plotly; per-sample Experimental_Share_%)
    • all_samples_experimental_alleles.html
    • hidden_allele_interpretation.txt

HIDDEN ALLELE LAYER (v21.1 addition):
    Adds Hidden_Allele_Flag, Hidden_Component_Pct (allele_level),
    Potential_Hidden_Allele, Hidden_Allele_Evidence_Score (locus_level),
    N_Loci_Potential_Hidden_Allele, Potential_Mixture_Flag, Mixture_Type
    (sample_level).  All Q_i formulas are UNCHANGED.

Version: 21.1.0
"""

import os
import logging
import argparse
from datetime import datetime

import numpy as np
import pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.formatting.rule import FormulaRule
from openpyxl.utils import get_column_letter

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    _PLOTLY_AVAILABLE = True
except ImportError:
    _PLOTLY_AVAILABLE = False

__version__ = "21.1.0"

# ================== CONFIGURABLE PARAMETERS ==================

# Minimum number of clean monoallelic loci required for the sample
# baseline (B_sample) to be flagged as reliable.
MIN_CLEAN_MONOALLELIC_FOR_SAMPLE = 5

# Minimum number of clean monoallelic observations required across the
# cohort for the cohort locus baseline (B_locus_cohort) to be flagged
# as reliable.
MIN_CLEAN_MONOALLELIC_FOR_LOCUS  = 5

# Minimum Raw_Count required to consider an allele detected (for
# inclusion at all in the allele_level table). v20 default is kept.
MIN_ALLELE_SUPPORT = 1

# Locus-level statuses that disqualify a locus from being considered
# "clean" (used both for clean monoallelic selection and for the
# unobservable/dropout detection in denominator completeness).
EXCLUDE_DROPOUT_STATUSES = {
    "PARTIALLY_UNOBSERVABLE_MIXTURE",
    "UNRESOLVED_DUE_TO_DETECTABILITY",
    "POTENTIAL_ALLELE_DROPOUT",
}

# Threshold on Hidden_Component_Pct = 100·(1 - Locus_Completeness_Obs)
# above which the locus denominator is considered incomplete.
HIDDEN_COMPONENT_THRESHOLD = 30.0

# Numeric rounding for output tables.
ROUND_Q    = 6
ROUND_PCT  = 4
# =============================================================

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────
# Locus name parsing  (kept from v20)
# ─────────────────────────────────────────────────────────────

def parse_locus_id(locus_id):
    parts = str(locus_id).split('_')
    return {
        'Locus_Name':                parts[0] if parts else locus_id,
        'Repeat_Descriptor':         parts[1] if len(parts) > 1 else None,
        'Reference_Fragment_Label':  parts[2] if len(parts) > 2 else None,
        'Reference_Allele_Label':    parts[3] if len(parts) > 3 else None,
    }


def detect_allele_blocks(df):
    """Return sorted list of allele indices found in wide-format columns."""
    cols = [c for c in df.columns if c.startswith('Allele ') and c[7:].isdigit()]
    return sorted({int(c.split()[1]) for c in cols})


# ─────────────────────────────────────────────────────────────
# Small numeric helpers
# ─────────────────────────────────────────────────────────────

def _is_missing(v):
    """True if v is None or NaN. Avoids pd.isna on scalars that break on strings."""
    if v is None:
        return True
    if isinstance(v, float) and np.isnan(v):
        return True
    return False


def _to_float(v):
    """Safe numeric cast. Returns None for missing/non-numeric."""
    if _is_missing(v):
        return None
    try:
        x = float(v)
    except (ValueError, TypeError):
        return None
    if np.isnan(x):
        return None
    return x


def _to_bool(v):
    """Safe boolean cast. Returns None for missing."""
    if _is_missing(v):
        return None
    if isinstance(v, str):
        s = v.strip().lower()
        if s in ('true', 't', '1', 'yes'):  return True
        if s in ('false', 'f', '0', 'no'):  return False
        return None
    try:
        return bool(v)
    except (ValueError, TypeError):
        return None


def _median_finite(values):
    """Median over finite values; returns None if none are finite."""
    arr = np.array([v for v in values if v is not None and np.isfinite(v)],
                   dtype=float)
    if arr.size == 0:
        return None
    return float(np.median(arr))


# ─────────────────────────────────────────────────────────────
# STEP 1.  Flatten wide-format input to long-format allele_level
# (raw attributes only — idealized/experimental fields filled later).
# ─────────────────────────────────────────────────────────────

def build_allele_level_raw(data, indices, max_alleles, sample_meta):
    """
    Convert the wide-format MLVA_v2 output into a long-format DataFrame
    with one row per (sample, locus, allele).

    Attached columns (raw/input only at this stage):
        Sample_ID, Access_Number, Locus_ID, Locus_Name, Repeat_Descriptor,
        Reference_Fragment_Label, Reference_Allele_Label, Allele_Rank_Input,
        Allele_Value, Raw_Count, Fragment_Length_bp, Mean_Read_Length_bp,
        p_star_input, Detectability_Class, Use_in_Idealized,
        Is_Potential_Dropout,
        Dropout_Adjusted_Status, Confidence_Score_Locus,
        Locus_Completeness_Obs, Incomplete_Denominator.
    """
    records = []

    for _, row in data.iterrows():
        acc        = str(row.get('Access_number', ''))
        sample_id  = row.get('Sample_ID')
        if _is_missing(sample_id) or not str(sample_id).strip():
            sample_id = acc
        locus_id   = str(row.get('Primer', ''))
        lp         = parse_locus_id(locus_id)

        meta       = sample_meta.get(acc, {})
        mean_rl    = meta.get('mean_read_length')

        dropout_status    = row.get('Dropout_Adjusted_Status')
        confidence_score  = _to_float(row.get('Confidence_Score'))
        locus_completeness= _to_float(row.get('Locus_Completeness_Obs'))
        incomplete_denom  = _to_bool(row.get('Incomplete_Denominator'))

        for i in indices:
            if i > max_alleles:
                continue

            allele_val = row.get(f'Allele {i}')
            if _is_missing(allele_val) or str(allele_val).strip() == '':
                continue

            raw_count   = _to_float(row.get(f'Raw_Count {i}'))
            frag_len    = _to_float(row.get(f'Loci_Size {i}'))
            p_star_in   = _to_float(row.get(f'Detectability_p {i}'))

            det_class_raw = row.get(f'Detectability_Class {i}')
            det_class     = (str(det_class_raw).strip()
                             if not _is_missing(det_class_raw) else None)

            use_in_ideal  = _to_bool(row.get(f'Use_in_Idealized {i}'))
            is_dropout    = _to_bool(row.get(f'Is_Potential_Dropout {i}'))

            if raw_count is None or raw_count < MIN_ALLELE_SUPPORT:
                continue

            records.append({
                # Identification
                'Sample_ID':                 sample_id,
                'Access_Number':             acc,
                'Locus_ID':                  locus_id,
                'Locus_Name':                lp['Locus_Name'],
                'Repeat_Descriptor':         lp['Repeat_Descriptor'],
                'Reference_Fragment_Label':  lp['Reference_Fragment_Label'],
                'Reference_Allele_Label':    lp['Reference_Allele_Label'],
                'Allele_Rank_Input':         i,
                'Allele_Value':              allele_val,
                # Raw data
                'Raw_Count':                 raw_count,
                'Fragment_Length_bp':        frag_len,
                'Mean_Read_Length_bp':       mean_rl,
                # Detectability (input)
                'p_star_input':              p_star_in,
                'Detectability_Class':       det_class,
                'Use_in_Idealized':          bool(use_in_ideal) if use_in_ideal is not None else False,
                'Is_Potential_Dropout':      bool(is_dropout)   if is_dropout   is not None else False,
                # Locus-level attributes (copied per allele row)
                'Dropout_Adjusted_Status':   dropout_status,
                'Confidence_Score_Locus':    confidence_score,
                'Locus_Completeness_Obs':    locus_completeness,
                'Incomplete_Denominator':    bool(incomplete_denom) if incomplete_denom is not None else False,
            })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
# STEP 2.  Idealized model (mini-genome coverage)
#
#   N_filt   = Σ Raw_Count       per sample
#   G_mini   = Σ Fragment_Length per sample
#   C_mini   = N_filt · R / G_mini      (R = Mean_Read_Length_bp)
#   p_i      = p_star_input (fallback)
#   E_ideal  = C_mini · p_i             ← NO length multiplier
#   Q_i      = Raw_Count / E_ideal
# ─────────────────────────────────────────────────────────────



# ─────────────────────────────────────────────────────────────
# STEP 1b.  build_allele_level_full — REPORTING-ONLY allele list
#
#  Same logic as build_allele_level_raw BUT:
#  • No MIN_ALLELE_SUPPORT filter → Raw_Count = 0 rows are KEPT
#  • Adds Allele_Detection_Status and Is_Reporting_Only columns
#
#  ❗ CRITICAL: this DataFrame is used ONLY for publication sheets
#     (reliable_alleles_by_sample, excluded_hidden_or_incomplete_loci).
#     It MUST NOT be passed to any calculation function
#     (apply_idealized_model, apply_experimental_model, etc.).
#
#  Hidden_Allele_Flag / Hidden_Component_Pct semantics (specification section 7):
#     These describe observability limits, NOT allele presence.
#     "Not_detected" means the allele was not observed; it does NOT
#     mean it is absent from the biological sample.
# ─────────────────────────────────────────────────────────────

def _compute_p_i_scalar(p_star_input, frag_len, mean_rl):
    """Return p_i for a single allele (same formula as apply_idealized_model)."""
    ps = _to_float(p_star_input)
    if ps is not None:
        return max(0.0, ps)
    rl = _to_float(mean_rl)
    fl = _to_float(frag_len)
    if rl is None or rl <= 0 or fl is None:
        return None
    if fl > rl:
        return 0.0
    return max(0.0, (rl - fl + 1.0) / rl)


def _allele_detection_status(raw_count, p_i, det_class):
    """
    Allele_Detection_Status per specification section 3.1.
    'Detected'                     : Raw_Count > 0
    'Not_detected'                 : Raw_Count == 0, p_i > 0
    'Not_evaluable_by_read_length' : p_i <= 0 or class Unobservable
    """
    if raw_count is not None and raw_count > 0:
        return "Detected"
    if det_class == "Unobservable" or (p_i is not None and p_i <= 0.0):
        return "Not_evaluable_by_read_length"
    return "Not_detected"


def build_allele_level_full(data, indices, max_alleles, sample_meta):
    """
    Build a REPORTING-ONLY long-format allele DataFrame.

    Differences from build_allele_level_raw:
      • Includes alleles with Raw_Count == 0 (no MIN_ALLELE_SUPPORT filter).
      • Adds Is_Reporting_Only  = True  for Raw_Count == 0 rows
                               = False for Raw_Count  > 0 rows
      • Adds Allele_Detection_Status column.
      • Adds p_i_full column (for downstream filter without running idealized model).

    This function is intentionally separate so that the calculation pipeline
    (build_allele_level_raw → apply_idealized_model → ...) remains unchanged.
    """
    records = []

    for _, row in data.iterrows():
        acc        = str(row.get('Access_number', ''))
        sample_id  = row.get('Sample_ID')
        if _is_missing(sample_id) or not str(sample_id).strip():
            sample_id = acc
        locus_id   = str(row.get('Primer', ''))
        lp         = parse_locus_id(locus_id)

        meta       = sample_meta.get(acc, {})
        mean_rl    = meta.get('mean_read_length')

        dropout_status    = row.get('Dropout_Adjusted_Status')
        confidence_score  = _to_float(row.get('Confidence_Score'))
        locus_completeness= _to_float(row.get('Locus_Completeness_Obs'))
        incomplete_denom  = _to_bool(row.get('Incomplete_Denominator'))

        for i in indices:
            if i > max_alleles:
                continue

            allele_val = row.get(f'Allele {i}')
            if _is_missing(allele_val) or str(allele_val).strip() == '':
                continue  # allele slot completely absent in input → skip

            raw_count   = _to_float(row.get(f'Raw_Count {i}'))
            frag_len    = _to_float(row.get(f'Loci_Size {i}'))
            p_star_in   = _to_float(row.get(f'Detectability_p {i}'))

            det_class_raw = row.get(f'Detectability_Class {i}')
            det_class     = (str(det_class_raw).strip()
                             if not _is_missing(det_class_raw) else None)

            use_in_ideal  = _to_bool(row.get(f'Use_in_Idealized {i}'))
            is_dropout    = _to_bool(row.get(f'Is_Potential_Dropout {i}'))

            # ── normalise raw_count: treat None / missing as 0 ──────────────
            if raw_count is None:
                raw_count = 0.0

            # ── compute p_i for this allele ──────────────────────────────────
            p_i = _compute_p_i_scalar(p_star_in, frag_len, mean_rl)

            # ── detection status ──────────────────────────────────────────────
            det_status = _allele_detection_status(raw_count, p_i, det_class)
            is_rep_only = (raw_count == 0)

            records.append({
                # Identification
                'Sample_ID':                 sample_id,
                'Access_Number':             acc,
                'Locus_ID':                  locus_id,
                'Locus_Name':                lp['Locus_Name'],
                'Allele_Rank_Input':         i,
                'Allele_Value':              allele_val,
                # Raw data
                'Raw_Count':                 raw_count,
                'Fragment_Length_bp':        frag_len,
                'Mean_Read_Length_bp':       mean_rl,
                # Detectability
                'p_i_full':                  p_i,
                'p_star_input':              p_star_in,
                'Detectability_Class':       det_class,
                'Use_in_Idealized':          bool(use_in_ideal) if use_in_ideal is not None else False,
                'Is_Potential_Dropout':      bool(is_dropout)   if is_dropout   is not None else False,
                # Locus-level (from input)
                'Dropout_Adjusted_Status':   dropout_status,
                'Confidence_Score_Locus':    confidence_score,
                'Locus_Completeness_Obs':    locus_completeness,
                'Incomplete_Denominator':    bool(incomplete_denom) if incomplete_denom is not None else False,
                # Reporting-only fields
                'Allele_Detection_Status':   det_status,
                'Is_Reporting_Only':         is_rep_only,
                # Locus-level flags from df_locus (merged later)
                'Denominator_Complete':      None,
                'Potential_Hidden_Allele':   None,
                'Hidden_Allele_Flag':        None,
                'Hidden_Component_Pct':      None,
                'Denominator_Issue_Type':    None,
            })

    return pd.DataFrame(records)

def apply_idealized_model(df):
    """
    Populate mini-genome block (N_filt, G_mini, C_mini), detectability
    (p_i, Is_Unobservable), and idealized ratios (E_ideal, Q_i) per allele.

    Fills the per-locus Idealized_Share_% after Q_i is known.
    """
    df = df.copy()

    # --- per-sample aggregates ---------------------------------
    # N_filt — sum of all Raw_Count of detected alleles in the sample.
    # G_mini — sum of Fragment_Length of all detected alleles of the sample.
    # Both sums include EVERY detected allele, regardless of Use_in_Idealized.
    sample_agg = (df
                  .groupby('Sample_ID')
                  .agg(
                      N_filt=('Raw_Count',          lambda s: float(pd.to_numeric(s, errors='coerce').dropna().sum())),
                      G_mini=('Fragment_Length_bp', lambda s: float(pd.to_numeric(s, errors='coerce').dropna().sum())),
                      R_mean=('Mean_Read_Length_bp', lambda s: (
                          float(pd.to_numeric(s, errors='coerce').dropna().iloc[0])
                          if pd.to_numeric(s, errors='coerce').dropna().size > 0 else np.nan))
                  )
                  .reset_index())

    def _c_mini(r):
        n, g, rl = r['N_filt'], r['G_mini'], r['R_mean']
        if any(_is_missing(x) for x in (n, g, rl)) or g <= 0 or rl <= 0:
            return np.nan
        return n * rl / g

    sample_agg['C_mini'] = sample_agg.apply(_c_mini, axis=1)

    df = df.merge(
        sample_agg[['Sample_ID', 'N_filt', 'G_mini', 'C_mini']],
        on='Sample_ID', how='left'
    )

    # --- per-allele p_i ----------------------------------------
    # Preference: p_star_input; fallback computed from R and L.
    def _p_i(row):
        ps = _to_float(row.get('p_star_input'))
        if ps is not None:
            return max(0.0, ps)
        rl = _to_float(row.get('Mean_Read_Length_bp'))
        fl = _to_float(row.get('Fragment_Length_bp'))
        if rl is None or rl <= 0 or fl is None:
            return None
        if fl > rl:
            return 0.0
        return max(0.0, (rl - fl + 1.0) / rl)

    df['p_i'] = df.apply(_p_i, axis=1)

    # Is_Unobservable: p_i <= 0 OR Detectability_Class == 'Unobservable'
    def _is_unobservable(row):
        pi = _to_float(row.get('p_i'))
        if pi is not None and pi <= 0:
            return True
        if row.get('Detectability_Class') == 'Unobservable':
            return True
        return False

    df['Is_Unobservable'] = df.apply(_is_unobservable, axis=1)

    # --- per-allele E_ideal and Q_i ----------------------------
    def _e_ideal(row):
        cov = _to_float(row.get('C_mini'))
        pi  = _to_float(row.get('p_i'))
        if cov is None or cov <= 0:
            return None
        if pi is None:
            return None
        if pi == 0.0:
            return 0.0
        return cov * pi

    df['E_ideal'] = df.apply(_e_ideal, axis=1)

    def _q_i(row):
        raw = _to_float(row.get('Raw_Count'))
        e   = _to_float(row.get('E_ideal'))
        if raw is None or e is None or e <= 0:
            return None
        return raw / e

    df['Q_i'] = df.apply(_q_i, axis=1)

    # --- Idealized_Share_% per (Sample_ID, Locus_ID) -----------
    # Share = Q_i / Σ(Q_j of same locus with valid Q_i) · 100
    df['Idealized_Share_%'] = np.nan
    for (sid, lid), idx in df.groupby(['Sample_ID', 'Locus_ID']).groups.items():
        q_vals = df.loc[idx, 'Q_i']
        valid  = q_vals.dropna()
        total  = float(valid.sum()) if len(valid) > 0 else None
        if total is not None and total > 0:
            df.loc[idx, 'Idealized_Share_%'] = (
                df.loc[idx, 'Q_i'].astype(float) / total * 100.0
            )

    return df


# ─────────────────────────────────────────────────────────────
# STEP 3.  Clean monoallelic locus flag
#
#   A locus (within a given sample) is clean monoallelic iff:
#     • exactly one detected allele in that locus,
#     • the allele has Raw_Count > 0,
#     • the allele has p_i > 0,
#     • Is_Unobservable == False for the allele,
#     • Dropout_Adjusted_Status ∉ EXCLUDE_DROPOUT_STATUSES.
#
#   This flag is attached per allele row (constant across the locus)
#   so it can be used both for baseline selection and for reporting.
# ─────────────────────────────────────────────────────────────

def apply_clean_monoallelic_flag(df):
    """Attach Is_Clean_Monoallelic_Locus (bool) on every allele row."""
    df = df.copy()
    df['Is_Clean_Monoallelic_Locus'] = False

    for (sid, lid), grp_idx in df.groupby(['Sample_ID', 'Locus_ID']).groups.items():
        grp = df.loc[grp_idx]

        if len(grp) != 1:
            # more than one allele detected at this locus → not monoallelic
            continue

        row = grp.iloc[0]
        raw = _to_float(row.get('Raw_Count'))
        pi  = _to_float(row.get('p_i'))
        ds  = row.get('Dropout_Adjusted_Status')

        if raw is None or raw <= 0:
            continue
        if pi is None or pi <= 0:
            continue
        if bool(row.get('Is_Unobservable', False)):
            continue
        if str(ds) in EXCLUDE_DROPOUT_STATUSES:
            continue

        df.loc[grp_idx, 'Is_Clean_Monoallelic_Locus'] = True

    return df


# ─────────────────────────────────────────────────────────────
# STEP 4.  Experimental model — two-step
#
#   STEP 4A — B_sample (per sample):
#       B_sample = median(Q_i over clean monoallelic loci of sample)
#       Sample_Baseline_Reliable = N_clean ≥ MIN_CLEAN_MONOALLELIC_FOR_SAMPLE
#       Q_i_sample_adjusted = Q_i / B_sample          (NA if B_sample None)
#
#   STEP 4B — B_locus_cohort (per locus, over cohort):
#       B_locus_cohort = median(Q_i_sample_adjusted over clean monoallelic
#                               observations of the same locus across samples)
#       Locus_Cohort_Baseline_Reliable = N_obs ≥ MIN_CLEAN_MONOALLELIC_FOR_LOCUS
#       Q_i_experimental = Q_i_sample_adjusted / B_locus_cohort
#       If cohort baseline unavailable: Q_i_experimental = Q_i_sample_adjusted
#                                       and Experimental_Cohort_Adjustment_Applied = False
# ─────────────────────────────────────────────────────────────

def apply_experimental_model(df):
    """
    Attach:
        B_sample, Sample_Baseline_Reliable, Q_i_sample_adjusted,
        B_locus_cohort, Locus_Cohort_Baseline_Reliable,
        Experimental_Cohort_Adjustment_Applied,
        Q_i_experimental, Experimental_Share_%.

    Returns the enriched DataFrame.
    """
    df = df.copy()

    # ---------------- STEP 4A — per-sample baseline ----------------
    b_sample_map  = {}  # sample_id → (B_sample, reliable)
    for sid, grp_idx in df.groupby('Sample_ID').groups.items():
        grp = df.loc[grp_idx]
        clean = grp[grp['Is_Clean_Monoallelic_Locus'] == True]
        # Guarantee: each clean monoallelic locus contributes exactly one
        # allele row, so clean['Q_i'] is already per-locus.
        q_vals = [_to_float(v) for v in clean['Q_i'].tolist()]
        q_vals = [v for v in q_vals if v is not None and v > 0]

        n_clean = len(q_vals)
        med     = _median_finite(q_vals) if n_clean > 0 else None
        reliable = (n_clean >= MIN_CLEAN_MONOALLELIC_FOR_SAMPLE
                    and med is not None and med > 0)

        b_sample_map[sid] = (med, bool(reliable), n_clean)

    df['B_sample'] = df['Sample_ID'].map(lambda s: b_sample_map.get(s, (None, False, 0))[0])
    df['Sample_Baseline_Reliable'] = df['Sample_ID'].map(
        lambda s: b_sample_map.get(s, (None, False, 0))[1])

    def _q_i_sample_adj(row):
        qi = _to_float(row.get('Q_i'))
        bs = _to_float(row.get('B_sample'))
        if qi is None:
            return None
        if bs is None or bs <= 0:
            return None
        return qi / bs

    df['Q_i_sample_adjusted'] = df.apply(_q_i_sample_adj, axis=1)

    # ---------------- STEP 4B — per-locus cohort baseline ----------
    # Build pool: clean monoallelic observations with valid Q_i_sample_adjusted,
    # grouped by Locus_ID (across all samples).
    b_locus_cohort_map = {}  # locus_id → (B_locus_cohort, reliable, n_obs)
    pool = df[(df['Is_Clean_Monoallelic_Locus'] == True) &
              (df['Q_i_sample_adjusted'].notna())]

    for lid, grp_idx in pool.groupby('Locus_ID').groups.items():
        vals = [_to_float(v) for v in pool.loc[grp_idx, 'Q_i_sample_adjusted'].tolist()]
        vals = [v for v in vals if v is not None and v > 0]

        n_obs = len(vals)
        med   = _median_finite(vals) if n_obs > 0 else None
        reliable = (n_obs >= MIN_CLEAN_MONOALLELIC_FOR_LOCUS
                    and med is not None and med > 0)
        b_locus_cohort_map[lid] = (med, bool(reliable), n_obs)

    df['B_locus_cohort'] = df['Locus_ID'].map(
        lambda l: b_locus_cohort_map.get(l, (None, False, 0))[0])
    df['Locus_Cohort_Baseline_Reliable'] = df['Locus_ID'].map(
        lambda l: b_locus_cohort_map.get(l, (None, False, 0))[1])

    # ---------------- Q_i_experimental and adjustment flag ----------
    def _q_i_exp(row):
        qsa = _to_float(row.get('Q_i_sample_adjusted'))
        blc = _to_float(row.get('B_locus_cohort'))
        if qsa is None:
            return None
        if blc is None or blc <= 0:
            # no cohort baseline → fallback to sample-only adjustment
            return qsa
        return qsa / blc

    def _adj_applied(row):
        qsa = _to_float(row.get('Q_i_sample_adjusted'))
        blc = _to_float(row.get('B_locus_cohort'))
        if qsa is None:
            return False
        return (blc is not None and blc > 0)

    df['Q_i_experimental'] = df.apply(_q_i_exp, axis=1)
    df['Experimental_Cohort_Adjustment_Applied'] = df.apply(_adj_applied, axis=1)

    # ---------------- Experimental_Share_% per locus ----------------
    df['Experimental_Share_%'] = np.nan
    for (sid, lid), idx in df.groupby(['Sample_ID', 'Locus_ID']).groups.items():
        e_vals = df.loc[idx, 'Q_i_experimental']
        valid  = e_vals.dropna()
        total  = float(valid.sum()) if len(valid) > 0 else None
        if total is not None and total > 0:
            df.loc[idx, 'Experimental_Share_%'] = (
                df.loc[idx, 'Q_i_experimental'].astype(float) / total * 100.0
            )

    return df


# ─────────────────────────────────────────────────────────────
# STEP 5.  Denominator completeness (per locus)
#
#   Issue types:
#     UNOBSERVABLE     — any allele with Is_Unobservable=True
#     LOW_COMPLETENESS — Incomplete_Denominator flag OR Hidden>threshold
#     MIXED            — both above
#     NONE             — clean
# ─────────────────────────────────────────────────────────────

def apply_denominator_completeness(df):
    """Attach Denominator_Complete (bool) and Denominator_Issue_Type (str)
    per allele row, constant across the locus."""
    df = df.copy()
    df['Denominator_Complete']   = True
    df['Denominator_Issue_Type'] = 'NONE'

    for (sid, lid), idx in df.groupby(['Sample_ID', 'Locus_ID']).groups.items():
        grp = df.loc[idx]

        has_unobs = bool(grp['Is_Unobservable'].any())

        incompl_vals = grp['Incomplete_Denominator'].dropna().astype(bool)
        incompl_flag = bool(incompl_vals.any()) if len(incompl_vals) > 0 else False

        lc_vals = grp['Locus_Completeness_Obs'].dropna()
        lc_val  = float(lc_vals.iloc[0]) if len(lc_vals) > 0 else None
        hidden_pct = (100.0 * (1.0 - lc_val)) if lc_val is not None else None
        low_completeness = (hidden_pct is not None
                            and hidden_pct > HIDDEN_COMPONENT_THRESHOLD)

        is_low = incompl_flag or low_completeness

        if not has_unobs and not is_low:
            it = 'NONE'
        elif has_unobs and not is_low:
            it = 'UNOBSERVABLE'
        elif not has_unobs and is_low:
            it = 'LOW_COMPLETENESS'
        else:
            it = 'MIXED'

        df.loc[idx, 'Denominator_Complete']   = (it == 'NONE')
        df.loc[idx, 'Denominator_Issue_Type'] = it

    return df


# ─────────────────────────────────────────────────────────────
# STEP 5b.  Hidden allele annotation (per allele row)
#
#   Hidden_Allele_Flag  — True when the allele is physically undetectable
#       because its amplicon is too long for the sequencing read (p_i ≤ 0)
#       or its Detectability_Class is explicitly 'Unobservable'.
#       Biologically: the allele may EXIST in the sample but is INVISIBLE.
#
#   Hidden_Component_Pct — fraction of allele signal lost due to detectability:
#       Hidden_Component_Pct = (1 − p_i) × 100   [clamped 0–100]
#       If p_i = 0   → 100 % (nothing observed)
#       If p_i = 1   →   0 % (fully observed)
#       If p_i = NaN → 100 % (worst-case assumption)
#
#   This is a QC / interpretation OVERLAY — it does NOT alter Q_i or any
#   other formula from STEP 2.
# ─────────────────────────────────────────────────────────────

def apply_hidden_allele_flags(df):
    """
    Attach Hidden_Allele_Flag (bool) and Hidden_Component_Pct (float) per
    allele row.  Must be called after apply_idealized_model (needs p_i).
    """
    df = df.copy()
    p_i_num = pd.to_numeric(df['p_i'], errors='coerce')

    # Hidden_Allele_Flag: p_i ≤ 0 (NaN counts as 0→False here) OR class == Unobservable
    p_zero   = (p_i_num <= 0).fillna(False)
    det_unob = df['Detectability_Class'].fillna('').eq('Unobservable')
    df['Hidden_Allele_Flag'] = (p_zero | det_unob)

    # Hidden_Component_Pct = (1 − p_i) × 100; NaN p_i → 100 (worst-case)
    p_clipped = p_i_num.clip(lower=0.0, upper=1.0)
    df['Hidden_Component_Pct'] = ((1.0 - p_clipped) * 100.0).fillna(100.0).round(2)

    return df


# ─────────────────────────────────────────────────────────────
# STEP 6.  Final allele_level column ordering + rounding
# ─────────────────────────────────────────────────────────────

ALLELE_LEVEL_COLUMNS = [
    # Identification
    'Sample_ID', 'Access_Number', 'Locus_ID', 'Locus_Name',
    'Repeat_Descriptor', 'Reference_Fragment_Label', 'Reference_Allele_Label',
    'Allele_Rank_Input', 'Allele_Value',
    # Raw data
    'Raw_Count', 'Fragment_Length_bp', 'Mean_Read_Length_bp',
    'p_i', 'Detectability_Class',
    # Mini-genome block
    'N_filt', 'G_mini', 'C_mini',
    # Idealized block
    'E_ideal', 'Q_i', 'Idealized_Share_%',
    # Sample baseline block
    'Is_Clean_Monoallelic_Locus', 'B_sample', 'Sample_Baseline_Reliable',
    'Q_i_sample_adjusted',
    # Cohort locus block
    'B_locus_cohort', 'Locus_Cohort_Baseline_Reliable',
    'Experimental_Cohort_Adjustment_Applied',
    'Q_i_experimental', 'Experimental_Share_%',
    # QC / completeness
    'Is_Unobservable', 'Incomplete_Denominator',
    'Denominator_Complete', 'Denominator_Issue_Type',
    'Dropout_Adjusted_Status',
    # Hidden allele interpretation layer
    'Hidden_Allele_Flag', 'Hidden_Component_Pct',
]


def finalize_allele_level(df):
    """Order columns per spec, round numerics, coerce booleans."""
    df = df.copy()

    # Round numerics
    for c in ('p_i', 'E_ideal', 'Q_i', 'B_sample', 'Q_i_sample_adjusted',
              'B_locus_cohort', 'Q_i_experimental', 'C_mini'):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce').round(ROUND_Q)

    for c in ('Idealized_Share_%', 'Experimental_Share_%'):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce').round(ROUND_PCT)

    for c in ('N_filt', 'G_mini'):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

    for c in ('Hidden_Component_Pct',):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce').round(2)

    # Boolean dtypes (nullable BooleanDtype so NA can propagate)
    bool_cols = [
        'Is_Clean_Monoallelic_Locus', 'Sample_Baseline_Reliable',
        'Locus_Cohort_Baseline_Reliable',
        'Experimental_Cohort_Adjustment_Applied',
        'Is_Unobservable', 'Incomplete_Denominator', 'Denominator_Complete',
        'Hidden_Allele_Flag',
    ]
    for c in bool_cols:
        if c in df.columns:
            df[c] = df[c].astype('boolean')

    # Column ordering — keep any extras at the end (for debug).
    final_cols = [c for c in ALLELE_LEVEL_COLUMNS if c in df.columns]
    extras     = [c for c in df.columns if c not in final_cols]
    df         = df[final_cols + extras]

    return df


# ─────────────────────────────────────────────────────────────
# STEP 7.  locus_level builder
# ─────────────────────────────────────────────────────────────

LOCUS_LEVEL_COLUMNS = [
    'Sample_ID', 'Access_Number', 'Locus_ID', 'Locus_Name',
    'N_Alleles_Detected', 'Is_Clean_Monoallelic_Locus',
    'Total_Raw_Count', 'Total_Q_i', 'Total_Q_i_experimental',
    'Dominant_Allele_By_Qi', 'Dominant_Allele_By_Experimental',
    'Dominant_Idealized_Share_%', 'Dominant_Experimental_Share_%',
    'Has_Unobservable_Alleles', 'Incomplete_Denominator',
    'Denominator_Complete', 'Denominator_Issue_Type',
    # Hidden allele interpretation layer
    'Mean_Hidden_Component_Pct', 'Max_Hidden_Component_Pct',
    'Potential_Hidden_Allele', 'Hidden_Allele_Evidence_Score',
]


def build_locus_level(df_allele):
    """Aggregate allele_level into one row per (sample, locus)."""
    rows = []
    for (sid, acc, lid, lname), grp in df_allele.groupby(
            ['Sample_ID', 'Access_Number', 'Locus_ID', 'Locus_Name'],
            sort=False):

        # N_Alleles_Detected: grp size (each row IS a detected allele here,
        # since build_allele_level_raw dropped Raw_Count < MIN_ALLELE_SUPPORT).
        n_alleles = len(grp)

        # Is_Clean_Monoallelic_Locus is constant across grp
        is_clean = bool(grp['Is_Clean_Monoallelic_Locus'].iloc[0])

        total_raw = float(pd.to_numeric(grp['Raw_Count'],
                                        errors='coerce').dropna().sum())

        q_valid   = pd.to_numeric(grp['Q_i'], errors='coerce').dropna()
        total_q   = float(q_valid.sum()) if len(q_valid) > 0 else None

        qe_valid  = pd.to_numeric(grp['Q_i_experimental'],
                                  errors='coerce').dropna()
        total_qe  = float(qe_valid.sum()) if len(qe_valid) > 0 else None

        # Dominant by Q_i
        if len(q_valid) > 0:
            idxmax    = q_valid.idxmax()
            dom_qi_al = grp.loc[idxmax, 'Allele_Value']
            dom_qi_sh = _to_float(grp.loc[idxmax, 'Idealized_Share_%'])
        else:
            dom_qi_al = None
            dom_qi_sh = None

        # Dominant by experimental
        if len(qe_valid) > 0:
            idxmax_e  = qe_valid.idxmax()
            dom_ex_al = grp.loc[idxmax_e, 'Allele_Value']
            dom_ex_sh = _to_float(grp.loc[idxmax_e, 'Experimental_Share_%'])
        else:
            dom_ex_al = None
            dom_ex_sh = None

        has_unobs = bool(grp['Is_Unobservable'].any())

        incompl_any = bool(grp['Incomplete_Denominator'].dropna().astype(bool).any()
                           if grp['Incomplete_Denominator'].notna().any() else False)

        denom_complete  = bool(grp['Denominator_Complete'].iloc[0])
        denom_issue     = str(grp['Denominator_Issue_Type'].iloc[0])

        # ── Hidden allele interpretation ───────────────────────────────────
        hcp_vals = pd.to_numeric(
            grp['Hidden_Component_Pct'] if 'Hidden_Component_Pct' in grp.columns
            else pd.Series(dtype=float), errors='coerce')
        mean_hcp = round(float(hcp_vals.mean()), 2) if hcp_vals.notna().any() else None
        max_hcp  = round(float(hcp_vals.max()),  2) if hcp_vals.notna().any() else None

        # Potential_Hidden_Allele: any of four conditions
        cond_high_hidden    = (hcp_vals > 50).any()  if hcp_vals.notna().any() else False
        cond_unobservable   = has_unobs
        cond_denom_issue    = (denom_issue != 'NONE')
        cond_incomplete     = incompl_any
        potential_hidden    = bool(cond_high_hidden or cond_unobservable
                                   or cond_denom_issue or cond_incomplete)

        # Hidden_Allele_Evidence_Score  ∈ [0, 1]
        #   0.5 × mean(Hidden_Component_Pct)/100
        #   + 0.25 if any allele is fully unobservable
        #   + 0.15 if denominator issue detected
        #   + 0.10 if incomplete denominator flag set
        score  = 0.5 * (float(hcp_vals.mean()) / 100.0) if hcp_vals.notna().any() else 0.0
        if cond_unobservable:  score += 0.25
        if cond_denom_issue:   score += 0.15
        if cond_incomplete:    score += 0.10
        score = round(min(1.0, score), 4)

        rows.append({
            'Sample_ID':                        sid,
            'Access_Number':                    acc,
            'Locus_ID':                         lid,
            'Locus_Name':                       lname,
            'N_Alleles_Detected':               n_alleles,
            'Is_Clean_Monoallelic_Locus':       is_clean,
            'Total_Raw_Count':                  round(total_raw, 4),
            'Total_Q_i':                        round(total_q,  ROUND_Q) if total_q is not None else None,
            'Total_Q_i_experimental':           round(total_qe, ROUND_Q) if total_qe is not None else None,
            'Dominant_Allele_By_Qi':            dom_qi_al,
            'Dominant_Allele_By_Experimental':  dom_ex_al,
            'Dominant_Idealized_Share_%':       round(dom_qi_sh, ROUND_PCT) if dom_qi_sh is not None else None,
            'Dominant_Experimental_Share_%':    round(dom_ex_sh, ROUND_PCT) if dom_ex_sh is not None else None,
            'Has_Unobservable_Alleles':         has_unobs,
            'Incomplete_Denominator':           incompl_any,
            'Denominator_Complete':             denom_complete,
            'Denominator_Issue_Type':           denom_issue,
            'Mean_Hidden_Component_Pct':        mean_hcp,
            'Max_Hidden_Component_Pct':         max_hcp,
            'Potential_Hidden_Allele':          potential_hidden,
            'Hidden_Allele_Evidence_Score':     score,
        })

    df_locus = pd.DataFrame(rows, columns=LOCUS_LEVEL_COLUMNS)

    # Boolean dtypes
    for c in ('Is_Clean_Monoallelic_Locus', 'Has_Unobservable_Alleles',
              'Incomplete_Denominator', 'Denominator_Complete',
              'Potential_Hidden_Allele'):
        if c in df_locus.columns:
            df_locus[c] = df_locus[c].astype('boolean')

    return df_locus


# ─────────────────────────────────────────────────────────────
# STEP 8.  sample_level builder
# ─────────────────────────────────────────────────────────────

SAMPLE_LEVEL_COLUMNS = [
    'Sample_ID', 'Access_Number',
    'N_filt', 'G_mini', 'C_mini',
    'N_Loci_Total', 'N_Clean_Monoallelic_Loci',
    'B_sample', 'Sample_Baseline_Reliable',
    'N_Loci_With_Cohort_Baseline', 'N_Loci_Without_Cohort_Baseline',
    'N_Loci_With_Unobservable_Alleles', 'N_Loci_Incomplete_Denominator',
    # Mixture / hidden-allele detection
    'N_Loci_Potential_Hidden_Allele',
    'Potential_Mixture_Flag',
    'Mixture_Type',
]


def build_sample_level(df_allele, df_locus):
    """Aggregate locus_level + allele_level into one row per sample."""
    rows = []

    # precompute sample-level scalars from allele_level
    allele_by_sample = df_allele.groupby('Sample_ID')

    for sid, grp in df_locus.groupby('Sample_ID', sort=False):
        acc = grp['Access_Number'].dropna().iloc[0] if grp['Access_Number'].notna().any() else sid

        # Pull one row per sample from allele_level for sample-level scalars
        sdf = allele_by_sample.get_group(sid)
        n_filt = _to_float(sdf['N_filt'].dropna().iloc[0]) if sdf['N_filt'].notna().any() else None
        g_mini = _to_float(sdf['G_mini'].dropna().iloc[0]) if sdf['G_mini'].notna().any() else None
        c_mini = _to_float(sdf['C_mini'].dropna().iloc[0]) if sdf['C_mini'].notna().any() else None
        b_samp = _to_float(sdf['B_sample'].dropna().iloc[0]) if sdf['B_sample'].notna().any() else None
        samp_reliable = bool(sdf['Sample_Baseline_Reliable'].dropna().astype(bool).iloc[0]) \
            if sdf['Sample_Baseline_Reliable'].notna().any() else False

        n_loci_total  = int(len(grp))
        n_clean_mono  = int(grp['Is_Clean_Monoallelic_Locus'].fillna(False).astype(bool).sum())

        # Loci with cohort baseline applied: at least one allele row for this
        # (sample, locus) has Experimental_Cohort_Adjustment_Applied == True.
        def _locus_cohort_applied(lid):
            sub = sdf[sdf['Locus_ID'] == lid]
            if sub.empty:
                return False
            v = sub['Experimental_Cohort_Adjustment_Applied'].dropna()
            return bool(v.astype(bool).any()) if len(v) else False

        cohort_applied_flags = grp['Locus_ID'].apply(_locus_cohort_applied)
        n_with_cohort    = int(cohort_applied_flags.sum())
        n_without_cohort = n_loci_total - n_with_cohort

        n_unobs   = int(grp['Has_Unobservable_Alleles'].fillna(False).astype(bool).sum())
        n_incompl = int(grp['Incomplete_Denominator'].fillna(False).astype(bool).sum())

        # ── Mixture / hidden-allele detection ─────────────────────────────
        n_hidden = 0
        if 'Potential_Hidden_Allele' in grp.columns:
            n_hidden = int(grp['Potential_Hidden_Allele'].fillna(False).astype(bool).sum())

        has_observed_multi = bool(grp['N_Alleles_Detected'].gt(1).any())
        has_hidden_loci    = (n_hidden > 0)
        potential_mix      = has_observed_multi or has_hidden_loci

        if has_observed_multi and has_hidden_loci:
            mix_type = 'AMBIGUOUS'
        elif has_observed_multi:
            mix_type = 'OBSERVED_MIX'
        elif has_hidden_loci:
            mix_type = 'HIDDEN_MIX'
        else:
            mix_type = 'CLEAN_SINGLE'

        rows.append({
            'Sample_ID':                        sid,
            'Access_Number':                    acc,
            'N_filt':                           n_filt,
            'G_mini':                           g_mini,
            'C_mini':                           round(c_mini, ROUND_Q) if c_mini is not None else None,
            'N_Loci_Total':                     n_loci_total,
            'N_Clean_Monoallelic_Loci':         n_clean_mono,
            'B_sample':                         round(b_samp, ROUND_Q) if b_samp is not None else None,
            'Sample_Baseline_Reliable':         samp_reliable,
            'N_Loci_With_Cohort_Baseline':      n_with_cohort,
            'N_Loci_Without_Cohort_Baseline':   n_without_cohort,
            'N_Loci_With_Unobservable_Alleles': n_unobs,
            'N_Loci_Incomplete_Denominator':    n_incompl,
            'N_Loci_Potential_Hidden_Allele':   n_hidden,
            'Potential_Mixture_Flag':           potential_mix,
            'Mixture_Type':                     mix_type,
        })

    df_sample = pd.DataFrame(rows, columns=SAMPLE_LEVEL_COLUMNS)
    for c in ('Sample_Baseline_Reliable', 'Potential_Mixture_Flag'):
        if c in df_sample.columns:
            df_sample[c] = df_sample[c].astype('boolean')

    return df_sample


# ─────────────────────────────────────────────────────────────
# STEP 9.  Publication tables  (v23 redesign)
#
#   Sheet 1: major_allele_qc_matrix
#       Sample × Locus pivot.
#       Cell = allele value (integer) if PASS, else status token:
#           ND          — locus not detected in this sample
#           MIXED       — N_Alleles_Detected > 1  (contamination visible)
#           HIDDEN      — Potential_Hidden_Allele=True (monoallelic but at-risk)
#           REJECTED    — Denominator_Complete=False (shares unreliable)
#
#   Sheet 2: major_allele_interpretation
#       One row per (Sample_ID × Locus_ID).  Dominant allele only.
#       Columns: allele value, QC status, raw count, idealized/experimental
#       shares, minor allele summary, contamination fraction + method,
#       denominator flags, hidden allele evidence, human-readable interpretation.
#
#   All calculation functions (STEPS 1–8) are UNCHANGED.
#   Only the output layer (this section) is modified.
# ─────────────────────────────────────────────────────────────


def _try_float(s):
    """Safe numeric cast for sort keys."""
    try:
        return float(s)
    except (ValueError, TypeError):
        return float('inf')


# ── QC status logic ────────────────────────────────────────────

def _final_qc_status(locus_row, n_detected):
    """
    Return one of: PASS | ND | MIXED | HIDDEN | REJECTED

    Priority order (highest to lowest):
        ND       — no allele detected
        MIXED    — more than one allele detected (visible contamination)
        HIDDEN   — monoallelic but Potential_Hidden_Allele = True
        REJECTED — Denominator_Complete = False (shares not trustworthy)
        PASS     — all checks clear
    """
    if locus_row is None:
        return 'ND'
    n = int(n_detected) if n_detected is not None else 0
    if n == 0:
        return 'ND'
    if n > 1:
        return 'MIXED'
    potential_hidden = bool(locus_row.get('Potential_Hidden_Allele') or False)
    denom_complete   = bool(locus_row.get('Denominator_Complete',  True) or False)
    if potential_hidden:
        return 'HIDDEN'
    if not denom_complete:
        return 'REJECTED'
    return 'PASS'


def _dom_allele_str(locus_row):
    """Return dominant allele as clean string (integer if numeric)."""
    val = locus_row.get('Dominant_Allele_By_Experimental')
    if val is None or (isinstance(val, float) and np.isnan(val)):
        val = locus_row.get('Dominant_Allele_By_Qi')
    if val is None:
        return ''
    try:
        return str(int(float(val)))
    except (ValueError, TypeError):
        return str(val)


# ── Sheet 1: major_allele_qc_matrix ───────────────────────────

def build_qc_matrix(df_locus, all_sample_locus_pairs):
    """
    Build Sample × Locus matrix for Sheet 1.

    Parameters
    ----------
    df_locus : DataFrame
        locus_level output (one row per sample × locus, calculated alleles only).
    all_sample_locus_pairs : set of (Sample_ID, Locus_Name)
        All (sample, locus_name) pairs present in input (including Raw_Count=0),
        used to distinguish genuine ND from a locus that was never run.

    Cell encoding
    -------------
    PASS     → dominant allele value (integer string, e.g. "8")
    ND       → "ND"
    MIXED    → "<value> | MIXED"   (dominant allele shown for context)
    HIDDEN   → "<value> | HIDDEN"
    REJECTED → "<value> | REJECTED"
    """
    if df_locus.empty:
        return pd.DataFrame()

    # Stable locus order: sort by locus name (numeric-aware)
    loci = sorted(
        df_locus['Locus_Name'].unique(),
        key=lambda x: (_try_float(str(x).replace('Ft-M', '').replace('Ft-m', '')), str(x))
    )
    samples = sorted(df_locus['Sample_ID'].unique())

    # Build lookup: (Sample_ID, Locus_Name) → locus row
    locus_lookup = {}
    for _, row in df_locus.iterrows():
        locus_lookup[(row['Sample_ID'], str(row['Locus_Name']))] = row

    matrix_rows = []
    for sid in samples:
        row_data = {'Sample_ID': sid}
        for lname in loci:
            key  = (sid, lname)
            lrow = locus_lookup.get(key)

            if lrow is None:
                # Locus was in input (possibly Raw_Count=0) or simply not run
                cell = 'ND'
            else:
                n_det  = lrow.get('N_Alleles_Detected', 0)
                status = _final_qc_status(lrow, n_det)
                dom    = _dom_allele_str(lrow)

                if status == 'PASS':
                    cell = dom if dom else 'ND'
                elif status == 'ND':
                    cell = 'ND'
                else:
                    # Show allele + flag so researcher sees both value and warning
                    cell = f'{dom} | {status}' if dom else status

            row_data[lname] = cell
        matrix_rows.append(row_data)

    return pd.DataFrame(matrix_rows, columns=['Sample_ID'] + loci)


# ── Sheet 2: major_allele_interpretation ──────────────────────

def _contamination_fraction(locus_row):
    """
    Contamination_Fraction_% = 100 − Dominant_Experimental_Share_%
    Fallback: 100 − Dominant_Idealized_Share_%

    Returns (fraction_pct: float|None, method: str)
    """
    exp = _to_float(locus_row.get('Dominant_Experimental_Share_%'))
    if exp is not None:
        return round(100.0 - exp, 4), 'Experimental'
    ideal = _to_float(locus_row.get('Dominant_Idealized_Share_%'))
    if ideal is not None:
        return round(100.0 - ideal, 4), 'Idealized (fallback)'
    return None, 'N/A'


def _contamination_note(method):
    """Human-readable formula note for the contamination column."""
    if method == 'Experimental':
        return ('100% \u2212 Dominant_Experimental_Share_% '
                '(two-step cohort-adjusted model; '
                'Q_i_experimental / \u03a3Q_j_experimental \u00d7 100)')
    if method == 'Idealized (fallback)':
        return ('100% \u2212 Dominant_Idealized_Share_% '
                '(fallback: experimental baseline unavailable; '
                'Q_i / \u03a3Q_j \u00d7 100)')
    return 'Not calculable \u2014 no share data available'


def _interpretation_text(status, contam_frac, contam_method,
                          potential_hidden, hidden_score, n_alleles):
    """
    Generate one concise English sentence per locus for the Interpretation column.
    """
    if status == 'ND':
        return 'Allele not detected in this sample at this locus.'

    parts = []

    if status == 'PASS':
        if contam_frac is not None and contam_frac <= 0.05:
            parts.append('Single dominant allele; locus is clean (contamination fraction \u22640.05%).')
        elif contam_frac is not None and contam_frac < 2.0:
            parts.append(f'Single dominant allele; trace minor signal ({contam_frac:.2f}%) '
                         f'within noise range.')
        elif contam_frac is not None:
            parts.append(f'Single dominant allele; minor signal {contam_frac:.2f}% '
                         f'detected \u2014 may indicate low-level admixture.')
        else:
            parts.append('Single dominant allele; denominator complete.')

    elif status == 'MIXED':
        lvl = ''
        if contam_frac is not None:
            if contam_frac < 1.0:
                lvl = f'trace admixture ({contam_frac:.2f}%)'
            elif contam_frac < 10.0:
                lvl = f'low-level admixture ({contam_frac:.2f}%)'
            else:
                lvl = f'significant admixture ({contam_frac:.1f}%)'
        else:
            lvl = 'admixture fraction not calculable'
        parts.append(f'Multiple alleles detected ({n_alleles} alleles): {lvl}. '
                     f'Possible contamination or mixed infection.')

    elif status == 'HIDDEN':
        parts.append(
            f'Single allele detected, but hidden-allele risk present '
            f'(evidence score {hidden_score:.3f} if hidden_score is not None else "N/A"). '
            f'A biologically present allele with amplicon > read length may be invisible. '
            f'Interpretation requires caution; consider long-read sequencing.')

    elif status == 'REJECTED':
        parts.append(
            'Locus REJECTED: denominator is incomplete '
            '(unobservable/low-completeness alleles present). '
            'Idealized and experimental share fractions are unreliable for this locus.')

    if potential_hidden and status not in ('HIDDEN', 'ND', 'REJECTED'):
        parts.append(
            f'Note: hidden-allele evidence also present '
            f'(score {hidden_score:.3f} if hidden_score is not None else "N/A"). '
            f'Verify with increased read length.')

    return ' '.join(parts) if parts else 'See QC columns for details.'


def build_interpretation_sheet(df_allele, df_locus):
    """
    Build major_allele_interpretation sheet (Sheet 2).

    One row per (Sample_ID, Locus_ID).
    Columns describe the dominant allele with all QC metrics and
    a human-readable Interpretation string.
    """
    if df_locus.empty:
        return pd.DataFrame()

    rows = []
    for _, lrow in df_locus.iterrows():
        sid   = lrow['Sample_ID']
        lid   = lrow['Locus_ID']
        lname = lrow.get('Locus_Name', '')
        n_det = int(lrow.get('N_Alleles_Detected', 0) or 0)

        status       = _final_qc_status(lrow, n_det)
        dom_val      = _dom_allele_str(lrow)
        contam_frac, contam_method = _contamination_fraction(lrow)
        potential_h  = bool(lrow.get('Potential_Hidden_Allele') or False)
        hidden_score = _to_float(lrow.get('Hidden_Allele_Evidence_Score'))
        max_hcp      = _to_float(lrow.get('Max_Hidden_Component_Pct'))
        denom_comp   = bool(lrow.get('Denominator_Complete', True) or False)
        denom_issue  = str(lrow.get('Denominator_Issue_Type', 'NONE'))
        dom_ideal    = _to_float(lrow.get('Dominant_Idealized_Share_%'))
        dom_exp      = _to_float(lrow.get('Dominant_Experimental_Share_%'))

        # Raw count for the dominant allele (from allele_level)
        dom_rows = df_allele[
            (df_allele['Sample_ID'] == sid) &
            (df_allele['Locus_ID']  == lid) &
            (df_allele['Allele_Value'].astype(str) == dom_val)
        ]
        dom_raw = (int(float(dom_rows['Raw_Count'].iloc[0]))
                   if not dom_rows.empty else None)

        # Per-allele idealized / experimental shares (from allele_level, more precise)
        if not dom_rows.empty:
            dom_ideal = _to_float(dom_rows['Idealized_Share_%'].iloc[0]) or dom_ideal
            dom_exp   = _to_float(dom_rows['Experimental_Share_%'].iloc[0]) or dom_exp

        # Minor allele summary string
        all_loc = df_allele[
            (df_allele['Sample_ID'] == sid) & (df_allele['Locus_ID'] == lid)
        ]
        minor_rows = all_loc[all_loc['Allele_Value'].astype(str) != dom_val]
        if not minor_rows.empty:
            minor_parts = []
            for _, mr in minor_rows.iterrows():
                mv  = mr['Allele_Value']
                msh = (_to_float(mr.get('Experimental_Share_%')) or
                       _to_float(mr.get('Idealized_Share_%')))
                msh_str = f'{msh:.2f}%' if msh is not None else 'N/A'
                minor_parts.append(f'{mv} ({msh_str})')
            minor_summary = '; '.join(minor_parts)
        else:
            minor_summary = 'None'

        # Contamination note (formula explanation)
        contam_note = _contamination_note(contam_method)

        # Human-readable interpretation
        interp = _interpretation_text(
            status, contam_frac, contam_method,
            potential_h, hidden_score, n_det
        )

        rows.append({
            'Sample_ID':                     sid,
            'Locus_ID':                      lid,
            'Locus_Name':                    lname,
            'Major_Allele':                  dom_val if dom_val else 'ND',
            'Final_QC_Status':               status,
            'Raw_Count_Major':               dom_raw,
            'N_Alleles_Detected':            n_det,
            'Idealized_Share_%':             (round(dom_ideal, ROUND_PCT)
                                              if dom_ideal is not None else None),
            'Experimental_Share_%':          (round(dom_exp, ROUND_PCT)
                                              if dom_exp is not None else None),
            'Minor_Alleles_Summary':         minor_summary,
            'Contamination_Fraction_%':      contam_frac,
            'Contamination_Calc_Method':     contam_note,
            'Denominator_Complete':          denom_comp,
            'Denominator_Issue_Type':        denom_issue,
            'Potential_Hidden_Allele':       potential_h,
            'Hidden_Allele_Evidence_Score':  (round(hidden_score, 4)
                                              if hidden_score is not None else None),
            'Max_Hidden_Component_%':        (round(max_hcp, 2)
                                              if max_hcp is not None else None),
            'Interpretation':                interp,
        })

    INTERP_COLS = [
        'Sample_ID', 'Locus_ID', 'Locus_Name',
        'Major_Allele', 'Final_QC_Status',
        'Raw_Count_Major', 'N_Alleles_Detected',
        'Idealized_Share_%', 'Experimental_Share_%',
        'Minor_Alleles_Summary',
        'Contamination_Fraction_%', 'Contamination_Calc_Method',
        'Denominator_Complete', 'Denominator_Issue_Type',
        'Potential_Hidden_Allele', 'Hidden_Allele_Evidence_Score',
        'Max_Hidden_Component_%',
        'Interpretation',
    ]
    return pd.DataFrame(rows, columns=INTERP_COLS)


# ── Entry point: build_publication_tables ─────────────────────

def build_publication_tables(df_allele, df_locus, output_path,
                              df_all_alleles_full=None):
    """
    Write publication_tables.xlsx  (v23 — two-sheet design).

    Sheet 1 — major_allele_qc_matrix
        Sample × Locus matrix. Each cell shows the dominant allele value
        (integer) or a status token: ND / MIXED / HIDDEN / REJECTED.
        Ready for copy-paste into a manuscript table.

    Sheet 2 — major_allele_interpretation
        One row per (Sample_ID × Locus_ID).  Contains all QC metrics,
        idealized and experimental shares, contamination fraction with
        formula annotation, hidden-allele evidence, and a plain-English
        Interpretation sentence.

    IMPORTANT: This function is the ONLY changed part of the pipeline.
    All calculation steps (build_allele_level_raw → apply_idealized_model
    → apply_experimental_model → denominator completeness → hidden allele
    flags → build_locus_level) are UNCHANGED.
    """
    # Collect all (sample, locus_name) pairs seen in the full allele list
    # (includes Raw_Count=0 rows) for robust ND detection.
    all_sample_locus_pairs = set()
    if df_all_alleles_full is not None and not df_all_alleles_full.empty:
        for _, r in df_all_alleles_full.iterrows():
            all_sample_locus_pairs.add((r['Sample_ID'], str(r.get('Locus_Name', ''))))
    for _, r in df_allele.iterrows():
        all_sample_locus_pairs.add((r['Sample_ID'], str(r.get('Locus_Name', ''))))

    df_matrix = build_qc_matrix(df_locus, all_sample_locus_pairs)
    df_interp = build_interpretation_sheet(df_allele, df_locus)

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        if not df_matrix.empty:
            df_matrix.to_excel(writer,
                               sheet_name='major_allele_qc_matrix',
                               index=False)
        if not df_interp.empty:
            df_interp.to_excel(writer,
                               sheet_name='major_allele_interpretation',
                               index=False)

    _beautify_pub_workbook(output_path)
    logger.info(f'Saved publication tables -> {output_path}')


# ── Beautification for publication workbook ────────────────────

# Status colour map (Sheet 1 and Sheet 2 Final_QC_Status column)
_QC_FILLS = {
    'PASS':     PatternFill('solid', fgColor='C6EFCE'),   # green
    'ND':       PatternFill('solid', fgColor='D9D9D9'),   # grey
    'MIXED':    PatternFill('solid', fgColor='F4B084'),   # orange-red
    'HIDDEN':   PatternFill('solid', fgColor='FFE699'),   # amber
    'REJECTED': PatternFill('solid', fgColor='FF6B6B'),   # red
}
_QC_FONTS = {
    'PASS':     Font(color='276221', bold=False),
    'ND':       Font(color='595959', bold=False),
    'MIXED':    Font(color='843C0C', bold=True),
    'HIDDEN':   Font(color='7F5A00', bold=True),
    'REJECTED': Font(color='FFFFFF', bold=True),
}


def _color_qc_matrix_sheet(ws):
    """
    Apply per-cell colouring to Sheet 1 (major_allele_qc_matrix).
    Cells containing status tokens are coloured; numeric cells get green.
    """
    thin = Side(style='thin', color='BFBFBF')
    border = Border(left=thin, right=thin, top=thin, bottom=thin)
    center = Alignment(horizontal='center', vertical='center')

    for row in ws.iter_rows(min_row=2):
        for cell in row:
            if cell.column == 1:   # Sample_ID column — skip
                continue
            val = str(cell.value) if cell.value is not None else ''
            cell.alignment = center
            cell.border = border

            # Detect status token (may be "8 | MIXED" or plain "MIXED")
            token = None
            for tok in ('MIXED', 'HIDDEN', 'REJECTED', 'ND', 'PASS'):
                if tok in val:
                    token = tok
                    break
            # Pure numeric → PASS
            if token is None:
                try:
                    float(val)
                    token = 'PASS'
                except (ValueError, TypeError):
                    pass

            if token and token in _QC_FILLS:
                cell.fill = _QC_FILLS[token]
                cell.font = _QC_FONTS[token]


def _color_interp_status_column(ws):
    """
    Apply cell colours to the Final_QC_Status column in Sheet 2.
    """
    headers = {str(c.value): c.column for c in ws[1] if c.value is not None}
    col = headers.get('Final_QC_Status')
    if col is None:
        return
    ltr = get_column_letter(col)
    thin = Side(style='thin', color='BFBFBF')
    border = Border(left=thin, right=thin, top=thin, bottom=thin)
    center = Alignment(horizontal='center', vertical='center')
    for row in ws.iter_rows(min_row=2, min_col=col, max_col=col):
        cell = row[0]
        val  = str(cell.value) if cell.value is not None else ''
        cell.alignment = center
        cell.border    = border
        if val in _QC_FILLS:
            cell.fill = _QC_FILLS[val]
            cell.font = _QC_FONTS[val]


def _beautify_pub_workbook(path):
    """
    Apply styling to publication_tables.xlsx:
    • Standard header + autofit on all sheets
    • Per-cell QC colour coding on Sheet 1
    • Status column colour on Sheet 2
    • Wrap text + wider Interpretation column on Sheet 2
    """
    from openpyxl import load_workbook
    wb = load_workbook(path)

    for ws in wb.worksheets:
        _style_header(ws, fill_color='1F4E78')
        _auto_fit_columns(ws, max_width=38)
        _apply_numeric_formats(ws)

        if ws.title == 'major_allele_qc_matrix':
            _color_qc_matrix_sheet(ws)
            # Freeze Sample_ID column + header row
            ws.freeze_panes = 'B2'

        elif ws.title == 'major_allele_interpretation':
            _style_issue_columns(ws)
            _color_interp_status_column(ws)
            # Widen Interpretation column and wrap text
            headers = {str(c.value): c.column for c in ws[1] if c.value is not None}
            interp_col = headers.get('Interpretation')
            if interp_col:
                ltr = get_column_letter(interp_col)
                ws.column_dimensions[ltr].width = 60
                for cell in ws[ltr]:
                    cell.alignment = Alignment(wrap_text=True,
                                               vertical='top')
            # Widen Contamination_Calc_Method column
            calc_col = headers.get('Contamination_Calc_Method')
            if calc_col:
                ws.column_dimensions[get_column_letter(calc_col)].width = 52
            ws.freeze_panes = 'A2'

    wb.save(path)


# ─────────────────────────────────────────────────────────────
# STEP 10.  Output dictionary (definition of every column)
# ─────────────────────────────────────────────────────────────

_DICT_COLS = ['Column_Name', 'Sheet', 'Definition', 'Formula']

_DICT_ROWS = [
    # Identification
    ('Sample_ID',               'allele,locus,sample', 'Sample identifier',                              'from input'),
    ('Access_Number',           'allele,locus,sample', 'FASTQ accession number',                         'from input'),
    ('Locus_ID',                'allele,locus',        'Full locus identifier',                          'from input (Primer column)'),
    ('Locus_Name',              'allele,locus',        'Short locus name',                               "split('_')[0] of Locus_ID"),
    ('Repeat_Descriptor',       'allele',              'Repeat unit descriptor',                         "split('_')[1] of Locus_ID"),
    ('Reference_Fragment_Label','allele',              'Reference fragment length label',                "split('_')[2] of Locus_ID"),
    ('Reference_Allele_Label',  'allele',              'Reference allele label',                         "split('_')[3] of Locus_ID"),
    ('Allele_Rank_Input',       'allele',              'Ordinal position of allele in wide-format input','direct'),
    ('Allele_Value',            'allele,locus',        'Allele value (repeat units)',                    'from input'),

    # Raw data
    ('Raw_Count',          'allele,locus', 'Number of reads assigned to this allele',                    'direct from input'),
    ('Fragment_Length_bp', 'allele',       'Amplicon length in bp',                                       'from input (Loci_Size)'),
    ('Mean_Read_Length_bp','allele',       'Mean sequencing read length for the sample in bp',           'from input (Mean_Read_Length)'),
    ('p_i',                'allele',       'Per-read detectability probability',                         'p_star if present else max(0, (R - L + 1)/R); 0 if L>R'),
    ('Detectability_Class','allele',       'Observable / Marginal / Effectively_Unobservable / Unobservable', 'from input'),

    # Mini-genome block
    ('N_filt', 'allele,sample', 'Sum of Raw_Count over all detected alleles of the sample',    'Σ Raw_Count per sample'),
    ('G_mini', 'allele,sample', 'Sum of Fragment_Length over all detected alleles',            'Σ Fragment_Length_bp per sample'),
    ('C_mini', 'allele,sample', 'Mini-genome coverage',                                        'C_mini = N_filt · R / G_mini'),

    # Idealized block
    ('E_ideal',            'allele', 'Idealized expected reads (no length multiplier)', 'E_ideal = C_mini · p_i'),
    ('Q_i',                'allele', 'Idealized ratio',                                 'Q_i = Raw_Count / E_ideal'),
    ('Idealized_Share_%',  'allele', 'Allele share within locus by idealized model',    'Q_i / Σ(Q_j in locus) · 100'),

    # Sample baseline block
    ('Is_Clean_Monoallelic_Locus','allele,locus', 'Locus has exactly 1 detected allele with Raw_Count>0, p_i>0, not unobservable, not in dropout set', 'combined flag'),
    ('B_sample',                  'allele,sample', 'Sample baseline',                                    'median(Q_i over clean monoallelic loci of sample)'),
    ('Sample_Baseline_Reliable',  'allele,sample', f'True if ≥{MIN_CLEAN_MONOALLELIC_FOR_SAMPLE} clean monoallelic loci contributed', f'N_clean ≥ {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE}'),
    ('Q_i_sample_adjusted',       'allele', 'Sample-adjusted idealized ratio',                           'Q_i / B_sample'),

    # Cohort locus block
    ('B_locus_cohort',                'allele', 'Cohort locus baseline',                                 'median(Q_i_sample_adjusted over clean monoallelic observations of the locus across samples)'),
    ('Locus_Cohort_Baseline_Reliable','allele', f'True if ≥{MIN_CLEAN_MONOALLELIC_FOR_LOCUS} observations', f'N_obs ≥ {MIN_CLEAN_MONOALLELIC_FOR_LOCUS}'),
    ('Experimental_Cohort_Adjustment_Applied','allele', 'True if cohort baseline was applied (not fallback)', 'B_locus_cohort present and positive'),
    ('Q_i_experimental',              'allele', 'Final experimental ratio',                              'Q_i_sample_adjusted / B_locus_cohort; fallback = Q_i_sample_adjusted'),
    ('Experimental_Share_%',          'allele', 'Allele share within locus by experimental model',       'Q_i_experimental / Σ(..) · 100'),

    # QC / completeness
    ('Is_Unobservable',         'allele',       'True if p_i≤0 or Detectability_Class=="Unobservable"',  'derived'),
    ('Incomplete_Denominator',  'allele,locus', 'True if hidden allele with significant expected signal', 'from input'),
    ('Denominator_Complete',    'allele,locus', 'True if no denominator issue (UNOBSERVABLE/LOW_COMPLETENESS)', 'issue_type == NONE'),
    ('Denominator_Issue_Type',  'allele,locus', 'NONE / UNOBSERVABLE / LOW_COMPLETENESS / MIXED',        'rule: unobs + (incompl OR hidden>threshold)'),
    ('Dropout_Adjusted_Status', 'allele',       'Locus dropout-aware status',                            'from input'),

    # locus_level extras
    ('N_Alleles_Detected',               'locus', 'Number of detected alleles at locus in sample', 'grp size'),
    ('Total_Raw_Count',                  'locus', 'Sum of Raw_Count across locus alleles',         'Σ Raw_Count'),
    ('Total_Q_i',                        'locus', 'Sum of Q_i across locus alleles',               'Σ Q_i'),
    ('Total_Q_i_experimental',           'locus', 'Sum of Q_i_experimental across locus alleles',  'Σ Q_i_experimental'),
    ('Dominant_Allele_By_Qi',            'locus', 'Allele with max Q_i',                           'argmax Q_i'),
    ('Dominant_Allele_By_Experimental',  'locus', 'Allele with max Q_i_experimental',              'argmax Q_i_experimental'),
    ('Dominant_Idealized_Share_%',       'locus', 'Idealized_Share_% of dominant-by-Q_i allele',   'lookup'),
    ('Dominant_Experimental_Share_%',    'locus', 'Experimental_Share_% of dominant-by-Exp allele','lookup'),
    ('Has_Unobservable_Alleles',         'locus', 'True if any allele Is_Unobservable',            'any(Is_Unobservable)'),

    # Hidden allele interpretation layer (allele_level)
    ('Hidden_Allele_Flag',         'allele', 'True when allele is physically undetectable (p_i ≤ 0 or Detectability_Class = Unobservable). Biologically: allele may exist but is invisible to sequencing reads due to amplicon length > read length.',
     'p_i ≤ 0 OR Detectability_Class == Unobservable'),
    ('Hidden_Component_Pct',       'allele', 'Percentage of allele signal lost due to detectability constraints: (1 − p_i) × 100. Range 0–100. 0 = fully observed; 100 = completely undetectable. NOT used in Q_i calculation — QC/interpretation only.',
     '(1 − p_i) × 100; NaN p_i → 100'),

    # Hidden allele interpretation layer (locus_level)
    ('Mean_Hidden_Component_Pct',       'locus', 'Mean Hidden_Component_Pct over all alleles in the locus', 'mean(Hidden_Component_Pct)'),
    ('Max_Hidden_Component_Pct',        'locus', 'Max Hidden_Component_Pct over all alleles in the locus',  'max(Hidden_Component_Pct)'),
    ('Potential_Hidden_Allele',         'locus', 'True if locus shows evidence of one or more hidden (undetectable) alleles. Triggered by: Max_Hidden_Component_Pct > 50, OR Has_Unobservable_Alleles, OR Denominator_Issue_Type != NONE, OR Incomplete_Denominator.',
     'cond_high_hidden OR cond_unobservable OR cond_denom_issue OR cond_incomplete'),
    ('Hidden_Allele_Evidence_Score',    'locus', 'Composite score 0–1 estimating probability that a hidden allele is present. Score = 0.5×mean(HCP/100) + 0.25 if unobservable + 0.15 if denom issue + 0.10 if incomplete denominator.',
     '0.5×mean_hcp/100 + 0.25×unobs + 0.15×denom_issue + 0.10×incomplete; capped at 1.0'),

    # Mixture detection (sample_level)
    ('N_Loci_Potential_Hidden_Allele',  'sample', 'Number of loci with Potential_Hidden_Allele = True', 'count(Potential_Hidden_Allele)'),
    ('Potential_Mixture_Flag',          'sample', 'True if the sample shows evidence of mixture — either observed multi-allelic loci OR loci with Potential_Hidden_Allele. Critical: False does NOT guarantee a pure sample.',
     'has_observed_multi OR N_Loci_Potential_Hidden_Allele > 0'),
    ('Mixture_Type',                    'sample', 'Classification of mixture evidence: CLEAN_SINGLE = no multi-allelic and no hidden allele evidence; OBSERVED_MIX = visible multi-allelic loci detected; HIDDEN_MIX = no visible multi-allelic but hidden allele evidence present; AMBIGUOUS = both observed and hidden signals present.',
     'classification rule on has_observed_multi + has_hidden_loci'),

    # sample_level extras
    ('N_Loci_Total',                    'sample', 'Total number of loci reported for the sample',  'count'),
    ('N_Clean_Monoallelic_Loci',        'sample', 'Clean monoallelic loci count',                  'Σ Is_Clean_Monoallelic_Locus'),
    ('N_Loci_With_Cohort_Baseline',     'sample', 'Loci where cohort baseline applied',            'Σ cohort_applied'),
    ('N_Loci_Without_Cohort_Baseline',  'sample', 'Loci with sample-only fallback',                'total − with_cohort'),
    ('N_Loci_With_Unobservable_Alleles','sample', 'Loci carrying any unobservable allele',          'Σ Has_Unobservable_Alleles'),
    ('N_Loci_Incomplete_Denominator',   'sample', 'Loci with Incomplete_Denominator flag',         'Σ Incomplete_Denominator'),
]


def build_output_dictionary():
    return pd.DataFrame(_DICT_ROWS, columns=_DICT_COLS)


# ─────────────────────────────────────────────────────────────
# STEP 11.  run_metadata.txt
# ─────────────────────────────────────────────────────────────

def write_run_metadata(output_path, input_file, args, run_ts, n_samples, n_loci,
                        n_clean_mono_samples, n_loci_cohort_ok):
    lines = [
        '=' * 64,
        'MLVA/VNTR Pipeline v21 — Run Metadata',
        '=' * 64,
        f'Pipeline version               : {__version__}',
        f'Run date/time                  : {run_ts}',
        f'Input file                     : {os.path.abspath(input_file)}',
        f'Output directory               : {os.path.abspath(args.output)}',
        f'Max alleles per locus          : {args.max_alleles}',
        f'Samples processed              : {n_samples}',
        f'Sample-locus records           : {n_loci}',
        f'Samples with reliable B_sample : {n_clean_mono_samples}',
        f'Loci with reliable B_locus_cohort : {n_loci_cohort_ok}',
        '',
        'Thresholds:',
        f'  MIN_ALLELE_SUPPORT                    = {MIN_ALLELE_SUPPORT}',
        f'  MIN_CLEAN_MONOALLELIC_FOR_SAMPLE      = {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE}',
        f'  MIN_CLEAN_MONOALLELIC_FOR_LOCUS       = {MIN_CLEAN_MONOALLELIC_FOR_LOCUS}',
        f'  HIDDEN_COMPONENT_THRESHOLD            = {HIDDEN_COMPONENT_THRESHOLD}',
        f'  EXCLUDE_DROPOUT_STATUSES              = {sorted(EXCLUDE_DROPOUT_STATUSES)}',
        '',
        'Formulas (v21):',
        '  N_filt   definition            = sum(Raw_Count over all detected alleles of sample)',
        '  G_mini   definition            = sum(Fragment_Length_bp over all detected alleles of sample)',
        '  C_mini   definition            = N_filt · R / G_mini, where R = Mean_Read_Length_bp',
        '  p_i      definition            = p_star if available else max(0, (R - L + 1)/R); 0 if L>R',
        '  E_ideal  formula               = C_mini · p_i      (NO length multiplier)',
        '  Q_i      formula               = Raw_Count / E_ideal',
        '  B_sample formula               = median(Q_i over clean monoallelic loci of sample)',
        '  Q_i_sample_adjusted formula    = Q_i / B_sample',
        '  B_locus_cohort formula         = median(Q_i_sample_adjusted over clean monoallelic '
        'observations of same locus across samples)',
        '  Q_i_experimental formula       = Q_i_sample_adjusted / B_locus_cohort'
        '  (fallback: Q_i_experimental = Q_i_sample_adjusted if B_locus_cohort missing)',
        '  Idealized_Share_%  formula     = Q_i / Σ(Q_j in locus) · 100',
        '  Experimental_Share_% formula   = Q_i_experimental / Σ(..) · 100',
        '',
        'Clean monoallelic locus definition:',
        '  Exactly 1 detected allele in the locus for the sample, AND',
        '  Raw_Count > 0 AND p_i > 0 AND Is_Unobservable == False AND',
        '  Dropout_Adjusted_Status NOT IN EXCLUDE_DROPOUT_STATUSES.',
        '',
        'Removed in this version vs. v20:',
        '  - Q_e (experimental ratio via old per-locus B_locus)',
        '  - Old B_locus (median Q_i within locus used as bias)',
        '  - Integrated_Mass, Integrated_Percent, Sample_Percent',
        '  - Major_Allele_Reliability_Pct, Minor_to_Major_Ratio',
        '  - Legacy Idealized_Support = Raw / p_star',
        '  - All plot-ready columns (Plot_*) and plots_experimental_idealized.html',
        '',
        'Output files:',
        '  master_allele_to_sample_statistics.xlsx',
        '  publication_tables.xlsx',
        '  output_dictionary.tsv',
        '  run_metadata.txt',
        '  interpretation_notes.txt',
        '  hidden_allele_interpretation.txt',
        '  plots_idealized_share.html          (per-sample idealized, requires plotly)',
        '  all_samples_idealized_alleles.html  (all-samples idealized, requires plotly)',
        '  plots_experimental_share.html       (per-sample experimental, requires plotly)',
        '  all_samples_experimental_alleles.html (all-samples experimental, requires plotly)',
        '=' * 64,
    ]
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    logger.info(f'Saved run metadata -> {output_path}')


# ─────────────────────────────────────────────────────────────
# STEP 12.  Excel beautification
# (lifted from v20 styling helpers, trimmed — status colors kept
#  for any remaining status columns we may add later.)
# ─────────────────────────────────────────────────────────────

def _auto_fit_columns(ws, max_width=42):
    for col_cells in ws.columns:
        letter = get_column_letter(col_cells[0].column)
        values = ["" if c.value is None else str(c.value) for c in col_cells]
        width  = min(max(len(v) for v in values) + 2, max_width) if values else 12
        ws.column_dimensions[letter].width = max(width, 10)


def _style_header(ws, fill_color='1F4E78'):
    fill   = PatternFill('solid', fgColor=fill_color)
    font   = Font(color='FFFFFF', bold=True)
    align  = Alignment(horizontal='center', vertical='center', wrap_text=True)
    thin   = Side(style='thin', color='D9E2F3')
    border = Border(left=thin, right=thin, top=thin, bottom=thin)
    for cell in ws[1]:
        cell.fill = fill; cell.font = font
        cell.alignment = align; cell.border = border
    ws.freeze_panes = 'A2'
    if ws.max_row >= 1 and ws.max_column >= 1:
        ws.auto_filter.ref = ws.dimensions


def _apply_numeric_formats(ws):
    percent_kw = ('Percent', 'Share_%', 'Fraction', 'Coverage')
    float_kw   = ('C_mini', 'E_ideal', 'Q_i', 'B_sample', 'B_locus_cohort',
                  'Q_i_sample_adjusted', 'Q_i_experimental', 'Total_Q_i',
                  'p_i', 'Confidence_Score')
    int_kw     = ('N_filt', 'G_mini', 'N_Alleles', 'N_Loci',
                  'N_Clean', 'Allele_Rank_Input', 'Fragment_Length_bp',
                  'Mean_Read_Length_bp', 'Raw_Count', 'Total_Raw_Count')
    headers = [c.value for c in ws[1]]
    for idx, header in enumerate(headers, start=1):
        if header is None:
            continue
        h = str(header)
        if any(k in h for k in percent_kw):
            fmt = '0.0000'
        elif any(k == h for k in float_kw) or any(h.startswith(k) for k in float_kw):
            fmt = '0.000000'
        elif any(k == h for k in int_kw) or any(h.startswith(k) for k in int_kw):
            fmt = '0'
        else:
            continue
        for row in ws.iter_rows(min_row=2, min_col=idx, max_col=idx):
            row[0].number_format = fmt


def _style_issue_columns(ws):
    """Color Denominator_Issue_Type cells."""
    headers = {str(c.value): c.column for c in ws[1] if c.value is not None}
    if 'Denominator_Issue_Type' not in headers:
        return
    col = headers['Denominator_Issue_Type']
    letter = get_column_letter(col)
    if ws.max_row is None or ws.max_row < 2:
        return
    rng = f'{letter}2:{letter}{ws.max_row}'
    fills = {
        'NONE':             PatternFill('solid', fgColor='C6E0B4'),
        'UNOBSERVABLE':     PatternFill('solid', fgColor='FFE699'),
        'LOW_COMPLETENESS': PatternFill('solid', fgColor='FFE699'),
        'MIXED':            PatternFill('solid', fgColor='F4B084'),
    }
    for text_val, fill in fills.items():
        try:
            ws.conditional_formatting.add(
                rng, FormulaRule(formula=[f'${letter}2="{text_val}"'],
                                 stopIfTrue=False, fill=fill))
        except (TypeError, ValueError):
            pass


def _beautify_sheet(ws, max_width=40):
    _style_header(ws)
    _auto_fit_columns(ws, max_width=max_width)
    _apply_numeric_formats(ws)
    _style_issue_columns(ws)


def _beautify_workbook(path):
    from openpyxl import load_workbook
    wb = load_workbook(path)
    for ws in wb.worksheets:
        _beautify_sheet(ws, max_width=40)
    wb.save(path)


# ─────────────────────────────────────────────────────────────
# STEP 13.  Master Excel
# ─────────────────────────────────────────────────────────────



# ─────────────────────────────────────────────────────────────
# STEP 13b.  Column description sheets for master workbook (specification section 8)
# ─────────────────────────────────────────────────────────────

_DESC_COLS = ['Column_Name', 'Definition', 'Formula', 'Important_for_Researcher']

_ALLELE_LEVEL_DESC = [
    ('Sample_ID',              'Sample identifier',                                      'from input file',                 ''),
    ('Access_Number',          'FASTQ file name or accession number',                     'from input file',                 ''),
    ('Locus_ID',               'Full locus identifier (name_repeat_length_reference)',    'from input file (Primer)',        ''),
    ('Locus_Name',             'Short locus name',                                        "split('_')[0] of Locus_ID",       ''),
    ('Repeat_Descriptor',      'Repeat unit size descriptor, for example 6bp',            "split('_')[1] of Locus_ID",       ''),
    ('Reference_Fragment_Label','Expected reference fragment length in bp',               "split('_')[2] of Locus_ID",       ''),
    ('Reference_Allele_Label', 'Repeat count in the reference allele',                    "split('_')[3] of Locus_ID",       ''),
    ('Allele_Rank_Input',      'Ordinal allele number in the original wide-format input', 'direct match to Allele {i}',      ''),
    ('Allele_Value',           'Allele value in repeat units (u); final genotype call',   'computed by MLVA_v2',             'Main genotyping result'),
    ('Raw_Count',              'Uncorrected number of reads assigned to this allele',     'direct counter',                  'Not corrected for fragment length'),
    ('Fragment_Length_bp',     'Amplicon length in bp; determines detectability',         'from MLVA_v2',                    'Values near Mean_Read_Length_bp can become Marginal or Unobservable'),
    ('Mean_Read_Length_bp',    'Mean read length for the sample in bp',                   'total_bases / total_reads',       'One value per sample'),
    ('p_i',                    'Per-read detection probability (0-1)',                    'p_star if present; otherwise (R-L+1)/R; 0 if L>R', 'Key model parameter'),
    ('Detectability_Class',    'Observable / Marginal / Effectively_Unobservable / Unobservable', 'from MLVA_v2',            'Only Observable and Marginal map to Use_in_Idealized=True'),
    ('N_filt',                 'Sum of Raw_Count across all alleles in the sample',       'Σ Raw_Count per sample',          ''),
    ('G_mini',                 'Sum of Fragment_Length across all alleles in the sample', 'Σ Fragment_Length_bp per sample', ''),
    ('C_mini',                 'Mini-genome coverage',                                    'C_mini = N_filt · R / G_mini',    ''),
    ('E_ideal',                'Expected reads under the idealized model',                'E_ideal = C_mini · p_i',          ''),
    ('Q_i',                    'Idealized ratio',                                         'Q_i = Raw_Count / E_ideal',       'Main normalized metric'),
    ('Idealized_Share_%',      'Allele share within the locus under the idealized model (%)', 'Q_i / Σ(Q_j) · 100',          'Not a reliability score; this is a share only'),
    ('Is_Clean_Monoallelic_Locus','True when exactly one allele is detected, p_i>0, not unobservable, and not dropout', 'combined flag', 'Most reliable loci for calibration'),
    ('B_sample',               'Sample baseline',                                         'median(Q_i over clean monoallelic loci)', ''),
    ('Sample_Baseline_Reliable',f'True if at least {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE} clean monoallelic loci contributed', f'N_clean ≥ {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE}', ''),
    ('Q_i_sample_adjusted',    'Q_i adjusted by the sample baseline',                     'Q_i / B_sample',                  ''),
    ('B_locus_cohort',         'Cohort locus baseline',                                   'median(Q_i_sample_adjusted) across the cohort', ''),
    ('Locus_Cohort_Baseline_Reliable',f'True if at least {MIN_CLEAN_MONOALLELIC_FOR_LOCUS} observations contributed', f'N_obs ≥ {MIN_CLEAN_MONOALLELIC_FOR_LOCUS}', ''),
    ('Experimental_Cohort_Adjustment_Applied','True when cohort calibration was applied', 'B_locus_cohort is available and > 0', ''),
    ('Q_i_experimental',       'Final experimental ratio',                                'Q_i_sample_adjusted / B_locus_cohort; fallback = Q_i_sample_adjusted', ''),
    ('Experimental_Share_%',   'Allele share within the locus under the experimental model (%)', 'Q_i_experimental / Σ(..) · 100', 'Not a reliability score; this is a share only'),
    ('Is_Unobservable',        'True if p_i ≤ 0 or Detectability_Class is Unobservable',  'derived',                         ''),
    ('Incomplete_Denominator', 'True when a physically unobservable allele has meaningful expected signal', 'from MLVA_v2', 'Strong warning; Confidence_Score is penalized by 0.5'),
    ('Denominator_Complete',   'True when no denominator problem was detected',           'Denominator_Issue_Type == NONE',  'Used by reliable_alleles_by_sample filters'),
    ('Denominator_Issue_Type', 'NONE / UNOBSERVABLE / LOW_COMPLETENESS / MIXED',          'rule: unobservable + (incomplete OR hidden>threshold)', ''),
    ('Dropout_Adjusted_Status','Locus status after detectability and dropout adjustment', 'from MLVA_v2',                    ''),
    ('Hidden_Allele_Flag',     'True when the allele is physically unobservable (p_i ≤ 0 or Unobservable). This is not proof that the allele is present.', 'p_i ≤ 0 OR Detectability_Class == Unobservable', 'Observability limitation only; not a genotype call'),
    ('Hidden_Component_Pct',   'Percentage of allele signal hidden by detectability constraints. Not used in Q_i.', '(1 − p_i) × 100; NaN p_i → 100', 'QC and interpretation layer only'),
]

_LOCUS_LEVEL_DESC = [
    ('Sample_ID',                    'Sample identifier',                                      'from allele_level', ''),
    ('Access_Number',                'FASTQ file name',                                        'from allele_level', ''),
    ('Locus_ID',                     'Full locus identifier',                                  'from allele_level', ''),
    ('Locus_Name',                   'Short locus name',                                       'from allele_level', ''),
    ('N_Alleles_Detected',           'Number of detected alleles at the locus',                'group size',        ''),
    ('Is_Clean_Monoallelic_Locus',   'True when the locus has exactly one clean allele',       'flag from allele_level', 'Most reliable loci'),
    ('Total_Raw_Count',              'Sum of Raw_Count across locus alleles',                  'Σ Raw_Count',       ''),
    ('Total_Q_i',                    'Sum of Q_i across locus alleles',                        'Σ Q_i',             ''),
    ('Total_Q_i_experimental',       'Sum of Q_i_experimental across locus alleles',           'Σ Q_i_experimental', ''),
    ('Dominant_Allele_By_Qi',        'Allele with the highest Q_i',                            'argmax Q_i',        ''),
    ('Dominant_Allele_By_Experimental','Allele with the highest Q_i_experimental',             'argmax Q_i_experimental', ''),
    ('Dominant_Idealized_Share_%',   'Idealized_Share_% for the dominant allele',              'lookup',            ''),
    ('Dominant_Experimental_Share_%','Experimental_Share_% for the dominant allele',           'lookup',            ''),
    ('Has_Unobservable_Alleles',     'True if any allele has Is_Unobservable=True',            'any(Is_Unobservable)', ''),
    ('Incomplete_Denominator',       'True if any allele has Incomplete_Denominator=True',     'any(Incomplete_Denominator)', ''),
    ('Denominator_Complete',         'True when no denominator problem was detected',          'Denominator_Issue_Type == NONE', 'Used by reliable_alleles_by_sample filters'),
    ('Denominator_Issue_Type',       'NONE / UNOBSERVABLE / LOW_COMPLETENESS / MIXED',         'classification rule', ''),
    ('Mean_Hidden_Component_Pct',    'Mean Hidden_Component_Pct across locus alleles. QC only.', 'mean(Hidden_Component_Pct)', 'Not proof that a hidden allele is present'),
    ('Max_Hidden_Component_Pct',     'Maximum Hidden_Component_Pct across locus alleles. QC only.', 'max(Hidden_Component_Pct)', 'Not proof that a hidden allele is present'),
    ('Potential_Hidden_Allele',      'True when the locus has evidence consistent with a hidden allele; not proof of presence.', 'Max_HCP>50 OR Has_Unobservable OR Denom_Issue OR Incomplete_Denom', 'Used by reliable_alleles_by_sample filters'),
    ('Hidden_Allele_Evidence_Score', 'Composite 0-1 score for hidden-allele evidence. QC only.', '0.5×mean_HCP/100 + 0.25×unobs + 0.15×denom_issue + 0.10×incompl; max=1.0', ''),
]

_SAMPLE_LEVEL_DESC = [
    ('Sample_ID',                        'Sample identifier',                                    'from locus_level', ''),
    ('Access_Number',                    'FASTQ file name',                                      'from locus_level', ''),
    ('N_filt',                           'Sum of Raw_Count across all alleles in the sample',     'Σ Raw_Count',      ''),
    ('G_mini',                           'Sum of Fragment_Length across all alleles in the sample', 'Σ Fragment_Length_bp', ''),
    ('C_mini',                           'Mini-genome coverage',                                  'N_filt · R / G_mini', ''),
    ('N_Loci_Total',                     'Total number of loci in the sample',                    'count(Locus_ID)',  ''),
    ('N_Clean_Monoallelic_Loci',         'Number of clean monoallelic loci',                      'Σ Is_Clean_Monoallelic_Locus', ''),
    ('B_sample',                         'Sample baseline',                                       'median(Q_i over clean monoallelic loci)', ''),
    ('Sample_Baseline_Reliable',         f'True if at least {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE} clean loci contributed', f'N_clean ≥ {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE}', ''),
    ('N_Loci_With_Cohort_Baseline',      'Loci with cohort calibration applied',                  'Σ cohort_applied', ''),
    ('N_Loci_Without_Cohort_Baseline',   'Loci using only sample-level adjustment',               'total − with_cohort', ''),
    ('N_Loci_With_Unobservable_Alleles', 'Loci containing Is_Unobservable alleles',               'Σ Has_Unobservable_Alleles', ''),
    ('N_Loci_Incomplete_Denominator',    'Loci with Incomplete_Denominator=True',                 'Σ Incomplete_Denominator', ''),
    ('N_Loci_Potential_Hidden_Allele',   'Loci with hidden-allele evidence; not proof of presence', 'count(Potential_Hidden_Allele)', 'QC only'),
    ('Potential_Mixture_Flag',           'True if multi-allelic loci or hidden-allele evidence are present', 'has_observed_multi OR N_Loci_Potential_Hidden_Allele > 0', 'False does not guarantee purity'),
    ('Mixture_Type',                     'CLEAN_SINGLE / OBSERVED_MIX / HIDDEN_MIX / AMBIGUOUS', 'classification rule', 'AMBIGUOUS means both observed and hidden signals are present'),
]


def build_allele_level_description():
    return pd.DataFrame(_ALLELE_LEVEL_DESC, columns=_DESC_COLS)


def build_locus_level_description():
    return pd.DataFrame(_LOCUS_LEVEL_DESC, columns=_DESC_COLS)


def build_sample_level_description():
    return pd.DataFrame(_SAMPLE_LEVEL_DESC, columns=_DESC_COLS)

def write_master_excel(output_path, df_allele, df_locus, df_sample, df_dict):
    """Write master workbook with data sheets + column description sheets."""
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        # ── Data sheets (unchanged) ───────────────────────────────────────────
        df_allele.to_excel(writer, sheet_name='allele_level',       index=False)
        df_locus.to_excel( writer, sheet_name='locus_level',        index=False)
        df_sample.to_excel(writer, sheet_name='sample_level',       index=False)
        df_dict.to_excel(  writer, sheet_name='output_dictionary',  index=False)
        # ── Column description sheets ────────────────────────────────────────
        build_allele_level_description().to_excel(
            writer, sheet_name='allele_level_description',  index=False)
        build_locus_level_description().to_excel(
            writer, sheet_name='locus_level_description',   index=False)
        build_sample_level_description().to_excel(
            writer, sheet_name='sample_level_description',  index=False)
    _beautify_workbook(output_path)
    logger.info(f'Saved master Excel -> {output_path}')


# ─────────────────────────────────────────────────────────────
# STEP 14.  Sample metadata extraction
# ─────────────────────────────────────────────────────────────

def _build_sample_meta(data):
    """Extract one-sample-per-row metadata (mean_read_length, etc.)."""
    sample_meta = {}
    meta_cols = {
        'Total_Bases_Sequenced': 'total_bases',
        'Read_Count':            'read_count',
        'Mean_Read_Length':      'mean_read_length',
    }
    for _, row in data.iterrows():
        acc = str(row.get('Access_number', ''))
        if acc not in sample_meta:
            m = {}
            for col, key in meta_cols.items():
                v = row.get(col)
                m[key] = _to_float(v)
            sample_meta[acc] = m
    return sample_meta


# ─────────────────────────────────────────────────────────────
# STEP 15.  QC checks
# ─────────────────────────────────────────────────────────────

def run_qc_checks(df_allele, df_locus):
    """Validate per-locus share summations (Idealized and Experimental)."""
    issues = 0
    for share_col in ('Idealized_Share_%', 'Experimental_Share_%'):
        valid = df_allele[df_allele[share_col].notna()]
        sums  = valid.groupby(['Sample_ID', 'Locus_ID'])[share_col].sum()
        bad   = sums[(sums - 100.0).abs() > 0.5]
        for idx, val in bad.items():
            logger.warning(f'QC: {share_col} sum = {val:.4f} '
                           f'(expected 100) sample={idx[0]} locus={idx[1]}')
            issues += 1
    if issues == 0:
        logger.info('QC checks passed.')
    else:
        logger.warning(f'QC: {issues} issue(s) found.')


# ─────────────────────────────────────────────────────────────
# STEP 16.  Interactive HTML plots  (Plotly — optional)
# ─────────────────────────────────────────────────────────────

def _check_plotly():
    if not _PLOTLY_AVAILABLE:
        logger.warning(
            'Plotly is not installed; skipping interactive plots. '
            'Install with: pip install plotly'
        )
    return _PLOTLY_AVAILABLE


def _allele_role_label(rank):
    return 'Dominant' if rank == 1 else 'Secondary'


def _build_ranked_df(df_allele, share_col):
    """Return allele_level copy with Allele_Rank_In_Locus and _share/_role cols."""
    df = df_allele.copy()
    df['_qi_num'] = pd.to_numeric(df['Q_i'], errors='coerce')
    df['Allele_Rank_In_Locus'] = (
        df.groupby(['Sample_ID', 'Locus_ID'])['_qi_num']
        .rank(ascending=False, method='first', na_option='bottom')
        .fillna(1)
        .astype(int)
    )
    df['_share'] = pd.to_numeric(df[share_col], errors='coerce').fillna(0.0)
    df['_role']  = df['Allele_Rank_In_Locus'].apply(_allele_role_label)
    return df


def _get_locus_colors(loci):
    """Return {locus_id: hex_color} for up to ~60 loci."""
    palette = (
        px.colors.qualitative.Plotly
        + px.colors.qualitative.D3
        + px.colors.qualitative.G10
        + px.colors.qualitative.T10
        + px.colors.qualitative.Alphabet
    )
    seen, uniq = set(), []
    for c in palette:
        if c not in seen:
            uniq.append(c); seen.add(c)
    return {lid: uniq[i % len(uniq)] for i, lid in enumerate(loci)}


_SEC_PATTERNS = ['/', '\\', 'x', '+']


def _make_customdata(r):
    """Build 9-element list for hovertemplate from a pandas Series row."""
    def fmt(key, decimals=None):
        v = r.get(key) if hasattr(r, 'get') else None
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return 'N/A'
        if decimals is not None:
            try:
                return round(float(v), decimals)
            except (ValueError, TypeError):
                return 'N/A'
        return str(v)
    return [
        fmt('Sample_ID'),
        fmt('Locus_ID'),
        fmt('Allele_Value'),
        fmt('_role'),
        fmt('Raw_Count'),
        fmt('Q_i', 6),
        fmt('Fragment_Length_bp'),
        fmt('Detectability_Class'),
        fmt('Denominator_Issue_Type'),
    ]


_HOVER_TMPL = (
    '<b>Sample:</b> %{customdata[0]}<br>'
    '<b>Locus:</b> %{customdata[1]}<br>'
    '<b>Allele:</b> %{customdata[2]}<br>'
    '<b>Role:</b> %{customdata[3]}<br>'
    '<b>Share_%:</b> %{y:.4f}<br>'
    '<b>Raw_Count:</b> %{customdata[4]}<br>'
    '<b>Q_i:</b> %{customdata[5]}<br>'
    '<b>Fragment_Length_bp:</b> %{customdata[6]}<br>'
    '<b>Detectability_Class:</b> %{customdata[7]}<br>'
    '<b>Denominator_Issue_Type:</b> %{customdata[8]}<br>'
    '<extra></extra>'
)


def build_idealized_plots_html(df_allele, output_dir):
    """
    Per-sample stacked bar subplots of Idealized_Share_%.
    Dominant allele (max Q_i) = solid color; secondary = striped pattern.
    Saves: plots_idealized_share.html
    """
    if not _check_plotly():
        return

    df = _build_ranked_df(df_allele, 'Idealized_Share_%')
    locus_colors = _get_locus_colors(df['Locus_ID'].unique())
    samples  = list(df['Sample_ID'].unique())
    n        = len(samples)
    max_rank = int(df['Allele_Rank_In_Locus'].max())

    v_spacing = 0.4 / (n - 1) if n > 1 else 0.0
    fig = make_subplots(
        rows=n, cols=1,
        subplot_titles=[str(s) for s in samples],
        vertical_spacing=v_spacing,
        shared_xaxes=False,
    )

    for row_idx, sample in enumerate(samples, start=1):
        sdf   = df[df['Sample_ID'] == sample]
        loci  = sorted(sdf['Locus_ID'].unique())

        for rank in range(1, max_rank + 1):
            rank_df = sdf[sdf['Allele_Rank_In_Locus'] == rank]
            x_vals, y_vals, colors_list, cd_list = [], [], [], []

            for locus in loci:
                sub   = rank_df[rank_df['Locus_ID'] == locus]
                color = locus_colors.get(locus, '#1f77b4')
                if len(sub) > 0:
                    r = sub.iloc[0]
                    x_vals.append(locus)
                    y_vals.append(float(r['_share']))
                    colors_list.append(color)
                    cd_list.append(_make_customdata(r))
                else:
                    x_vals.append(locus)
                    y_vals.append(0.0)
                    colors_list.append(color)
                    cd_list.append([sample, locus, '', '', '', '', '', '', ''])

            is_sec   = (rank > 1)
            pat_shp  = _SEC_PATTERNS[(rank - 2) % len(_SEC_PATTERNS)] if is_sec else ''
            pat_kw   = dict(shape=pat_shp, size=8, solidity=0.4) if is_sec else dict(shape='')

            fig.add_trace(go.Bar(
                x=x_vals,
                y=y_vals,
                name='Dominant' if rank == 1 else f'Secondary #{rank - 1}',
                legendgroup=f'rank_{rank}',
                showlegend=(row_idx == 1),
                marker=dict(
                    color=colors_list,
                    opacity=1.0 if not is_sec else 0.65,
                    pattern=pat_kw,
                ),
                customdata=cd_list,
                hovertemplate=_HOVER_TMPL,
            ), row=row_idx, col=1)

        fig.update_yaxes(range=[0, 108], title_text='Idealized_Share_%',
                         row=row_idx, col=1)

    fig.update_layout(
        barmode='stack',
        title=dict(text='MLVA Idealized Share % — Per Sample (alleles within locus)',
                   font=dict(size=16)),
        height=max(300, 310 * n),
        template='plotly_white',
        legend_title_text='Allele Role',
        hoverlabel=dict(bgcolor='white', font_size=12),
    )

    out_path = os.path.join(output_dir, 'plots_idealized_share.html')
    fig.write_html(out_path, include_plotlyjs='cdn')
    logger.info(f'Saved idealized share plots    -> {out_path}')


def build_all_samples_idealized_plot(df_allele, output_dir):
    """
    One wide stacked bar chart of Idealized_Share_% across all samples & loci.
    X = "Sample · Locus", Y = Idealized_Share_%.
    Saves: all_samples_idealized_alleles.html
    """
    if not _check_plotly():
        return

    df        = _build_ranked_df(df_allele, 'Idealized_Share_%')
    df['_xl'] = df['Sample_ID'].astype(str) + ' · ' + df['Locus_ID'].astype(str)
    df_sorted = df.sort_values(['Sample_ID', 'Locus_ID', 'Allele_Rank_In_Locus'])
    x_order   = list(dict.fromkeys(df_sorted['_xl'].tolist()))  # stable unique order
    max_rank  = int(df['Allele_Rank_In_Locus'].max())

    dom_color  = '#2196F3'
    sec_colors = ['#FF9800', '#E91E63', '#9C27B0']

    fig = go.Figure()

    for rank in range(1, max_rank + 1):
        rank_df = df_sorted[df_sorted['Allele_Rank_In_Locus'] == rank]
        x_map   = {row['_xl']: row for _, row in rank_df.iterrows()}

        x_vals, y_vals, cd_list = [], [], []
        for xl in x_order:
            if xl in x_map:
                r = x_map[xl]
                x_vals.append(xl)
                y_vals.append(float(r['_share']))
                cd_list.append(_make_customdata(r))
            else:
                x_vals.append(xl)
                y_vals.append(0.0)
                cd_list.append(['', '', '', '', '', '', '', '', ''])

        is_sec  = (rank > 1)
        color   = dom_color if not is_sec else sec_colors[(rank - 2) % len(sec_colors)]
        pat_shp = _SEC_PATTERNS[(rank - 2) % len(_SEC_PATTERNS)] if is_sec else ''

        fig.add_trace(go.Bar(
            x=x_vals,
            y=y_vals,
            name='Dominant' if rank == 1 else f'Secondary #{rank - 1}',
            marker=dict(
                color=color,
                opacity=1.0 if not is_sec else 0.70,
                pattern=dict(shape=pat_shp, size=8, solidity=0.4) if is_sec else dict(shape=''),
            ),
            customdata=cd_list,
            hovertemplate=_HOVER_TMPL,
        ))

    chart_width = max(1400, len(x_order) * 36)
    fig.update_layout(
        barmode='stack',
        title=dict(text='MLVA Idealized Share % — All Samples & Loci',
                   font=dict(size=16)),
        xaxis=dict(title='Sample · Locus', tickangle=-60,
                   tickfont=dict(size=9), automargin=True),
        yaxis=dict(title='Idealized_Share_%', range=[0, 108]),
        height=620,
        width=chart_width,
        template='plotly_white',
        legend_title_text='Allele Role',
        hoverlabel=dict(bgcolor='white', font_size=12),
        margin=dict(b=200),
    )

    html = fig.to_html(
        include_plotlyjs='cdn', full_html=True,
        config={'scrollZoom': True, 'displayModeBar': True},
    )
    html = html.replace('</head>',
                        '<style>body{overflow-x:auto;}</style>\n</head>')

    out_path = os.path.join(output_dir, 'all_samples_idealized_alleles.html')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(html)
    logger.info(f'Saved all-samples idealized    -> {out_path}')


def build_experimental_plots_html(df_allele, output_dir):
    """
    Per-sample and all-samples stacked bar plots of Experimental_Share_%.
    Saves: plots_experimental_share.html, all_samples_experimental_alleles.html
    """
    if not _check_plotly():
        return

    share_col = 'Experimental_Share_%'
    if df_allele[share_col].isna().all():
        logger.info('No Experimental_Share_% data; skipping experimental plots.')
        return

    df = _build_ranked_df(df_allele, share_col)
    locus_colors = _get_locus_colors(df['Locus_ID'].unique())
    samples  = list(df['Sample_ID'].unique())
    n        = len(samples)
    max_rank = int(df['Allele_Rank_In_Locus'].max())

    # ---- per-sample subplots ----
    v_spacing = 0.4 / (n - 1) if n > 1 else 0.0
    fig = make_subplots(
        rows=n, cols=1,
        subplot_titles=[str(s) for s in samples],
        vertical_spacing=v_spacing,
        shared_xaxes=False,
    )

    for row_idx, sample in enumerate(samples, start=1):
        sdf  = df[df['Sample_ID'] == sample]
        loci = sorted(sdf['Locus_ID'].unique())

        for rank in range(1, max_rank + 1):
            rank_df = sdf[sdf['Allele_Rank_In_Locus'] == rank]
            x_vals, y_vals, colors_list, cd_list = [], [], [], []

            for locus in loci:
                sub   = rank_df[rank_df['Locus_ID'] == locus]
                color = locus_colors.get(locus, '#1f77b4')
                if len(sub) > 0:
                    r = sub.iloc[0]
                    x_vals.append(locus)
                    y_vals.append(float(r['_share']))
                    colors_list.append(color)
                    cd_list.append(_make_customdata(r))
                else:
                    x_vals.append(locus)
                    y_vals.append(0.0)
                    colors_list.append(color)
                    cd_list.append([sample, locus, '', '', '', '', '', '', ''])

            is_sec  = (rank > 1)
            pat_shp = _SEC_PATTERNS[(rank - 2) % len(_SEC_PATTERNS)] if is_sec else ''

            fig.add_trace(go.Bar(
                x=x_vals, y=y_vals,
                name='Dominant' if rank == 1 else f'Secondary #{rank - 1}',
                legendgroup=f'rank_{rank}',
                showlegend=(row_idx == 1),
                marker=dict(
                    color=colors_list,
                    opacity=1.0 if not is_sec else 0.65,
                    pattern=dict(shape=pat_shp, size=8, solidity=0.4) if is_sec
                            else dict(shape=''),
                ),
                customdata=cd_list,
                hovertemplate=_HOVER_TMPL,
            ), row=row_idx, col=1)

        fig.update_yaxes(range=[0, 108], title_text='Experimental_Share_%',
                         row=row_idx, col=1)

    fig.update_layout(
        barmode='stack',
        title=dict(text='MLVA Experimental Share % — Per Sample',
                   font=dict(size=16)),
        height=max(300, 310 * n),
        template='plotly_white',
        legend_title_text='Allele Role',
        hoverlabel=dict(bgcolor='white', font_size=12),
    )

    out1 = os.path.join(output_dir, 'plots_experimental_share.html')
    fig.write_html(out1, include_plotlyjs='cdn')
    logger.info(f'Saved experimental share plots -> {out1}')

    # ---- all-samples chart ----
    df['_xl'] = df['Sample_ID'].astype(str) + ' · ' + df['Locus_ID'].astype(str)
    df_sorted = df.sort_values(['Sample_ID', 'Locus_ID', 'Allele_Rank_In_Locus'])
    x_order   = list(dict.fromkeys(df_sorted['_xl'].tolist()))

    dom_color  = '#4CAF50'
    sec_colors = ['#FF9800', '#E91E63', '#9C27B0']

    fig2 = go.Figure()
    for rank in range(1, max_rank + 1):
        rank_df = df_sorted[df_sorted['Allele_Rank_In_Locus'] == rank]
        x_map   = {row['_xl']: row for _, row in rank_df.iterrows()}

        x_vals, y_vals, cd_list = [], [], []
        for xl in x_order:
            if xl in x_map:
                r = x_map[xl]
                x_vals.append(xl); y_vals.append(float(r['_share']))
                cd_list.append(_make_customdata(r))
            else:
                x_vals.append(xl); y_vals.append(0.0)
                cd_list.append(['', '', '', '', '', '', '', '', ''])

        is_sec  = (rank > 1)
        color   = dom_color if not is_sec else sec_colors[(rank - 2) % len(sec_colors)]
        pat_shp = _SEC_PATTERNS[(rank - 2) % len(_SEC_PATTERNS)] if is_sec else ''

        fig2.add_trace(go.Bar(
            x=x_vals, y=y_vals,
            name='Dominant' if rank == 1 else f'Secondary #{rank - 1}',
            marker=dict(
                color=color,
                opacity=1.0 if not is_sec else 0.70,
                pattern=dict(shape=pat_shp, size=8, solidity=0.4) if is_sec
                        else dict(shape=''),
            ),
            customdata=cd_list,
            hovertemplate=_HOVER_TMPL,
        ))

    chart_width = max(1400, len(x_order) * 36)
    fig2.update_layout(
        barmode='stack',
        title=dict(text='MLVA Experimental Share % — All Samples & Loci',
                   font=dict(size=16)),
        xaxis=dict(title='Sample · Locus', tickangle=-60,
                   tickfont=dict(size=9), automargin=True),
        yaxis=dict(title='Experimental_Share_%', range=[0, 108]),
        height=620, width=chart_width,
        template='plotly_white',
        legend_title_text='Allele Role',
        hoverlabel=dict(bgcolor='white', font_size=12),
        margin=dict(b=200),
    )

    html2 = fig2.to_html(
        include_plotlyjs='cdn', full_html=True,
        config={'scrollZoom': True, 'displayModeBar': True},
    )
    html2 = html2.replace('</head>',
                           '<style>body{overflow-x:auto;}</style>\n</head>')

    out2 = os.path.join(output_dir, 'all_samples_experimental_alleles.html')
    with open(out2, 'w', encoding='utf-8') as f:
        f.write(html2)
    logger.info(f'Saved all-samples experimental -> {out2}')


def write_hidden_allele_notes(output_dir, df_sample=None):
    """
    Write hidden_allele_interpretation.txt.
    Optionally appends a per-sample mixture summary table.
    """
    sep  = '─' * 72
    sep2 = '=' * 72
    lines = [
        sep2,
        'MLVA/VNTR Pipeline v21 — Hidden Allele Interpretation Guide',
        sep2,
        '',
        '─── CORE CONCEPT ──────────────────────────────────────────────────────',
        '',
        'In NGS-based MLVA/VNTR analysis, an allele can be ABSENT FROM THE DATA',
        'without being absent from the biological sample.',
        '',
        'This happens when the PCR amplicon for that allele is LONGER than the',
        'sequencing read length (e.g. amplicon = 350 bp, read = 150 bp).',
        '',
        'In that case:',
        '  p_i  ≈ 0   →   detectability probability approaches zero',
        '  E_ideal = C_mini × p_i  ≈ 0   →   expected reads ≈ 0',
        '  Raw_Count observed  ≈ 0   →   allele appears absent in output',
        '',
        'CRITICAL DISTINCTION:',
        '  ❌  "Allele not in data"  ≠  "Allele not in sample"',
        '  ✅  "Allele not in data with p_i ≈ 0"  =  "Allele likely present',
        '       but physically undetectable with current read length"',
        '',
        '─── FORMULAS ──────────────────────────────────────────────────────────',
        '',
        'Hidden_Allele_Flag:',
        '  True  if  p_i ≤ 0  OR  Detectability_Class == "Unobservable"',
        '  Signals: this allele exists in the marker profile but CANNOT be',
        '  reliably sequenced.',
        '',
        'Hidden_Component_Pct:',
        '  (1 − p_i) × 100   [clamped 0–100%]',
        '',
        '  Meaning:',
        '    0 %   → allele is FULLY OBSERVABLE (short amplicon, high p_i)',
        '   30 %   → 30 % of allele signal is "hidden" by detectability loss',
        '  100 %   → allele is COMPLETELY UNDETECTABLE (p_i = 0)',
        '',
        '  Note: Hidden_Component_Pct does NOT enter into any Q_i or',
        '  share calculation — it is a DIAGNOSTIC METRIC ONLY.',
        '',
        '─── LOCUS-LEVEL FLAGS ─────────────────────────────────────────────────',
        '',
        'Potential_Hidden_Allele (locus_level):',
        '  True when any of these conditions are met:',
        '  A)  Max_Hidden_Component_Pct > 50 %',
        '      → at least one allele in the locus is mostly undetectable',
        '  B)  Has_Unobservable_Alleles == True',
        '      → at least one allele explicitly flagged Is_Unobservable',
        '  C)  Denominator_Issue_Type != NONE',
        '      → UNOBSERVABLE / LOW_COMPLETENESS / MIXED issue present',
        '  D)  Incomplete_Denominator == True',
        '      → upstream pipeline flagged the denominator as incomplete',
        '',
        'Hidden_Allele_Evidence_Score (0–1):',
        '  Composite score estimating the probability of a hidden allele.',
        '  Higher = stronger evidence.',
        '  Formula:',
        '    score = 0.50 × mean(Hidden_Component_Pct) / 100',
        '          + 0.25  [if any allele Is_Unobservable]',
        '          + 0.15  [if Denominator_Issue_Type != NONE]',
        '          + 0.10  [if Incomplete_Denominator == True]',
        '    score = min(1.0, score)',
        '',
        '─── SAMPLE-LEVEL MIXTURE CLASSIFICATION ───────────────────────────────',
        '',
        'Mixture_Type values:',
        '',
        '  CLEAN_SINGLE   — No multi-allelic loci AND no Potential_Hidden_Allele.',
        '                   Consistent with a pure single-strain sample.',
        '                   ⚠ Does not PROVE purity — hidden alleles with',
        '                   p_i = 0 leave no trace at all.',
        '',
        '  OBSERVED_MIX   — At least one locus has N_Alleles_Detected > 1.',
        '                   Mixture is directly visible in the data.',
        '',
        '  HIDDEN_MIX     — No visible multi-allelic loci, BUT at least one',
        '                   locus has Potential_Hidden_Allele = True.',
        '                   The mix is biologically plausible but the',
        '                   secondary allele cannot be observed directly.',
        '                   Example: ft18 + ft29 mixture where the ft29-',
        '                   specific m3 allele is too long to sequence.',
        '',
        '  AMBIGUOUS      — Both observed multi-allelic loci AND hidden allele',
        '                   evidence present. Complex mixture — requires',
        '                   manual review.',
        '',
        '─── BIOLOGICAL EXAMPLE: ft18 vs ft29 ─────────────────────────────────',
        '',
        'Consider a mixed infection of two F. tularensis strains:',
        '  • Strain ft18: does NOT carry allele m3',
        '  • Strain ft29: DOES carry allele m3, but the m3 amplicon is',
        '                 ~ 350 bp while the sequencing read is 150 bp.',
        '',
        'In a ft18+ft29 mixture:',
        '  1. m3 allele exists in the sample (from ft29)',
        '  2. p_i(m3) ≈ 0  →  Expected_Reads ≈ 0',
        '  3. Sequencer produces 0 reads for m3',
        '  4. Naïve script: "m3 absent → sample looks like ft18 only"',
        '  5. This pipeline: Potential_Hidden_Allele = True for locus m3',
        '                    Mixture_Type = HIDDEN_MIX',
        '',
        'The Mixture_Type = HIDDEN_MIX flag triggers additional scrutiny,',
        'preventing misidentification of the mixed sample as pure ft18.',
        '',
        '─── HOW TO USE THESE RESULTS ──────────────────────────────────────────',
        '',
        '  1. Filter locus_level WHERE Potential_Hidden_Allele = True',
        '     → identify loci that warrant caution in interpretation.',
        '',
        '  2. Check sample_level Mixture_Type:',
        '     CLEAN_SINGLE → likely pure, but verify Hidden_Allele_Evidence_Score',
        '     OBSERVED_MIX → mixture confirmed by visible data',
        '     HIDDEN_MIX   → potential hidden mixture, increase read length',
        '                    or use long-read sequencing to resolve',
        '     AMBIGUOUS    → complex case, manual expert review required',
        '',
        '  3. To resolve HIDDEN_MIX: increase sequencing read length so that',
        '     Fragment_Length_bp ≤ Mean_Read_Length_bp for all alleles,',
        '     thereby raising p_i above 0.',
        '',
        sep2,
    ]

    # Optional per-sample summary
    if df_sample is not None and 'Mixture_Type' in df_sample.columns:
        lines += [
            '',
            '─── PER-SAMPLE MIXTURE SUMMARY ────────────────────────────────────',
            '',
            f'  {"Sample_ID":<40} {"Mixture_Type":<15} {"N_Hidden_Loci":<15} {"Pot_Mix"}',
            f'  {"─"*40} {"─"*14} {"─"*14} {"─"*8}',
        ]
        for _, row in df_sample.iterrows():
            sid      = str(row.get('Sample_ID', ''))[:39]
            mtype    = str(row.get('Mixture_Type', ''))
            n_hidden = str(row.get('N_Loci_Potential_Hidden_Allele', ''))
            pot_mix  = str(row.get('Potential_Mixture_Flag', ''))
            lines.append(f'  {sid:<40} {mtype:<15} {n_hidden:<15} {pot_mix}')
        lines.append('')

    out_path = os.path.join(output_dir, 'hidden_allele_interpretation.txt')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    logger.info(f'Saved hidden allele notes      -> {out_path}')


def write_interpretation_notes(output_dir):
    """Write interpretation_notes.txt explaining Idealized_Share_% vs confidence."""
    sep  = '─' * 72
    sep2 = '=' * 72
    lines = [
        sep2,
        'MLVA/VNTR Pipeline v21 — Interpretation Notes',
        sep2,
        '',
        'This file explains the meaning of key output columns, in particular',
        'why Idealized_Share_% and Experimental_Share_% are NOT confidence scores.',
        '',
        sep,
        'Q1.  Is Idealized_Share_% a "% of reliability / confidence"',
        '     for the idealized model?',
        sep,
        '',
        'SHORT ANSWER: No.',
        '',
        'Idealized_Share_% is the RELATIVE SHARE (proportion) of this allele',
        'within a specific locus for a specific sample, computed under the',
        'idealized (theoretical coverage) model.',
        '',
        'Formula:',
        '    Idealized_Share_% = Q_i / Σ(Q_j for all alleles j in the same locus) × 100',
        '    where Q_i = Raw_Count / E_ideal  (E_ideal = C_mini · p_i)',
        '',
        'What it tells you:',
        '  • Idealized_Share_% ≈ 100 %  →  single meaningful allele (dominant or mono).',
        '  • Idealized_Share_% < 100 %  →  secondary alleles present; value = their',
        '    proportional contribution (proxy for contamination / mixture fraction).',
        '  • For secondary alleles: how large their relative signal is compared to',
        '    the dominant allele, after correcting for detectability differences.',
        '',
        'What it does NOT tell you:',
        '  • Whether the idealized model fit is trustworthy.',
        '  • Whether Q_i was computed from a complete allele denominator.',
        '  • Whether the locus has dropout or unobservable alleles hidden from the sum.',
        '',
        sep,
        'Q2.  Is Experimental_Share_% a "% of reliability / confidence"',
        '     for the experimental model?',
        sep,
        '',
        'SHORT ANSWER: No.',
        '',
        'Experimental_Share_% is the RELATIVE SHARE of this allele within a locus',
        'after a two-step empirical normalization.',
        '',
        'Formula:',
        '    Q_i_sample_adjusted = Q_i / B_sample',
        '    Q_i_experimental    = Q_i_sample_adjusted / B_locus_cohort',
        '                          (fallback = Q_i_sample_adjusted if no cohort baseline)',
        '    Experimental_Share_% = Q_i_experimental / Σ(Q_j_experimental in locus) × 100',
        '',
        'What it tells you:',
        '  • An empirically bias-corrected view of allele proportion.',
        '  • Accounts for sample-level coverage bias (B_sample) and cohort-level',
        '    locus bias (B_locus_cohort).',
        '  • Secondary alleles with significant Experimental_Share_% after both',
        '    correction steps are stronger evidence of true mixture / contamination.',
        '',
        'What it does NOT tell you:',
        '  • If Sample_Baseline_Reliable = False or Locus_Cohort_Baseline_Reliable',
        '    = False, the normalization is approximate — interpret with caution.',
        '',
        sep,
        'Key conceptual difference between the two metrics',
        sep,
        '',
        'Idealized_Share_%     — proportion after THEORETICAL correction',
        '                        (mini-genome coverage model; no empirical calibration).',
        '',
        'Experimental_Share_%  — proportion after EMPIRICAL two-step normalization:',
        '                        (1) sample-level baseline B_sample,',
        '                        (2) cohort locus baseline B_locus_cohort.',
        '',
        'Neither metric is a confidence / reliability score.',
        'Both are REPRESENTATIONAL metrics: what fraction of the locus signal',
        'comes from this allele under a given normalization model.',
        '',
        sep,
        'Columns that DO reflect reliability / quality of the result',
        sep,
        '',
        f'  Denominator_Complete           True if no unobservable/hidden alleles',
        f'                                 inflate the share denominator.',
        f'  Denominator_Issue_Type         NONE / UNOBSERVABLE / LOW_COMPLETENESS / MIXED.',
        f'  Confidence_Score_Locus         Upstream per-locus confidence score.',
        f'  Sample_Baseline_Reliable       True if B_sample was built from',
        f'                                 ≥ {MIN_CLEAN_MONOALLELIC_FOR_SAMPLE} clean monoallelic loci.',
        f'  Locus_Cohort_Baseline_Reliable True if B_locus_cohort was built from',
        f'                                 ≥ {MIN_CLEAN_MONOALLELIC_FOR_LOCUS} cohort observations.',
        f'  Is_Clean_Monoallelic_Locus     True = safest, unambiguous single allele.',
        '',
        sep2,
    ]
    out_path = os.path.join(output_dir, 'interpretation_notes.txt')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    logger.info(f'Saved interpretation notes     -> {out_path}')


# ─────────────────────────────────────────────────────────────
# STEP 17.  Main pipeline
# ─────────────────────────────────────────────────────────────

def run_pipeline(data, input_file, args):
    run_ts     = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # -- Numeric conversion for relevant input columns -------------
    numeric_markers = (
        'Raw_Count', 'Corrected_Count', 'Loci_Size',
        'Detectability_p', 'Expected_Reads', 'Observation_Ratio',
        'Locus_Factor', 'Z_Score', 'Baseline_Sample', 'Confidence_Score',
        'Mean_Read_Length', 'Genome_Size', 'Total_Bases',
        'Locus_Completeness',
    )
    for col in data.columns:
        if any(k in col for k in numeric_markers):
            data[col] = pd.to_numeric(data[col], errors='coerce')

    indices = detect_allele_blocks(data)
    logger.info(f'Detected allele indices: {indices}')

    sample_meta = _build_sample_meta(data)

    # -- STEP 1 : flatten to allele_level (raw only) ---------------
    logger.info('Flattening wide-format input → long-format allele_level ...')
    df_allele = build_allele_level_raw(data, indices, args.max_alleles, sample_meta)
    df_allele = df_allele.sort_values(
        ['Sample_ID', 'Locus_ID', 'Allele_Rank_Input']
    ).reset_index(drop=True)
    logger.info(f'  allele_level rows : {len(df_allele)}')

    # -- STEP 1b : full allele list for reporting (Raw_Count=0 included) ---
    # ❗ REPORTING ONLY — never passed to calculation steps.
    #    Raw_Count=0 rows do NOT affect N_filt, G_mini, C_mini, Q_i,
    #    B_sample, B_locus_cohort, or Q_i_experimental.
    logger.info('Building full allele list (including Raw_Count=0 for reporting) ...')
    df_all_alleles_full = build_allele_level_full(
        data, indices, args.max_alleles, sample_meta)
    logger.info(f'  allele_level_full rows: {len(df_all_alleles_full)} '
                f'(includes {(df_all_alleles_full["Raw_Count"] == 0).sum()} zero-count rows)')

    # -- STEP 2 : idealized model (mini-genome → Q_i) --------------
    logger.info('Computing idealized model (mini-genome coverage → Q_i) ...')
    df_allele = apply_idealized_model(df_allele)

    # -- STEP 3 : clean monoallelic flag ---------------------------
    logger.info('Marking clean monoallelic loci ...')
    df_allele = apply_clean_monoallelic_flag(df_allele)

    # -- STEP 4 : experimental model (B_sample → B_locus_cohort) ---
    logger.info('Applying two-step experimental adjustment ...')
    df_allele = apply_experimental_model(df_allele)

    # -- STEP 5 : denominator completeness -------------------------
    logger.info('Computing denominator completeness ...')
    df_allele = apply_denominator_completeness(df_allele)

    # -- STEP 5b : hidden allele annotation -----------------------
    logger.info('Annotating hidden alleles (Hidden_Allele_Flag, Hidden_Component_Pct) ...')
    df_allele = apply_hidden_allele_flags(df_allele)

    # -- STEP 6 : finalize allele_level ----------------------------
    df_allele = finalize_allele_level(df_allele)

    # -- STEP 7 : locus_level --------------------------------------
    logger.info('Building locus_level ...')
    df_locus = build_locus_level(df_allele)

    # -- STEP 7b : merge locus-level flags into df_all_alleles_full ----
    if not df_all_alleles_full.empty:
        locus_flags = df_locus[[
            'Sample_ID', 'Locus_ID',
            'Denominator_Complete', 'Potential_Hidden_Allele',
            'Denominator_Issue_Type',
        ]].drop_duplicates(subset=['Sample_ID', 'Locus_ID'])

        # Also bring Hidden_Allele_Flag / Hidden_Component_Pct from allele_level
        allele_hidden = df_allele[[
            'Sample_ID', 'Locus_ID', 'Allele_Value',
            'Hidden_Allele_Flag', 'Hidden_Component_Pct',
        ]].copy()
        allele_hidden['Allele_Value'] = allele_hidden['Allele_Value'].astype(str)

        # Drop old placeholder columns from df_all_alleles_full before merge
        for col in ['Denominator_Complete', 'Potential_Hidden_Allele',
                    'Denominator_Issue_Type']:
            if col in df_all_alleles_full.columns:
                df_all_alleles_full = df_all_alleles_full.drop(columns=[col])

        df_all_alleles_full = df_all_alleles_full.merge(
            locus_flags, on=['Sample_ID', 'Locus_ID'], how='left')

        df_all_alleles_full['Allele_Value_str'] = (
            df_all_alleles_full['Allele_Value'].astype(str))
        allele_hidden = allele_hidden.rename(
            columns={'Allele_Value': 'Allele_Value_str'})
        for col in ['Hidden_Allele_Flag', 'Hidden_Component_Pct']:
            if col in df_all_alleles_full.columns:
                df_all_alleles_full = df_all_alleles_full.drop(columns=[col])
        df_all_alleles_full = df_all_alleles_full.merge(
            allele_hidden,
            on=['Sample_ID', 'Locus_ID', 'Allele_Value_str'],
            how='left'
        ).drop(columns=['Allele_Value_str'])

    # -- STEP 8 : sample_level -------------------------------------
    logger.info('Building sample_level ...')
    df_sample = build_sample_level(df_allele, df_locus)

    # -- STEP 9 : QC checks ----------------------------------------
    run_qc_checks(df_allele, df_locus)

    # -- STEP 10 : output dictionary -------------------------------
    df_dict = build_output_dictionary()

    # -- STEP 11 : write master Excel ------------------------------
    master_path = os.path.join(output_dir, 'master_allele_to_sample_statistics.xlsx')
    write_master_excel(master_path, df_allele, df_locus, df_sample, df_dict)

    # -- STEP 12 : publication_tables.xlsx -------------------------
    pub_path = os.path.join(output_dir, 'publication_tables.xlsx')
    build_publication_tables(df_allele, df_locus, pub_path,
                             df_all_alleles_full=df_all_alleles_full)

    # -- STEP 13 : output_dictionary.tsv ---------------------------
    dict_path = os.path.join(output_dir, 'output_dictionary.tsv')
    df_dict.to_csv(dict_path, sep='\t', index=False)
    logger.info(f'Saved output dictionary -> {dict_path}')

    # -- STEP 14 : run_metadata.txt --------------------------------
    n_samples = df_allele['Sample_ID'].nunique()
    n_records = len(df_locus)
    n_clean_mono_samples = int(df_sample['Sample_Baseline_Reliable']
                               .fillna(False).astype(bool).sum()) \
                           if not df_sample.empty else 0
    # Number of loci (distinct Locus_ID) with cohort baseline reliable
    cohort_ok_loci = df_allele.loc[
        df_allele['Locus_Cohort_Baseline_Reliable'].fillna(False).astype(bool),
        'Locus_ID'
    ].unique()
    n_loci_cohort_ok = int(len(cohort_ok_loci))

    meta_path = os.path.join(output_dir, 'run_metadata.txt')
    write_run_metadata(meta_path, input_file, args, run_ts,
                        n_samples, n_records,
                        n_clean_mono_samples, n_loci_cohort_ok)

    # -- STEP 15 : interactive HTML plots -------------------------
    logger.info('Building interactive HTML plots ...')
    build_idealized_plots_html(df_allele, output_dir)
    build_all_samples_idealized_plot(df_allele, output_dir)
    build_experimental_plots_html(df_allele, output_dir)

    # -- STEP 15b : interpretation notes --------------------------
    write_interpretation_notes(output_dir)
    write_hidden_allele_notes(output_dir, df_sample=df_sample)

    # -- STEP 16 : debug outputs -----------------------------------
    if args.debug_outputs:
        debug_dir = os.path.join(output_dir, 'debug')
        os.makedirs(debug_dir, exist_ok=True)
        df_allele.to_csv(os.path.join(debug_dir, 'allele_level_full.csv'), index=False)
        df_locus .to_csv(os.path.join(debug_dir, 'locus_level.csv'),       index=False)
        df_sample.to_csv(os.path.join(debug_dir, 'sample_level.csv'),      index=False)
        logger.info(f'Saved debug CSVs to {debug_dir}')

    logger.info(f'Done. Outputs in: {output_dir}')
    for f in ('master_allele_to_sample_statistics.xlsx',
              'publication_tables.xlsx',
              'output_dictionary.tsv',
              'run_metadata.txt',
              'interpretation_notes.txt',
              'hidden_allele_interpretation.txt',
              'plots_idealized_share.html',
              'all_samples_idealized_alleles.html',
              'plots_experimental_share.html',
              'all_samples_experimental_alleles.html'):
        logger.info('  ' + f)

    return {
        'allele_level':      df_allele,
        'locus_level':       df_locus,
        'sample_level':      df_sample,
        'output_dictionary': df_dict,
        'run_metadata': {
            'run_ts':              run_ts,
            'n_samples':           n_samples,
            'n_records':           n_records,
            'n_clean_mono_samples':n_clean_mono_samples,
            'n_loci_cohort_ok':    n_loci_cohort_ok,
        },
    }


# ─────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='MLVA/VNTR contamination analysis — output module v21'
    )
    parser.add_argument('-i', '--input',       required=True,
                        help='Path to MLVA_v2 output file (.csv or .xlsx)')
    parser.add_argument('-o', '--output',      required=True,
                        help='Output directory')
    parser.add_argument('-m', '--max-alleles', type=int, default=4,
                        help='Maximum alleles per locus (default: 4)')
    parser.add_argument('--debug_outputs',     action='store_true',
                        help='Save debug CSV files to debug/')
    args = parser.parse_args()

    if args.input.endswith('.csv'):
        data = pd.read_csv(args.input)
    elif args.input.endswith(('.xlsx', '.xls')):
        data = pd.read_excel(args.input)
    else:
        raise ValueError('Input must be .csv or .xlsx')

    run_pipeline(data, args.input, args)
