#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# HYBRID MLVA PIPELINE — Idealized Model v2.4
# =============================================================================
#
# Detectability model (uniform start-site approximation):
#   p*(L, R) = (R - L + 1) / R   for L <= R   [per-read R, not mean_read_len]
#   p*(L, R) = 0                  for L >  R
#
# Detectability classes (configurable via CLI):
#   Observable              : p*(L) >= obs_threshold  [default 0.5]
#   Marginal                : p_min <= p*(L) < obs_threshold
#   Effectively_Unobservable: 0 < p*(L) < p_min       [default 0.05]
#   Unobservable            : p*(L) = 0   (L > R)
#
# Changes vs v2.3:
#   STRICT-STAGE REWRITE — full-length primer validation replaces seed-automaton:
#   - BBDuk (18 nt 3'-end seeds, hdist=2) kept as coarse pre-filter only.
#   - New collect_full_length_primer_hits(): searches all 4 primer forms
#     (FWD, REV, RC_FWD, RC_REV) at full primer length with ≤ nbmismatch Hamming.
#   - New build_candidate_pairs(): enumerates both valid orientation configs
#     (Config A: FWD+RC_REV, Config B: REV+RC_FWD) and all hit combinations.
#   - New evaluate_candidate_pair(): computes fragment size, split/non-split,
#     insert sequence, and validates 1 ≤ size ≤ max_distance.
#   - New select_best_candidate(): deterministic tie-breaking:
#     non-split > fewer mismatches > smaller size > smaller left_start.
#   - process_chunk() rewritten to use the new pipeline; all downstream
#     calculations (sizeU, p_star, idealized model) unchanged.
#   - New per-read audit fields in result_entry:
#     Selected_Config, Selected_Total_Mismatches, Selected_Is_Split,
#     Valid_Candidate_Count.
#   - Per-chunk counters for BBDuk→strict funnel diagnostics:
#     reads_after_bbduk, reads_passing_strict, reads_rejected_no_pair,
#     reads_rejected_invalid_size, reads_with_multiple_pairs.
#   - Old build_primer_automaton_exact_seeds / collect_verified_positions
#     removed from the hot path (kept only for BBDuk coarse filter step).
#
# Changes vs v2.1 (patch-list applied):
#   PATCH 1 : proxy_expected in compute_locus_completeness now uses
#             allele-specific p_eff = max(p_star, P_MIN) × locus_factor,
#             not the flat P_MIN constant. locus_factor added as argument.
#   PATCH 2 : (already correct in v2.1 — mean_p_star from per-read accum.)
#   PATCH 3 : (already correct in v2.1 — mean_p_star stored in agg maps.)
#   PATCH 4 : is_allele_dropout adds third guard: expected_reads ≥ baseline*0.3.
#             Prevents dropout flags when expected signal is inherently tiny.
#   PATCH 5 : hard_unobservable in assign_dropout_status now explicitly
#             includes Effectively_Unobservable (was already there; kept explicit).
#   PATCH 6 : Incomplete_Denominator threshold changed from fixed BASELINE_MIN_RAW
#             to sample-relative baseline_sample*0.3 (adaptive per sample).
#   PATCH 7 : Z_Score = (raw_count − expected_reads) / sqrt(expected_reads)
#             added per allele (Poisson deviation; None when expected unavailable).
#   PATCH 8 : Confidence_Score added per locus (multiplicative penalty for
#             Incomplete_Denominator, low completeness, unresolved status).
#   PATCH 9 : MAX_WEIGHT = 1/P_MIN already set at module level (no change needed).
#   PATCH 10: existing_status now derived from lc["Base_Status"] (computed
#             inside compute_locus_completeness) instead of hardcoded "UNKNOWN".
#             Enables assign_dropout_status to integrate with contamination
#             classification rather than operating as a detached layer.
#
# New output columns vs v2.1 (allele-level):
#   Z_Score {i}         — Poisson deviation from expected
# New output columns vs v2.1 (locus/sample-level):
#   Confidence_Score    — multiplicative quality flag (1.0 = no penalties)
#
# REQUIRED user argument: --genome-size INT  (genome size in bp)
# All other parameters derived from input files.
# CLI options: --p-min FLOAT  --obs-threshold FLOAT
# =============================================================================

from collections import defaultdict
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import re, os.path, getopt, csv, itertools
import time
from Bio import SeqIO
import os
import subprocess
import shutil
import pandas as pd
from tqdm import tqdm
import sys
import gzip
import math
import traceback
import psutil
import logging
import queue
import threading
import uuid
import numpy as np

# Note: pyahocorasick is no longer required in v2.4.
# The strict-stage now uses full-length Hamming search (collect_full_length_primer_hits).
# BBDuk (external tool) handles the coarse filter before reads reach this code.

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# BBDuk executable resolution.
# BBDuk is part of BBMap/BBTools and is not installed by pip. The path can be
# supplied with --bbduk-path, the BBDUK_PATH environment variable, or PATH.
BBDUK_EXECUTABLE_NAMES = ("bbduk.sh", "bbduk", "bbduk.bat")


def resolve_bbduk_path(explicit_path=None):
    """Resolve the BBDuk executable path from CLI, environment, or PATH."""
    candidates = []
    if explicit_path:
        candidates.append(explicit_path)
    env_path = os.environ.get("BBDUK_PATH")
    if env_path:
        candidates.append(env_path)
    candidates.extend(BBDUK_EXECUTABLE_NAMES)

    for candidate in candidates:
        expanded = os.path.expanduser(candidate)
        has_path_separator = (
            os.path.sep in expanded
            or (os.path.altsep is not None and os.path.altsep in expanded)
        )
        if has_path_separator and os.path.isfile(expanded):
            return expanded

        found = shutil.which(expanded)
        if found:
            return found

    raise FileNotFoundError(
        "BBDuk executable was not found. Use the recommended installation "
        "(`conda env create -f environment.yml`), which installs BBDuk "
        "automatically. If you already created the conda environment, run "
        "`conda activate readtrail` before launching inference."
    )

# Global constants
SENTINEL          = "SENTINEL"
CHUNK_SIZE        = 50000
MEMORY_THRESHOLD  = 0.85
CPU_THRESHOLD     = 0.90
MAX_PROCESSES     = 24
THREADS_PER_PROCESS = 2
MAX_READS_PER_TASK  = 1000000

# Idealized Model v2 constants
P_MIN             = 0.05   # Effectively_Unobservable threshold (configurable via --p-min)
OBS_THRESHOLD     = 0.5    # Observable threshold (configurable via --obs-threshold)
MAX_WEIGHT        = 1.0 / P_MIN   # = 20.0  — logical cap: 1/p_min (not arbitrary)
BASELINE_MIN_LOCI = 5      # minimum loci for sample-level baseline
BASELINE_MIN_RAW  = 10     # minimum raw count per locus for baseline candidacy
BASELINE_MAX_L_FRAC = 0.8  # L < this fraction of mean_read_len for baseline candidate
# Fix 1: relative minor-allele threshold for baseline candidacy
# A locus is excluded from baseline only if the minor allele has >= 2 raw reads
# AND its fraction >= BASELINE_MINOR_FRAC_THRESHOLD. This prevents over-exclusion
# at high coverage (e.g. 2 PCR-noise reads out of 200 = 1% → keep as candidate).
BASELINE_MINOR_FRAC_THRESHOLD = 0.05   # 5% relative threshold
# Fix 3: warning threshold for low-p_star Marginal alleles
LOW_P_MARGINAL_THRESHOLD = 0.15   # p_star < this → Has_Low_P_Marginal warning

# IUPAC dictionaries (unchanged)
dico_comp = {
    'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M",
    "V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"
}
dico_degenerate = {
    'R':'AG', 'Y':'CT', 'S':'CG', 'W':'AT', 'K':'GT', 'M':'AC',
    'B':'CGT', "D":'AGT', 'H':'ACT', "V":"ACG", "N":"ACGT", ".":""
}

# =============================================================================
# Idealized Model v2 — core functions
# =============================================================================

def compute_detectability(fragment_len, read_len, p_min=P_MIN):
    """
    Compute detectability for a SINGLE read of length read_len.

    p*(L, R) = (R - L + 1) / R   for L <= R   (uniform start-site approximation)
    p*(L, R) = 0                  for L >  R

    This must be called with the ACTUAL read_len of each individual read
    (not mean_read_len) to avoid systematic bias in per-read corrections.
    mean_read_len is only used for baseline computation and output summaries.

    Classes (p_min configurable via --p-min, default 0.05):
        Observable              : p*(L) >= OBS_THRESHOLD (default 0.5)
        Marginal                : p_min <= p*(L) < OBS_THRESHOLD
        Effectively_Unobservable: 0 < p*(L) < p_min
        Unobservable            : p*(L) = 0  (L > R)

    Returns (p_star, detectability_class, use_in_idealized).
    """
    L = fragment_len
    R = read_len

    if R <= 0:
        return 0.0, "Unobservable", False
    if L > R:
        return 0.0, "Unobservable", False

    p = (R - L + 1) / R

    if p < p_min:
        return p, "Effectively_Unobservable", False
    elif p < OBS_THRESHOLD:
        return p, "Marginal", True
    else:
        return p, "Observable", True




# =============================================================================
# Helper: Allele_Detection_Status  (ТЗ §3.1 addition)
# =============================================================================
#
# ❗ IMPORTANT SEMANTIC NOTE (ТЗ §7 — strict interpretation of Hidden):
#    Hidden_Allele_Flag and Hidden_Component_Pct in the downstream Statistik
#    script reflect OBSERVABILITY CONSTRAINTS only.  They do NOT prove that an
#    allele is present in the sample.  Allele_Detection_Status likewise
#    describes what the sequencer could see, not biological truth.
#
#    'Not_detected'                 → allele could have been observed (p_i > 0)
#                                     but produced zero reads in this sample.
#    'Not_evaluable_by_read_length' → allele is physically invisible because
#                                     the amplicon exceeds the read length.
#                                     Absence of reads carries no genotypic
#                                     information for this allele.
#    'Detected'                     → Raw_Count > 0; allele is confirmed present.

def get_allele_detection_status(raw_count, p_star, det_class):
    """
    Classify allele observability per ТЗ §3.1.
    Called for EVERY allele in the wide output, including Raw_Count = 0 rows.
    Does NOT alter any calculation (Q_i, E_ideal, etc.).
    """
    if raw_count > 0:
        return "Detected"
    if det_class == "Unobservable" or (p_star is not None and p_star <= 0.0):
        return "Not_evaluable_by_read_length"
    return "Not_detected"

def compute_idealized_support(raw_count, p_star, use_in_idealized,
                              min_raw_for_weight=5):
    """
    Length-bias corrected support: Raw_Count / p*(L).
    Cap = 1 / P_MIN  (logical maximum: correcting to the hardest-observable allele).
    Returns None if allele is Unobservable / Effectively_Unobservable.

    Fix 2: low-count guard (min_raw_for_weight=5).
    If raw_count < min_raw_for_weight the signal is statistically too weak to
    amplify reliably via 1/p_star. The raw count is returned unchanged so the
    allele remains visible in the data but cannot dominate the normalisation.
    Returns a tuple (support_value, weight_applied, weight_suppressed) to allow
    the caller to record the audit flags.
    """
    if not use_in_idealized or p_star <= 0:
        return None, False, False
    if raw_count < min_raw_for_weight:
        # Low-count: return raw, do NOT amplify
        return float(raw_count), False, True
    cap = 1.0 / P_MIN          # logical cap: MAX_WEIGHT = 1/p_min
    support = min(raw_count / p_star, raw_count * cap)
    return support, True, False


def compute_expected_reads(baseline_sample, p_star, use_in_idealized,
                           locus_factor=1.0):
    """
    Expected reads for one allele:
        E_i = Baseline_sample * p*(L_i) * locus_factor

    locus_factor accounts for locus-specific amplification efficiency.
    Default = 1.0 (no correction). Set from cohort locus calibration when available.

    Returns None if allele is not observable or baseline unavailable.
    """
    if not use_in_idealized or p_star <= 0 or baseline_sample is None:
        return None
    return baseline_sample * p_star * locus_factor


def compute_observation_ratio(raw_count, expected_reads):
    """
    Observation_Ratio = Raw_Count / Expected_Reads.
    Returns None if Expected_Reads is None or zero.
    """
    if expected_reads is None or expected_reads <= 0:
        return None
    return raw_count / expected_reads


def compute_locus_factor(locus, cohort_locus_medians, global_median):
    """
    Locus-specific amplification factor relative to cohort global median.

    locus_factor = median(total_raw_count for locus across cohort) / global_median

    This corrects for GC-bias, primer efficiency, and other locus-specific effects.
    Returns 1.0 if data unavailable (neutral correction).
    """
    locus_med = cohort_locus_medians.get(locus)
    if locus_med is None or global_median is None or global_median <= 0:
        return 1.0
    return locus_med / global_median


def compute_baseline(allele_counts_by_locus, locus_sizes, mean_read_len,
                     locus_statuses=None, p_min=P_MIN):
    """
    Compute Baseline_sample = median(Total_Raw_Count) over high-confidence loci.

    Candidate locus criteria (ALL must be met):
        1. Fragment length L < BASELINE_MAX_L_FRAC * mean_read_len   (short enough)
        2. Total raw count > BASELINE_MIN_RAW                        (well covered)
        3. p*(L, mean_read_len) >= OBS_THRESHOLD                     (Observable class)
        4. No problematic existing status flags
        5. Monoallelic OR minor allele raw support < 2               (clean signal)

    Returns:
        baseline_value (float or None),
        baseline_loci_count (int),
        baseline_method: 'sample_median' | 'cohort_median' | 'unavailable'

    NOTE: mean_read_len is used here as a proxy for per-read R since we are
    characterising the locus globally, not per-read. Per-read corrections
    use actual read_len in process_chunk().
    """
    if mean_read_len <= 0:
        return None, 0, "unavailable"

    bad_statuses = {
        "LOW_SUPPORT", "AMBIGUOUS_DOMINANCE",
        "POSSIBLE_CONTAMINATION", "STRONG_MIXTURE"
    }

    candidate_totals = []
    for locus, allele_dict in allele_counts_by_locus.items():
        if not allele_dict:
            continue
        if locus_statuses and locus_statuses.get(locus) in bad_statuses:
            continue

        frag_len = locus_sizes.get(locus)
        if frag_len is None:
            continue
        if frag_len >= BASELINE_MAX_L_FRAC * mean_read_len:
            continue

        p_star, det_class, _ = compute_detectability(frag_len, mean_read_len, p_min)
        if det_class != "Observable":
            continue

        total_raw = sum(allele_dict.values())
        if total_raw <= BASELINE_MIN_RAW:
            continue

        sorted_counts = sorted(allele_dict.values(), reverse=True)
        # Fix 1: hybrid minor-allele filter — require BOTH absolute and relative threshold.
        # Old logic (minor >= 2 absolute) was too strict at high coverage:
        #   e.g. 2 PCR-noise reads out of 200 (1%) excluded a perfectly clean locus.
        # New logic: exclude only when minor signal is both >= 2 reads AND >= 5% of total.
        if len(sorted_counts) > 1:
            minor_raw      = sorted_counts[1]
            minor_fraction = minor_raw / total_raw if total_raw > 0 else 0.0
            if minor_raw >= 2 and minor_fraction >= BASELINE_MINOR_FRAC_THRESHOLD:
                continue   # genuinely multi-allelic → skip

        candidate_totals.append(total_raw)

    n = len(candidate_totals)
    if n >= BASELINE_MIN_LOCI:
        candidate_totals.sort()
        mid = n // 2
        baseline = (candidate_totals[mid] if n % 2 == 1
                    else (candidate_totals[mid - 1] + candidate_totals[mid]) / 2)
        return baseline, n, "sample_median"

    return None, n, "cohort_median" if n > 0 else "unavailable"


def compute_locus_completeness(allele_data, baseline_sample, locus_factor=1.0):
    """
    Locus-level completeness metrics.

    allele_data: list of dicts with keys:
        raw_count, p_star, use_in_idealized, expected_reads, det_class

    locus_factor: cohort locus amplification factor (default 1.0).
        Used in proxy_expected so the critical-unobservable threshold is
        allele-specific rather than a flat P_MIN proxy.

    Incomplete_Denominator logic (PATCH 1 / PATCH 6):
        Only set True when there is a CRITICAL unobservable allele —
        i.e. one with det_class in {Unobservable, Effectively_Unobservable}
        AND its proxy expected reads exceed 30 % of the sample baseline.
        proxy_expected uses the allele's own p_star (or P_MIN as floor)
        so long alleles are judged individually, not with a flat threshold.
        A Marginal allele with p close to p_min does NOT trigger this flag.
    """
    observed_total = sum(d['raw_count'] for d in allele_data)
    expected_obs   = sum(
        d['expected_reads'] for d in allele_data
        if d['expected_reads'] is not None
    )

    # Critical unobservable: det_class is hard-unobservable AND the allele
    # would have been expected to carry reads if it were shorter.
    # PATCH 1: proxy uses allele-specific p_star (floor P_MIN) × locus_factor.
    # PATCH 6: threshold is 30 % of baseline_sample, not the fixed BASELINE_MIN_RAW.
    has_critical_unobservable = False
    baseline_s = baseline_sample or 0
    for d in allele_data:
        if d['det_class'] in ("Unobservable", "Effectively_Unobservable"):
            p_eff = max(d.get('p_star', 0.0), P_MIN)
            proxy_expected = baseline_s * p_eff * locus_factor
            if proxy_expected > baseline_s * 0.3:
                has_critical_unobservable = True
                break

    obs_completeness = (observed_total / expected_obs) if expected_obs > 0 else None

    # Fix 3: soft warning for near-unobservable Marginal alleles.
    # These don't trigger Incomplete_Denominator (that would be over-strict),
    # but they pose an interpretability risk because the denominator is not
    # fully reliable. A confidence penalty is applied downstream.
    has_low_p_marginal   = False
    low_p_marginal_count = 0
    for d in allele_data:
        if (d.get('det_class') == "Marginal"
                and d.get('p_star') is not None
                and d['p_star'] < LOW_P_MARGINAL_THRESHOLD):
            has_low_p_marginal    = True
            low_p_marginal_count += 1

    # PATCH 10: derive Base_Status so assign_dropout_status can integrate it.
    # Simple heuristic pre-classification before dropout layer is applied:
    #   HIGH_CONFIDENCE_CLEAN  : completeness ≥ 0.8, no unobservable alleles
    #   LOW_COMPLETENESS       : completeness < 0.5
    #   INCOMPLETE_LOCUS       : Incomplete_Denominator = True
    #   UNKNOWN                : insufficient data
    if obs_completeness is None:
        base_status = "UNKNOWN"
    elif has_critical_unobservable:
        base_status = "INCOMPLETE_LOCUS"
    elif obs_completeness >= 0.8 and not any(not d['use_in_idealized'] for d in allele_data):
        base_status = "HIGH_CONFIDENCE_CLEAN"
    elif obs_completeness < 0.5:
        base_status = "LOW_COMPLETENESS"
    else:
        base_status = "UNKNOWN"

    return {
        "Observed_Total_Raw":          observed_total,
        "Observable_Expected_Total":   round(expected_obs, 4) if expected_obs else 0,
        "Locus_Completeness_Obs":      round(obs_completeness, 4) if obs_completeness is not None else None,
        "Incomplete_Denominator":      has_critical_unobservable,
        "Completeness_Is_Lower_Bound": has_critical_unobservable,
        "Has_Unobservable_Expected":   any(not d['use_in_idealized'] for d in allele_data),
        "Base_Status":                 base_status,   # PATCH 10
        # Fix 3: soft Marginal warning
        "Has_Low_P_Marginal":          has_low_p_marginal,
        "Low_P_Marginal_Count":        low_p_marginal_count,
    }


def is_allele_dropout(raw_count, expected_reads, baseline_sample, obs_ratio_threshold=0.1):
    """
    Per-allele dropout flag (ПРОБЛЕМА №4 fix).

    An allele is considered a dropout candidate when:
        obs_ratio < threshold  AND  raw_count < baseline * 0.2

    Both conditions must hold simultaneously to avoid false positives
    in low-coverage samples where variance is high.
    Returns (is_dropout, obs_ratio).
    """
    if expected_reads is None or expected_reads <= 0:
        return False, None
    obs_ratio = raw_count / expected_reads
    raw_threshold = (baseline_sample * 0.2) if baseline_sample else 0
    # PATCH 4: third guard — expected signal must itself be meaningful
    # (≥ 30 % of baseline). Prevents dropout flags when the allele was
    # barely expected in the first place (e.g. very low locus_factor).
    expected_significant = expected_reads >= (baseline_sample or 0) * 0.3
    dropout = (obs_ratio < obs_ratio_threshold) and (raw_count < raw_threshold) and expected_significant
    return dropout, obs_ratio


def assign_dropout_status(locus_completeness_dict, allele_data, existing_status,
                          baseline_sample=None):
    """
    Assign Dropout_Adjusted_Status (ПРОБЛЕМА №5 fix).

    Rules (priority order):
    1. PARTIALLY_UNOBSERVABLE_MIXTURE:
         - Observed monoallelic (n_observed == 1)
         - AND at least 1 allele with det_class Unobservable/Effectively_Unobservable
           (i.e. a physically hidden expected allele exists in this locus)
         This is the Ft03 case: one allele visible, one hidden due to length.

    2. POTENTIAL_ALLELE_DROPOUT:
         - Observed monoallelic
         - The visible allele is Observable (good detectability)
         - Locus_Completeness_Obs < 0.5  (we expected more reads than we got)

    3. UNRESOLVED_DUE_TO_DETECTABILITY:
         - Incomplete_Denominator = True AND existing status is not already clean

    4. Fall-through: return existing_status.
       HIGH_CONFIDENCE_CLEAN is blocked if Incomplete_Denominator = True.
    """
    completeness = locus_completeness_dict.get("Locus_Completeness_Obs")
    n_observed   = sum(1 for d in allele_data if d.get('raw_count', 0) > 0)

    # Partition alleles by detectability
    # PATCH 5: hard_unobservable covers both Unobservable AND Effectively_Unobservable
    hard_unobservable = [
        d for d in allele_data
        if d.get('det_class') in ("Unobservable", "Effectively_Unobservable")
    ]
    observable = [d for d in allele_data if d.get('use_in_idealized')]

    # Rule 1 — Ft03 case: one allele seen, but at least one is physically hidden
    if n_observed == 1 and hard_unobservable:
        return "PARTIALLY_UNOBSERVABLE_MIXTURE"

    # Rule 2 — dropout by completeness: allele present but under-represented
    if (n_observed == 1
            and observable
            and completeness is not None
            and completeness < 0.5):
        return "POTENTIAL_ALLELE_DROPOUT"

    # Rule 3 — incomplete denominator blocks clean call
    if locus_completeness_dict.get("Incomplete_Denominator"):
        if existing_status == "HIGH_CONFIDENCE_CLEAN":
            return "UNRESOLVED_DUE_TO_DETECTABILITY"
        if existing_status not in ("PARTIALLY_UNOBSERVABLE_MIXTURE",
                                   "POTENTIAL_ALLELE_DROPOUT"):
            return "UNRESOLVED_DUE_TO_DETECTABILITY"

    return existing_status


# =============================================================================
# Helper functions (unchanged from original)
# =============================================================================

def timethis(func):
    def foo(*args, **kwargs):
        t1 = time.time()
        vals = func(*args, **kwargs)
        t2 = time.time()
        logger.info(f"Function='{func.__name__}' took {t2-t1:.2f} seconds")
        return vals
    return foo

def get_system_resources():
    mem = psutil.virtual_memory()
    cpu_usage = psutil.cpu_percent(interval=0.1)
    available_cores = psutil.cpu_count(logical=True)
    return {
        'available_memory': mem.available / (1024 ** 3),
        'total_memory': mem.total / (1024 ** 3),
        'memory_percent': mem.percent,
        'cpu_usage': cpu_usage,
        'available_cores': available_cores
    }

def adjust_resources(tasks_count, threads_per_task):
    resources = get_system_resources()
    logger.info(f"System resources: {resources}")
    max_processes = min(
        MAX_PROCESSES,
        resources['available_cores'] // max(1, threads_per_task),
        tasks_count
    ) or 1
    mem_per_process = resources['available_memory'] / max_processes
    if mem_per_process < 4:
        max_processes = max(1, int(resources['available_memory'] // 4))
    if resources['memory_percent'] > MEMORY_THRESHOLD * 100:
        logger.warning("Memory usage exceeds threshold. Reducing processes.")
        max_processes = max(1, max_processes // 2)
    if resources['cpu_usage'] > CPU_THRESHOLD * 100:
        logger.warning("CPU usage exceeds threshold. Reducing processes.")
        max_processes = max(1, max_processes // 2)
    logger.info(f"Adjusted to {max_processes} processes with {threads_per_task} threads each")
    return max_processes, threads_per_task

def count_reads(file_path: str) -> int:
    count = 0
    try:
        opener = gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')
        with opener as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    count += 1
    except Exception as e:
        logger.error(f"Error counting reads in {file_path}: {e}")
        return 0
    logger.info(f"File {file_path} contains {count} reads")
    return count

def split_file_tasks(file_paths, max_reads_per_task):
    tasks = []
    for file_path in file_paths:
        total_reads = count_reads(file_path)
        if total_reads <= max_reads_per_task:
            tasks.append((file_path, 0, total_reads))
        else:
            num_chunks = math.ceil(total_reads / max_reads_per_task)
            for i in range(num_chunks):
                tasks.append((file_path, i * max_reads_per_task,
                              min(max_reads_per_task, total_reads - i * max_reads_per_task)))
    return tasks

def build_dictionnary(bin_file_path):
    logger.info(f"Building binning dictionary from {bin_file_path}")
    try:
        with open(bin_file_path, "r", encoding='utf-8') as f:
            bin_content = f.read().replace("\t",";").replace(",",";").replace(" ",";").replace("\r","").split("\n")
            bin_content = [line for line in bin_content if line.strip()]
            bin_file_data = [primer.split(";") for primer in bin_content]
        dico_bin = {}
        for primer in bin_file_data:
            if len(primer) < 3:
                logger.warning(f"Skipping malformed binning line: {';'.join(primer)}")
                continue
            locus, size_range, allele_val = map(str.strip, primer[:3])
            if "-" in size_range:
                try:
                    start_str, end_str = size_range.split("-")
                    start, end = int(start_str.strip()), int(end_str.strip())
                    if start > end:
                        start, end = end, start
                        logger.warning(f"Corrected reversed range '{size_range}' for locus {locus}.")
                    tmp = [(str(e), allele_val) for e in range(start, end + 1)]
                except ValueError:
                    logger.warning(f"Invalid range '{size_range}' for locus {locus}. Skipping.")
                    continue
                if locus in dico_bin:
                    dico_bin[locus].extend(tmp)
                else:
                    dico_bin[locus] = tmp
            else:
                try:
                    size_val = str(int(size_range))
                    if locus in dico_bin:
                        dico_bin[locus].append((size_val, allele_val))
                    else:
                        dico_bin[locus] = [(size_val, allele_val)]
                except ValueError:
                    logger.warning(f"Invalid size '{size_range}' for locus {locus}. Skipping.")
                    continue
        logger.info(f"Binning dictionary built with {len(dico_bin)} loci")
        return dico_bin
    except Exception as e:
        logger.error(f"Error reading binning file {bin_file_path}: {e}")
        sys.exit(2)

def clean_primers(primers_list):
    logger.info("Cleaning primers")
    tmp = []
    allowed_chars = set("ACGTMRWSYKVHDBN.")
    for primers_line in primers_list:
        primers = [p.strip() for p in re.split(r'[;, \t]', primers_line) if p.strip()]
        if len(primers) >= 3:
            if all(c in allowed_chars for c in primers[1].upper()) and all(c in allowed_chars for c in primers[2].upper()):
                primers = [primers[0]] + [primer.upper() for primer in primers[1:]]
                tmp.append(primers)
            else:
                logger.warning(f"Skipping primer line with invalid characters: {primers_line}")
        else:
            logger.warning(f"Skipping malformed primer line: {primers_line}")
    logger.info(f"Cleaned {len(tmp)} primers")
    return tmp

_COMP_TABLE = str.maketrans(
    "ACGTMRWSYKVHDBNacgtmrwsykvhdbn",
    "TGCAKYWSRMBDHVNtgcakywsrmbdhvn"
)

def inverComp(seq: str) -> str:
    seq = seq.upper()
    try:
        return seq.translate(_COMP_TABLE)[::-1]
    except Exception:
        unknown_chars = {nuc for nuc in seq if nuc not in dico_comp}
        logger.warning(f"Unknown nucleotide(s) '{','.join(unknown_chars)}' in {seq}. Replacing with 'X'.")
        return "".join([dico_comp.get(nuc, 'X') for nuc in seq[::-1]])

def degenerated_primers(primer):
    degenerated_nucs = []
    positions = []
    primers = set()
    base_primer_list = list(primer)
    has_degenerate = False
    for i, nuc in enumerate(primer):
        if nuc in dico_degenerate:
            has_degenerate = True
            positions.append(i)
            degenerated_nucs.append(dico_degenerate[nuc])
    if not has_degenerate:
        return [primer]
    combinations = list(itertools.product(*degenerated_nucs))
    for comb in combinations:
        current_primer_list = list(base_primer_list)
        for pos, nuc in zip(positions, comb):
            current_primer_list[pos] = nuc
        primers.add("".join(current_primer_list))
    return list(primers)

def binning_correction(primer_locus, size, sizeU, dico_bin):
    if primer_locus in dico_bin:
        correction_found = False
        for bin_size_str, corrected_u_str in dico_bin[primer_locus]:
            try:
                bin_size = int(bin_size_str)
                if size == bin_size:
                    sizeU = float(corrected_u_str.replace('u','').strip())
                    correction_found = True
                    break
            except ValueError:
                logger.warning(f"Invalid size ('{bin_size_str}') or allele ('{corrected_u_str}') for locus {primer_locus}")
                continue
        if not correction_found:
            logger.debug(f"No correction found for locus {primer_locus}, size {size}")
    else:
        logger.debug(f"Locus {primer_locus} not found in binning dictionary")
    return sizeU

def get_flanking(seq, primers, pos1, pos2, splitted, flanking_len, config="A"):
    """
    Extract flanking sequences around the primer pair.
    config: "A" = FWD+RC_REV (left=pos1=FWD start), "B" = REV+RC_FWD (left=pos1=REV start).
    In both cases pos1 is the left primer start and pos2 is the right primer start,
    so the orientation logic is identical; we keep the old structure.
    """
    if flanking_len <= 0:
        return ["", ""]
    p1_len = len(primers[1])
    p2_len = len(primers[2])
    seq_len = len(seq)
    f1, f2 = "", ""
    try:
        if splitted:
            # left primer near end of read, right primer at start
            f1_end   = pos1
            f1_start = max(0, f1_end - flanking_len)
            f1       = seq[f1_start:f1_end]
            f2_start = pos2 + p2_len
            f2_end   = min(seq_len, f2_start + flanking_len)
            f2       = seq[f2_start:f2_end]
        else:
            first_primer_pos  = min(pos1, pos2)
            second_primer_pos = max(pos1, pos2)
            second_primer_len = p1_len if second_primer_pos == pos1 else p2_len
            f1_end   = first_primer_pos
            f1_start = max(0, f1_end - flanking_len)
            f1       = seq[f1_start:f1_end]
            f2_start = second_primer_pos + second_primer_len
            f2_end   = min(seq_len, f2_start + flanking_len)
            f2       = seq[f2_start:f2_end]
    except Exception as e:
        logger.warning(f"Error calculating flanking sequence: {e}")
        return ["FLANK_ERR", "FLANK_ERR"]
    return [f1, f2]

def run_bbduk(file_path, output_dir, primers, nbmismatch, threads,
              bbduk_path, max_distance=1000):
    logger.info(f"Running BBDuk for {file_path}")
    name = os.path.basename(file_path)
    name, ext = os.path.splitext(name)
    if ext == '.gz':
        name, _ = os.path.splitext(name)
    matched_files = {}
    for primer in primers:
        primer_name = primer[0].replace('_', '-')
        fwd = primer[1][-18:]
        rev = primer[2][-18:]
        primer_string = f"{fwd},{rev}"
        out_matched = os.path.join(output_dir, f"{name}_{primer_name}_MATCHED_{uuid.uuid4().hex[:8]}.fastq")
        cmd = [
            bbduk_path,
            f"in={file_path}",
            f"outm={out_matched}",
            f"literal={primer_string}",
            "mm=f",
            f"k={min(18, len(primer[1]))}",
            "minlen=1",
            f"hdist={nbmismatch}",
            "forbidn=t",
            f"threads={threads}",
            "overwrite=t",
            "Xmx=4g",
        ]
        logger.info(f"Running BBDuk command for {primer_name}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"BBDuk command failed for {primer_name}: {result.stderr[:500]}")
            continue
        if os.path.exists(out_matched):
            file_size = os.path.getsize(out_matched)
            if file_size == 0:
                logger.warning(f"Created 0-byte file {out_matched}. Removing.")
                os.remove(out_matched)
            else:
                matched_files[primer_name] = out_matched
                logger.debug(f"Matched file created: {out_matched} ({file_size} bytes)")
        else:
            logger.warning(f"Matched file {out_matched} was not created for {primer_name}.")
    return matched_files

def count_nucleotides(file_path: str):
    total_nt   = 0
    read_count = 0
    try:
        opener = gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')
        with opener as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total_nt   += len(line.rstrip('\n'))
                    read_count += 1
    except Exception as e:
        logger.error(f"Error counting nucleotides in {file_path}: {e}")
        return 0, 0
    logger.info(f"Total nucleotides: {total_nt}, Read count: {read_count}")
    return total_nt, read_count

def hamming_mismatches(seq_a: str, seq_b: str) -> int:
    if len(seq_a) != len(seq_b):
        return max(len(seq_a), len(seq_b))
    a = np.frombuffer(seq_a.encode(), dtype=np.uint8)
    b = np.frombuffer(seq_b.encode(), dtype=np.uint8)
    return int(np.count_nonzero(a != b))

# =============================================================================
# Strict-stage full-length primer validation (v2.4)
# Replaces seed-automaton as the basis for accept/reject decisions.
# BBDuk remains the coarse filter; these functions handle exact validation.
# =============================================================================

def collect_full_length_primer_hits(seq: str, fwd_primer: str, rev_primer: str,
                                    nbmismatch: int,
                                    rc_fwd: str = None, rc_rev: str = None) -> dict:
    """
    Search all 4 primer forms in *seq* at full primer length with <= nbmismatch
    Hamming distance. Returns a dict of hit lists keyed by form name.

    rc_fwd / rc_rev can be pre-computed outside the read loop for speed.
    If not provided, they are computed here (backwards-compatible).
    """
    if rc_fwd is None:
        rc_fwd = inverComp(fwd_primer)
    if rc_rev is None:
        rc_rev = inverComp(rev_primer)
    seq_len = len(seq)

    hits = {"FWD": [], "REV": [], "RC_FWD": [], "RC_REV": []}

    for form_name, query in (("FWD",    fwd_primer),
                              ("REV",    rev_primer),
                              ("RC_FWD", rc_fwd),
                              ("RC_REV", rc_rev)):
        qlen = len(query)
        if qlen == 0 or seq_len < qlen:
            continue

        # preflight: guaranteed-matching prefix (pigeonhole principle)
        guaranteed = max(1, qlen - nbmismatch)
        seed = query[:min(4, guaranteed)]
        if seed and seed not in seq:
            continue

        q_arr = np.frombuffer(query.encode(), dtype=np.uint8)
        max_start = seq_len - qlen + 1

        for start in range(max_start):
            window = seq[start:start + qlen]
            w_arr = np.frombuffer(window.encode(), dtype=np.uint8)
            mm = int(np.count_nonzero(q_arr != w_arr))
            if mm <= nbmismatch:
                hits[form_name].append({"start": start, "length": qlen,
                                        "mismatches": mm})
                if mm == 0 and nbmismatch == 0:
                    break  # exact match found — no need to continue

    return hits


def build_candidate_pairs(hits, fwd_primer, rev_primer):
    """
    Enumerate all orientation-compatible primer pairs from the hit lists.

    Config A: FWD (left anchor) + RC_REV (right anchor)  — standard orientation
    Config B: REV (left anchor) + RC_FWD (right anchor)  — reverse orientation

    For each config every combination of left-hit × right-hit is a candidate.
    Returns a list of candidate dicts.
    """
    candidates = []

    # Config A: FWD + RC_REV
    for lh in hits["FWD"]:
        for rh in hits["RC_REV"]:
            candidates.append({
                "config":           "A",
                "left_type":        "FWD",
                "right_type":       "RC_REV",
                "left_start":       lh["start"],
                "right_start":      rh["start"],
                "left_len":         lh["length"],
                "right_len":        rh["length"],
                "left_mismatches":  lh["mismatches"],
                "right_mismatches": rh["mismatches"],
                "total_mismatches": lh["mismatches"] + rh["mismatches"],
            })

    # Config B: REV + RC_FWD
    for lh in hits["REV"]:
        for rh in hits["RC_FWD"]:
            candidates.append({
                "config":           "B",
                "left_type":        "REV",
                "right_type":       "RC_FWD",
                "left_start":       lh["start"],
                "right_start":      rh["start"],
                "left_len":         lh["length"],
                "right_len":        rh["length"],
                "left_mismatches":  lh["mismatches"],
                "right_mismatches": rh["mismatches"],
                "total_mismatches": lh["mismatches"] + rh["mismatches"],
            })

    return candidates


def evaluate_candidate_pair(candidate, seq, contig, max_distance):
    """
    Compute fragment geometry for one candidate pair.

    Fragment size = end_of_right_primer − start_of_left_primer
    (i.e. full amplicon length inclusive of both primers).

    Split case: right primer wraps around the read end (circular / split-read).
    Allowed only when contig == False.

    Returns a result dict; "is_valid" is False when the pair must be discarded.
    """
    seq_len   = len(seq)
    ls        = candidate["left_start"]
    rs        = candidate["right_start"]
    ll        = candidate["left_len"]
    rl        = candidate["right_len"]
    splitted  = False
    size      = -1

    if ls < rs:
        # Normal: left primer starts before right primer
        size = (rs + rl) - ls
    elif ls > rs and not contig:
        # Wrap-around: left primer is near the end of the read,
        # right primer appears at the beginning (split/circular)
        size     = rl + (seq_len - ls)
        splitted = True
    # ls == rs → same position → invalid

    if size < 1 or size > max_distance:
        return {"is_valid": False, "config": candidate["config"],
                "size": size, "splitted": splitted}

    # Build insert sequence (full amplicon, always in FWD orientation)
    if splitted:
        raw = seq[ls:] + seq[:(rs + rl)]
        # For Config B the amplicon is on the antisense strand — reverse complement
        insert = raw if candidate["config"] == "A" else inverComp(raw)
    else:
        raw = seq[ls:(rs + rl)]
        insert = raw if candidate["config"] == "A" else inverComp(raw)

    return {
        "is_valid":         True,
        "config":           candidate["config"],
        "left_type":        candidate["left_type"],
        "right_type":       candidate["right_type"],
        "size":             size,
        "splitted":         splitted,
        "insert":           insert,
        "p1":               ls,            # left primer start (for flanking, sorting)
        "p2":               rs,            # right primer start
        "left_len":         ll,
        "right_len":        rl,
        "left_mismatches":  candidate["left_mismatches"],
        "right_mismatches": candidate["right_mismatches"],
        "total_mismatches": candidate["total_mismatches"],
    }


def select_best_candidate(valid_candidates):
    """
    Deterministic selection from multiple valid candidates.

    Priority (ascending sort key, so first element is best):
      1. non-split preferred over split        (splitted: False < True)
      2. fewer total mismatches                (total_mismatches)
      3. smaller fragment size                 (size)
      4. smaller left primer start position    (p1)
    """
    return min(
        valid_candidates,
        key=lambda c: (c["splitted"], c["total_mismatches"], c["size"], c["p1"])
    )

def fastq_to_fasta(fastq_file, fasta_file):
    logger.info(f"Converting {fastq_file} to {fasta_file}")
    try:
        with open(fastq_file, "r") as handle_in, open(fasta_file, "w") as handle_out:
            count = SeqIO.convert(handle_in, "fastq", handle_out, "fasta")
        logger.debug(f"Converted {count} records to {fasta_file}")
        return True
    except Exception as e:
        logger.error(f'Error converting to fasta: {e}')
        return False

def save_matched_reads(matched_files, matched_reads, output_path):
    fastq_output_folder = os.path.join(output_path, "matched_fastq")
    fasta_output_folder = os.path.join(output_path, "matched_fasta")
    os.makedirs(fastq_output_folder, exist_ok=True)
    os.makedirs(fasta_output_folder, exist_ok=True)
    for primer_name in matched_files:
        fastq_output_path = matched_files[primer_name]
        fastq_copy_path = os.path.join(fastq_output_folder, f"{os.path.basename(fastq_output_path)}_{uuid.uuid4().hex[:8]}.fastq")
        try:
            shutil.copy(fastq_output_path, fastq_copy_path)
            file_size = os.path.getsize(fastq_copy_path)
            if file_size == 0:
                logger.warning(f"Copied 0-byte FASTQ file {fastq_copy_path}. Removing.")
                os.remove(fastq_copy_path)
            else:
                logger.debug(f"Copied FASTQ to {fastq_copy_path} ({file_size} bytes)")
                fasta_output_path = os.path.join(fasta_output_folder,
                    f"{os.path.basename(fastq_copy_path).replace('.fastq', '.fasta')}")
                if fastq_to_fasta(fastq_copy_path, fasta_output_path):
                    fasta_size = os.path.getsize(fasta_output_path)
                    if fasta_size == 0:
                        logger.warning(f"Created 0-byte FASTA file {fasta_output_path}. Removing.")
                        os.remove(fasta_output_path)
                    else:
                        logger.debug(f"Converted matched reads to {fasta_output_path} ({fasta_size} bytes)")
        except Exception as e:
            logger.error(f"Error saving matched reads for {primer_name}: {e}")


# =============================================================================
# Core processing — read-level loop (unchanged logic, extended output)
# =============================================================================

def process_chunk(records, primers, nbmismatch, contig, binning, sequence_prefix,
                  flanking, output_path, round_val, max_distance, dico_bin,
                  chunk_index, primer_name, chunk_uuid):
    """
    Process a chunk of reads for one primer pair.

    v2.4: Two-stage validation:
      Stage 1 — BBDuk (upstream, not here): coarse 18-nt 3'-end filter, hdist=2.
      Stage 2 — Full-length strict validation (this function):
        collect_full_length_primer_hits → build_candidate_pairs →
        evaluate_candidate_pair → select_best_candidate.

    The seed-automaton is no longer used as the decision basis.
    All downstream calculations (sizeU, detectability, idealized model) unchanged.
    """
    logger.debug(f"Processing chunk {chunk_index} for {primer_name} with {len(records)} records")
    all_sequence_results      = defaultdict(list)
    allele_counts             = defaultdict(lambda: defaultdict(int))
    corrected_allele_counts   = defaultdict(lambda: defaultdict(float))
    matched_reads             = defaultdict(list)
    locus_sizes               = defaultdict(dict)
    p_star_sums               = defaultdict(lambda: defaultdict(float))
    p_star_counts_acc         = defaultdict(lambda: defaultdict(int))
    det_class_votes           = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    # v2.4 funnel counters (diagnostic)
    reads_after_bbduk             = len(records)
    reads_passing_strict          = 0
    reads_rejected_no_pair        = 0
    reads_rejected_invalid_size   = 0
    reads_with_multiple_pairs     = 0

    filtered_fastq_folder = os.path.join(output_path, "filtered_both_primers_fastq")
    filtered_fasta_folder = os.path.join(output_path, "filtered_both_primers_fasta")
    os.makedirs(filtered_fastq_folder, exist_ok=True)
    os.makedirs(filtered_fasta_folder, exist_ok=True)

    primer = next((p for p in primers if p[0].replace('_', '-') == primer_name), None)
    if not primer:
        logger.error(f"Primer {primer_name} not found")
        return (all_sequence_results, allele_counts, corrected_allele_counts,
                matched_reads, locus_sizes, p_star_sums, det_class_votes)

    primer_full_name  = primer[0]
    primer_short_name = primer_full_name.split('_')[0]
    primer_file_name  = primer_name

    filtered_fastq_path = os.path.join(filtered_fastq_folder,
        f"{primer_file_name}_chunk_{chunk_index}_{chunk_uuid}_BOTH_PRIMERS.fastq")
    filtered_fasta_path = os.path.join(filtered_fasta_folder,
        f"{primer_file_name}_chunk_{chunk_index}_{chunk_uuid}_BOTH_PRIMERS.fasta")
    filtered_records = []

    fwd_primer = primer[1]
    rev_primer = primer[2]
    rc_fwd_cached = inverComp(fwd_primer)   # computed once per chunk
    rc_rev_cached = inverComp(rev_primer)   # computed once per chunk

    for s, record in enumerate(records):
        seq      = str(record.seq).upper()
        seq_id   = record.id
        read_len = len(seq)
        current_sequence_name = f"{sequence_prefix}{chunk_index}_{chunk_uuid}_{s+1}"

        matched_reads[primer_file_name].append(record)

        # ── Stage 2: full-length strict validation ────────────────────────────
        try:
            full_hits = collect_full_length_primer_hits(
                seq, fwd_primer, rev_primer, nbmismatch,
                rc_fwd=rc_fwd_cached, rc_rev=rc_rev_cached)
        except Exception as e:
            logger.warning(f"Error in full-length primer search for "
                           f"{primer_file_name} read {seq_id}: {e}")
            reads_rejected_no_pair += 1
            continue

        candidates = build_candidate_pairs(full_hits, fwd_primer, rev_primer)
        if not candidates:
            reads_rejected_no_pair += 1
            continue

        evaluated = []
        for cand in candidates:
            ev = evaluate_candidate_pair(cand, seq, contig, max_distance)
            if ev["is_valid"]:
                evaluated.append(ev)

        if not evaluated:
            reads_rejected_invalid_size += 1
            continue

        if len(evaluated) > 1:
            reads_with_multiple_pairs += 1

        best = select_best_candidate(evaluated)
        reads_passing_strict += 1
        filtered_records.append(record)
        # ── End strict validation ─────────────────────────────────────────────

        size     = best["size"]
        splitted = best["splitted"]
        insert   = best["insert"]
        pos1     = best["p1"]    # left primer start
        pos2     = best["p2"]    # right primer start

        # ── sizeU calculation (unchanged) ─────────────────────────────────────
        try:
            repeat_bp   = float(primer[0].split('_')[1].lower().replace("bp", ""))
            expected_bp = float(primer[0].split('_')[2].lower().replace("bp", ""))
            ref_u       = float(primer[0].split('_')[3].upper().replace("U",  ""))
            sizeU = ref_u + ((size - expected_bp) / repeat_bp) if repeat_bp > 0 else 'DivByZero'
            if binning:
                sizeU = binning_correction(primer_short_name, size, sizeU, dico_bin)
            if round_val > 0:
                if isinstance(sizeU, (int, float)):
                    if sizeU >= math.floor(sizeU) and sizeU < (math.floor(sizeU) + round_val):
                        sizeU = math.floor(sizeU)
                    elif sizeU <= math.ceil(sizeU) and sizeU > (math.ceil(sizeU) - round_val):
                        sizeU = math.ceil(sizeU)
                    else:
                        sizeU = math.floor(sizeU) + 0.5
                    if str(sizeU)[-2:] == '.0':
                        sizeU = int(sizeU)
        except Exception as e:
            logger.warning(f"Error calculating sizeU for {primer_file_name}, {seq_id}: {e}")
            sizeU = 'CalcError'

        if not isinstance(sizeU, (int, float)) or (isinstance(sizeU, float) and math.isnan(sizeU)):
            continue
        # ── End sizeU ────────────────────────────────────────────────────────

        # ── Idealized Model v2 — detectability & corrected weight (unchanged) ─
        p_star, det_class, use_in_idealized = compute_detectability(size, read_len)

        corrected_weight = 0.0
        if use_in_idealized and p_star > 0:
            corrected_weight = min(1.0 / p_star, MAX_WEIGHT)
        elif not use_in_idealized and size > read_len:
            logger.warning(
                f"Fragment longer than read (Unobservable): "
                f"size={size}, read_len={read_len}, read={seq_id}"
            )
        # ── End idealized model ───────────────────────────────────────────────

        result_entry = [
            primer_full_name, pos1, pos2, size, sizeU, current_sequence_name,
            nbmismatch, fwd_primer, "", rev_primer, "", insert,
            # v2 fields:
            p_star, det_class, use_in_idealized,
            # v2.4 per-read audit fields:
            best["config"],                  # Selected_Config
            best["total_mismatches"],        # Selected_Total_Mismatches
            best["splitted"],                # Selected_Is_Split
            len(evaluated),                  # Valid_Candidate_Count
        ]
        if flanking:
            result_entry.extend(
                get_flanking(seq, primer, pos1, pos2, splitted, flanking_len,
                             config=best["config"]))

        all_sequence_results[primer_short_name].append(result_entry)
        allele_counts[primer_short_name][sizeU] += 1

        if use_in_idealized:
            corrected_allele_counts[primer_short_name][sizeU] += corrected_weight

        if sizeU not in locus_sizes[primer_short_name]:
            locus_sizes[primer_short_name][sizeU] = size

        p_star_sums[primer_short_name][sizeU]       += p_star
        p_star_counts_acc[primer_short_name][sizeU] += 1
        det_class_votes[primer_short_name][sizeU][det_class] += 1

    # Log funnel stats for this chunk
    logger.debug(
        f"Chunk {chunk_index} [{primer_name}] funnel: "
        f"after_bbduk={reads_after_bbduk}, "
        f"passing_strict={reads_passing_strict}, "
        f"rejected_no_pair={reads_rejected_no_pair}, "
        f"rejected_bad_size={reads_rejected_invalid_size}, "
        f"multiple_valid_pairs={reads_with_multiple_pairs}"
    )

    # Write filtered reads
    if filtered_records:
        try:
            with open(filtered_fastq_path, 'w') as handle:
                SeqIO.write(filtered_records, handle, 'fastq')
            file_size = os.path.getsize(filtered_fastq_path)
            if file_size == 0:
                logger.warning(f"Created 0-byte FASTQ {filtered_fastq_path}. Removing.")
                os.remove(filtered_fastq_path)
            else:
                logger.debug(f"Saved {len(filtered_records)} reads to {filtered_fastq_path}")
                if fastq_to_fasta(filtered_fastq_path, filtered_fasta_path):
                    fasta_size = os.path.getsize(filtered_fasta_path)
                    if fasta_size == 0:
                        logger.warning(f"Created 0-byte FASTA {filtered_fasta_path}. Removing.")
                        os.remove(filtered_fasta_path)
        except Exception as e:
            logger.error(f"Error writing filtered files for {primer_file_name}: {e}")

    # Return RAW accumulators so process_file can aggregate correctly across chunks.
    # BUG FIX v2.4.1: the previous version pre-computed mean_p_star and dominant_det_class
    # here and returned them. But process_file / _merge_chunk expect RAW structures:
    #   c_ps = p_star_sums        -> dict {locus: {allele: float_sum}}
    #   c_dc = det_class_votes    -> dict {locus: {allele: {class_str: count}}}
    # Returning pre-computed means caused:
    #   c_dc[locus][allele] = "Observable"  (a string, not a dict)
    #   process_file then called sum("Observable".values()) -> AttributeError
    #   which was caught by the broad except and silently returned an empty result,
    #   making ALL chunks appear to have zero alleles.
    return (all_sequence_results, allele_counts, corrected_allele_counts,
            matched_reads, locus_sizes, p_star_sums, det_class_votes)


def process_bbduk_results(matched_files, primers, nbmismatch, contig, binning,
                          sequence_prefix, flanking, output_path, round_val=0.25,
                          max_distance=1000, dico_bin=None):
    """Process all matched files for a sample. Returns extended stats including locus_sizes."""
    logger.info("Starting process_bbduk_results")
    all_sequence_results    = defaultdict(list)
    allele_counts           = defaultdict(lambda: defaultdict(int))
    corrected_allele_counts = defaultdict(lambda: defaultdict(float))
    matched_reads           = defaultdict(list)
    locus_sizes_agg         = defaultdict(dict)
    # ПРОБЛЕМА 1 fix: aggregate mean_p_star and dominant_det_class across chunks
    p_star_sums_agg         = defaultdict(lambda: defaultdict(float))
    p_star_counts_agg       = defaultdict(lambda: defaultdict(int))
    det_class_votes_agg     = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    total_reads = 0
    total_bases = 0
    min_len     = float('inf')
    max_len     = 0

    def _merge_chunk(cr, ca, cc, cm, cs, c_ps, c_pc, c_dc,
                     local_ar, local_ac, local_cc2, local_mr, local_ls,
                     local_ps, local_pc, local_dc):
        """Merge one chunk's results into local accumulators."""
        for k, v in cr.items():  local_ar[k].extend(v)
        for k, v in ca.items():
            for allele, cnt in v.items(): local_ac[k][allele] += cnt
        for k, v in cc.items():
            for allele, wt  in v.items(): local_cc2[k][allele] += wt
        for k, v in cm.items():  local_mr[k].extend(v)
        for locus, sizes in cs.items():
            for allele, sz in sizes.items():
                if allele not in local_ls[locus]:
                    local_ls[locus][allele] = sz
        for locus, ad in c_ps.items():
            for allele, s in ad.items():
                local_ps[locus][allele] += s
        for locus, ad in c_pc.items():
            for allele, n in ad.items():
                local_pc[locus][allele] += n
        for locus, ad in c_dc.items():
            for allele, cls_dict in ad.items():
                for cls, n in cls_dict.items():
                    local_dc[locus][allele][cls] += n

    def process_file(primer_name, matched_file):
        if not os.path.exists(matched_file):
            logger.warning(f"Matched file for {primer_name} not found at {matched_file}.")
            empty = (defaultdict(list), defaultdict(lambda: defaultdict(int)),
                     defaultdict(lambda: defaultdict(float)), defaultdict(list),
                     defaultdict(dict),
                     defaultdict(lambda: defaultdict(float)),
                     defaultdict(lambda: defaultdict(int)),
                     defaultdict(lambda: defaultdict(lambda: defaultdict(int))),
                     (0, 0, 0, 0))
            return empty

        chunk_records = []
        chunk_index   = 0
        total_reads_p = 0
        total_bases_p = 0
        min_len_p     = float('inf')
        max_len_p     = 0

        local_ar  = defaultdict(list)
        local_ac  = defaultdict(lambda: defaultdict(int))
        local_cc2 = defaultdict(lambda: defaultdict(float))
        local_mr  = defaultdict(list)
        local_ls  = defaultdict(dict)
        local_ps  = defaultdict(lambda: defaultdict(float))
        local_pc  = defaultdict(lambda: defaultdict(int))
        local_dc  = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

        try:
            with (gzip.open(matched_file, 'rt') if matched_file.endswith('.gz') else open(matched_file, 'r')) as f:
                for record in SeqIO.parse(f, 'fastq'):
                    rl = len(record.seq)
                    total_bases_p += rl
                    total_reads_p += 1
                    if rl < min_len_p: min_len_p = rl
                    if rl > max_len_p: max_len_p = rl
                    chunk_records.append(record)

                    if len(chunk_records) >= CHUNK_SIZE:
                        chunk_uuid = uuid.uuid4().hex[:8]
                        cr, ca, cc, cm, cs, c_ps, c_dc = process_chunk(
                            chunk_records, primers, nbmismatch, contig, binning,
                            sequence_prefix, flanking, output_path, round_val,
                            max_distance, dico_bin, chunk_index, primer_name, chunk_uuid)
                        # rebuild c_pc from c_ps (counts equal to p_star_counts_acc in chunk)
                        c_pc = defaultdict(lambda: defaultdict(int))
                        for locus, ad in c_ps.items():
                            for allele in ad:
                                c_pc[locus][allele] = sum(c_dc[locus][allele].values())
                        _merge_chunk(cr, ca, cc, cm, cs, c_ps, c_pc, c_dc,
                                     local_ar, local_ac, local_cc2, local_mr, local_ls,
                                     local_ps, local_pc, local_dc)
                        chunk_records = []
                        chunk_index  += 1
                        if get_system_resources()['memory_percent'] > MEMORY_THRESHOLD * 100:
                            logger.warning("Memory threshold exceeded. Pausing.")
                            time.sleep(5)

                if chunk_records:
                    chunk_uuid = uuid.uuid4().hex[:8]
                    cr, ca, cc, cm, cs, c_ps, c_dc = process_chunk(
                        chunk_records, primers, nbmismatch, contig, binning,
                        sequence_prefix, flanking, output_path, round_val,
                        max_distance, dico_bin, chunk_index, primer_name, chunk_uuid)
                    c_pc = defaultdict(lambda: defaultdict(int))
                    for locus, ad in c_ps.items():
                        for allele in ad:
                            c_pc[locus][allele] = sum(c_dc[locus][allele].values())
                    _merge_chunk(cr, ca, cc, cm, cs, c_ps, c_pc, c_dc,
                                 local_ar, local_ac, local_cc2, local_mr, local_ls,
                                 local_ps, local_pc, local_dc)

            logger.info(f"Processed {total_reads_p} reads for {primer_name}")
        except Exception as e:
            logger.error(f"Error processing {matched_file}: {e}")
            traceback.print_exc()   # v2.4.1: always print full traceback so silent fails are visible
            empty = (defaultdict(list), defaultdict(lambda: defaultdict(int)),
                     defaultdict(lambda: defaultdict(float)), defaultdict(list),
                     defaultdict(dict),
                     defaultdict(lambda: defaultdict(float)),
                     defaultdict(lambda: defaultdict(int)),
                     defaultdict(lambda: defaultdict(lambda: defaultdict(int))),
                     (0, 0, 0, 0))
            return empty

        stats = (total_reads_p, total_bases_p, min_len_p, max_len_p)
        return (local_ar, local_ac, local_cc2, local_mr, local_ls,
                local_ps, local_pc, local_dc, stats)

    with ThreadPoolExecutor(max_workers=THREADS_PER_PROCESS * len(matched_files)) as executor:
        futures = [executor.submit(process_file, pn, mf) for pn, mf in matched_files.items()]
        for future in futures:
            result = future.result()
            (local_ar, local_ac, local_cc, local_mr, local_ls,
             local_ps, local_pc, local_dc, stats) = result
            total_reads += stats[0]
            total_bases += stats[1]
            if stats[2] < min_len: min_len = stats[2]
            if stats[3] > max_len: max_len = stats[3]
            for k, v in local_ar.items():
                all_sequence_results[k].extend(v)
            for k, v in local_ac.items():
                for allele, cnt in v.items():
                    allele_counts[k][allele] += cnt
            for k, v in local_cc.items():
                for allele, wt in v.items():
                    corrected_allele_counts[k][allele] += wt
            for k, v in local_mr.items():
                matched_reads[k].extend(v)
            for locus, sizes in local_ls.items():
                for allele, sz in sizes.items():
                    if allele not in locus_sizes_agg[locus]:
                        locus_sizes_agg[locus][allele] = sz
            # Merge p_star accumulators (ПРОБЛЕМА 1 fix)
            for locus, ad in local_ps.items():
                for allele, s in ad.items():
                    p_star_sums_agg[locus][allele] += s
            for locus, ad in local_pc.items():
                for allele, n in ad.items():
                    p_star_counts_agg[locus][allele] += n
            for locus, ad in local_dc.items():
                for allele, cls_dict in ad.items():
                    for cls, n in cls_dict.items():
                        det_class_votes_agg[locus][allele][cls] += n

    mean_len = total_bases / total_reads if total_reads > 0 else 0
    if min_len == float('inf'): min_len = 0

    stats_dict = {
        'total_reads': total_reads,
        'total_bases': total_bases,
        'mean_len':    mean_len,
        'min_len':     min_len,
        'max_len':     max_len,
    }

    # Derive representative fragment size per locus (median across alleles)
    locus_repr_size = {}
    for locus, allele_sizes in locus_sizes_agg.items():
        sizes_list = sorted(allele_sizes.values())
        mid = len(sizes_list) // 2
        locus_repr_size[locus] = (sizes_list[mid] if len(sizes_list) % 2 == 1
                                  else (sizes_list[mid-1] + sizes_list[mid]) / 2)

    # Build mean_p_star and dominant_det_class per allele from aggregated per-read data
    # ПРОБЛЕМА 1 fix: these replace any post-hoc recalculation from mean_read_len
    mean_p_star_agg     = {}   # locus → {allele → mean_p_star across reads}
    dominant_det_cls    = {}   # locus → {allele → most common det_class across reads}
    for locus in set(list(p_star_sums_agg.keys()) + list(allele_counts.keys())):
        mean_p_star_agg[locus]  = {}
        dominant_det_cls[locus] = {}
        for allele in allele_counts.get(locus, {}):
            n   = p_star_counts_agg[locus].get(allele, 0)
            s   = p_star_sums_agg[locus].get(allele, 0.0)
            mean_p_star_agg[locus][allele] = s / n if n > 0 else 0.0
            votes = det_class_votes_agg[locus].get(allele, {})
            dominant_det_cls[locus][allele] = (max(votes, key=votes.get)
                                                if votes else "Unobservable")

    return (
        {k: v for k, v in all_sequence_results.items()},
        {k: dict(v) for k, v in allele_counts.items()},
        {k: dict(v) for k, v in corrected_allele_counts.items()},
        matched_reads,
        stats_dict,
        dict(locus_sizes_agg),
        locus_repr_size,
        mean_p_star_agg,       # NEW: per-allele mean p_star from per-read accumulation
        dominant_det_cls,      # NEW: per-allele dominant detectability class
    )


# =============================================================================
# Per-sample runner
# =============================================================================

def run(sample_name, file_paths, primers, round_val, nbmismatch_max, threads,
        output_queue, contig_setting, binning_setting, sequence_prefix,
        flanking_setting, output_dir, dico_bin, bbduk_path):
    logger.info(f"Processing sample: {sample_name} with files: {file_paths}")
    start_time = time.time()
    try:
        total_nt = 0
        read_count = 0
        matched_files_combined = {}

        file_tasks = split_file_tasks(file_paths, MAX_READS_PER_TASK)
        for file_path, start_read, num_reads in file_tasks:
            nt, rc = count_nucleotides(file_path)
            total_nt   += nt
            read_count += rc
            matched_files = run_bbduk(
                file_path, output_dir, primers, nbmismatch_max, threads,
                bbduk_path
            )
            for primer_name, matched_file in matched_files.items():
                if primer_name not in matched_files_combined:
                    matched_files_combined[primer_name] = []
                matched_files_combined[primer_name].append(matched_file)

        combined_matched_files = {}
        suffix = "_R1_R2" if len(file_paths) == 2 else "_R1" if "R1" in file_paths[0] else "_R2"
        for primer_name, files in matched_files_combined.items():
            combined_file = os.path.join(output_dir,
                f"combined_{sample_name}{suffix}_{primer_name}_MATCHED_{uuid.uuid4().hex[:8]}.fastq")
            try:
                with open(combined_file, 'w') as outfile:
                    for file in files:
                        with open(file, 'r') as infile:
                            shutil.copyfileobj(infile, outfile)
                        os.remove(file)
                file_size = os.path.getsize(combined_file)
                if file_size == 0:
                    logger.warning(f"Created 0-byte combined file {combined_file}. Removing.")
                    os.remove(combined_file)
                else:
                    combined_matched_files[primer_name] = combined_file
                    logger.info(f"Created combined MATCHED file: {combined_file} ({file_size} bytes)")
            except Exception as e:
                logger.error(f"Error combining files for {primer_name}: {e}")

        matched_reads_temp = defaultdict(list)
        for primer_name in combined_matched_files:
            try:
                with open(combined_matched_files[primer_name], 'r') as f:
                    for record in SeqIO.parse(f, 'fastq'):
                        matched_reads_temp[primer_name].append(record)
            except Exception as e:
                logger.error(f"Error reading combined file {primer_name}: {e}")
        save_matched_reads(combined_matched_files, matched_reads_temp, output_dir)

        (all_results, allele_counts, corrected_allele_counts,
         matched_reads, stats, locus_sizes_agg, locus_repr_size,
         mean_p_star_agg, dominant_det_cls) = process_bbduk_results(
            combined_matched_files, primers, nbmismatch_max, contig_setting,
            binning_setting, sequence_prefix, flanking_setting, output_dir,
            round_val, dico_bin=dico_bin
        )

        for combined_file in combined_matched_files.values():
            try:
                if os.path.exists(combined_file):
                    os.remove(combined_file)
            except Exception as e:
                logger.warning(f"Could not delete combined file {combined_file}: {e}")

        layout = "PE" if len(file_paths) == 2 else "SE"
        end_time = time.time()
        logger.info(f"Sample {sample_name} processing took {end_time - start_time:.2f} seconds")

        output_queue.put((
            {}, all_results, allele_counts, corrected_allele_counts,
            sample_name, total_nt, read_count, combined_matched_files,
            matched_reads, stats, layout,
            locus_sizes_agg, locus_repr_size,
            mean_p_star_agg, dominant_det_cls   # NEW: per-read aggregated detectability
        ))

    except Exception as e:
        logger.error(f"Error processing {sample_name}: {e}")
        traceback.print_exc()
        output_queue.put((None, None, None, None, sample_name, 0, 0, {}, {}, {}, "SE",
                          {}, {}, {}, {}))
    finally:
        output_queue.put(SENTINEL)


# =============================================================================
# Usage / argument parsing
# =============================================================================

def usage():
    print("Usage: python ReadTRAIL-Inference.py -i <input_dir> -o <output_dir> -p <primers_file> --genome-size <INT>")
    print()
    print("Required:")
    print("  -i, --input         Directory with FASTQ (.fastq / .fq / .fastq.gz) files.")
    print("  -o, --output        Directory to save output files.")
    print("  -p, --primer        Primer file (CSV/TSV: Locus_RepBP_ExpBP_RefU;FwdSeq;RevSeq).")
    print("  --genome-size INT   Genome size in bp (e.g. 4800000 for E. coli).")
    print()
    print("Optional:")
    print("  -m, --mismatch INT  Max mismatches allowed [default: 2].")
    print("  -t, --threads  INT  Threads per sample [default: 4].")
    print("  -b, --binning  FILE Path to binning file.")
    print("  --flanking-seq INT  Flanking sequence length [default: 0].")
    print("  --p-min         FLOAT Detectability threshold p_min [default: 0.05].")
    print("  --obs-threshold FLOAT Observable class threshold [default: 0.5].")
    print("  -h, --help          Show this help message.")
    print()
    print("Notes:")
    print("  Read count, mean read length, and total bases sequenced are derived")
    print("  automatically from the input FASTQ files.")
    print("  Genome size must be specified by the user as it cannot be inferred")
    print("  from short-read data without a reference genome.")
    print("  BBDuk is installed automatically by environment.yml.")


# =============================================================================
# Main
# =============================================================================

@timethis
def main():
    global dico_bin, flanking_len, P_MIN, OBS_THRESHOLD, MAX_WEIGHT

    dico_bin         = {}
    fasta_path       = None
    output_path      = None
    primer_file_path = None
    nb_mismatch      = 2
    threads          = 4
    binning_file     = None
    binning_setting  = False
    flanking_setting = False
    flanking_len     = 0
    round_val        = 0.25
    genome_size      = None   # REQUIRED — must be provided by user
    bbduk_path       = None

    logger.info("Parsing command-line arguments")
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "hi:o:p:m:t:b:",
            ["help", "input=", "output=", "primer=",
             "mismatch=", "threads=", "binning=",
             "flanking-seq=", "genome-size=", "p-min=", "obs-threshold=",
             "bbduk-path="])
    except getopt.GetoptError as err:
        logger.error(f"Argument error: {err}")
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            fasta_path = arg
        elif opt in ("-o", "--output"):
            output_path = arg
        elif opt in ("-p", "--primer"):
            primer_file_path = arg
        elif opt in ("-m", "--mismatch"):
            nb_mismatch = int(arg)
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-b", "--binning"):
            binning_setting = True
            binning_file = arg
        elif opt == "--flanking-seq":
            flanking_len     = int(arg)
            flanking_setting = flanking_len > 0
        elif opt == "--genome-size":
            genome_size = int(arg)
            logger.info(f"Genome size set to {genome_size} bp")
        elif opt == "--p-min":
            P_MIN = float(arg)
            MAX_WEIGHT = 1.0 / P_MIN   # keep cap consistent with p_min
            logger.info(f"p_min set to {P_MIN}, MAX_WEIGHT updated to {MAX_WEIGHT:.2f}")
        elif opt == "--obs-threshold":
            OBS_THRESHOLD = float(arg)
            logger.info(f"obs_threshold set to {OBS_THRESHOLD}")
        elif opt == "--bbduk-path":
            bbduk_path = arg

    # Validate required arguments
    if not all([fasta_path, output_path, primer_file_path]):
        logger.error("Missing required arguments: -i, -o, -p are all required.")
        usage()
        sys.exit(2)

    if genome_size is None:
        logger.error("--genome-size is required. Example: --genome-size 4800000")
        usage()
        sys.exit(2)

    if genome_size <= 0:
        logger.error(f"Invalid genome size: {genome_size}. Must be a positive integer.")
        sys.exit(2)

    if not os.path.isdir(fasta_path):
        logger.error(f"Input path '{fasta_path}' is not a directory.")
        sys.exit(2)
    fasta_path = os.path.join(fasta_path, '')

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_path = os.path.join(output_path, '')

    if not os.path.isfile(primer_file_path):
        logger.error(f"Primer file '{primer_file_path}' not found.")
        sys.exit(2)

    try:
        bbduk_path = resolve_bbduk_path(bbduk_path)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(2)
    logger.info(f"Using BBDuk executable: {bbduk_path}")

    # Load primers
    logger.info(f"Reading primers from {primer_file_path}")
    with open(primer_file_path, "r", encoding='utf-8') as f:
        primer_lines = f.read().replace("\r", "").split("\n")
        primer_lines = [line for line in primer_lines if line.strip()]
    Primers = clean_primers(primer_lines)
    if not Primers:
        logger.error("No valid primers found.")
        sys.exit(2)

    if binning_setting:
        dico_bin = build_dictionnary(binning_file)

    def parse_fastq_sample_read(base_name):
        """Return (sample_name, read_type) for supported PE FASTQ names.

        Supported examples:
          Sample_R1_001.fastq.gz  -> (Sample, R1)
          Sample_R2_001.fastq.gz  -> (Sample, R2)
          Sample_R1.fastq.gz      -> (Sample, R1)
          Sample_R2.fastq.gz      -> (Sample, R2)
          Sample_1.fastq.gz       -> (Sample, R1)
          Sample_2.fastq.gz       -> (Sample, R2)
          Sample_1_001.fastq.gz   -> (Sample, R1)
          Sample_2_001.fastq.gz   -> (Sample, R2)
        """
        patterns = (
            r'^(.+?)_(R[12])(?:_001)?$',
            r'^(.+?)_([12])(?:_001)?$',
        )
        for pattern in patterns:
            match = re.match(pattern, base_name)
            if not match:
                continue
            sample_name, read_token = match.groups()
            read_type = read_token if read_token.startswith('R') else f"R{read_token}"
            return sample_name, read_type
        return None, None

    def infer_read_suffix(file_paths):
        read_types = []
        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            base_name, ext = os.path.splitext(file_name)
            if ext.lower() == '.gz':
                base_name, ext2 = os.path.splitext(base_name)
                ext = ext2 + ext
            _, read_type = parse_fastq_sample_read(base_name)
            if read_type:
                read_types.append(read_type)
        read_types = sorted(set(read_types))
        if read_types == ['R1', 'R2']:
            return '_R1_R2'
        if read_types == ['R1']:
            return '_R1'
        if read_types == ['R2']:
            return '_R2'
        return '_R1_R2' if len(file_paths) == 2 else '_SE'

    # Discover samples
    logger.info(f"Scanning input directory {fasta_path}")
    sample_files = defaultdict(lambda: {'R1': None, 'R2': None})
    for file in sorted(os.listdir(fasta_path)):
        file_path = os.path.join(fasta_path, file)
        if os.path.isdir(file_path) or file.startswith('.'):
            continue
        base_name, ext = os.path.splitext(file)
        if ext.lower() == '.gz':
            base_name, ext2 = os.path.splitext(base_name)
            ext = ext2 + ext
        if ext.lower() in ['.fastq', '.fq', '.fastq.gz']:
            sample_name, read_type = parse_fastq_sample_read(base_name)
            if sample_name and read_type:
                sample_files[sample_name][read_type] = file_path
            else:
                logger.warning(
                    f"File {file} does not match expected naming pattern "
                    f"(SampleName_R1_001.fastq, SampleName_R1.fastq, SampleName_1.fastq)."
                )

    tasks = []
    for sample_name, files in sample_files.items():
        if files['R1'] and files['R2']:
            tasks.append((sample_name, [files['R1'], files['R2']]))
        elif files['R1']:
            tasks.append((sample_name, [files['R1']]))
        elif files['R2']:
            tasks.append((sample_name, [files['R2']]))

    if not tasks:
        logger.error("No valid FASTQ files found in input directory.")
        sys.exit(2)

    logger.info(f"Tasks to process: {len(tasks)} samples")

    max_processes, threads = adjust_resources(len(tasks), threads)
    output_queue = mp.Manager().Queue()
    all_outputs  = []

    logger.info(f"Starting processing with {max_processes} processes, {threads} threads each")

    def _collect_output_queue(non_blocking=False):
        collected = 0
        while True:
            try:
                if non_blocking:
                    queue_val = output_queue.get_nowait()
                else:
                    queue_val = output_queue.get(timeout=300)
                if queue_val != SENTINEL:
                    all_outputs.append(queue_val)
                    collected += 1
                    logger.info(f"Received data for sample: {queue_val[4]}")
                else:
                    logger.info(f"Received SENTINEL, processed {len(all_outputs)}/{expected_tasks} tasks")
            except queue.Empty:
                break
            except Exception as e:
                logger.error(f"Error retrieving from queue: {e}")
                traceback.print_exc()
                break
        return collected

    with ProcessPoolExecutor(max_workers=max_processes) as executor:
        futures = []
        for sample_name, file_paths in tasks:
            future = executor.submit(
                run, sample_name, file_paths, Primers, round_val, nb_mismatch, threads,
                output_queue, False, binning_setting, "seq", flanking_setting,
                output_path, dico_bin, bbduk_path
            )
            futures.append(future)

        expected_tasks = len(tasks)
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error(f"Future error: {e}")
                traceback.print_exc()
            _collect_output_queue(non_blocking=True)

        # Drain anything written by late-finishing workers after the last future
        # completed. This avoids losing samples when a later worker batch runs
        # longer than the old queue timeout.
        _collect_output_queue(non_blocking=True)

    if len(all_outputs) != expected_tasks:
        logger.warning(f"Expected {expected_tasks} results, got {len(all_outputs)}")

    if not all_outputs:
        logger.error("No results received.")
        sys.exit(1)

    # =========================================================================
    # Cohort-level baseline fallback
    # Collect all sample baselines; if a sample has too few loci, use cohort median.
    # =========================================================================
    cohort_baselines = []
    sample_baseline_info = {}  # sample_name → (baseline_val, loci_count, method)

    for entry in all_outputs:
        if entry[1] is None:
            continue
        (_, all_results, allele_counts, corrected_allele_counts,
         sample_name, total_nt, read_count, combined_matched_files,
         matched_reads, stats, layout,
         locus_sizes_agg, locus_repr_size,
         mean_p_star_agg, dominant_det_cls) = entry

        mean_read_len = stats['mean_len']
        baseline_val, n_loci, method = compute_baseline(
            allele_counts, locus_repr_size, mean_read_len)

        if method == "sample_median" and baseline_val is not None:
            cohort_baselines.append(baseline_val)
        sample_baseline_info[sample_name] = (baseline_val, n_loci, method)

    cohort_median = None
    if cohort_baselines:
        cohort_baselines.sort()
        mid = len(cohort_baselines) // 2
        cohort_median = (cohort_baselines[mid] if len(cohort_baselines) % 2 == 1
                         else (cohort_baselines[mid-1] + cohort_baselines[mid]) / 2)
        logger.info(f"Cohort baseline median: {cohort_median:.2f} from {len(cohort_baselines)} samples")

    # Resolve unavailable baselines
    for sample_name, (bval, n_loci, method) in sample_baseline_info.items():
        if bval is None and cohort_median is not None:
            sample_baseline_info[sample_name] = (cohort_median, n_loci, "cohort_median")
            logger.info(f"Sample {sample_name}: using cohort baseline {cohort_median:.2f}")
        elif bval is None:
            logger.warning(f"Sample {sample_name}: no baseline available (unavailable)")
            sample_baseline_info[sample_name] = (None, n_loci, "unavailable")

    # =========================================================================
    # ПРОБЛЕМА 2 fix: Cohort locus_factor
    # locus_factor_i = median(total_raw_count for locus i across cohort) / cohort_global_median
    # This corrects for locus-specific amplification efficiency, GC-bias, primer efficiency.
    # =========================================================================
    cohort_locus_raw = defaultdict(list)   # locus → [total_raw_count per sample]
    for entry in all_outputs:
        if entry[1] is None:
            continue
        allele_counts_entry = entry[2]
        if allele_counts_entry is None:
            continue
        for locus, allele_dict in allele_counts_entry.items():
            total = sum(allele_dict.values())
            if total > 0:
                cohort_locus_raw[locus].append(total)

    # Median per locus across cohort
    cohort_locus_medians = {}
    for locus, vals in cohort_locus_raw.items():
        vals_sorted = sorted(vals)
        mid = len(vals_sorted) // 2
        cohort_locus_medians[locus] = (vals_sorted[mid] if len(vals_sorted) % 2 == 1
                                       else (vals_sorted[mid-1] + vals_sorted[mid]) / 2)

    # Global cohort median (across all loci medians)
    all_locus_meds = sorted(cohort_locus_medians.values())
    global_cohort_median = None
    if all_locus_meds:
        mid = len(all_locus_meds) // 2
        global_cohort_median = (all_locus_meds[mid] if len(all_locus_meds) % 2 == 1
                                else (all_locus_meds[mid-1] + all_locus_meds[mid]) / 2)
        logger.info(f"Cohort global median for locus_factor: {global_cohort_median:.2f} "
                    f"across {len(all_locus_meds)} loci")

    # =========================================================================
    # Build output table  (two-pass: pre-compute → cohort registry → augment)
    # =========================================================================
    #
    # PASS 1: compute allele_v2_data per (sample_name, primer_full) and store.
    # After pass 1: build cohort_allele_info (all alleles observed anywhere).
    # PASS 2: augment each sample with zero-count alleles from cohort, write rows.
    #
    # ❗ Zero-count rows are written for REPORTING ONLY (reliable_alleles matrix
    #    and excluded_hidden_or_incomplete_loci in the downstream Statistik script).
    #    They MUST NOT enter any calculation:  N_filt, G_mini, C_mini, Q_i,
    #    B_sample, B_locus_cohort, Q_i_experimental are unchanged.
    # =========================================================================

    output_base_name    = os.path.basename(os.path.normpath(fasta_path)) or "mlva_input"
    analysis_excel_path = os.path.join(output_path, f"MLVA_analysis_{output_base_name}.xlsx")
    analysis_csv_path   = os.path.join(output_path, f"MLVA_analysis_{output_base_name}.csv")

    def _allele_sort_key(x):
        allele, raw_count, corr_count = x
        try:
            allele_num = float(allele)
        except Exception:
            allele_num = float("-inf")
        return (corr_count, raw_count, allele_num)

    # ── PASS 1: pre-compute per (sample, primer) ──────────────────────────────
    pre_computed = {}   # (sample_name, primer_full) → dict

    for entry in sorted(all_outputs, key=lambda x: x[4]):
        if entry[1] is None:
            continue
        (_, all_results, allele_counts, corrected_allele_counts,
         sample_name, total_nt, read_count, combined_matched_files,
         matched_reads, stats, layout,
         locus_sizes_agg, locus_repr_size,
         mean_p_star_agg, dominant_det_cls) = entry

        mean_read_len = stats['mean_len'] if stats['mean_len'] > 0 else 0
        total_bases   = stats['total_bases']

        baseline_val, baseline_n, baseline_method = sample_baseline_info.get(
            sample_name, (None, 0, "unavailable"))

        sample_file_paths = [t for t in tasks if t[0] == sample_name][0][1]
        suffix    = infer_read_suffix(sample_file_paths)
        accession = f"combined_{sample_name}{suffix}.fastq"

        for primer in Primers:
            primer_full  = primer[0]
            primer_short = primer_full.split('_')[0]
            key = (sample_name, primer_full)

            raw_counts  = allele_counts.get(primer_short, {})
            corr_counts = corrected_allele_counts.get(primer_short, {})

            if not raw_counts:
                pre_computed[key] = {
                    'empty': True,
                    'accession': accession,
                    'sample_name': sample_name,
                    'read_count': read_count,
                    'mean_read_len': mean_read_len,
                    'total_bases': total_bases,
                    'layout': layout,
                    'baseline_val': baseline_val,
                    'baseline_n': baseline_n,
                    'baseline_method': baseline_method,
                    'primer_short': primer_short,
                }
                continue

            total_raw_hits  = sum(raw_counts.values())
            total_corr_hits = sum(corr_counts.values())

            all_alleles = set(raw_counts.keys()) | set(corr_counts.keys())
            allele_list = sorted(
                [(a, raw_counts.get(a, 0), corr_counts.get(a, 0.0)) for a in all_alleles],
                key=_allele_sort_key, reverse=True
            )

            locus_factor = compute_locus_factor(
                primer_short, cohort_locus_medians, global_cohort_median)

            allele_v2_data = []
            for a, r, c in allele_list:
                frag_size = locus_sizes_agg.get(primer_short, {}).get(a)
                if frag_size is None:
                    for rec in all_results.get(primer_short, []):
                        if rec[4] == a:
                            frag_size = rec[3]
                            break
                if frag_size is None:
                    frag_size = 0

                mean_pstar = mean_p_star_agg.get(primer_short, {}).get(a)
                dom_class  = dominant_det_cls.get(primer_short, {}).get(a)

                if mean_pstar is None:
                    mean_pstar, dom_class, _ = compute_detectability(
                        frag_size, mean_read_len, P_MIN)
                    logger.debug(f"No per-read p_star for {primer_short}/{a}, "
                                 f"falling back to mean_read_len estimate")

                use_in_idealized = (dom_class not in
                                    ("Unobservable", "Effectively_Unobservable"))

                _is_val, _weight_applied, _weight_suppressed = compute_idealized_support(
                    r, mean_pstar, use_in_idealized)

                _is_low_p_marginal = (
                    dom_class == "Marginal"
                    and mean_pstar is not None
                    and mean_pstar < LOW_P_MARGINAL_THRESHOLD
                )
                _is_reliable = (
                    r >= 10
                    and mean_pstar is not None
                    and mean_pstar >= 0.10
                    and use_in_idealized
                )

                expected_reads = compute_expected_reads(
                    baseline_val, mean_pstar, use_in_idealized, locus_factor)

                is_dropout, obs_ratio = is_allele_dropout(
                    r, expected_reads, baseline_val)

                if expected_reads is not None and expected_reads > 0:
                    z_score = (r - expected_reads) / math.sqrt(expected_reads)
                else:
                    z_score = None

                is_eff_unobs = dom_class == "Effectively_Unobservable"

                allele_v2_data.append({
                    'allele':                       a,
                    'raw_count':                    r,
                    'corr_count':                   c,
                    'frag_size':                    frag_size,
                    'p_star':                       mean_pstar,
                    'det_class':                    dom_class,
                    'use_in_idealized':             use_in_idealized,
                    'expected_reads':               expected_reads,
                    'obs_ratio':                    obs_ratio,
                    'is_eff_unobs':                 is_eff_unobs,
                    'is_dropout':                   is_dropout,
                    'locus_factor':                 round(locus_factor, 4),
                    'z_score':                      z_score,
                    'idealized_support':            _is_val,
                    'idealized_weight_applied':     _weight_applied,
                    'idealized_weight_suppressed':  _weight_suppressed,
                    'idealized_support_reliable':   _is_reliable,
                    'is_low_p_marginal':            _is_low_p_marginal,
                    'is_reporting_only':            False,   # detected allele
                })

            lc = compute_locus_completeness(allele_v2_data, baseline_val,
                                            locus_factor=locus_factor)
            existing_status = lc.get("Base_Status", "UNKNOWN")
            dropout_status  = assign_dropout_status(
                lc, allele_v2_data, existing_status, baseline_val)

            confidence = 1.0
            if lc["Incomplete_Denominator"]:
                confidence *= 0.5
            if lc["Locus_Completeness_Obs"] is not None and lc["Locus_Completeness_Obs"] < 0.5:
                confidence *= 0.5
            if dropout_status in ("PARTIALLY_UNOBSERVABLE_MIXTURE",
                                  "UNRESOLVED_DUE_TO_DETECTABILITY"):
                confidence *= 0.7
            if lc.get("Has_Low_P_Marginal", False):
                confidence *= 0.85

            pre_computed[key] = {
                'empty':            False,
                'accession':        accession,
                'sample_name':      sample_name,
                'read_count':       read_count,
                'mean_read_len':    mean_read_len,
                'total_bases':      total_bases,
                'layout':           layout,
                'baseline_val':     baseline_val,
                'baseline_n':       baseline_n,
                'baseline_method':  baseline_method,
                'total_raw_hits':   total_raw_hits,
                'total_corr_hits':  total_corr_hits,
                'allele_v2_data':   allele_v2_data,
                'lc':               lc,
                'dropout_status':   dropout_status,
                'confidence':       confidence,
                'primer_short':     primer_short,
                'locus_factor':     locus_factor,
            }

    # ── Build cohort allele registry ──────────────────────────────────────────
    # For every locus, collect all allele values observed across ANY sample.
    # Zero-count rows are added when an allele was seen elsewhere but not here.
    # Reporting-only: these rows MUST NOT enter Q_i / E_ideal calculations.
    cohort_allele_info = defaultdict(dict)  # primer_short → {allele: {frag, p, cls}}
    for (sn, pf), d in pre_computed.items():
        if d.get('empty', True):
            continue
        ps = d['primer_short']
        for ad in d['allele_v2_data']:
            a = ad['allele']
            if a not in cohort_allele_info[ps]:
                cohort_allele_info[ps][a] = {
                    'frag_size': ad['frag_size'],
                    'p_star':    ad['p_star'],
                    'det_class': ad['det_class'],
                }

    # ── Compute final max_alleles (detected + cohort zero-count combined) ─────
    max_alleles_detected = 0
    for entry in all_outputs:
        if entry[2] is None:
            continue
        for counts in entry[2].values():
            max_alleles_detected = max(max_alleles_detected, len(counts))
        for counts in entry[3].values():
            max_alleles_detected = max(max_alleles_detected, len(counts))

    max_alleles_cohort = max(
        (len(v) for v in cohort_allele_info.values()), default=0)
    max_alleles = max(max_alleles_detected, max_alleles_cohort)

    # ── Build header (20 cols per allele — adds Allele_Detection_Status) ──────
    allele_columns = []
    for i in range(1, max_alleles + 1):
        allele_columns.extend([
            f"Allele {i}",
            f"Raw_Count {i}",
            f"Corrected_Count {i}",
            f"% of Hits {i}",
            f"% Corrected Hits {i}",
            f"Loci_Size {i}",
            f"Detectability_p {i}",
            f"Detectability_Class {i}",
            f"Use_in_Idealized {i}",
            f"Is_Effectively_Unobservable {i}",
            f"Expected_Reads {i}",
            f"Observation_Ratio {i}",
            f"Is_Potential_Dropout {i}",
            f"Locus_Factor {i}",
            f"Z_Score {i}",
            f"Idealized_Weight_Applied {i}",
            f"Idealized_Weight_Suppressed_LowRaw {i}",
            f"Idealized_Support_Reliable {i}",
            f"Is_Low_P_Marginal {i}",
            f"Allele_Detection_Status {i}",   # NEW — ТЗ §3.1
        ])

    analysis_header = [
        "Access_number",
        "Sample_ID",
        "Primer",
        "Read_Count",
        "Mean_Read_Length",
        "Total_Bases_Sequenced",
        "Layout",
        "Baseline_Sample",
        "Baseline_Loci_Count",
        "Baseline_Method",
        "Observed_Total_Raw",
        "Observable_Expected_Total",
        "Locus_Completeness_Obs",
        "Incomplete_Denominator",
        "Completeness_Is_Lower_Bound",
        "Has_Unobservable_Expected",
        "Dropout_Adjusted_Status",
        "Confidence_Score",
        "Has_Low_P_Marginal",
        "Low_P_Marginal_Count",
    ] + allele_columns

    # ── PASS 2: augment + write rows ──────────────────────────────────────────
    analysis_data = []
    logger.info("Processing results for output (pass 2: augmenting with zero-count alleles)")

    for entry in sorted(all_outputs, key=lambda x: x[4]):
        if entry[1] is None:
            continue
        sample_name = entry[4]

        for primer in Primers:
            primer_full  = primer[0]
            primer_short = primer_full.split('_')[0]
            key = (sample_name, primer_full)
            d   = pre_computed.get(key)

            # ── Empty locus ───────────────────────────────────────────────────
            if d is None or d.get('empty', True):
                acc_e  = d['accession']   if d else f"combined_{sample_name}.fastq"
                bl_e   = d['baseline_val'] if d else None
                bln_e  = d['baseline_n']   if d else 0
                blm_e  = d['baseline_method'] if d else "unavailable"
                rc_e   = d['read_count']   if d else 0
                mrl_e  = d['mean_read_len'] if d else 0
                tb_e   = d['total_bases']  if d else 0
                lay_e  = d['layout']       if d else ""
                row = [acc_e, sample_name, primer_full, rc_e,
                       round(mrl_e, 2), tb_e, lay_e,
                       bl_e if bl_e is not None else "",
                       bln_e, blm_e,
                       "", "", "", "", "", "", "", "", "", ""]
                row.extend([""] * (20 * max_alleles))
                analysis_data.append(row)
                continue

            # ── Copy pre-computed data ────────────────────────────────────────
            allele_v2_data  = list(d['allele_v2_data'])  # shallow copy
            detected_alleles = {ad['allele'] for ad in allele_v2_data}
            total_raw_hits   = d['total_raw_hits']
            total_corr_hits  = d['total_corr_hits']
            lc               = d['lc']
            dropout_status   = d['dropout_status']
            confidence       = d['confidence']
            accession        = d['accession']
            mean_read_len    = d['mean_read_len']
            total_bases      = d['total_bases']
            layout           = d['layout']
            baseline_val     = d['baseline_val']
            baseline_n       = d['baseline_n']
            baseline_method  = d['baseline_method']
            read_count       = d['read_count']
            locus_factor     = d['locus_factor']

            # ── Augment with zero-count alleles from cohort ───────────────────
            # These rows are REPORTING-ONLY; they do not change lc/confidence
            # because lc was already computed from detected alleles only (above).
            for a, info in sorted(
                cohort_allele_info[primer_short].items(),
                key=lambda x: (
                    float(x[0]) if str(x[0]).replace('.', '').replace('-', '').isdigit()
                    else float('inf'),
                    str(x[0])
                )
            ):
                if a in detected_alleles:
                    continue
                ps_zero = info['p_star']
                dc_zero = info['det_class']
                uid_zero = dc_zero not in ("Unobservable", "Effectively_Unobservable")
                er_zero  = compute_expected_reads(
                    baseline_val, ps_zero, uid_zero, locus_factor)

                allele_v2_data.append({
                    'allele':                       a,
                    'raw_count':                    0,
                    'corr_count':                   0.0,
                    'frag_size':                    info['frag_size'],
                    'p_star':                       ps_zero,
                    'det_class':                    dc_zero,
                    'use_in_idealized':             uid_zero,
                    'expected_reads':               er_zero,
                    'obs_ratio':                    None,
                    'is_eff_unobs':                 dc_zero == "Effectively_Unobservable",
                    'is_dropout':                   False,
                    'locus_factor':                 round(locus_factor, 4),
                    'z_score':                      None,
                    'idealized_support':            None,
                    'idealized_weight_applied':     False,
                    'idealized_weight_suppressed':  False,
                    'idealized_support_reliable':   False,
                    'is_low_p_marginal':            False,
                    'is_reporting_only':            True,   # ← NEVER used in formulas
                })

            # ── Compute Allele_Detection_Status for all alleles ───────────────
            for ad in allele_v2_data:
                ad['detection_status'] = get_allele_detection_status(
                    ad['raw_count'], ad['p_star'], ad['det_class'])

            # ── Build row ─────────────────────────────────────────────────────
            row = [
                accession, sample_name, primer_full,
                read_count, round(mean_read_len, 2), total_bases, layout,
                round(baseline_val, 4) if baseline_val is not None else "",
                baseline_n, baseline_method,
                lc["Observed_Total_Raw"],
                lc["Observable_Expected_Total"],
                lc["Locus_Completeness_Obs"] if lc["Locus_Completeness_Obs"] is not None else "",
                lc["Incomplete_Denominator"],
                lc["Completeness_Is_Lower_Bound"],
                lc["Has_Unobservable_Expected"],
                dropout_status,
                round(confidence, 4),
                lc.get("Has_Low_P_Marginal", False),
                lc.get("Low_P_Marginal_Count", 0),
            ]

            for ad in allele_v2_data:
                raw_perc  = (ad['raw_count'] / total_raw_hits * 100) if total_raw_hits else 0
                corr_perc = (ad['corr_count'] / total_corr_hits * 100) if total_corr_hits else 0
                allele_display = (int(float(ad['allele']))
                                  if str(ad['allele']).endswith('.0')
                                  else ad['allele'])
                row.extend([
                    allele_display,
                    ad['raw_count'],
                    round(ad['corr_count'], 4),
                    f"{raw_perc:.2f}%",
                    f"{corr_perc:.2f}%",
                    ad['frag_size'],
                    round(ad['p_star'], 4) if ad['p_star'] is not None else "",
                    ad['det_class'],
                    ad['use_in_idealized'],
                    ad['is_eff_unobs'],
                    round(ad['expected_reads'], 2) if ad['expected_reads'] is not None else "",
                    round(ad['obs_ratio'],    4) if ad['obs_ratio']      is not None else "",
                    ad['is_dropout'],
                    ad.get('locus_factor', ""),
                    round(ad['z_score'], 4) if ad.get('z_score') is not None else "",
                    ad.get('idealized_weight_applied',    ""),
                    ad.get('idealized_weight_suppressed', ""),
                    ad.get('idealized_support_reliable',  ""),
                    ad.get('is_low_p_marginal', False),
                    ad.get('detection_status', 'Detected'),  # NEW column
                ])

            # Pad to max_alleles  (20 cols per allele)
            row.extend([""] * (20 * (max_alleles - len(allele_v2_data))))
            analysis_data.append(row)

    # Validate row lengths
    logger.info("Validating analysis data")
    expected_columns = len(analysis_header)
    for i, row in enumerate(analysis_data):
        if len(row) != expected_columns:
            logger.error(f"Row {i} has {len(row)} columns, expected {expected_columns}")
            sys.exit(1)
        for j, val in enumerate(row):
            if val is None or (isinstance(val, float) and math.isnan(val)):
                analysis_data[i][j] = ""

    # Save CSV
    logger.info(f"Saving CSV to {analysis_csv_path}")
    try:
        csv_df = pd.DataFrame(analysis_data, columns=analysis_header)
        csv_df.to_csv(analysis_csv_path, index=False)
        logger.info(f"CSV saved: {analysis_csv_path}")
    except Exception as e:
        logger.error(f"Error saving CSV: {e}")
        traceback.print_exc()
        sys.exit(1)

    # Save Excel
    logger.info(f"Saving Excel to {analysis_excel_path}")
    try:
        analysis_df = pd.DataFrame(analysis_data, columns=analysis_header)
        with pd.ExcelWriter(analysis_excel_path, engine='openpyxl') as writer:
            analysis_df.to_excel(writer, index=False, sheet_name='Analysis')
        logger.info(f"Excel saved: {analysis_excel_path}")
    except Exception as e:
        logger.error(f"Error saving Excel: {e}")
        traceback.print_exc()
        sys.exit(1)

    logger.info("MLVA v2 analysis complete.")


if __name__ == "__main__":
    main()
