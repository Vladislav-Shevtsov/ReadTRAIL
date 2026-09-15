#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# ReadTRAIL-Inference — Idealized Model v2.5
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
# =============================================================================
# CHANGES IN v2.5  (audit round, see ReadTRAIL_audit_2026-09-10)
# =============================================================================
# P0 — correctness
#   [1] Stage II: 4-nt exact preflight REMOVED.  It discarded ~35.5 % of reads
#       whose primer mismatch fell inside the first 4 nt (measured).  Replaced
#       by an exhaustive, vectorised sliding Hamming search over the full
#       primer length (numpy sliding_window_view), which is also ~12x faster.
#   [2] Stage II is now IUPAC-aware: a degenerate primer base matches any base
#       in its set at zero mismatch cost.  An 'N' in the READ always costs one
#       mismatch (conservative), on both stages.
#   [3] split / wrap-around (circular) geometry DISABLED for FASTQ reads.
#       Only the two linear orientations remain: Config A (FWD + RC_REV) and
#       Config B (REV + RC_FWD).  The old branch also produced size != len(raw).
#   [4] Stage I (BBDuk): a SINGLE anchor = 3'-terminal 18 nt of the FORWARD
#       primer, with rcomp=t, in one pass.  REV/RC_REV anchors are no longer
#       used: every read able to pass Stage II must contain FWD or RC_FWD.
#       forbidn=f so the coarse filter can never be stricter than Stage II.
#   [5] Samples with zero BBDuk candidates now return a valid empty result
#       (Processing_Status = NO_CANDIDATE_READS) instead of raising
#       ValueError: max_workers must be greater than 0.
#   [6] FASTQ discovery rewritten: .fastq/.fq/.fastq.gz/.fq.gz, multi-part
#       (_001/_002/...) and multi-lane (L001/L002/...) files merged into one
#       sample, single-end accepted, and any FASTQ that cannot be assigned
#       unambiguously aborts the run instead of being skipped silently.
#   [7] Input sequencing metrics are computed from the INPUT FASTQ and are
#       reported separately from BBDuk-candidate metrics:
#         Input_Read_Count / Input_Total_Bases / Input_Mean_Read_Length /
#         Input_Min_Read_Length / Input_Max_Read_Length
#         BBDuk_Candidate_Reads / BBDuk_Candidate_Bases /
#         BBDuk_Candidate_Mean_Read_Length
#       The legacy columns Read_Count / Mean_Read_Length /
#       Total_Bases_Sequenced mixed the two populations (the latter two were
#       BBDuk-derived and summed over every locus).
#   [8] FASTQ parsing uses a record-aware iterator; a malformed record aborts
#       the run with an explicit error instead of silently shifting the frame.
#   [9] reads_passing_strict split into Reads_Passing_Primer_Pair /
#       Reads_With_Valid_Fragment_Size / Reads_With_Valid_Allele_Call so the
#       audit funnel balances.
#
# P1 — model and performance
#   [10] Dead chunking infrastructure removed: count_reads(),
#        split_file_tasks(), MAX_READS_PER_TASK.  CHUNK_SIZE = 50000 (the
#        working record-level chunking) is unchanged.
#   [11] matched_reads is no longer accumulated as SeqRecord objects, merged,
#        or shipped through the result queue.  It was never read.  FASTQ audit
#        artefacts are unaffected: they are written from disk files.
#   [12] The extra combined R1+R2 FASTQ copy is no longer created.  The
#        calculation still aggregates R1 and R2 into one sample accumulator
#        and the Access_number label is unchanged, because every accumulator
#        is an additive per-read statistic.
#   [13] Manager().Queue() + SENTINEL replaced by ProcessPoolExecutor return
#        values; samples can no longer be lost by a silent queue drain.
#   [14] Input read-length HISTOGRAM is collected and carried to the output so
#        that the population detectability
#             P(L) = E[max(0, R - L + 1)] / E[R]
#        can be computed over the input read-length distribution instead of
#        the mean p* over reads that already contain the fragment (which is a
#        length-biased sample: measured Spearman(L, effective read length)
#        = 0.432, p = 1.75e-17).  Reported as an additional column; the
#        published per-read p* is unchanged.
#   [15] locus_sizes is now resolved deterministically (median of observed
#        fragment sizes) instead of 'first record wins'.
#
# P2 — interface / documentation
#   [16] --max-size is a real CLI option and reaches evaluate_candidate_pair().
#   [17] --genome-size is OPTIONAL.  When supplied it yields diagnostic
#        columns Genome_Coverage_Estimate and (per allele) Locus_Efficiency.
#   [18] Audit funnel written per sample to audit_funnel.tsv, and the invoking
#        command line + all parameters written to run_metadata.txt.
#   [19] R1 and R2 are analysed as independent read-level observations without
#        paired-read merging (unchanged behaviour, now documented).
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
CHUNK_SIZE        = 50000
MEMORY_THRESHOLD  = 0.85
CPU_THRESHOLD     = 0.90
MAX_PROCESSES     = 24
THREADS_PER_PROCESS = 2
SAVE_FILTERED_READS = False  # audit artifacts; never needed for calculations
ANCHOR_LEN        = 18     # Stage-I anchor: 3'-terminal k-mer of the FWD primer
DEFAULT_MAX_SIZE  = 1000   # L_max; overridable with --max-size
BBDUK_XMX         = "4g"   # JVM heap handed to BBDuk; overridable with --bbduk-xmx.
# v2.5: this was hardcoded as Xmx=4g. With several samples in flight that is
# 4 GB per BBDuk process, which fails outright on any machine with less RAM
# than max_processes x 4 GB.
READ_LEN_BIN      = 1      # bin width (bp) of the input read-length histogram.
# 1 = exact read lengths.  A wider bin makes P(L) inherit the bin's rounding
# error (with a 10-bp bin a 301-bp library is represented as 304.5 bp), and the
# histogram is tiny anyway: at most one entry per distinct read length.
# NOTE: MAX_READS_PER_TASK / count_reads() / split_file_tasks() were removed in
# v2.5.  They defined logical read ranges that BBDuk cannot consume, so every
# logical chunk re-ran BBDuk on the whole FASTQ and duplicated matched reads.

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
# v2.5 — IUPAC-aware primer matching
# =============================================================================
# A primer base is a SET of acceptable read bases.  A read base matches at zero
# cost when it belongs to that set; otherwise it costs one mismatch.  An 'N' in
# the READ belongs to no set and therefore always costs one mismatch — the same
# rule is applied on Stage I (forbidn=f keeps BBDuk from being stricter).
IUPAC_SETS = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'T',
    'R':'AG', 'Y':'CT', 'S':'CG', 'W':'AT', 'K':'GT', 'M':'AC',
    'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG', 'N':'ACGT',
    '.':'ACGT',
}

def primer_has_degenerate(primer: str) -> bool:
    """True when the primer contains any non-ACGT IUPAC code."""
    return any(c not in 'ACGT' for c in str(primer).upper())


def build_iupac_table(query: str) -> np.ndarray:
    """
    Boolean lookup of shape (len(query), 256).

    table[j, b] is True when read byte b satisfies query position j.
    Unknown query characters accept only themselves, so behaviour for plain
    A/C/G/T primers is byte-identical to the previous exact comparison.
    """
    q = str(query).upper()
    table = np.zeros((len(q), 256), dtype=bool)
    for j, ch in enumerate(q):
        for base in IUPAC_SETS.get(ch, ch):
            table[j, ord(base)] = True
    return table


def hamming_mismatches(seq_a: str, seq_b: str) -> int:
    """Plain positional mismatch count (kept for backwards compatibility)."""
    if len(seq_a) != len(seq_b):
        return max(len(seq_a), len(seq_b))
    a = np.frombuffer(seq_a.encode(), dtype=np.uint8)
    b = np.frombuffer(seq_b.encode(), dtype=np.uint8)
    return int(np.count_nonzero(a != b))


# =============================================================================
# v2.5 — record-aware FASTQ reading and input read-length statistics
# =============================================================================

class FastqFormatError(RuntimeError):
    """Raised when an input FASTQ record is malformed."""


def open_fastq(file_path):
    """Open a plain or gzipped FASTQ for text reading."""
    return (gzip.open(file_path, 'rt') if str(file_path).endswith('.gz')
            else open(file_path, 'r'))


def iter_fastq(file_path):
    """
    Yield (title, sequence, quality) for every record.

    Unlike the previous `i % 4` counting this validates record structure: a
    truncated or corrupted file raises FastqFormatError instead of silently
    shifting the reading frame and mis-counting every subsequent read.
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    with open_fastq(file_path) as handle:
        try:
            for title, seq, qual in FastqGeneralIterator(handle):
                yield title, seq, qual
        except ValueError as exc:
            raise FastqFormatError(f"Malformed FASTQ record in {file_path}: {exc}") from exc


def scan_input_fastq(file_path, bin_width=READ_LEN_BIN):
    """
    Single full pass over one INPUT FASTQ.

    Returns a dict with the true sequencing metrics of that file plus a binned
    read-length histogram.  These are properties of the input library and must
    not be confused with the BBDuk-candidate subset, which is enriched for the
    loci being queried and is summed over every locus.
    """
    total_nt = 0
    read_count = 0
    min_len = None
    max_len = 0
    hist = defaultdict(int)
    for _title, seq, _qual in iter_fastq(file_path):
        n = len(seq)
        total_nt += n
        read_count += 1
        if min_len is None or n < min_len:
            min_len = n
        if n > max_len:
            max_len = n
        hist[(n // bin_width) * bin_width] += 1
    logger.info(f"Input {os.path.basename(file_path)}: {read_count} reads, {total_nt} bases")
    return {
        'total_nt': total_nt,
        'read_count': read_count,
        'min_len': min_len or 0,
        'max_len': max_len,
        'hist': dict(hist),
    }


def merge_length_hist(target: dict, source: dict) -> dict:
    """Accumulate one read-length histogram into another."""
    for k, v in source.items():
        target[k] = target.get(k, 0) + v
    return target


def population_detectability(fragment_len, length_hist, bin_width=READ_LEN_BIN):
    """
    Detectability of a fragment of length L over the INPUT read-length
    distribution f(R):

        P(L) = E_f[max(0, R - L + 1)] / E_f[R]

    This is the quantity that the number of primer-delimited observations is
    actually proportional to.  It is NOT the same as either
      * p*(L, mean R)  — Jensen's inequality: max(0, .) is convex, so using the
        mean read length UNDER-estimates P(L); nor
      * mean p* over the reads that were observed — those reads were sampled
        with probability proportional to (R - L + 1), so their mean read length
        is inflated and P(L) is OVER-estimated, increasingly so for long L.
    Returns None when the histogram is unavailable.
    """
    if not length_hist:
        return None
    num = 0.0
    den = 0.0
    for start, count in length_hist.items():
        if count <= 0:
            continue
        # bin representative: midpoint of [start, start + bin_width)
        r = start + (bin_width - 1) / 2.0
        den += count * r
        span = r - fragment_len + 1.0
        if span > 0:
            num += count * span
    if den <= 0:
        return None
    return max(0.0, num / den)

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

# v2.5: count_reads() and split_file_tasks() were removed — see the change log
# at the top of this file.  BBDuk consumes a file path, not a record range, so
# the logical "chunks" they produced made BBDuk re-scan the whole FASTQ once per
# chunk and duplicate every matched read.  Record-level chunking is done in
# process_bbduk_results() with CHUNK_SIZE and is unaffected.

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
              bbduk_path, max_distance=DEFAULT_MAX_SIZE):
    """
    Stage I — coarse BBDuk filter (v2.5).

    A single anchor is used: the 3'-terminal ANCHOR_LEN nt of the FORWARD
    primer, searched together with its reverse complement via rcomp=t in ONE
    pass over the FASTQ.

    Completeness.  Any read that can pass Stage II carries a full-length hit of
    FWD (Config A) or of RC_FWD (Config B) with <= nbmismatch mismatches.  By
    the pigeonhole principle such a hit implies <= nbmismatch mismatches over
    the 3'-terminal ANCHOR_LEN window, which is exactly what BBDuk searches
    with hdist=nbmismatch.  The REV / RC_REV anchors used before v2.5 were
    therefore redundant and only enlarged the candidate set.

    forbidn=f: an N in a read costs one mismatch here exactly as it does in
    Stage II.  With forbidn=t the coarse filter was STRICTER than the strict
    stage and dropped reads that Stage II would have accepted, silently and
    without any counter.
    """
    logger.info(f"Running BBDuk for {file_path}")
    name = os.path.basename(file_path)
    name, ext = os.path.splitext(name)
    if ext == '.gz':
        name, _ = os.path.splitext(name)
    matched_files = {}
    for primer in primers:
        primer_name = primer[0].replace('_', '-')
        anchor = primer[1][-ANCHOR_LEN:]
        out_matched = os.path.join(output_dir, f"{name}_{primer_name}_MATCHED_{uuid.uuid4().hex[:8]}.fastq")
        cmd = [
            bbduk_path,
            f"in={file_path}",
            f"outm={out_matched}",
            f"literal={anchor}",
            "mm=f",
            "rcomp=t",
            f"k={len(anchor)}",
            "minlen=1",
            f"hdist={nbmismatch}",
            "forbidn=f",
            f"threads={threads}",
            "overwrite=t",
            f"Xmx={BBDUK_XMX}",
        ]
        if primer_has_degenerate(anchor):
            # expand IUPAC codes of the anchor into concrete k-mers
            cmd.append("copyundefined=t")
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

# v2.5: count_nucleotides() removed — replaced by scan_input_fastq(), which
# parses complete FASTQ records (so a truncated file raises instead of
# silently shifting the frame) and also returns min/max length and a binned
# read-length histogram.

# =============================================================================
# Strict-stage full-length primer validation (v2.4)
# Replaces seed-automaton as the basis for accept/reject decisions.
# BBDuk remains the coarse filter; these functions handle exact validation.
# =============================================================================

def build_primer_forms(fwd_primer: str, rev_primer: str):
    """
    Pre-compute the four primer forms and their IUPAC lookup tables once per
    locus, outside the per-read loop.
    """
    rc_fwd = inverComp(fwd_primer)
    rc_rev = inverComp(rev_primer)
    forms = []
    for form_name, query in (("FWD", fwd_primer), ("REV", rev_primer),
                             ("RC_FWD", rc_fwd), ("RC_REV", rc_rev)):
        q = str(query).upper()
        forms.append((form_name, q, len(q), build_iupac_table(q)))
    return forms, rc_fwd, rc_rev


def collect_full_length_primer_hits(seq: str, fwd_primer: str, rev_primer: str,
                                    nbmismatch: int,
                                    rc_fwd: str = None, rc_rev: str = None,
                                    forms=None) -> dict:
    """
    Exhaustive full-length search of all 4 primer forms in *seq* with
    <= nbmismatch mismatches.  Returns a dict of hit lists keyed by form name.

    v2.5 changes:
      * The 4-nt exact prefix preflight is GONE.  It required the first
        min(4, |q| - m) bases of the primer to occur verbatim somewhere in the
        read, so a permitted mismatch inside that prefix could reject a valid
        read before it was ever tested.  Measured on simulated 301-bp reads
        carrying one mismatch in the first 4 nt: 35.5 % of otherwise valid
        primer hits were discarded.  It also invalidated the completeness
        claim of Proposition 3.1.
      * The per-position Python loop is replaced by a vectorised sliding window
        (numpy.lib.stride_tricks.sliding_window_view), ~12x faster, so removing
        the preflight costs nothing in runtime.
      * Matching is IUPAC-aware: a degenerate primer base matches any base of
        its set at zero cost; an 'N' in the read matches nothing and costs one
        mismatch.

    `forms` may be supplied by build_primer_forms() to hoist table construction
    out of the read loop; rc_fwd / rc_rev are accepted for compatibility.
    """
    if forms is None:
        forms, rc_fwd, rc_rev = build_primer_forms(fwd_primer, rev_primer)

    seq_len = len(seq)
    hits = {"FWD": [], "REV": [], "RC_FWD": [], "RC_REV": []}
    if seq_len == 0:
        return hits
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)

    for form_name, query, qlen, table in forms:
        if qlen == 0 or seq_len < qlen:
            continue
        windows = np.lib.stride_tricks.sliding_window_view(arr, qlen)
        # matches[i, j] — read base at window i, offset j satisfies query[j]
        matches = table[np.arange(qlen)[None, :], windows]
        mismatches = qlen - matches.sum(axis=1)
        for start in np.flatnonzero(mismatches <= nbmismatch):
            hits[form_name].append({"start": int(start), "length": qlen,
                                    "mismatches": int(mismatches[start])})
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


def evaluate_candidate_pair(candidate, seq, contig, max_distance,
                            allow_split=False):
    """
    Compute fragment geometry for one candidate pair.

    Fragment size = end_of_right_primer - start_of_left_primer
    (full amplicon length, inclusive of both primers).

    v2.5: the split / wrap-around branch is DISABLED for FASTQ input
    (allow_split defaults to False).  An Illumina read is a linear molecule, so
    a primer pair is only meaningful when the two primers occur in physical
    order within the read; the two biologically valid orientations are already
    covered by Config A (FWD + RC_REV) and Config B (REV + RC_FWD).  The old
    branch additionally computed size = rl + (seq_len - ls) while extracting
    raw = seq[ls:] + seq[:rs + rl], so the reported fragment length disagreed
    with the length of the extracted insert whenever rs > 0.  Circular contigs,
    if ever needed, belong in a separate mode.
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
    elif ls > rs and allow_split and not contig:
        # Retained only for an explicit circular mode; never used for FASTQ.
        size     = rl + (seq_len - ls)
        splitted = True
    # ls == rs -> same position -> invalid
    # ls  > rs without allow_split -> discarded (was: wrap-around)

    if size < 1 or size > max_distance:
        return {"is_valid": False, "config": candidate["config"],
                "size": size, "splitted": splitted}

    # Build insert sequence (full amplicon, always in FWD orientation)
    if splitted:
        raw = seq[ls:] + seq[:(rs + rl)]
        insert = raw if candidate["config"] == "A" else inverComp(raw)
    else:
        raw = seq[ls:(rs + rl)]
        # For Config B the amplicon is on the antisense strand
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

def save_matched_audit_copies(matched_files, output_path):
    """
    Copy the BBDuk candidate FASTQ files into matched_fastq/ and a FASTA
    rendition into matched_fasta/, for user-side validation.

    v2.5: renamed from save_matched_reads().  A boolean parameter of the same
    name shadowed this function inside run(), so passing --save-matched-reads
    raised "TypeError: 'bool' object is not callable"; the broad except then
    discarded the ENTIRE sample.  The unused `matched_reads` argument is also
    gone — the copies are made from the files on disk and never needed the
    in-memory records.

    matched_files maps a locus name to a path or to a list of paths.
    """
    fastq_output_folder = os.path.join(output_path, "matched_fastq")
    fasta_output_folder = os.path.join(output_path, "matched_fasta")
    os.makedirs(fastq_output_folder, exist_ok=True)
    os.makedirs(fasta_output_folder, exist_ok=True)
    for primer_name in matched_files:
        entry = matched_files[primer_name]
        paths = entry if isinstance(entry, (list, tuple)) else [entry]
        for fastq_output_path in paths:
            _copy_one_audit_artifact(primer_name, fastq_output_path,
                                     fastq_output_folder, fasta_output_folder)


def _copy_one_audit_artifact(primer_name, fastq_output_path,
                             fastq_output_folder, fasta_output_folder):
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
    size_counts               = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    p_star_sums               = defaultdict(lambda: defaultdict(float))
    det_class_votes           = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    # v2.5 audit funnel.  The previous single counter `reads_passing_strict`
    # was incremented before sizeU had been computed and validated, so reads
    # that produced no allele were still counted as having passed the strict
    # stage and the funnel did not balance.  Three explicit stages instead.
    #
    # The stages are NESTED, not parallel, so the invariant is a chain rather
    # than one flat sum (a read that finds a primer pair and then fails the
    # size check belongs to both Reads_Passing_Primer_Pair and
    # Rejected_Invalid_Size, so adding all three to the candidate total would
    # double-count it):
    #     BBDuk_Candidate_Reads          = Reads_Passing_Primer_Pair
    #                                      + Rejected_No_Pair
    #     Reads_Passing_Primer_Pair      = Reads_With_Valid_Fragment_Size
    #                                      + Rejected_Invalid_Size
    #     Reads_With_Valid_Fragment_Size = Reads_With_Valid_Allele_Call
    #                                      + Rejected_Invalid_Allele
    funnel = {
        "BBDuk_Candidate_Reads":          len(records),
        "Reads_Passing_Primer_Pair":      0,
        "Reads_With_Valid_Fragment_Size": 0,
        "Reads_With_Valid_Allele_Call":   0,
        "Rejected_No_Pair":               0,
        "Rejected_Invalid_Size":          0,
        "Rejected_Invalid_Allele":        0,
        "Reads_With_Multiple_Valid_Pairs": 0,
    }

    filtered_fastq_folder = os.path.join(output_path, "filtered_both_primers_fastq")
    filtered_fasta_folder = os.path.join(output_path, "filtered_both_primers_fasta")
    if SAVE_FILTERED_READS:
        os.makedirs(filtered_fastq_folder, exist_ok=True)
        os.makedirs(filtered_fasta_folder, exist_ok=True)

    primer = next((p for p in primers if p[0].replace('_', '-') == primer_name), None)
    if not primer:
        logger.error(f"Primer {primer_name} not found")
        return (all_sequence_results, allele_counts, corrected_allele_counts,
                size_counts, p_star_sums, det_class_votes, funnel)

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
    # IUPAC lookup tables and reverse complements built once per chunk
    primer_forms, rc_fwd_cached, rc_rev_cached = build_primer_forms(fwd_primer, rev_primer)

    for s, record in enumerate(records):
        seq      = str(record.seq).upper()
        seq_id   = record.id
        read_len = len(seq)
        current_sequence_name = f"{sequence_prefix}{chunk_index}_{chunk_uuid}_{s+1}"

        # v2.5: the BBDuk candidate is NOT retained as a SeqRecord.  The old
        # matched_reads accumulator held every candidate of every locus in RAM,
        # was merged across chunks, was pickled through the result queue and
        # was then never read in main().  Read-level FASTQ artefacts are
        # unaffected: BBDuk candidates are copied from the files on disk by
        # save_matched_audit_copies(), and reads passing full primer-pair
        # validation are written from `filtered_records` below.

        # -- Stage 2: full-length strict validation ---------------------------
        try:
            full_hits = collect_full_length_primer_hits(
                seq, fwd_primer, rev_primer, nbmismatch,
                rc_fwd=rc_fwd_cached, rc_rev=rc_rev_cached,
                forms=primer_forms)
        except Exception as e:
            logger.warning(f"Error in full-length primer search for "
                           f"{primer_file_name} read {seq_id}: {e}")
            funnel["Rejected_No_Pair"] += 1
            continue

        candidates = build_candidate_pairs(full_hits, fwd_primer, rev_primer)
        if not candidates:
            funnel["Rejected_No_Pair"] += 1
            continue
        funnel["Reads_Passing_Primer_Pair"] += 1

        evaluated = []
        for cand in candidates:
            ev = evaluate_candidate_pair(cand, seq, contig, max_distance)
            if ev["is_valid"]:
                evaluated.append(ev)

        if not evaluated:
            funnel["Rejected_Invalid_Size"] += 1
            continue
        funnel["Reads_With_Valid_Fragment_Size"] += 1

        if len(evaluated) > 1:
            funnel["Reads_With_Multiple_Valid_Pairs"] += 1

        best = select_best_candidate(evaluated)
        if SAVE_FILTERED_READS:
            filtered_records.append(record)
        # -- End strict validation --------------------------------------------

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
            funnel["Rejected_Invalid_Allele"] += 1
            continue
        funnel["Reads_With_Valid_Allele_Call"] += 1
        # -- End sizeU --------------------------------------------------------

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

        # v2.5: record the full distribution of observed fragment sizes for
        # this allele instead of keeping whichever size happened to be seen
        # first.  With binning correction several sizes can map to the same
        # allele, so "first record wins" made Loci_Size depend on chunk order
        # and undermined the bit-reproducibility claim.
        size_counts[primer_short_name][sizeU][size] += 1

        p_star_sums[primer_short_name][sizeU]       += p_star
        det_class_votes[primer_short_name][sizeU][det_class] += 1

    # Log funnel stats for this chunk
    logger.debug(
        f"Chunk {chunk_index} [{primer_name}] funnel: " +
        ", ".join(f"{k}={v}" for k, v in funnel.items())
    )

    # Write filtered reads
    if SAVE_FILTERED_READS and filtered_records:
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
    # BUG FIX v2.4.1 (kept): pre-computing mean_p_star / dominant_det_class here
    # produced a string where _merge_chunk expected a dict, which the broad
    # except swallowed, making every chunk look empty.  Only raw structures are
    # returned:
    #   p_star_sums     -> {locus: {allele: sum of p*}}
    #   det_class_votes -> {locus: {allele: {class: count}}}
    #   size_counts     -> {locus: {allele: {fragment_size: count}}}
    return (all_sequence_results, allele_counts, corrected_allele_counts,
            size_counts, p_star_sums, det_class_votes, funnel)


def process_bbduk_results(matched_files, primers, nbmismatch, contig, binning,
                          sequence_prefix, flanking, output_path, round_val=0.25,
                          max_distance=DEFAULT_MAX_SIZE, dico_bin=None):
    """
    Process every BBDuk-matched file of one sample.

    v2.5 notes
    ----------
    * The returned `stats` describe the BBDuk CANDIDATE population only, and
      they are summed over loci: a read matching k loci is counted k times.
      They are therefore reported as BBDuk_Candidate_* and must not be used as
      sequencing metrics of the library — those come from scan_input_fastq().
    * Fragment size per allele is resolved deterministically (median over all
      observed sizes) rather than "first record wins".
    * A sample with no BBDuk hits returns a valid empty result instead of
      constructing ThreadPoolExecutor(max_workers=0), which raised ValueError
      and lost the whole sample through the broad except in run().
    """
    logger.info("Starting process_bbduk_results")
    all_sequence_results    = defaultdict(list)
    allele_counts           = defaultdict(lambda: defaultdict(int))
    corrected_allele_counts = defaultdict(lambda: defaultdict(float))
    size_counts_agg         = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    p_star_sums_agg         = defaultdict(lambda: defaultdict(float))
    p_star_counts_agg       = defaultdict(lambda: defaultdict(int))
    det_class_votes_agg     = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    funnel_agg              = defaultdict(int)

    total_reads = 0
    total_bases = 0
    min_len     = float('inf')
    max_len     = 0

    def _empty_result():
        return (defaultdict(list),
                defaultdict(lambda: defaultdict(int)),
                defaultdict(lambda: defaultdict(float)),
                defaultdict(lambda: defaultdict(lambda: defaultdict(int))),
                defaultdict(lambda: defaultdict(float)),
                defaultdict(lambda: defaultdict(lambda: defaultdict(int))),
                defaultdict(int),
                (0, 0, 0, 0))

    def _merge_chunk(cr, ca, cc, csz, c_ps, c_dc, c_fn,
                     local_ar, local_ac, local_cc2, local_sz,
                     local_ps, local_pc, local_dc, local_fn):
        """Merge one chunk's raw accumulators into the per-file accumulators."""
        for k, v in cr.items():
            local_ar[k].extend(v)
        for k, v in ca.items():
            for allele, cnt in v.items():
                local_ac[k][allele] += cnt
        for k, v in cc.items():
            for allele, wt in v.items():
                local_cc2[k][allele] += wt
        for locus, alleles in csz.items():
            for allele, sizes in alleles.items():
                for sz, cnt in sizes.items():
                    local_sz[locus][allele][sz] += cnt
        for locus, ad in c_ps.items():
            for allele, val in ad.items():
                local_ps[locus][allele] += val
        for locus, ad in c_dc.items():
            for allele, cls_dict in ad.items():
                for cls, cnt in cls_dict.items():
                    local_dc[locus][allele][cls] += cnt
        # p_star_counts are the number of contributing reads, which is exactly
        # the number of detectability-class votes for that allele.
        for locus, ad in c_dc.items():
            for allele, cls_dict in ad.items():
                local_pc[locus][allele] += sum(cls_dict.values())
        for key, val in c_fn.items():
            local_fn[key] += val

    def process_file(primer_name, matched_entry):
        """
        v2.5: `matched_entry` is a LIST of BBDuk output files for this locus —
        normally one per input FASTQ (R1 and R2).  They are streamed one after
        another into the SAME per-locus accumulators, which is arithmetically
        identical to concatenating them into one combined FASTQ first, because
        every accumulator is an additive per-read statistic with no cross-read
        term.  R1 and R2 therefore remain independent read-level observations
        aggregated into a single sample, exactly as before, but without writing
        and re-reading a full extra copy of every matched read.
        """
        matched_paths = ([matched_entry] if isinstance(matched_entry, str)
                         else list(matched_entry))
        matched_paths = [p for p in matched_paths if p]
        existing = [p for p in matched_paths if os.path.exists(p)]
        for missing in [p for p in matched_paths if p not in existing]:
            logger.warning(f"Matched file for {primer_name} not found at {missing}.")
        if not existing:
            return _empty_result()

        chunk_records = []
        chunk_index   = 0
        total_reads_p = 0
        total_bases_p = 0
        min_len_p     = float('inf')
        max_len_p     = 0

        local_ar  = defaultdict(list)
        local_ac  = defaultdict(lambda: defaultdict(int))
        local_cc2 = defaultdict(lambda: defaultdict(float))
        local_sz  = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        local_ps  = defaultdict(lambda: defaultdict(float))
        local_pc  = defaultdict(lambda: defaultdict(int))
        local_dc  = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        local_fn  = defaultdict(int)

        def _run_chunk(recs, idx):
            chunk_uuid = uuid.uuid4().hex[:8]
            cr, ca, cc, csz, c_ps, c_dc, c_fn = process_chunk(
                recs, primers, nbmismatch, contig, binning,
                sequence_prefix, flanking, output_path, round_val,
                max_distance, dico_bin, idx, primer_name, chunk_uuid)
            _merge_chunk(cr, ca, cc, csz, c_ps, c_dc, c_fn,
                         local_ar, local_ac, local_cc2, local_sz,
                         local_ps, local_pc, local_dc, local_fn)

        try:
            for matched_file in existing:
                with (gzip.open(matched_file, 'rt') if matched_file.endswith('.gz')
                      else open(matched_file, 'r')) as f:
                    for record in SeqIO.parse(f, 'fastq'):
                        rl = len(record.seq)
                        total_bases_p += rl
                        total_reads_p += 1
                        if rl < min_len_p: min_len_p = rl
                        if rl > max_len_p: max_len_p = rl
                        chunk_records.append(record)

                        if len(chunk_records) >= CHUNK_SIZE:
                            _run_chunk(chunk_records, chunk_index)
                            chunk_records = []
                            chunk_index  += 1
                            if get_system_resources()['memory_percent'] > MEMORY_THRESHOLD * 100:
                                logger.warning("Memory threshold exceeded. Pausing.")
                                time.sleep(5)

            if chunk_records:
                _run_chunk(chunk_records, chunk_index)

            logger.info(f"Processed {total_reads_p} candidate reads for {primer_name}")
        except FastqFormatError:
            raise
        except Exception as e:
            logger.error(f"Error processing candidates for {primer_name}: {e}")
            traceback.print_exc()
            return _empty_result()

        stats = (total_reads_p, total_bases_p, min_len_p, max_len_p)
        return (local_ar, local_ac, local_cc2, local_sz,
                local_ps, local_pc, local_dc, local_fn, stats)

    # v2.5: a sample with zero BBDuk hits is a legitimate outcome, not an error.
    if not matched_files:
        logger.warning("No BBDuk candidate files for this sample — returning an "
                       "empty but valid result (NO_CANDIDATE_READS).")
        return ({}, {}, {},
                {'total_reads': 0, 'total_bases': 0, 'mean_len': 0,
                 'min_len': 0, 'max_len': 0},
                {}, {}, {}, {}, dict(funnel_agg))

    max_workers = max(1, THREADS_PER_PROCESS * len(matched_files))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, pn, mf) for pn, mf in matched_files.items()]
        for future in futures:
            (local_ar, local_ac, local_cc, local_sz,
             local_ps, local_pc, local_dc, local_fn, stats) = future.result()
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
            for locus, alleles in local_sz.items():
                for allele, sizes in alleles.items():
                    for sz, cnt in sizes.items():
                        size_counts_agg[locus][allele][sz] += cnt
            for locus, ad in local_ps.items():
                for allele, val in ad.items():
                    p_star_sums_agg[locus][allele] += val
            for locus, ad in local_pc.items():
                for allele, cnt in ad.items():
                    p_star_counts_agg[locus][allele] += cnt
            for locus, ad in local_dc.items():
                for allele, cls_dict in ad.items():
                    for cls, cnt in cls_dict.items():
                        det_class_votes_agg[locus][allele][cls] += cnt
            for key, val in local_fn.items():
                funnel_agg[key] += val

    mean_len = total_bases / total_reads if total_reads > 0 else 0
    if min_len == float('inf'): min_len = 0

    stats_dict = {
        # BBDuk candidate population, summed over loci — NOT library metrics
        'total_reads': total_reads,
        'total_bases': total_bases,
        'mean_len':    mean_len,
        'min_len':     min_len,
        'max_len':     max_len,
    }

    # Deterministic fragment size per allele: median of all observed sizes.
    def _weighted_median(size_hist):
        items = sorted(size_hist.items())
        total = sum(c for _, c in items)
        if total == 0:
            return None
        half, run = total / 2.0, 0
        for sz, cnt in items:
            run += cnt
            if run >= half:
                return sz
        return items[-1][0]

    locus_sizes_agg = defaultdict(dict)
    for locus, alleles in size_counts_agg.items():
        for allele, size_hist in alleles.items():
            med = _weighted_median(size_hist)
            if med is not None:
                locus_sizes_agg[locus][allele] = med

    # Representative fragment size per locus (median across its alleles)
    locus_repr_size = {}
    for locus, allele_sizes in locus_sizes_agg.items():
        sizes_list = sorted(allele_sizes.values())
        if not sizes_list:
            continue
        mid = len(sizes_list) // 2
        locus_repr_size[locus] = (sizes_list[mid] if len(sizes_list) % 2 == 1
                                  else (sizes_list[mid - 1] + sizes_list[mid]) / 2)

    # Per-allele mean p_star and dominant detectability class from per-read data
    mean_p_star_agg  = {}
    dominant_det_cls = {}
    for locus in set(list(p_star_sums_agg.keys()) + list(allele_counts.keys())):
        mean_p_star_agg[locus]  = {}
        dominant_det_cls[locus] = {}
        for allele in allele_counts.get(locus, {}):
            n = p_star_counts_agg[locus].get(allele, 0)
            sm = p_star_sums_agg[locus].get(allele, 0.0)
            mean_p_star_agg[locus][allele] = sm / n if n > 0 else 0.0
            votes = det_class_votes_agg[locus].get(allele, {})
            dominant_det_cls[locus][allele] = (max(votes, key=votes.get)
                                               if votes else "Unobservable")

    return (
        {k: v for k, v in all_sequence_results.items()},
        {k: dict(v) for k, v in allele_counts.items()},
        {k: dict(v) for k, v in corrected_allele_counts.items()},
        stats_dict,
        dict(locus_sizes_agg),
        locus_repr_size,
        mean_p_star_agg,
        dominant_det_cls,
        dict(funnel_agg),
    )


# =============================================================================
# Per-sample runner
# =============================================================================

def run(sample_name, file_paths, primers, round_val, nbmismatch_max, threads,
        contig_setting, binning_setting, sequence_prefix,
        flanking_setting, output_dir, dico_bin, bbduk_path,
        save_matched_audit=False, max_size=DEFAULT_MAX_SIZE, read_types=None):
    """
    Process one sample end to end and RETURN its result.

    v2.5 changes
    ------------
    * Returns a dict instead of pushing tuples onto a Manager().Queue() with a
      SENTINEL.  The queue carried the (never-read) matched_reads structure
      through a socket, and its drain loop swallowed exceptions, so a sample
      could disappear without any error.  ProcessPoolExecutor already delivers
      exactly one result per submitted task.
    * Input sequencing metrics come from a full record-aware pass over the
      INPUT FASTQ (scan_input_fastq), together with a binned read-length
      histogram.  They are reported separately from BBDuk-candidate metrics.
    * The BBDuk outputs of R1 and R2 are no longer concatenated into an extra
      combined FASTQ; they are streamed into the same per-locus accumulators.
    * `save_matched_audit` no longer shadows the audit-writing function.
    """
    logger.info(f"Processing sample: {sample_name} with files: {file_paths}")
    start_time = time.time()
    try:
        # -- Input library metrics (full pass over each input FASTQ) ----------
        input_total_nt = 0
        input_read_count = 0
        input_min_len = None
        input_max_len = 0
        length_hist = {}
        matched_files_combined = defaultdict(list)

        # BBDuk consumes a file path, so each physical input is scanned once.
        # (Before v2.5 a logical chunk loop re-ran BBDuk on the whole file once
        # per chunk, duplicating every matched read for files > 1e6 reads.)
        for file_path in file_paths:
            info = scan_input_fastq(file_path)
            input_total_nt += info['total_nt']
            input_read_count += info['read_count']
            if info['read_count']:
                if input_min_len is None or info['min_len'] < input_min_len:
                    input_min_len = info['min_len']
                if info['max_len'] > input_max_len:
                    input_max_len = info['max_len']
            merge_length_hist(length_hist, info['hist'])

            matched_files = run_bbduk(
                file_path, output_dir, primers, nbmismatch_max, threads,
                bbduk_path, max_distance=max_size
            )
            for primer_name, matched_file in matched_files.items():
                matched_files_combined[primer_name].append(matched_file)

        matched_files_combined = dict(matched_files_combined)
        input_mean_len = (input_total_nt / input_read_count) if input_read_count else 0.0
        input_min_len = input_min_len or 0

        if save_matched_audit:
            save_matched_audit_copies(matched_files_combined, output_dir)
        else:
            logger.info("Skipping matched_fastq/matched_fasta audit copies (default).")

        (all_results, allele_counts, corrected_allele_counts,
         stats, locus_sizes_agg, locus_repr_size,
         mean_p_star_agg, dominant_det_cls, funnel) = process_bbduk_results(
            matched_files_combined, primers, nbmismatch_max, contig_setting,
            binning_setting, sequence_prefix, flanking_setting, output_dir,
            round_val, max_distance=max_size, dico_bin=dico_bin
        )

        for paths in matched_files_combined.values():
            for tmp in paths:
                try:
                    if os.path.exists(tmp):
                        os.remove(tmp)
                except Exception as e:
                    logger.warning(f"Could not delete temporary file {tmp}: {e}")

        processing_status = "OK" if matched_files_combined else "NO_CANDIDATE_READS"
        if not matched_files_combined:
            logger.warning(f"Sample {sample_name}: no BBDuk candidate reads for any "
                           f"locus; reported as NO_CANDIDATE_READS.")

        # v2.5: layout is decided by which read directions are present, not by
        # the number of files.  A paired sample split across two lanes has four
        # files and was previously mislabelled SE; a single-end sample split
        # into two parts was mislabelled PE.
        rts = sorted(set(read_types)) if read_types else []
        if rts:
            layout = "PE" if rts == ["R1", "R2"] else "SE"
        else:
            layout = "PE" if len(file_paths) == 2 else "SE"
        logger.info(f"Sample {sample_name} processing took {time.time() - start_time:.2f} seconds")

        funnel = dict(funnel)
        funnel.setdefault("BBDuk_Candidate_Reads", 0)
        return {
            'ok':                 True,
            'sample_name':        sample_name,
            'all_results':        all_results,
            'allele_counts':      allele_counts,
            'corrected_counts':   corrected_allele_counts,
            'input_read_count':   input_read_count,
            'input_total_bases':  input_total_nt,
            'input_mean_len':     input_mean_len,
            'input_min_len':      input_min_len,
            'input_max_len':      input_max_len,
            'length_hist':        length_hist,
            'candidate_stats':    stats,
            'layout':             layout,
            'locus_sizes':        locus_sizes_agg,
            'locus_repr_size':    locus_repr_size,
            'mean_p_star':        mean_p_star_agg,
            'dominant_det_cls':   dominant_det_cls,
            'funnel':             funnel,
            'input_files':        list(file_paths),
            'processing_status':  processing_status,
        }

    except Exception as e:
        logger.error(f"Error processing {sample_name}: {e}")
        traceback.print_exc()
        return {
            'ok':                False,
            'sample_name':       sample_name,
            'error':             f"{type(e).__name__}: {e}",
            'input_files':       list(file_paths),
            'processing_status': "FAILED",
        }


# =============================================================================
# Usage / argument parsing
# =============================================================================

def usage():
    print("Usage: python ReadTRAIL-Inference.py -i <input_dir> -o <output_dir> -p <primers_file>")
    print()
    print("Required:")
    print("  -i, --input         Directory with FASTQ (.fastq / .fq / .fastq.gz / .fq.gz) files.")
    print("  -o, --output        Directory to save output files.")
    print("  -p, --primer        Primer file (CSV/TSV: Locus_RepBP_ExpBP_RefU;FwdSeq;RevSeq).")
    print()
    print("Optional:")
    print("  -m, --mismatch INT  Max mismatches allowed per primer [default: 2].")
    print("  -t, --threads  INT  Threads per sample [default: 4].")
    print("  -b, --binning  FILE Path to binning file.")
    print("  --flanking-seq INT  Flanking sequence length [default: 0].")
    print(f"  --max-size     INT  Upper bound L_max on reconstructed fragment length [default: {DEFAULT_MAX_SIZE}].")
    print("  --round-val   FLOAT Allele discretization tolerance delta [default: 0.25].")
    print("  --p-min       FLOAT Detectability threshold p_min [default: 0.05].")
    print("  --obs-threshold FLOAT Observable class threshold [default: 0.5].")
    print("  --genome-size  INT  Genome size in bp. OPTIONAL; when given it adds the")
    print("                      diagnostic columns Genome_Coverage_Estimate and")
    print("                      Locus_Efficiency. It does not enter any classification.")
    print("  --bbduk-path  PATH  Path to bbduk.sh.")
    print(f"  --bbduk-xmx    STR  JVM heap per BBDuk process [default: {BBDUK_XMX}]. Total\n"
          "                      demand is this value times the number of parallel samples.")
    print("  --save-matched-reads  Save BBDuk candidate FASTQ/FASTA audit copies [default: off].")
    print("  --save-filtered-reads Save reads passing full primer-pair validation [default: off].")
    print("  -h, --help          Show this help message.")
    print()
    print("Notes:")
    print("  Input read count, mean/min/max read length and total bases are measured")
    print("  by a full record-aware pass over the input FASTQ files and are reported")
    print("  separately from the BBDuk candidate metrics.")
    print("  R1 and R2 are analysed as independent read-level observations; reads are")
    print("  not merged into pairs.")
    print("  Enabling either --save-* option must not change the analysis output.")
    print("  BBDuk is installed automatically by environment.yml.")


# =============================================================================
# Main
# =============================================================================

@timethis
def main():
    global dico_bin, flanking_len, P_MIN, OBS_THRESHOLD, MAX_WEIGHT, SAVE_FILTERED_READS, BBDUK_XMX

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
    genome_size      = None   # OPTIONAL in v2.5 — diagnostics only
    bbduk_path       = None
    save_matched     = False
    max_size         = DEFAULT_MAX_SIZE

    logger.info("Parsing command-line arguments")
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "hi:o:p:m:t:b:",
            ["help", "input=", "output=", "primer=",
             "mismatch=", "threads=", "binning=",
             "flanking-seq=", "genome-size=", "p-min=", "obs-threshold=",
             "bbduk-path=", "save-matched-reads", "save-filtered-reads",
             "max-size=", "round-val=", "bbduk-xmx="])
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
            logger.info(f"Genome size set to {genome_size} bp (diagnostics only)")
        elif opt == "--bbduk-xmx":
            BBDUK_XMX = str(arg)
            logger.info(f"BBDuk JVM heap set to {BBDUK_XMX}")
        elif opt == "--max-size":
            max_size = int(arg)
            logger.info(f"max fragment size L_max set to {max_size} bp")
        elif opt == "--round-val":
            round_val = float(arg)
            logger.info(f"allele discretization tolerance delta set to {round_val}")
        elif opt == "--p-min":
            P_MIN = float(arg)
            MAX_WEIGHT = 1.0 / P_MIN   # keep cap consistent with p_min
            logger.info(f"p_min set to {P_MIN}, MAX_WEIGHT updated to {MAX_WEIGHT:.2f}")
        elif opt == "--obs-threshold":
            OBS_THRESHOLD = float(arg)
            logger.info(f"obs_threshold set to {OBS_THRESHOLD}")
        elif opt == "--bbduk-path":
            bbduk_path = arg
        elif opt == "--save-matched-reads":
            save_matched = True
        elif opt == "--save-filtered-reads":
            SAVE_FILTERED_READS = True

    # Validate required arguments
    if not all([fasta_path, output_path, primer_file_path]):
        logger.error("Missing required arguments: -i, -o, -p are all required.")
        usage()
        sys.exit(2)

    # v2.5: --genome-size is optional.  It was validated as mandatory but never
    # used in any computation.  When supplied it now produces the diagnostic
    # columns Genome_Coverage_Estimate and Locus_Efficiency.
    if genome_size is not None and genome_size <= 0:
        logger.error(f"Invalid genome size: {genome_size}. Must be a positive integer.")
        sys.exit(2)

    if max_size <= 0:
        logger.error(f"Invalid --max-size: {max_size}. Must be a positive integer.")
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

    # =====================================================================
    # FASTQ discovery (v2.5)
    # ---------------------------------------------------------------------
    # The previous implementation lost input files in three different ways:
    #   * `.fq` / `.fq.gz` were not in the accepted extension list and were
    #     skipped without even a warning;
    #   * the name pattern only accepted an optional `_001` suffix, so
    #     `Sample_R1_002.fastq` matched nothing and was dropped;
    #   * a lane token was treated as part of the sample name, so
    #     `S_L001_R1_001` and `S_L002_R1_001` became two different samples;
    #   * `sample_files[sample][read] = path` silently overwrote a previous
    #     file mapped to the same slot.
    # Now every part/lane of a sample is collected into a list, and any FASTQ
    # that cannot be assigned unambiguously aborts the run.
    # =====================================================================
    FASTQ_SUFFIXES = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

    def strip_fastq_suffix(file_name):
        low = file_name.lower()
        for suf in sorted(FASTQ_SUFFIXES, key=len, reverse=True):
            if low.endswith(suf):
                return file_name[:-len(suf)], suf
        return None, None

    def parse_fastq_sample_read(base_name):
        """
        Return (sample_name, read_type, part_key) for a FASTQ base name.

        Recognised read tokens: _R1/_R2 and _1/_2, optionally followed by a
        part number (_001, _002, ...).  A lane token (_L001, _L002, ...)
        immediately preceding the read token belongs to the RUN, not to the
        sample, and is stripped so that all lanes merge into one sample.

        part_key orders the parts of one sample deterministically.
        """
        patterns = (
            r'^(?P<sample>.+?)_(?P<read>R[12])(?:_(?P<part>\d+))?$',
            r'^(?P<sample>.+?)_(?P<read>[12])(?:_(?P<part>\d+))?$',
        )
        for pattern in patterns:
            match = re.match(pattern, base_name)
            if not match:
                continue
            sample = match.group('sample')
            token = match.group('read')
            part = match.group('part') or '000'
            read_type = token if token.startswith('R') else f"R{token}"
            lane = ''
            lane_match = re.match(r'^(?P<stem>.+?)_(?P<lane>L\d{3})$', sample)
            if lane_match:
                sample = lane_match.group('stem')
                lane = lane_match.group('lane')
            return sample, read_type, (lane, part)
        return None, None, None

    def infer_read_suffix(read_types):
        rts = sorted(set(read_types))
        if rts == ['R1', 'R2']:
            return '_R1_R2'
        if rts == ['R1']:
            return '_R1'
        if rts == ['R2']:
            return '_R2'
        return '_SE'

    logger.info(f"Scanning input directory {fasta_path}")
    sample_files = defaultdict(lambda: defaultdict(list))
    unassigned = []
    non_fastq = []
    for file in sorted(os.listdir(fasta_path)):
        file_path = os.path.join(fasta_path, file)
        if os.path.isdir(file_path) or file.startswith('.'):
            continue
        base_name, suffix = strip_fastq_suffix(file)
        if suffix is None:
            non_fastq.append(file)
            continue
        sample_name, read_type, part_key = parse_fastq_sample_read(base_name)
        if sample_name and read_type:
            sample_files[sample_name][read_type].append((part_key, file_path))
        else:
            unassigned.append(file)

    if unassigned:
        logger.error(
            "The following FASTQ files could not be assigned unambiguously to a "
            "sample and read direction. Refusing to continue, because silently "
            "skipping them would break the guarantee that every input read is "
            "examined. Expected forms: Sample_R1_001.fastq(.gz), Sample_R1.fastq(.gz), "
            "Sample_1.fastq(.gz), optionally with a lane token such as "
            "Sample_L001_R1_001.fastq.gz:")
        for f in unassigned:
            logger.error(f"    {f}")
        sys.exit(2)
    if non_fastq:
        logger.info(f"Ignoring {len(non_fastq)} non-FASTQ file(s) in the input directory.")

    tasks = []
    sample_read_types = {}
    for sample_name in sorted(sample_files):
        by_read = sample_files[sample_name]
        ordered_paths = []
        for read_type in ('R1', 'R2'):
            for _key, path in sorted(by_read.get(read_type, [])):
                ordered_paths.append(path)
        if not ordered_paths:
            continue
        sample_read_types[sample_name] = sorted(by_read.keys())
        n_parts = len(ordered_paths)
        if n_parts > 2:
            logger.info(f"Sample {sample_name}: merging {n_parts} FASTQ parts/lanes "
                        f"({', '.join(sorted(by_read.keys()))}).")
        tasks.append((sample_name, ordered_paths))

    if not tasks:
        logger.error("No valid FASTQ files found in input directory.")
        sys.exit(2)

    total_input_files = sum(len(p) for _s, p in tasks)
    logger.info(f"Tasks to process: {len(tasks)} samples, {total_input_files} FASTQ files")
    for sample_name, paths in tasks:
        logger.info(f"  {sample_name}: {len(paths)} file(s) -> "
                    f"{', '.join(os.path.basename(p) for p in paths)}")

    max_processes, threads = adjust_resources(len(tasks), threads)
    logger.info(f"Starting processing with {max_processes} processes, {threads} threads each")

    # v2.5: results are collected from the futures themselves.  The previous
    # Manager().Queue() + SENTINEL drain loop shipped the unused matched_reads
    # structure through a socket and swallowed exceptions with a bare
    # `except: break`, so a sample could vanish without an error.
    all_outputs = []
    failed_samples = []
    with ProcessPoolExecutor(max_workers=max_processes) as executor:
        future_to_sample = {
            executor.submit(
                run, sample_name, file_paths, Primers, round_val, nb_mismatch, threads,
                False, binning_setting, "seq", flanking_setting,
                output_path, dico_bin, bbduk_path, save_matched, max_size,
                sample_read_types.get(sample_name, [])
            ): sample_name
            for sample_name, file_paths in tasks
        }
        for future in as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                result = future.result()
            except Exception as e:
                logger.error(f"Sample {sample_name} raised: {e}")
                traceback.print_exc()
                failed_samples.append(sample_name)
                continue
            if not result or not result.get('ok'):
                logger.error(f"Sample {sample_name} failed: "
                             f"{(result or {}).get('error', 'unknown error')}")
                failed_samples.append(sample_name)
                continue
            logger.info(f"Received data for sample: {sample_name}")
            all_outputs.append(result)

    if len(all_outputs) != len(tasks):
        logger.error(f"Expected {len(tasks)} results, got {len(all_outputs)}. "
                     f"Failed samples: {failed_samples}")

    if not all_outputs:
        logger.error("No results received.")
        sys.exit(1)

    all_outputs.sort(key=lambda r: r['sample_name'])

    # =====================================================================
    # Audit funnel (v2.5)
    # ---------------------------------------------------------------------
    # Written per sample so that read losses are explicit at every stage.  The
    # stages are nested, so the invariant is a chain of three equalities:
    #     BBDuk_Candidate_Reads          = Reads_Passing_Primer_Pair
    #                                      + Rejected_No_Pair
    #     Reads_Passing_Primer_Pair      = Reads_With_Valid_Fragment_Size
    #                                      + Rejected_Invalid_Size
    #     Reads_With_Valid_Fragment_Size = Reads_With_Valid_Allele_Call
    #                                      + Rejected_Invalid_Allele
    # Funnel_Balanced is True only when all three hold.  BBDuk candidate counts
    # are per (read x locus) observations, so a read matching several loci
    # contributes to several of them.
    # =====================================================================
    funnel_path = os.path.join(output_path, "audit_funnel.tsv")
    funnel_cols = ["Sample_ID", "Processing_Status", "Input_FASTQ_Files",
                   "Input_Read_Count", "Input_Total_Bases", "Input_Mean_Read_Length",
                   "Input_Min_Read_Length", "Input_Max_Read_Length",
                   "BBDuk_Candidate_Reads", "BBDuk_Candidate_Bases",
                   "BBDuk_Candidate_Mean_Read_Length",
                   "Reads_Passing_Primer_Pair", "Reads_With_Valid_Fragment_Size",
                   "Reads_With_Valid_Allele_Call", "Rejected_No_Pair",
                   "Rejected_Invalid_Size", "Rejected_Invalid_Allele",
                   "Reads_With_Multiple_Valid_Pairs", "Funnel_Balanced"]
    try:
        with open(funnel_path, 'w', newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(funnel_cols)
            for res in all_outputs:
                fn = res.get('funnel', {}) or {}
                cand = fn.get("BBDuk_Candidate_Reads", 0)
                pair = fn.get("Reads_Passing_Primer_Pair", 0)
                size = fn.get("Reads_With_Valid_Fragment_Size", 0)
                call = fn.get("Reads_With_Valid_Allele_Call", 0)
                balanced = (
                    cand == pair + fn.get("Rejected_No_Pair", 0)
                    and pair == size + fn.get("Rejected_Invalid_Size", 0)
                    and size == call + fn.get("Rejected_Invalid_Allele", 0)
                )
                cs = res.get('candidate_stats', {}) or {}
                writer.writerow([
                    res['sample_name'], res.get('processing_status', ''),
                    ";".join(os.path.basename(p) for p in res.get('input_files', [])),
                    res.get('input_read_count', 0), res.get('input_total_bases', 0),
                    round(res.get('input_mean_len', 0.0), 2),
                    res.get('input_min_len', 0), res.get('input_max_len', 0),
                    cand, cs.get('total_bases', 0), round(cs.get('mean_len', 0.0), 2),
                    fn.get("Reads_Passing_Primer_Pair", 0),
                    fn.get("Reads_With_Valid_Fragment_Size", 0),
                    fn.get("Reads_With_Valid_Allele_Call", 0),
                    fn.get("Rejected_No_Pair", 0), fn.get("Rejected_Invalid_Size", 0),
                    fn.get("Rejected_Invalid_Allele", 0),
                    fn.get("Reads_With_Multiple_Valid_Pairs", 0),
                    balanced,
                ])
                if not balanced:
                    logger.warning(f"Audit funnel does not balance for {res['sample_name']}: {fn}")
        logger.info(f"Audit funnel written: {funnel_path}")
    except Exception as e:
        logger.error(f"Could not write audit funnel: {e}")

    # =====================================================================
    # run_metadata.txt — exact command line and every effective parameter
    # =====================================================================
    meta_path = os.path.join(output_path, "run_metadata.txt")
    try:
        with open(meta_path, 'w', encoding='utf-8') as fh:
            fh.write("ReadTRAIL-Inference run metadata\n")
            fh.write("=" * 60 + "\n")
            fh.write(f"script                 : {os.path.abspath(__file__)}\n")
            fh.write(f"version                : Idealized Model v2.5\n")
            fh.write(f"timestamp              : {time.strftime('%Y-%m-%d %H:%M:%S %z')}\n")
            fh.write(f"command_line           : {' '.join(sys.argv)}\n")
            fh.write(f"working_directory      : {os.getcwd()}\n")
            fh.write(f"python                 : {sys.version.split()[0]}\n")
            fh.write(f"numpy                  : {np.__version__}\n")
            fh.write("-" * 60 + "\n")
            fh.write(f"input_dir              : {fasta_path}\n")
            fh.write(f"output_dir             : {output_path}\n")
            fh.write(f"primer_file            : {primer_file_path}\n")
            fh.write(f"binning_file           : {binning_file}\n")
            fh.write(f"bbduk_path             : {bbduk_path}\n")
            fh.write(f"bbduk_xmx              : {BBDUK_XMX}\n")
            fh.write(f"mismatch (m)           : {nb_mismatch}\n")
            fh.write(f"anchor_len             : {ANCHOR_LEN}\n")
            fh.write(f"max_size (L_max)       : {max_size}\n")
            fh.write(f"round_val (delta)      : {round_val}\n")
            fh.write(f"p_min                  : {P_MIN}\n")
            fh.write(f"obs_threshold          : {OBS_THRESHOLD}\n")
            fh.write(f"max_weight (1/p_min)   : {MAX_WEIGHT}\n")
            fh.write(f"flanking_len           : {flanking_len}\n")
            fh.write(f"genome_size            : {genome_size}\n")
            fh.write(f"chunk_size             : {CHUNK_SIZE}\n")
            fh.write(f"read_len_bin           : {READ_LEN_BIN}\n")
            fh.write(f"threads_per_sample     : {threads}\n")
            fh.write(f"max_processes          : {max_processes}\n")
            fh.write(f"save_matched_reads     : {save_matched}\n")
            fh.write(f"save_filtered_reads    : {SAVE_FILTERED_READS}\n")
            fh.write(f"split_wraparound       : disabled (linear reads only)\n")
            fh.write(f"paired_read_merging    : none (R1 and R2 are independent observations)\n")
            fh.write("-" * 60 + "\n")
            fh.write(f"samples_expected       : {len(tasks)}\n")
            fh.write(f"samples_succeeded      : {len(all_outputs)}\n")
            fh.write(f"samples_failed         : {failed_samples}\n")
            fh.write(f"input_fastq_files      : {total_input_files}\n")
            for sample_name, paths in tasks:
                fh.write(f"  {sample_name}: {', '.join(os.path.basename(p) for p in paths)}\n")
        logger.info(f"Run metadata written: {meta_path}")
    except Exception as e:
        logger.error(f"Could not write run metadata: {e}")


    # =========================================================================
    # Cohort-level baseline fallback
    # Collect all sample baselines; if a sample has too few loci, use cohort median.
    # =========================================================================
    cohort_baselines = []
    sample_baseline_info = {}  # sample_name → (baseline_val, loci_count, method)

    # v2.5: the baseline candidate criteria BC-1 (L < 0.8 * R) and BC-2
    # (p*(L, R) >= obs_threshold) are evaluated against the INPUT mean read
    # length.  Before v2.5 they used the mean over the BBDuk candidate subset,
    # which is enriched for the queried loci and was summed over every locus,
    # so it was not a property of the library at all.
    for res in all_outputs:
        sample_name   = res['sample_name']
        allele_counts = res['allele_counts'] or {}
        mean_read_len = res['input_mean_len']
        baseline_val, n_loci, method = compute_baseline(
            allele_counts, res['locus_repr_size'] or {}, mean_read_len)

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
    for res in all_outputs:
        allele_counts_entry = res['allele_counts']
        if not allele_counts_entry:
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
    #    They MUST NOT enter any calculation: N_filt, G_mini,
    #    Allele_Bases_Covered, C_mini, Q_i,
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

    for res in all_outputs:
        sample_name             = res['sample_name']
        all_results             = res['all_results']
        allele_counts           = res['allele_counts'] or {}
        corrected_allele_counts = res['corrected_counts'] or {}
        locus_sizes_agg         = res['locus_sizes'] or {}
        locus_repr_size         = res['locus_repr_size'] or {}
        mean_p_star_agg         = res['mean_p_star'] or {}
        dominant_det_cls        = res['dominant_det_cls'] or {}
        layout                  = res['layout']
        cand_stats              = res['candidate_stats'] or {}

        # Library metrics (input FASTQ) vs BBDuk candidate metrics
        mean_read_len   = res['input_mean_len'] if res['input_mean_len'] > 0 else 0
        input_read_cnt  = res['input_read_count']
        input_bases     = res['input_total_bases']
        input_min_len   = res['input_min_len']
        input_max_len   = res['input_max_len']
        length_hist     = res['length_hist'] or {}
        cand_reads      = cand_stats.get('total_reads', 0)
        cand_bases      = cand_stats.get('total_bases', 0)
        cand_mean_len   = cand_stats.get('mean_len', 0.0)
        genome_cov      = (input_bases / genome_size) if genome_size else None

        baseline_val, baseline_n, baseline_method = sample_baseline_info.get(
            sample_name, (None, 0, "unavailable"))

        suffix    = infer_read_suffix(sample_read_types.get(sample_name, []))
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
                    'input_read_count': input_read_cnt,
                    'input_total_bases': input_bases,
                    'mean_read_len': mean_read_len,
                    'input_min_len': input_min_len,
                    'input_max_len': input_max_len,
                    'length_hist': length_hist,
                    'cand_reads': cand_reads,
                    'cand_bases': cand_bases,
                    'cand_mean_len': cand_mean_len,
                    'genome_cov': genome_cov,
                    'processing_status': res.get('processing_status', 'OK'),
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
                'input_read_count': input_read_cnt,
                'input_total_bases': input_bases,
                'mean_read_len':    mean_read_len,
                'input_min_len':    input_min_len,
                'input_max_len':    input_max_len,
                'length_hist':      length_hist,
                'cand_reads':       cand_reads,
                'cand_bases':       cand_bases,
                'cand_mean_len':    cand_mean_len,
                'genome_cov':       genome_cov,
                'processing_status': res.get('processing_status', 'OK'),
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
    for res in all_outputs:
        for counts in (res['allele_counts'] or {}).values():
            max_alleles_detected = max(max_alleles_detected, len(counts))
        for counts in (res['corrected_counts'] or {}).values():
            max_alleles_detected = max(max_alleles_detected, len(counts))

    max_alleles_cohort = max(
        (len(v) for v in cohort_allele_info.values()), default=0)
    max_alleles = max(max_alleles_detected, max_alleles_cohort)

    # -- Build header (22 cols per allele) ------------------------------------
    # v2.5 adds two diagnostic per-allele columns:
    #   Population_Detectability_p — P(L) computed over the INPUT read-length
    #       histogram, i.e. E[max(0, R-L+1)] / E[R].  Detectability_p remains
    #       the published quantity (mean p* over the reads supporting the
    #       allele); the two differ because those reads are a length-biased
    #       sample, increasingly so for long fragments.  Both are reported so
    #       the effect can be quantified before any model is changed.
    #   Locus_Efficiency — Raw_Count / (Genome_Coverage_Estimate * P(L)),
    #       written only when --genome-size was supplied.  It is an
    #       assay-efficiency diagnostic independent of the sample's own loci,
    #       unlike B_sample which is estimated from the same data it normalises.
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
            f"Allele_Detection_Status {i}",
            f"Population_Detectability_p {i}",   # NEW v2.5
            f"Locus_Efficiency {i}",             # NEW v2.5 (needs --genome-size)
        ])

    # v2.5 column change.  The legacy trio Read_Count / Mean_Read_Length /
    # Total_Bases_Sequenced mixed two different populations: Read_Count came
    # from the input FASTQ, while the other two were computed over the BBDuk
    # candidate subset summed across every locus (2.2-12.0 % of Read_Count in
    # the validation cohort).  They are now reported separately and named for
    # what they are.  ReadTRAIL-Statistics reads Input_Mean_Read_Length and
    # falls back to Mean_Read_Length for files produced before v2.5.
    analysis_header = [
        "Access_number",
        "Sample_ID",
        "Primer",
        "Processing_Status",
        "Input_Read_Count",
        "Input_Total_Bases",
        "Input_Mean_Read_Length",
        "Input_Min_Read_Length",
        "Input_Max_Read_Length",
        "BBDuk_Candidate_Reads",
        "BBDuk_Candidate_Bases",
        "BBDuk_Candidate_Mean_Read_Length",
        "Genome_Coverage_Estimate",
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

    for res in all_outputs:
        sample_name = res['sample_name']

        for primer in Primers:
            primer_full  = primer[0]
            primer_short = primer_full.split('_')[0]
            key = (sample_name, primer_full)
            d   = pre_computed.get(key)

            # ── Empty locus ───────────────────────────────────────────────────
            if d is None or d.get('empty', True):
                # v2.5: a locus with no candidate reads is still reported, with
                # the sample's real input metrics and an explicit
                # Processing_Status, so that the QC denominator stays complete
                # (loci that produced nothing used to vanish from the
                # publication table entirely).
                acc_e  = d['accession']       if d else f"combined_{sample_name}.fastq"
                bl_e   = d['baseline_val']    if d else None
                bln_e  = d['baseline_n']      if d else 0
                blm_e  = d['baseline_method'] if d else "unavailable"
                lay_e  = d['layout']          if d else ""
                ps_e   = d.get('processing_status', 'OK') if d else 'NO_CANDIDATE_READS'
                gcov_e = d['genome_cov']      if d else None
                row = [acc_e, sample_name, primer_full, ps_e,
                       d['input_read_count']  if d else 0,
                       d['input_total_bases'] if d else 0,
                       round(d['mean_read_len'], 2) if d else 0,
                       d['input_min_len']     if d else 0,
                       d['input_max_len']     if d else 0,
                       d['cand_reads']        if d else 0,
                       d['cand_bases']        if d else 0,
                       round(d['cand_mean_len'], 2) if d else 0,
                       round(gcov_e, 4) if gcov_e is not None else "",
                       lay_e,
                       bl_e if bl_e is not None else "",
                       bln_e, blm_e,
                       "", "", "", "", "", "", "", "", "", ""]
                row.extend([""] * (22 * max_alleles))
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
            layout           = d['layout']
            baseline_val     = d['baseline_val']
            baseline_n       = d['baseline_n']
            baseline_method  = d['baseline_method']
            locus_factor     = d['locus_factor']
            input_read_cnt   = d['input_read_count']
            input_bases      = d['input_total_bases']
            input_min_len    = d['input_min_len']
            input_max_len    = d['input_max_len']
            length_hist      = d.get('length_hist') or {}
            cand_reads       = d['cand_reads']
            cand_bases       = d['cand_bases']
            cand_mean_len    = d['cand_mean_len']
            genome_cov       = d['genome_cov']
            proc_status      = d.get('processing_status', 'OK')

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
                proc_status,
                input_read_cnt, input_bases, round(mean_read_len, 2),
                input_min_len, input_max_len,
                cand_reads, cand_bases, round(cand_mean_len, 2),
                round(genome_cov, 4) if genome_cov is not None else "",
                layout,
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
                # v2.5 diagnostics
                pop_p = (population_detectability(ad['frag_size'], length_hist)
                         if ad.get('frag_size') else None)
                locus_eff = None
                if genome_cov and pop_p and pop_p > 0:
                    locus_eff = ad['raw_count'] / (genome_cov * pop_p)
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
                    ad.get('detection_status', 'Detected'),
                    round(pop_p, 6) if pop_p is not None else "",
                    round(locus_eff, 4) if locus_eff is not None else "",
                ])

            # Pad to max_alleles  (22 cols per allele)
            row.extend([""] * (22 * (max_alleles - len(allele_v2_data))))
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
