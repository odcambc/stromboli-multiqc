"""Pure parsing of STROMBOLI per-sample QC summary JSON.

Kept free of any MultiQC import so it can be unit-tested on its own. STROMBOLI writes
one summary file per sample (e.g. results/qc/{sample}.stromboli_qc.json); this module
turns the JSON text into a validated dict.

The expected shape (schema_version 5) is documented in examples/sample.stromboli_qc.json.
Single maintainer, so the "contract" is just this documented shape plus the version int:
bump SCHEMA_VERSION here and in STROMBOLI together when the format changes.
"""
import json

# v5: per-call read support — call_depth_counts (exact {cluster_depth: n_calls}, binned by
# the plugin like cluster sizes) and a median_call_depth scalar, from the cluster_depth column
# STROMBOLI now writes on each variant call. Calls below ~10 reads are noise-dominated, so the
# depth distribution is the FDR-relevant view of a run.
# v4: barcode/variant coverage analytics — ORF positional coverage profile + evenness
# (frac_orf_covered, coverage_gini), variant_types (SNV/MNV/indel), barcodes_per_variant_counts
# (redundancy), and barcode composition (length, GC).
# v3: add metrics latent in the producer's per-barcode tables — variant_consequences
# (DMS class breakdown), variants_per_barcode_counts, cluster_purity_histogram, and
# scalars (n_positions_mutated, n_indels, n_impure_clusters, flagged confident/ambiguous).
# v2: cluster sizes are emitted as exact {size: n_clusters} counts (cluster_size_counts)
# rather than producer-side coarse bins, so the plugin can choose data-driven bins across
# the whole cohort. Allele fraction stays pre-binned (bounded [0,1], not sample-dependent).
SCHEMA_VERSION = 5

# Every scalar metric in the schema. The plugin styles a headline subset of these in
# the general-statistics table (see SCALAR_HEADERS in stromboli.py); all are kept in the
# exported per-sample data file.
SCALAR_FIELDS = [
    "reads_total",
    "reads_with_barcode",
    "n_clusters",
    "n_clusters_passing",
    "median_cluster_size",
    "n_variants",
    "n_distinct_variants",
    "median_call_depth",
    "n_positions_mutated",
    "n_codons_covered",
    "orf_codons",
    "frac_orf_covered",
    "coverage_gini",
    "n_indels",
    "mean_barcode_gc",
    "n_flagged_merged",
    "n_flagged_mixed",
    "n_impure_clusters",
    "n_flagged_confident",
    "n_flagged_ambiguous",
]

# Distribution metrics rendered as plots. cluster_size_counts, call_depth_counts,
# variants_per_barcode_counts, barcodes_per_variant_counts and barcode_length_counts are exact
# {key: n} (the plugin bins the first four; length is plotted at integer resolution).
# allele_fraction_histogram, cluster_purity_histogram and barcode_gc_histogram are pre-binned.
# variant_consequences and variant_types are categorical {class: n}. orf_coverage_profile is
# {codon_window: n}.
HISTOGRAM_FIELDS = [
    "cluster_size_counts",
    "call_depth_counts",
    "allele_fraction_histogram",
    "variants_per_barcode_counts",
    "barcodes_per_variant_counts",
    "cluster_purity_histogram",
    "barcode_length_counts",
    "barcode_gc_histogram",
    "variant_consequences",
    "variant_types",
    "orf_coverage_profile",
]

# 1-2-5 x 10^k "nice" edges give ~3 bins per order of magnitude — enough resolution to
# read a skewed cluster-size distribution without committing to a fixed, sample-blind
# bin scheme. Swap this sequence (e.g. powers of 2, or plain decades) to retune the bins.
_NICE_STEPS = (1, 2, 5)


def _nice_edges(max_size):
    """Ascending 1-2-5 x 10^k edges, returned once one strictly exceeds max_size.

    Always at least two edges, so there is always at least one bin.
    """
    edges = []
    k = 0
    while True:
        for step in _NICE_STEPS:
            edge = step * 10 ** k
            edges.append(edge)
            if edge > max_size:
                return edges
        k += 1


def bin_cluster_sizes(counts_by_sample):
    """Bin per-sample exact cluster-size counts into shared, data-driven bins.

    counts_by_sample: {sample: {size(int): n_clusters(int)}}.

    Returns (labels, {sample: {label: n_clusters}}). The bins are derived from the single
    largest cluster seen across *all* samples, so every sample's series shares one ordered
    x-axis (required for a meaningful overlay) yet the axis still scales to the data: a run
    topping out at 12 reads and one topping out at 3000 get appropriately sized bins. The
    final bin is open-ended ("500+") to absorb the long tail.
    """
    max_size = max(
        (size for counts in counts_by_sample.values() for size in counts),
        default=1,
    )
    edges = _nice_edges(max(max_size, 1))

    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i == len(edges) - 2:  # last bin: open-ended tail
            labels.append("{}+".format(lo))
        elif lo == hi - 1:  # single integer
            labels.append(str(lo))
        else:
            labels.append("{}-{}".format(lo, hi - 1))

    def _label_for(size):
        for i in range(len(edges) - 1):
            if i == len(edges) - 2 or size < edges[i + 1]:
                return labels[i]
        return labels[-1]

    binned = {}
    for sample, counts in counts_by_sample.items():
        row = {label: 0 for label in labels}
        for size, n in counts.items():
            row[_label_for(size)] += n
        binned[sample] = row
    return labels, binned


class SchemaVersionWarning(UserWarning):
    """Raised (as a warning) when a summary file's schema_version is unexpected."""


def parse_stromboli_qc(text):
    """Parse summary JSON text into a dict. Returns (sample_name, data).

    Raises ValueError on malformed JSON or a missing sample name. An unexpected
    schema_version is surfaced via a warning, not an error (forward-compatible-ish).
    """
    import warnings

    data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("STROMBOLI QC summary must be a JSON object")

    version = data.get("schema_version")
    if version != SCHEMA_VERSION:
        warnings.warn(
            "STROMBOLI QC schema_version {} != expected {}; parsing best-effort".format(
                version, SCHEMA_VERSION
            ),
            SchemaVersionWarning,
        )

    sample = data.get("sample")
    if not sample:
        raise ValueError("STROMBOLI QC summary is missing the 'sample' field")

    return sample, data
