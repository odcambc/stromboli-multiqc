"""Pure parsing of STROMBOLI per-sample QC summary JSON.

Kept free of any MultiQC import so it can be unit-tested on its own. STROMBOLI writes
one summary file per sample (e.g. results/qc/{sample}.stromboli_qc.json); this module
turns the JSON text into a validated dict.

The expected shape (schema_version 1) is documented in examples/sample.stromboli_qc.json.
Single maintainer, so the "contract" is just this documented shape plus the version int:
bump SCHEMA_VERSION here and in STROMBOLI together when the format changes.
"""
import json

SCHEMA_VERSION = 1

# Scalar metrics surfaced in the MultiQC general-statistics table.
SCALAR_FIELDS = [
    "reads_total",
    "reads_with_barcode",
    "n_clusters",
    "n_clusters_passing",
    "median_cluster_size",
    "n_variants",
    "n_flagged_merged",
    "n_flagged_mixed",
]

# Histogram metrics (ordered-bin dicts) rendered as plots.
HISTOGRAM_FIELDS = ["cluster_size_histogram", "allele_fraction_histogram"]


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
