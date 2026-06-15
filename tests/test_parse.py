"""Unit tests for the pure parser (no MultiQC needed). Run with: uv run pytest"""
import json
import os
import warnings

import pytest

from multiqc_stromboli.parse import (
    SCHEMA_VERSION,
    SchemaVersionWarning,
    bin_cluster_sizes,
    parse_stromboli_qc,
)

EXAMPLE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "examples", "sample.stromboli_qc.json"
)


def test_parses_example():
    sample, data = parse_stromboli_qc(open(EXAMPLE).read())
    assert sample == "EXP-NBD114_barcode23.subsample"
    assert data["schema_version"] == SCHEMA_VERSION
    assert data["n_flagged_mixed"] == 3


def test_missing_sample_raises():
    with pytest.raises(ValueError):
        parse_stromboli_qc(json.dumps({"schema_version": SCHEMA_VERSION}))


def test_unexpected_version_warns():
    with pytest.warns(SchemaVersionWarning):
        parse_stromboli_qc(json.dumps({"schema_version": 999, "sample": "s"}))


def test_bins_scale_to_largest_cluster():
    # A small-max cohort gets few bins; a thousands-max cohort gets more, reaching up.
    small_labels, _ = bin_cluster_sizes({"a": {1: 5, 2: 3, 9: 1}})
    big_labels, _ = bin_cluster_sizes({"a": {1: 5, 3000: 1}})
    assert len(big_labels) > len(small_labels)
    assert big_labels[-1].endswith("+")  # open-ended tail bin
    assert "2400" not in "".join(small_labels)  # small cohort never reaches the thousands


def test_bins_are_shared_and_total_preserving():
    # Two samples with different maxima share one ordered label list, and binning loses
    # no clusters (every count lands in exactly one bin).
    counts = {"a": {1: 10, 4: 5, 500: 2}, "b": {2: 7, 50: 3}}
    labels, binned = bin_cluster_sizes(counts)
    assert set(binned) == {"a", "b"}
    for sample, row in binned.items():
        assert list(row) == labels  # same bins, same order, for every sample
        assert sum(row.values()) == sum(counts[sample].values())


def test_size_lands_in_expected_bin():
    labels, binned = bin_cluster_sizes({"a": {1: 1, 4: 1, 12: 1}})
    # 1-2-5-10-20 edges -> labels "1", "2-4", "5-9", "10+"
    assert labels == ["1", "2-4", "5-9", "10+"]
    assert binned["a"] == {"1": 1, "2-4": 1, "5-9": 0, "10+": 1}
