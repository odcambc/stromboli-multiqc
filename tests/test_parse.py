"""Unit tests for the pure parser (no MultiQC needed). Run with: uv run pytest"""
import json
import os
import warnings

import pytest

from multiqc_stromboli.parse import (
    SCHEMA_VERSION,
    SchemaVersionWarning,
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
