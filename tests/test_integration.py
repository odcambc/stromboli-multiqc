"""End-to-end smoke test: run the real `multiqc` CLI over example QC JSONs and confirm
the STROMBOLI module discovers them and writes its data into the report.

Unlike test_parse.py (pure parser, no MultiQC), this exercises the whole plugin path —
the before_config search-pattern hook, file discovery, and the module's data export. It
is skipped automatically when the `multiqc` binary is not on PATH, so the parser tests
still run in a bare environment.
"""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "examples" / "sample.stromboli_qc.json"

pytestmark = pytest.mark.skipif(
    shutil.which("multiqc") is None, reason="multiqc not installed"
)


# Pin the report name so the data dir is a deterministic "<name>_data".
REPORT_NAME = "report"
DATA_DIR = REPORT_NAME + "_data"

# Every plot section the module renders for a schema_version 5 summary. Asserting the full
# set guards against a section silently dropping out (e.g. a renamed field or a plot that
# errors and is skipped). Keep in step with MultiqcModule's add_section() ids.
EXPECTED_PLOT_IDS = [
    "stromboli_consequences",
    "stromboli_variant_types",
    "stromboli_orf_coverage",
    "stromboli_clashes",
    "stromboli_cluster_sizes",
    "stromboli_call_depth",
    "stromboli_allele_fractions",
    "stromboli_variants_per_barcode",
    "stromboli_barcodes_per_variant",
    "stromboli_cluster_purity",
    "stromboli_barcode_length",
    "stromboli_barcode_gc",
]


def _run_multiqc(input_dir, out_dir):
    subprocess.run(
        ["multiqc", str(input_dir), "-o", str(out_dir),
         "-n", REPORT_NAME, "--force", "--no-ansi"],
        check=True,
        capture_output=True,
        text=True,
    )


def test_renders_single_sample(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    shutil.copy(EXAMPLE, in_dir)
    out_dir = tmp_path / "out"

    _run_multiqc(in_dir, out_dir)

    exported = (out_dir / DATA_DIR / "multiqc_stromboli.txt").read_text()
    assert "EXP-NBD114_barcode23.subsample" in exported
    assert "8210" in exported  # reads_with_barcode made it through
    # A v3-only field reached the exported data (the new metrics aren't dropped).
    assert "n_positions_mutated" in exported

    # The headline scalar reached the shared general-stats table.
    general = (out_dir / DATA_DIR / "multiqc_general_stats.txt").read_text()
    assert "reads_with_barcode" in general

    # Every v3 plot section was rendered into the report.
    report = (out_dir / (REPORT_NAME + ".html")).read_text()
    missing = [pid for pid in EXPECTED_PLOT_IDS if pid not in report]
    assert not missing, "missing plot sections: {}".format(missing)


def test_renders_multiple_samples(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    shutil.copy(EXAMPLE, in_dir)

    # Synthesize a second sample so the multi-sample discovery path is covered.
    second = json.loads(EXAMPLE.read_text())
    second["sample"] = "EXP-NBD114_barcode24.subsample"
    second["reads_with_barcode"] = 9100
    (in_dir / "barcode24.stromboli_qc.json").write_text(json.dumps(second))

    out_dir = tmp_path / "out"
    _run_multiqc(in_dir, out_dir)

    exported = (out_dir / DATA_DIR / "multiqc_stromboli.txt").read_text()
    assert "EXP-NBD114_barcode23.subsample" in exported
    assert "EXP-NBD114_barcode24.subsample" in exported
